"""The open run follows saved NCBI jobs without duplicates or stale navigation callbacks."""

import shutil
import subprocess
from pathlib import Path

import pytest


def test_waiting_jobs_are_collected_automatically_with_safe_timer_lifecycle() -> None:
    node = shutil.which("node")
    if not node:
        pytest.skip("Node is required for frontend timer checks")
    script = r'''
const fs=require('fs'), vm=require('vm');
const elements=new Map(), timers=new Map();let nextId=0;
const context={console, frozen:null,state:{result:null,checked:new Set()},
  escapeHTML:String, signed:String,fmt:String,error:()=>{},
  setTimeout:(callback,delay)=>{const id=++nextId;timers.set(id,{callback,delay});return id;},
  clearTimeout:id=>timers.delete(id), document:{addEventListener:()=>{}},
  $:id=>{
    if(!elements.has(id)) elements.set(id,{addEventListener:()=>{},textContent:'',innerHTML:''});
    return elements.get(id);
  },timers};
vm.createContext(context);vm.runInContext(fs.readFileSync(process.argv[1],'utf8'),context);
vm.runInContext(`
(async()=>{
  const assert=(test,message)=>{if(!test)throw Error(message);};
  const reset=()=>{
    timers.clear();AUTO_SCREEN_CHECKS.clear();
    state.result={app:{run_id:'run-1'},candidates:[{id:'c1',rank:1,target_site_dna:'ACGT'}]};
    state.screenings={};state.screeningJobs={transcript:{status:'waiting',id:'job-1'}};
  };
  reset();let posts=0;
  api=async(path,options)=>{
    if(options){assert(path==='/api/screenings/job-1/resume','Resubmitted a new search');posts++;
      assert(options.headers['Content-Type']==='application/json','Missing JSON media type');
      assert(options.body==='{}','Missing JSON request body');
      return {id:'job-1'};}
    return {id:'job-1',scope:'transcript',status:'complete',result:{specificity:[]}};
  };
  renderScreeningStatus();renderScreeningStatus();
  assert(timers.size===1,'Repeated renders created duplicate timers');
  const pending=[...timers.values()][0];
  assert(pending.delay===60000,'Automatic checks are too frequent');
  pending.callback();
  for(let i=0;i<10;i++)await Promise.resolve();
  assert(posts===1,'Waiting job was not collected once');
  assert(state.screeningJobs.transcript.status==='complete','Completion not displayed');
  reset();scheduleScreeningCheck('transcript');
  const staleRun=[...timers.values()][0];
  state.result.app.run_id='run-2';staleRun.callback();
  for(let i=0;i<3;i++)await Promise.resolve();
  assert(posts===1,'An old timer followed navigation to a new run');
  reset();scheduleScreeningCheck('transcript');
  const staleJob=[...timers.values()][0];
  state.screeningJobs.transcript.id='job-2';staleJob.callback();
  assert(posts===1,'An old timer collected a replaced job');
  reset();state.screeningJobs.transcript.retrievalIssue=true;
  renderScreeningStatus();
  assert(timers.size===0,'Unreadable reports entered an endless automatic retry loop');
  reset();scheduleScreeningCheck('transcript');
  await resumeScreening('transcript');
  assert(timers.size===0,'Manual collection left an automatic retry queued');
  assert(posts===2,'Manual collection did not use the saved search');
  reset();api=async(path,options)=>{
    if(options) throw Object.assign(new Error('Already collecting'),
      {code:'screening-already-collecting'});
    return {id:'job-1',scope:'transcript',status:'complete',result:{specificity:[]}};
  };
  await resumeScreening('transcript');
  assert(state.screeningJobs.transcript.status==='complete','Another tab looked like a failure');
  assert(timers.size===0,'Completed search kept an automatic retry');
  reset();api=async(path,options)=>{
    if(options) throw Object.assign(new Error('Already finished'),
      {code:'screening-already-complete'});
    return {id:'job-1',scope:'transcript',status:'complete',result:{specificity:[]}};
  };
  await resumeScreening('transcript');
  assert(state.screeningJobs.transcript.status==='complete','Another tab completion looked failed');
  // Migration: a previous NCBI queue and completed NCBI records cannot block local work.
  reset();timers.clear();AUTO_SCREEN_CHECKS.clear();
  $('screening-backend').value='local';
  REFERENCE_STATUS={ready:true,installed:{release:'R1'},catalog_release:'R1'};
  state.result.candidates=[{id:'c1',rank:1},{id:'c2',rank:2}];
  state.screeningJobs={genomic:{status:'waiting',id:'old-ncbi',backend:'ncbi'}};
  const row=id=>({candidate:id,hit_count:0,perfect_hits:0,strong_threshold:16,
    wang_homology_hits:0,near_perfect_hits:0});
  const old={id:'old-record',scope:'genomic',result:{backend:'ncbi',
    searched_at:'2026-01-01',specificity:[row('c1'),row('c2')]}};
  mergeScreening(old);
  assert(alreadyScreened('genomic').size===0,'NCBI findings counted as local coverage');
  assert(currentScreeningJob('genomic')===null,'Old NCBI queue blocks local work');
  renderScreeningStatus();
  assert(!$('screen-genomic').disabled,'Local search disabled by old NCBI queue');
  assert(timers.size===0,'Local mode resumed an old NCBI request');
  assert(state.screeningJobs.genomic.id==='old-ncbi','Switching methods deleted the old ID');
  mergeScreening({id:'local-record',scope:'genomic',result:{backend:'local',
    reference:{release:'R1'},searched_at:'2026-02-01',specificity:[row('c1')]}});
  assert(alreadyScreened('genomic').size===1,'Completed local coverage missing');
  mergeScreening({...old,id:'later-ncbi',result:{...old.result,searched_at:'2026-03-01'}});
  assert(alreadyScreened('genomic').size===1,'Later NCBI result erased local coverage');
  REFERENCE_STATUS.installed.release='R2';
  assert(alreadyScreened('genomic').size===0,'Earlier reference counted as current coverage');
  $('screening-backend').value='ncbi';
  assert(alreadyScreened('genomic').size===2,'Saved NCBI coverage was lost');
  assert(localScreeningDelay({result:{backend:'local'},poll_after_seconds:1})===1000,
    'Local result collection is unnecessarily slow');
  assert(localScreeningDelay({result:{backend:'local'},poll_after_seconds:10})===10000,
    'Local polling ignored server advice');
  // An in-flight old poll must not overwrite a new method's visible job.
  let releaseOld;
  api=async()=>new Promise(resolve=>{releaseOld=resolve;});
  const oldPoll=pollScreening('stale-ncbi','run-1','genomic','old',0,1,'ncbi');
  $('screening-backend').value='local';
  state.screeningJobs.genomic={backend:'local',status:'running',id:'new-local'};
  releaseOld({id:'stale-ncbi',scope:'genomic',status:'complete',result:old.result});
  assert(await oldPoll==='changed-backend','An old poll followed an engine switch');
  assert(state.screeningJobs.genomic.id==='new-local','Old poll overwrote local job');

})()`,context).then(()=>console.log('PASS')).catch(e=>{console.error(e);process.exitCode=1;});
'''
    source = Path("src/sarna/app/static/offtarget.js").resolve()
    result = subprocess.run([node, "-e", script, str(source)], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert "PASS" in result.stdout
