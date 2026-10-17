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
})()`,context).then(()=>console.log('PASS')).catch(e=>{console.error(e);process.exitCode=1;});
'''
    source = Path("src/sarna/app/static/offtarget.js").resolve()
    result = subprocess.run([node, "-e", script, str(source)], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert "PASS" in result.stdout
