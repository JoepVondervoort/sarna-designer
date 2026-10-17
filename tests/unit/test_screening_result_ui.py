"""Display units, filtering and efficiency must not change scientific findings."""

import shutil
import subprocess
from pathlib import Path

import pytest


def test_result_units_filters_limits_caches_and_annotation_deduplication() -> None:
    node = shutil.which("node")
    if not node:
        pytest.skip("Node is required for frontend state checks")
    script = r'''
const fs=require('fs'),vm=require('vm'),elements=new Map();
const context={console,frozen:null,state:{result:null,checked:new Set()},
  escapeHTML:String,fmt:String,signed:String,error:()=>{},setTimeout,clearTimeout,
  document:{hidden:false,addEventListener:()=>{}},
  $:id=>{if(!elements.has(id))elements.set(id,{
    textContent:'',innerHTML:'',addEventListener:()=>{}});return elements.get(id);}};
vm.createContext(context);vm.runInContext(fs.readFileSync(process.argv[1],'utf8'),context);
vm.runInContext(`(async()=>{
  const assert=(ok,message)=>{if(!ok)throw Error(message);};
  state.result={app:{run_id:'run'},candidates:Array.from({length:14},(_,i)=>
    ({id:'candidate-'+i,rank:i+1,target_site_dna:'ACGTACGTACGTACGTACG'}))};
  state.screeningScope='transcript';state.screeningOpen='candidate-0';
  const hits=Array.from({length:100},(_,i)=>({candidate:'candidate-0',accession:'NM_'+i,
    title:i===0?'SGPL1':'AOPEP',identities:i===0?16:15,alignment_length:19,gaps:0,
    subject_start:i+1,subject_end:i+19,subject_strand:'plus',query_start:1,query_end:19,
    evalue:1}));
  const entry={row:{candidate:'candidate-0',hit_count:100,wang_homology_hits:1,
    strong_threshold:16,perfect_hits:0,alignments_capped:true},hits,placements:[],
    search:{id:'search',parameters:{max_hits:100,max_alignments_per_candidate:200},
    annotations:{},searched_at:'2026-10-17',notes:[]}};
  const held={byCandidate:{'candidate-0':entry},searches:[entry.search]};
  state.screenings={transcript:held};
  const original=JSON.stringify(entry.hits);
  const summary=summaryLine(held,[entry],19,false);
  assert(summary.includes('saRNA candidates searched'),'Candidate unit missing');
  assert(summary.includes('alignments returned'),'Returned-alignment unit missing');
  assert(summary.includes('alignments with ≥16/19'),'Threshold-alignment unit missing');
  assert(screeningLimit(entry).includes('maximum 100 database sequences'),'Wrong server cap');
  const groups=groupHits(entry);assert(groups===groupHits(entry),'Groups were recomputed');
  state.screeningMatchFilter='threshold';
  assert(visibleHitGroups(entry,19).length===1,'Threshold filter counted weaker alignments');
  const details=alignmentDetail(entry,19,false);
  assert(details.includes('SGPL1')&&!details.includes('AOPEP'),'Weak matches leak into filter');
  assert(details.includes('99 weaker alignments'),'Weaker matches were hidden without explanation');
  state.screeningMatchFilter='all';
  assert(visibleHitGroups(entry,19).length===100,'All-returned view discarded matches');
  assert(JSON.stringify(entry.hits)===original,'Display filter changed scientific hits');
  assert(entry.row.wang_homology_hits===1,'Display filter changed threshold count');
  assert(localScreeningDelay({progress:{phase:'waiting'}})===15000,'Waiting polls too fast');
  assert(localScreeningDelay({progress:{phase:'placing'}})===5000,'Local progress too slow');
  document.hidden=true;
  assert(localScreeningDelay({progress:{phase:'placing'}})===30000,'Hidden page polls too fast');
  document.hidden=false;
  let writes=0;const target={set innerHTML(value){writes++;}};
  screeningHTML(target,'same');screeningHTML(target,'same');
  assert(writes===1,'Unchanged status destroyed existing DOM');
  let finish,calls=0;
  api=async()=>{calls++;return new Promise(resolve=>{finish=resolve;});};
  const first=loadHitAnnotations(entry,['NM_0']);
  await loadHitAnnotations(entry,['NM_0']);
  assert(calls===1,'Overlapping expansion fetched the same annotation twice');
  finish({'NM_0':{status:'available',chromosome:'10'}});await first;
  assert(SCREENING_ANNOTATIONS.size===0,'Completed annotation reservation leaked');
  state.screeningExpandAll=true;
  const many=loadHitAnnotations(entry,hits.map(h=>h.accession));
  assert(calls===2,'Visible next batch did not start');
  state.screeningOpen=null;finish({});await many;
  assert(calls===2,'Closed candidate continued fetching invisible annotations');
  state.screeningJobs={};
  let fetches=0,complete;
  api=()=>{fetches++;return new Promise(resolve=>{complete=resolve;});};
  const pollOne=pollScreening('job','run','transcript','test',0,1);
  const pollTwo=pollScreening('job','run','transcript','test',0,1);
  assert(pollOne===pollTwo&&fetches===1,'Same job started duplicate polling loops');
  complete({id:'job',scope:'transcript',status:'waiting',stage:'Report retrieval needs retry'});
  await pollOne;
  assert(SCREENING_POLLS.size===0,'Finished polling reservation leaked');
  assert(SCREENING_ANNOTATIONS.size===0,'Cancelled annotation reservation leaked');
  const clientCapped={...entry,row:{...entry.row,hit_count:200},
    hits:Array.from({length:200},()=>({...hits[0]}))};
  assert(screeningLimit(clientCapped).includes('app retained at most 200'),
    'Client truncation was described as an NCBI sequence limit');

})()`,context).then(()=>console.log('PASS')).catch(e=>{console.error(e);process.exitCode=1;});
'''
    source = Path("src/sarna/app/static/offtarget.js").resolve()
    result = subprocess.run([node, "-e", script, str(source)], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert "PASS" in result.stdout
