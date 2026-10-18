"""Regression: widespread sequence matches must not erase the design-ranked top list."""

import shutil
import subprocess
from pathlib import Path

import pytest


def test_advisory_findings_selection_and_plot_coverage() -> None:
    node = shutil.which("node")
    if not node:
        pytest.skip("Node is required for frontend state checks")
    script = r'''
const fs=require('fs'),vm=require('vm'),assert=require('assert');
const elements=new Map(),listeners=new Map();
const context={console,frozen:null,state:{result:null,checked:new Set()},
 escapeHTML:String,fmt:String,signed:String,error:()=>{},setTimeout,clearTimeout,
 renderTable:()=>{},document:{hidden:false,addEventListener:()=>{}},
 $:id=>{if(!elements.has(id))elements.set(id,{innerHTML:'',
 addEventListener:(kind,fn)=>listeners.set(id+':'+kind,fn)});return elements.get(id);}};
vm.createContext(context);vm.runInContext(fs.readFileSync(process.argv[1],'utf8'),context);
vm.runInContext(`
state.result={app:{run_id:'run'},candidates:Array.from({length:14},(_,i)=>
 ({id:'c'+i,rank:i+1,rating:{met:6,evaluated:7}}))};
state.screenings={};
for(const scope of ['genomic','transcript']){
 const entries=Object.fromEntries(state.result.candidates.map(c=>[c.id,{
 row:{candidate:c.id,identity_counts:{'19':0,'18':2,'17':30,'16':500},
 enumeration_complete:true,alignments_capped:false,
 selection_assessment:{status:'excluded',policy:{id:'conservative-full-length-v1'},
 reasons:['2 full-length matches at 18/19.']}},search:{id:scope}}]));
 state.screenings[scope]={byEvidence:{'local:':entries}};
}
renderScreenedShortlist();
`,context);
assert(context.$('screened-shortlist').innerHTML.includes('Review findings'));
assert(!context.$('screened-shortlist').innerHTML.includes('No candidates meet'));
listeners.get('screened-shortlist:click')({target:{id:'screened-use',closest:()=>null}});
assert.deepEqual([...context.state.checked],Array.from({length:10},(_,i)=>'c'+i));
vm.runInContext(`delete state.screenings.genomic.byEvidence['local:'].c0;
renderScreenedShortlist();validateScreenedExport();`,context);
assert(context.state.checked.has('c0'),'Missing evidence removed an export selection');
assert(context.$('screened-shortlist').innerHTML.includes('Not available'));
assert(context.$('top-offtarget-summary').innerHTML.includes('evidence incomplete'));
vm.runInContext(`state.screeningScope='genomic';state.screenings={};`,context);
const html=vm.runInContext(`tssPlot([{row:{hit_count:30243,details_retained:2800,candidate:'c0'},
 placements:Array.from({length:744},()=>({distance_to_tss:-100,
 hit:{identities:18,alignment_length:19},contig:'1',position:1000}))}])`,context);
assert(html.includes('744 of 30243 reported genomic alignments (2.5%)'));
assert(html.includes('not the fraction of all alignments near a TSS'));
assert(html.includes('may belong to another gene'));
assert(html.startsWith('<details'),'Plot should be available as a disclosure');
const empty=vm.runInContext(`tssPlot([{row:{hit_count:100,details_retained:20},placements:[]}])`,context);
assert(empty.includes('(0%)'));
assert(empty.includes('Other alignments may also lie near a start site'));
console.log('PASS');
'''
    result = subprocess.run(
        [node, "-e", script, str(Path("src/sarna/app/static/offtarget.js").resolve())],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
