"""Plots must preserve missing evidence, query orientation and saved rule outcomes."""

import shutil
import subprocess
from pathlib import Path

import pytest


def test_comparisons_use_saved_counts_verdicts_and_valid_mismatch_positions() -> None:
    node = shutil.which("node")
    if not node:
        pytest.skip("Node is required for frontend state checks")
    script = r'''
const fs=require('fs'),vm=require('vm'),assert=require('assert');
const context={console,frozen:null,state:{},escapeHTML:String,fmt:String,signed:String,
 setTimeout,clearTimeout,error:()=>{},
 document:{addEventListener:()=>{}},$:()=>({addEventListener:()=>{}})};
vm.createContext(context);vm.runInContext(fs.readFileSync(process.argv[1],'utf8'),context);
vm.runInContext(`
state.result={configuration:{rules:[{id:'p1',name:'GC preference',kind:'preference'},
 {id:'p2',name:'Flank preference',kind:'preference'},{id:'p3',name:'Missing',kind:'preference'}]}};
const candidate={id:'c1',rank:1,target_site_dna:'ACGTACGTACGTACGTACG',
 rating:{percent:100},rules:[{id:'p1',kind:'preference',outcome:'pass',detail:'Measured pass'},
 {id:'p2',kind:'preference',outcome:'fail',detail:'Measured failure'}]};
const valid={candidate:'c1',accession:'NC_1',identities:18,alignment_length:19,gaps:0,
 query_start:1,query_end:19,subject_start:100,subject_end:118,subject_strand:'minus',
 has_alignment_detail:true,mismatch_positions:[2]};
const entry={row:{hit_count:1000,enumeration_complete:true,alignments_capped:false,
 identity_counts:{'19':0,'18':2,'17':10}},hits:[valid,
 {...valid,query_start:2}, {...valid,gaps:1}, {...valid,has_alignment_detail:false},
 {...valid,mismatch_positions:[]},{...valid,mismatch_positions:[20]},
 {...valid,identities:17,mismatch_positions:[2,2]}]};
const decision={candidate,evidence:[{scope:'genomic',entry},{scope:'transcript'}]};
state.testDecision=decision;state.testEntry=entry;
`,context);
const counts=vm.runInContext('matchCountPlot([state.testDecision])',context);
assert(counts.includes('linear scale from 0 to 12'));
assert(counts.includes('0 / 2 / 10'));
assert(counts.includes('No identity counts available'));
assert(!counts.includes('NaN'));
const heatmap=vm.runInContext('preferenceHeatmap([state.testDecision])',context);
assert.equal((heatmap.match(/criterion-cell pass/g)||[]).length,1);
assert.equal((heatmap.match(/criterion-cell fail/g)||[]).length,1);
assert.equal((heatmap.match(/criterion-cell unknown/g)||[]).length,1);
assert(heatmap.includes('Measured failure'));
assert(heatmap.includes('Not evaluated'));
vm.runInContext(`state.testDecision.candidate.preference_verdicts=state.testDecision.candidate.rules;
delete state.testDecision.candidate.rules;`,context);
assert.equal(vm.runInContext('preferenceHeatmap([state.testDecision])',context),heatmap,
 'Compact API verdicts must give the same heatmap as full archival verdicts');
assert.equal(vm.runInContext('mismatchMapData(state.testEntry,19).usable.length',context),1);
let map=vm.runInContext('mismatchPositionPlot([state.testDecision])',context);
assert(map.includes('1 alignments shown out of 7 retained details and 1000 reported'));
assert(map.includes('target position 2: mismatch'));
assert(map.includes('target position 18: match'),'Minus strand must not reverse query positions');
assert(map.includes('not the antisense guide orientation'));
vm.runInContext(`state.mismatchPlotScope='transcript'`,context);
map=vm.runInContext('mismatchPositionPlot([state.testDecision])',context);
assert(map.includes('No screening evidence is available'));
assert(map.includes('Missing details do not mean zero matches'));
vm.runInContext(`state.testEntry.row.identity_counts={'19':0,'18':0,'17':0};
state.testEntry.row.enumeration_complete=false;`,context);
const zero=vm.runInContext('matchCountPlot([state.testDecision])',context);
assert(zero.includes('0 / 0 / 0 · incomplete'));
assert(!zero.includes('NaN'));
vm.runInContext(`state.mismatchPlotScope='genomic';state.testEntry.hits=Array.from({length:55},(_,i)=>
 ({...state.testEntry.hits[0],identities:i<5?19:18,mismatch_positions:i<5?[]:[2],subject_start:100+i}));`,context);
map=vm.runInContext('mismatchPositionPlot([state.testDecision])',context);
assert(map.includes('40 alignments shown out of 55 retained details'));
assert(map.includes('Search evidence is incomplete'));
assert(map.indexOf('19/19 bases') < map.indexOf('18/19 bases'));
'''
    result = subprocess.run(
        [node, "-e", script, str(Path("src/sarna/app/static/offtarget.js").resolve())],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
