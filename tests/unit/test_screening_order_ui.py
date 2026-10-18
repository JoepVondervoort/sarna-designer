"""Screening preserves design context and allows explicit fewest-match comparisons."""

import shutil
import subprocess
from pathlib import Path

import pytest


def test_screening_order_scores_zero_counts_and_incomplete_evidence() -> None:
    node = shutil.which("node")
    if not node:
        pytest.skip("Node is required for frontend state checks")
    script = r'''
const fs=require('fs'),vm=require('vm'),assert=require('assert');
const handlers=new Map(),elements=new Map();
const context={console,frozen:null,state:{},escapeHTML:String,fmt:String,signed:String,
 setTimeout,clearTimeout,error:()=>{},
 document:{addEventListener:(name,fn)=>handlers.set(name,fn)},
 $:id=>{if(!elements.has(id))elements.set(id,{addEventListener:()=>{}});return elements.get(id);}};
vm.createContext(context);vm.runInContext(fs.readFileSync(process.argv[1],'utf8'),context);
vm.runInContext(`
state.result={candidates:[{id:'best',rank:1,rating:{met:6,evaluated:7}},
 {id:'second',rank:2,rating:{met:5,evaluated:7}},
 {id:'zero',rank:3,rating:{met:4,evaluated:7}},
 {id:'partial',rank:4,rating:{met:3,evaluated:7}}]};
state.screeningScope='genomic';
const entry=(id,n,extra={})=>({row:{candidate:id,hit_count:n,perfect_hits:0,
 wang_homology_hits:n,enumeration_complete:true,strong_threshold:16,...extra},
 hits:[],placements:[],search:{parameters:{database:'reference'},annotations:{},backend:'local'}});
state.entries=[entry('second',100),entry('zero',0),entry('best',614),
 entry('partial',0,{alignments_capped:true}),entry('unknown',null)];
`,context);
const ids=expression=>Array.from(vm.runInContext(expression,context));
assert.deepEqual(ids('screeningTableOrder(state.entries).map(e=>e.row.candidate)'),
 ['best','second','zero','partial','unknown']);
assert.deepEqual(ids('screeningTableOrder(state.entries,"hit_count").map(e=>e.row.candidate)'),
 ['zero','second','best','unknown','partial']);
assert.deepEqual(ids('screeningTableOrder(state.entries,"hit_count",true).map(e=>e.row.candidate)'),
 ['best','second','zero','unknown','partial']);
let html=vm.runInContext('resultsTable({},state.entries,19,true)',context);
const order=[...html.matchAll(/data-screenrow="([^"]+)"/g)].map(x=>x[1]);
assert.deepEqual(order,['best','second','zero','partial','unknown']);
assert(html.includes('6 / 7'));
assert(html.includes('data-label="Design rank">1'));
assert(html.includes('No reported alignments; specificity untested'));
assert(html.includes('data-label="Design score"><strong>Not available'));
assert(html.includes('same default order as Candidates'));
vm.runInContext('state.screeningSort="hit_count"',context);
html=vm.runInContext('resultsTable({},state.entries,19,true)',context);
assert.equal([...html.matchAll(/data-screenrow="([^"]+)"/g)][0][1],'zero');
assert(html.includes('aria-sort="ascending"'));
const original=vm.runInContext('JSON.stringify(state.entries)',context);
vm.runInContext('screeningTableOrder(state.entries,"perfect_hits");',context);
assert.equal(vm.runInContext('JSON.stringify(state.entries)',context),original);
vm.runInContext('renderOffTarget=()=>{};',context);
handlers.get('click')({target:{closest:selector=>selector==='[data-screensort]'
 ? {dataset:{screensort:'hit_count'}}:null}});
assert.equal(context.state.screeningSortDescending,true);
handlers.get('click')({target:{closest:selector=>selector==='[data-screensort]'
 ? {dataset:{screensort:'design_rank'}}:null}});
assert.equal(context.state.screeningSort,'design_rank');
assert.equal(context.state.screeningSortDescending,false);
vm.runInContext('state.screeningOpen="zero"',context);
html=vm.runInContext('candidateRow(state.entries[1],19,true)',context);
assert(html.includes('colspan="10"'));
html=vm.runInContext('candidateRow(state.entries[1],19,false)',context);
assert(html.includes('colspan="9"'));
'''
    result = subprocess.run(
        [node, "-e", script, str(Path("src/sarna/app/static/offtarget.js").resolve())],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
