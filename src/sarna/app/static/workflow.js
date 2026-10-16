/* Researcher workspace: navigation, durable notes, and saved-run comparison. No design scoring. */
"use strict";
const workspaces = new Map();
const historyState = {offset:0, selected:new Set(), rows:[], request:0};
function activeWorkspace() { return state.result ? workspaces.get(state.result.app.run_id) : null; }
function workflowRun(result) {
  const id = result.app.run_id;
  const saved = result.workspace?.shortlist || {revision:0, entries:[]};
  const existing = workspaces.get(id);
  if (!existing || (!existing.dirty && !existing.saving && saved.revision>existing.revision)) {
    workspaces.set(id, {revision:saved.revision, entries:new Map(saved.entries.map(e=>[e.candidate,e.reason])), dirty:false, saving:null, problem:"", conflict:false});
  }
  document.body.classList.add("workflow-review");
  renderShortlist(); workflowStatus();
  if (state.selected) workflowCandidate(result.candidates.find(c=>c.id===state.selected));
}
function workflowStatus() {
  const run = state.result;
  for (const button of document.querySelectorAll("[data-stage]")) {
    button.disabled = !run && ["candidates","screening","export"].includes(button.dataset.stage);
  }
  $("edit-setup").hidden = !run || !!frozen;
  $("edit-setup").textContent = "Go to design setup";
  $("workflow-status").textContent = run
    ? `Reviewing saved run ${run.app.run_id.slice(0,8)}. ${run.counts.passing} candidates pass available gates; inspect missing evidence before selecting. ${frozen?"This is a read-only offline record.":"Settings edits create a new run."}`
    : "Choose your target and reference, then calculate candidates. Later stages become available with a saved result.";
}
function workflowCandidate(candidate) {
  if (!candidate || !activeWorkspace()) return;
  $("shortlist-candidate")?.remove();
  const button = document.createElement("button");
  button.id = "shortlist-candidate"; button.type = "button"; button.className = "secondary";
  button.textContent = activeWorkspace().entries.has(candidate.id) ? "Remove from shortlist" : "Add to research shortlist";
  button.disabled = !!frozen;
  button.addEventListener("click", () => {
    const workspace = activeWorkspace();
    workspace.entries.has(candidate.id) ? workspace.entries.delete(candidate.id) : workspace.entries.set(candidate.id, "");
    editWorkspace(workspace); renderShortlist(); workflowCandidate(candidate);
  });
  $("candidate-detail").querySelector(".detail-heading")?.after(button);
}
function shortlistStatus() {
  const workspace = activeWorkspace(); if (!workspace) return;
  $("shortlist-status").textContent = frozen ? "Archived researcher notes, read-only. Export selection remains available."
    : workspace.problem || (workspace.saving ? "Saving shortlist…" : workspace.dirty ? "Notes awaiting save…" : `Saved locally · revision ${workspace.revision}. Add a reason for each selection.`);
  $("shortlist-retry").hidden = !workspace.problem || workspace.conflict || !!frozen;
  $("shortlist-reload").hidden = !workspace.conflict || !!frozen;
}
function renderShortlist() {
  const workspace = activeWorkspace(); if (!workspace) return;
  $("shortlist-count").textContent = `· ${workspace.entries.size}`;
  $("shortlist-select").disabled = !workspace.entries.size;
  $("shortlist-entries").innerHTML = [...workspace.entries].map(([id,reason]) => {
    const candidate = state.result.candidates.find(c=>c.id===id);
    if (!candidate) return "";
    return `<article class="shortlist-entry" data-shortlist="${escapeHTML(id)}"><div><button type="button" class="text-button" data-shortlist-inspect="${escapeHTML(id)}">${escapeHTML(id)}</button><span class="tag">${candidate.rank===null?"EXCLUDED BY GATES":`CHECKLIST RANK ${candidate.rank}`}</span><code>${escapeHTML(candidate.target_site_dna)}</code><button type="button" class="text-button" data-shortlist-remove="${escapeHTML(id)}" ${frozen?"disabled":""}>Remove</button></div><label>Selection reason <textarea rows="2" maxlength="2000" data-shortlist-reason="${escapeHTML(id)}" ${frozen?"readonly":""} placeholder="Record your reasoning, evidence gaps, or planned follow-up.">${escapeHTML(reason)}</textarea></label>${candidate.rank===null?'<p class="hint">Shortlisting retains this excluded candidate for review; it does not override the failed gates.</p>':""}</article>`;
  }).join("") || '<p class="hint shortlist-empty">Inspect a candidate and choose “Add to research shortlist”. Table checkboxes are a temporary export/screening selection; the shortlist is saved with reasons.</p>';
  shortlistStatus();
}
function editWorkspace(workspace) {
  if (frozen) return;
  workspace.dirty = true; if (!workspace.conflict) workspace.problem = "";
  clearTimeout(workspace.timer);
  const id = state.result.app.run_id;
  workspace.timer = setTimeout(()=>saveWorkspace(id,workspace), 400);
  shortlistStatus();
}
async function saveWorkspace(id, workspace) {
  if (workspace.saving) return workspace.saving;
  if (!workspace.dirty || frozen) return;
  if (workspace.conflict) return;
  workspace.saving = (async()=>{
    while (workspace.dirty) {
      workspace.dirty = false;
      const entries = [...workspace.entries].map(([candidate,reason])=>({candidate,reason}));
      try {
        const saved = await api(`/api/runs/${id}/shortlist`, {method:"POST", headers:{"Content-Type":"application/json"}, body:JSON.stringify({revision:workspace.revision,entries})});
        workspace.revision = saved.revision;
        if (state.result?.app.run_id===id) {
          state.result.workspace = {shortlist:saved};
        }
      } catch (failure) {
        workspace.dirty = true; workspace.conflict = failure.code === "shortlist-conflict";
        workspace.problem = `Notes were not saved: ${failure.message} Your edits remain in this window.${workspace.conflict?" Reload saved notes replaces these local edits; download them first if needed.":" Retry saving before closing."}`;
        break;
      }
    }
  })();
  shortlistStatus();
  try { await workspace.saving; }
  finally { workspace.saving = null; if (state.result?.app.run_id===id) shortlistStatus(); }
}
async function flushWorkspace(id) {
  const workspace = workspaces.get(id); if (!workspace || frozen) return;
  clearTimeout(workspace.timer); await saveWorkspace(id,workspace);
  if (workspace.dirty || workspace.problem) throw new Error("Save or resolve shortlist notes before exporting the workspace. Download notes separately to keep unsaved edits.");
}
$("shortlist-entries").addEventListener("input",event=>{
  const id = event.target.dataset.shortlistReason; if (!id || frozen) return;
  const workspace=activeWorkspace(); workspace.entries.set(id,event.target.value); editWorkspace(workspace);
});
$("shortlist-entries").addEventListener("click",event=>{
  const inspect=event.target.closest("[data-shortlist-inspect]");
  if (inspect) { selectCandidate(inspect.dataset.shortlistInspect); $("candidate-layout").scrollIntoView({block:"start"}); }
  const remove=event.target.closest("[data-shortlist-remove]");
  if (remove && !frozen) { const workspace=activeWorkspace();workspace.entries.delete(remove.dataset.shortlistRemove);editWorkspace(workspace);renderShortlist();workflowCandidate(state.result.candidates.find(c=>c.id===state.selected)); }
});
$("shortlist-select").addEventListener("click",()=>{
  state.checked = new Set(activeWorkspace().entries.keys()); renderTable();
  $("shortlist-status").textContent = `${state.checked.size} shortlisted candidates selected. Screening uses passing candidates only. Explicit order selections can include excluded candidates; inspect their failed gates.`;
});
$("shortlist-retry").addEventListener("click",()=>saveWorkspace(state.result.app.run_id,activeWorkspace()));
$("shortlist-reload").addEventListener("click",async()=>{
  const id=state.result.app.run_id;
  try {
    const saved=await api(`/api/runs/${id}/shortlist`);
    workspaces.delete(id);
    if (state.result.app.run_id!==id) return;
    state.result.workspace={shortlist:saved};workflowRun(state.result);
  } catch(failure) { error(failure); }
});
const downloadNotes = document.createElement("button");
downloadNotes.type="button";downloadNotes.className="secondary";downloadNotes.textContent="Download shortlist notes";downloadNotes.id="shortlist-download";
downloadNotes.addEventListener("click",()=>{
  const workspace=activeWorkspace();if(!workspace)return;
  downloadBlob(new Blob([JSON.stringify({run_id:state.result.app.run_id,revision:workspace.revision,unsaved:workspace.dirty,entries:[...workspace.entries].map(([candidate,reason])=>({candidate,reason}))},null,2)],{type:"application/json"}),`shortlist-${state.result.app.run_id.slice(0,8)}.json`);
});
$("shortlist-card").append(downloadNotes);
window.addEventListener("beforeunload",event=>{
  if ([...workspaces.values()].some(workspace=>workspace.dirty || workspace.saving)) { event.preventDefault();event.returnValue=""; }
});
function stageTarget(stage) {
  if (stage==="target" || stage==="reference") {
    return frozen ? $("result-reference") : stage==="reference" && state.source==="gene" ? $("reference") : $("setup");
  }
  return $(stage==="candidates"?"candidate-layout":stage==="screening"?"offtarget-card":"export-workspace");
}
$("workflow-nav").addEventListener("click",event=>{
  const button=event.target.closest("[data-stage]");if(!button)return;
  document.querySelectorAll("[data-stage]").forEach(node=>node.removeAttribute("aria-current"));
  button.setAttribute("aria-current","step");
  const target=stageTarget(button.dataset.stage);target?.scrollIntoView({block:"start",behavior:"smooth"});
  if(target){target.setAttribute("tabindex","-1");target.focus({preventScroll:true});}
  workflowStatus();
});
$("edit-setup").addEventListener("click",()=>{
  $("setup").scrollIntoView({block:"start",behavior:"smooth"});
  $("setup").setAttribute("tabindex","-1");$("setup").focus({preventScroll:true});
});

/* History searches all saved runs in pages; comparison never modifies either snapshot. */
async function loadHistory(reset=false) {
  if (reset) historyState.offset=0;
  const request=++historyState.request;
  $("history-message").textContent="Searching saved runs…";
  try {
    const query=new URLSearchParams({q:$("history-search").value,status:$("history-status").value,offset:historyState.offset,limit:50});
    const [rows,groups]=await Promise.all([api(`/api/runs?${query}`),api("/api/batches")]);if(request!==historyState.request)return;
    const text=$("history-search").value.toLowerCase(), status=$("history-status").value;
    $("history-groups").innerHTML=groups.filter(group=>(!text || `${group.label} ${group.symbol} ${group.id}`.toLowerCase().includes(text)) && (!status || group.status===status)).map(group=>`<div class="history-row"><div><strong>${escapeHTML(group.label)}</strong><p class="hint">Promoter group · ${escapeHTML(group.status)}</p></div><button type="button" data-history-group="${escapeHTML(group.id)}" class="secondary">Open group</button></div>`).join("");
    historyState.rows=rows;
    $("history-rows").innerHTML=rows.map(row=>`<article class="history-row"><input type="checkbox" data-history-choose="${escapeHTML(row.id)}" aria-label="Compare ${escapeHTML(row.label)}" ${historyState.selected.has(row.id)?"checked":""} ${row.status!=="complete"?"disabled":""}><div><strong>${escapeHTML(row.label)}</strong><p class="hint">${escapeHTML(new Date(row.created).toLocaleString())} · ${escapeHTML(row.status)} · ${escapeHTML(row.id.slice(0,8))}</p><p class="hint">${escapeHTML(row.summary?.configuration?.preset_name || row.error || row.stage)}</p></div><button type="button" class="secondary" data-history-open="${escapeHTML(row.id)}" ${row.status!=="complete"?"disabled":""}>Open</button></article>`).join("");
    $("history-message").textContent=rows.length?`Showing ${historyState.offset+1}–${historyState.offset+rows.length}. ${historyState.selected.size} selected for comparison.`:"No saved runs match this search.";
    $("history-previous").disabled=!historyState.offset;$("history-next").disabled=rows.length<50;
  } catch(failure) {if(request===historyState.request)$("history-message").textContent=failure.message;}
}
$("open-history").addEventListener("click",()=>{$("history-dialog").showModal();loadHistory();});
$("close-history").addEventListener("click",()=>$("history-dialog").close());
let historyTimer;
$("history-search").addEventListener("input",()=>{clearTimeout(historyTimer);historyTimer=setTimeout(()=>loadHistory(true),200);});
$("history-status").addEventListener("change",()=>loadHistory(true));
$("history-next").addEventListener("click",()=>{historyState.offset+=50;loadHistory();});
$("history-previous").addEventListener("click",()=>{historyState.offset=Math.max(0,historyState.offset-50);loadHistory();});
$("history-rows").addEventListener("change",event=>{
  const id=event.target.dataset.historyChoose;if(!id)return;
  if(event.target.checked && historyState.selected.size===2) {event.target.checked=false;$("history-message").textContent="Select at most two runs. Uncheck one to change the comparison.";return;}
  event.target.checked?historyState.selected.add(id):historyState.selected.delete(id);
  $("history-comparison").innerHTML="";
  $("history-compare").disabled=historyState.selected.size!==2;
  $("history-message").textContent=`${historyState.selected.size} selected for comparison, including selections on other pages.`;
});
$("history-rows").addEventListener("click",async event=>{
  const button=event.target.closest("[data-history-open]");if(!button || state.busy)return;
  try {const run=await api(`/api/runs/${button.dataset.historyOpen}`);if(run.status!=="complete")throw new Error(run.error || run.stage);state.batch=null;$("batch-comparison").hidden=true;restoreDraft(run.result);showResult(run.result);$("saved-runs").value=run.id;$("history-dialog").close();}catch(failure){$("history-message").textContent=failure.message;}
});
$("history-groups").addEventListener("click",async event=>{
  const button=event.target.closest("[data-history-group]");if(!button || state.busy)return;
  try{await openSavedBatch(button.dataset.historyGroup);$("history-dialog").close();}catch(failure){$("history-message").textContent=failure.message;}
});
function compareRuns(left,right) {
  const sameInput=JSON.stringify(ordered(left.promoter))===JSON.stringify(ordered(right.promoter)) && JSON.stringify(ordered(left.app.anchor))===JSON.stringify(ordered(right.app.anchor));
  const describe=run=>({"Run":run.app.run_id,"Gene":run.gene?.symbol || "Pasted DNA","Reference":run.app.anchor?.label || "User-supplied sequence","Assembly":run.promoter.assembly || "Not assigned","Cell context":run.app.cell_line,"Promoter span":`${run.promoter.length} bp`,"Configuration":run.configuration.preset.name,"Code fingerprint":run.reproducibility?.software?.package_sha256 || "Not recorded","Passing / generated":`${run.counts.passing} / ${run.counts.generated}`});
  const a=describe(left),b=describe(right);
  const settings=run=>Object.fromEntries([...run.configuration.rules,...(run.configuration.disabled_rules || [])].map(rule=>[rule.id,{role:rule.kind,values:Object.fromEntries(rule.parameters.map(parameter=>[parameter.name,parameter.value])),option:rule.selected_option?.id || null}]));
  const before=settings(left),after=settings(right);
  const changes=[...new Set([...Object.keys(before),...Object.keys(after)])].filter(id=>JSON.stringify(ordered(before[id]))!==JSON.stringify(ordered(after[id]))).map(id=>`<li><b>${escapeHTML(id)}</b><br>First: <code>${escapeHTML(JSON.stringify(before[id] || "omitted"))}</code><br>Second: <code>${escapeHTML(JSON.stringify(after[id] || "omitted"))}</code></li>`);
  let rankNote="Different input or reference: candidate ranks are not matched across these runs. Passing counts describe each run independently; they do not establish better activity.";
  if(sameInput) {
    const old=new Map(left.candidates.filter(c=>c.rank!==null).map(c=>[`${c.offset}:${c.target_site_dna}`,c.rank]));
    const now=new Map(right.candidates.filter(c=>c.rank!==null).map(c=>[`${c.offset}:${c.target_site_dna}`,c.rank]));
    rankNote=`Same saved promoter and reference: ${[...now.keys()].filter(k=>!old.has(k)).length} newly passing, ${[...old.keys()].filter(k=>!now.has(k)).length} no longer passing, ${[...now].filter(([k,rank])=>old.has(k)&&old.get(k)!==rank).length} shared candidates changed rank. These are computational differences, not evidence of improved activity.`;
  }
  return `<h3>Run comparison</h3><p class="hint">${escapeHTML(rankNote)}</p><div class="table-scroll"><table><thead><tr><th>Recorded field</th><th>First run</th><th>Second run</th></tr></thead><tbody>${Object.keys(a).map(key=>`<tr><th>${escapeHTML(key)}</th><td>${escapeHTML(a[key])}</td><td>${escapeHTML(b[key])}</td></tr>`).join("")}</tbody></table></div><h4>Rule and parameter changes</h4><ul>${changes.join("") || "<li>No differences in effective rule roles, values or named options.</li>"}</ul><p class="hint">The same cell name does not establish equivalent cell-specific evidence. Inspect each run's sources and missing evidence. Comparisons do not update stored scientific results.</p>`;
}
$("history-compare").addEventListener("click",async()=>{
  const ids=[...historyState.selected];if(ids.length!==2)return;
  $("history-compare").disabled=true;$("history-message").textContent="Loading the two saved records…";
  try {const runs=await Promise.all(ids.map(id=>api(`/api/runs/${id}`)));if(runs.some(run=>run.status!=="complete"))throw new Error("Both runs must be complete.");$("history-comparison").innerHTML=compareRuns(runs[0].result,runs[1].result);$("history-message").textContent="Comparison loaded. First and second follow your selection order.";}catch(failure){$("history-message").textContent=failure.message;}finally{$("history-compare").disabled=historyState.selected.size!==2;}
});
if (frozen) { $("open-history").hidden=true; }
if (state.result) workflowRun(state.result); else workflowStatus();
