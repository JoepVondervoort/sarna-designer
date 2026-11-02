/* Researcher workspace: navigation, durable notes, and saved-run comparison. No design scoring. */
"use strict";
const workspaces = new Map();
const historyState = {offset:0, selected:new Set(), rows:[], request:0};
const WORKSPACE_TABS = ["design", "candidates", "screening", "export"];
function workspaceAvailable(tab) {
  return tab === "design" || tab === "screening" || !!state.result || !!state.batch;
}
function setWorkspaceTab(tab, {focus=false, scroll=true} = {}) {
  if (!WORKSPACE_TABS.includes(tab) || !workspaceAvailable(tab)) return;
  state.workspaceTab = tab;
  for (const name of WORKSPACE_TABS) {
    const active = name === tab, button = $(`tab-${name}`), panel = $(`workspace-${name}`);
    panel.hidden = !active;
    button.setAttribute("aria-selected", String(active));
    button.tabIndex = active ? 0 : -1;
  }
  document.body.dataset.workspace = tab;
  if (tab === "candidates" && state.result) requestAnimationFrame(()=>renderMap());
  if (tab === "export") renderExportSelection();
  workflowStatus();
  if (scroll) $("workflow-nav").scrollIntoView({block:"start"});
  if (focus) $(`workspace-${tab}`).focus({preventScroll:true});
}
function workflowBatch(batch) {
  state.workspaceBatchId = batch.id;
  workflowStatus();
  if (!state.workspaceTab) setWorkspaceTab(frozen ? "candidates" : "design", {scroll:false});
  if (batch.status !== "running" && !frozen) {
    $("design-ready").hidden = false;
    $("design-ready-message").textContent = "Promoter comparison finished. Review the completed runs when you are ready.";
  }
}
function activeSelection() {
  if (!state.result) return [];
  const managedEmpty=state.screenedManagedRun===state.result.app.run_id && !state.checked.size;
  return state.result.candidates.filter(c=>state.checked.size || managedEmpty ? state.checked.has(c.id) : c.rank!==null);
}
function renderExportSelection() {
  const host=$("export-selection"), preview=$("export-preview-host");
  if(!state.result){host.hidden=true;preview.hidden=true;return;}
  host.hidden=false;preview.hidden=false;
  const selected=activeSelection().sort((a,b)=>(a.rank ?? Infinity)-(b.rank ?? Infinity));
  const title=!selected.length?"No candidates selected":state.checked.size?`${selected.length} selected candidates`:`All ${selected.length} design-passing candidates`;
  screeningHTML(host,`<h3>${title}</h3><p class="hint">${!selected.length?"Check candidates in Candidates before downloading selected sequences.":state.checked.size?"Shared selection with screening. Change it in Candidates.":"No rows checked: downloads use all passing candidates."}</p>`);
  const evidenceById=new Map(screenedDecisions().map(d=>[d.candidate.id,d]));
  const previewOpen=!!preview.querySelector("#export-preview[open]");
  screeningHTML(preview,selected.length?`<details id="export-preview" ${previewOpen?"open":""}><summary>Preview selected sequences and screening evidence</summary><div class="table-scroll"><table><caption>Duplex preview · ordered RNA sequences, 5′ → 3′</caption><thead><tr><th>Candidate</th><th>Sense</th><th>Antisense</th><th>Evidence</th></tr></thead><tbody>${selected.slice(0,100).map(c=>{
    const order=orderRowFor(state.result,c.id);
    return `<tr><th><button type="button" class="text-button" data-export-inspect="${escapeHTML(c.id)}">${escapeHTML(c.id)}</button><small>${c.rank===null?"Excluded by design gates":`Design rank ${c.rank}`}</small></th><td><code>${escapeHTML(order?.sense_order_5to3 || "Not available")}</code></td><td><code>${escapeHTML(order?.antisense_order_5to3 || "Not available")}</code></td><td><div class="candidate-screening-note">${evidenceById.has(c.id)?candidateScreeningSummary(evidenceById.get(c.id)):escapeHTML(screeningSentence(c.id))}</div></td></tr>`;
  }).join("")}</tbody></table></div>${selected.length>100?`<p class="hint">Previewing the first 100 of ${selected.length}; downloads include the complete selection.</p>`:""}</details>`:"");
}
function viewScreeningCandidate(id, requestedScope=null) {
  if (!state.result?.candidates.some(c=>c.id===id)) return;
  state.screeningOpen = id;
  state.mismatchPlotCandidate = id;
  const preferred = requestedScope || state.screeningScope || "genomic";
  const scope = comparisonEntry(preferred,id) ? preferred
    : comparisonEntry("genomic",id) ? "genomic" : comparisonEntry("transcript",id) ? "transcript" : preferred;
  state.mismatchPlotScope = scope;
  renderOffTarget(scope);
  renderScreenedShortlist();
  setWorkspaceTab("screening", {focus:true});
  const details = $("offtarget-results").querySelector(".screening-detail-table");
  if (details) details.open = true;
  const row = [...$("offtarget-results").querySelectorAll("[data-screenrow]")].find(node=>node.dataset.screenrow===id);
  (row || $("offtarget-card")).scrollIntoView({block:"start"});
}

function activeWorkspace() { return state.result ? workspaces.get(state.result.app.run_id) : null; }
function workflowRun(result) {
  investigationRun(result);
  const id = result.app.run_id;
  const changedRun = state.workspaceRunId !== id;
  state.workspaceRunId = id;
  const saved = result.workspace?.shortlist || {revision:0, entries:[]};
  const existing = workspaces.get(id);
  if (!existing || (!existing.dirty && !existing.saving && saved.revision>existing.revision)) {
    workspaces.set(id, {revision:saved.revision, entries:new Map(saved.entries.map(e=>[e.candidate,e.reason])), dirty:false, saving:null, problem:"", conflict:false});
  }
  document.body.classList.add("workflow-review");
  renderShortlist(); workflowStatus();
  if (changedRun || !state.workspaceTab) {
    setWorkspaceTab(state.workspaceTab || (frozen ? "candidates" : "design"), {scroll:false});
    $("design-ready").hidden = !!frozen;
    $("design-ready-message").textContent = `Design ready: ${result.counts.passing} of ${result.counts.generated} targets pass the current gates. Your results are saved. Open Candidates when you are ready.`;
  }
  if (state.selected) workflowCandidate(result.candidates.find(c=>c.id===state.selected));
}
function workflowStatus() {
  const run = state.result;
  for (const button of document.querySelectorAll("[data-stage]")) {
    button.disabled = !workspaceAvailable(button.dataset.stage);
  }
  for (const button of document.querySelectorAll("[data-workspace-link]")) button.disabled = !workspaceAvailable(button.dataset.workspaceLink);
  $("run-context-details").hidden=!run && !state.batch;
  $("edit-setup").hidden = !run || !!frozen || state.workspaceTab === "design";
  $("edit-setup").textContent = "Edit design";
  const name = run?.gene?.symbol || state.gene?.symbol || (run ? "Pasted promoter" : "No target selected");
  const reference = run?.app.request?.custom_reference?.name || run?.transcript?.ensembl_id || run?.app.anchor?.label || "";
  $("workspace-context-label").textContent = `${name}${reference ? ` · ${reference}` : ""}${run ? "" : " · Design setup"}`;
  $("workspace-selection-count").textContent = run ? state.checked.size?`${state.checked.size} selected · screening & export`:state.screenedManagedRun===run.app.run_id?"No candidates selected":`${activeSelection().length} active · all passing` : "";
  renderCalmSummary();
  const jobs = Object.values(state.screeningJobs || {});
  const running = jobs.some(j=>j.status === "running"), waiting = jobs.some(j=>j.status === "waiting");
  const badge = $("screening-tab-status");
  badge.textContent = running ? "Running" : waiting ? "Waiting" : jobs.some(j=>j.status === "error") ? "Needs review" : Object.keys(state.screenings || {}).length ? "Results" : "";
  badge.classList.toggle("is-loading", running);
  badge.setAttribute("aria-label", running ? "Screening running in the background" : `Screening ${badge.textContent || "not started"}`);
  $("screening-before-run").hidden = !!run;
  $("offtarget-card").hidden = !run;
  $("archived-design").hidden = !frozen;
  $("top-duplex-disclosure").hidden = !run || !run.counts.passing;
  $("group-export").hidden = !state.batch;
  $("export-workspace").hidden = !run;
  $("reproducibility-panel").hidden = !run?.reproducibility;
  if ($("candidate-run")) $("candidate-run").hidden = !run;
  $("result-area").hidden = !run || state.workspaceTab === "design";
  const descriptions = {
    design: frozen ? "Saved reference and design settings · read-only offline record." : "Choose your target and context. Advanced settings are optional; edits create a new run.",
    candidates: "Compare design rank, inspect candidates and keep your selection. Screening summaries update here as searches finish.",
    screening: "Review sequence matches. Searches continue in the background while you use other tabs.",
    export: "Review the selected duplexes and download sequences or the complete scientific record."
  };
  const description = descriptions[state.workspaceTab || "design"];
  if ($("workflow-status").textContent !== description) $("workflow-status").textContent = description;
}
function workflowCandidate(candidate) {
  if (!candidate || !activeWorkspace()) return;
  $("shortlist-candidate")?.remove();
  const button = document.createElement("button");
  button.id = "shortlist-candidate"; button.type = "button"; button.className = "secondary";
  button.textContent = activeWorkspace().entries.has(candidate.id) ? "Remove saved candidate" : "Save to shortlist";
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
  if (workspace.problem) $("shortlist-card").open=true;
  $("shortlist-status").textContent = frozen ? "Archived researcher notes, read-only. Export selection remains available."
    : workspace.problem || (workspace.saving ? "Saving shortlist…" : workspace.dirty ? "Notes awaiting save…" : `Saved locally`);
  $("shortlist-retry").hidden = !workspace.problem || workspace.conflict || !!frozen;
  $("shortlist-reload").hidden = !workspace.conflict || !!frozen;
}
function renderShortlist() {
  const workspace = activeWorkspace(); if (!workspace) return;
  $("shortlist-count").textContent = `· ${workspace.entries.size}`;
  $("saved-shortlist-count").textContent = workspace.entries.size ? `(${workspace.entries.size})` : "";
  $("shortlist-select").disabled = !workspace.entries.size;
  $("shortlist-entries").innerHTML = [...workspace.entries].map(([id,reason]) => {
    const candidate = state.result.candidates.find(c=>c.id===id);
    if (!candidate) return "";
    return `<article class="shortlist-entry" data-shortlist="${escapeHTML(id)}"><div><button type="button" class="text-button" data-shortlist-inspect="${escapeHTML(id)}">${escapeHTML(id)}</button><span class="tag">${candidate.rank===null?"EXCLUDED BY GATES":`CHECKLIST RANK ${candidate.rank}`}</span><code>${escapeHTML(candidate.target_site_dna)}</code><button type="button" class="text-button" data-shortlist-remove="${escapeHTML(id)}" ${frozen?"disabled":""}>Remove</button></div><label>Selection reason <textarea rows="2" maxlength="2000" data-shortlist-reason="${escapeHTML(id)}" ${frozen?"readonly":""} placeholder="Record your reasoning, evidence gaps, or planned follow-up.">${escapeHTML(reason)}</textarea></label>${candidate.rank===null?'<p class="hint">Shortlisting retains this excluded candidate for review; it does not override the failed gates.</p>':""}</article>`;
  }).join("") || '<p class="hint shortlist-empty">Save candidates with notes from their inspector. “Use saved shortlist” copies them into the active screening and export selection.</p>';
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
          state.result.workspace = {...state.result.workspace,shortlist:saved};
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
  if (inspect) { setWorkspaceTab("candidates"); selectCandidate(inspect.dataset.shortlistInspect); $("candidate-layout").scrollIntoView({block:"start"}); }
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
    state.result.workspace={...state.result.workspace,shortlist:saved};workflowRun(state.result);
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
$("workflow-nav").addEventListener("click", event=>{
  const button = event.target.closest("[data-stage]");
  if (button && !button.disabled) setWorkspaceTab(button.dataset.stage, {focus:true});
});
$("workflow-nav").addEventListener("keydown", event=>{
  if (!["ArrowLeft","ArrowRight","Home","End"].includes(event.key)) return;
  const enabled = [...document.querySelectorAll("[data-stage]")].filter(b=>!b.disabled);
  const index = enabled.indexOf(event.target);
  if (index<0) return;
  event.preventDefault();
  const next = event.key === "Home" ? 0 : event.key === "End" ? enabled.length-1
    : (index + (event.key === "ArrowRight" ? 1 : -1) + enabled.length) % enabled.length;
  setWorkspaceTab(enabled[next].dataset.stage, {scroll:false});
  enabled[next].focus();
});
$("edit-setup").addEventListener("click",()=>setWorkspaceTab("design", {focus:true}));
document.addEventListener("click", event=>{
  const link = event.target.closest("[data-workspace-link]");
  if (link) setWorkspaceTab(link.dataset.workspaceLink, {focus:true});
  const evidence = event.target.closest("[data-open-run-evidence]");
  if (evidence) renderEvidence();
  const inspect = event.target.closest("[data-export-inspect]");
  if (inspect) { setWorkspaceTab("candidates"); selectCandidate(inspect.dataset.exportInspect); $("candidate-detail").scrollIntoView({block:"start"}); }
  const screening = event.target.closest("[data-screening-candidate]");
  if (screening) viewScreeningCandidate(screening.dataset.screeningCandidate,screening.dataset.screeningScope);
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
// Map investigation is researcher metadata, independent of the frozen design results.
function investigationRun(result) {
  if (state.investigationRunId !== result.app.run_id) {
    state.investigationRunId = result.app.run_id;
    state.investigation = clone(result.workspace?.investigation || {revision:0,annotations:[],motifs:null});
    state.mapSearchMatches = [];
    state.investigationConflict = false;
    $("map-search-results").innerHTML = "";
    $("map-search-status").textContent = "";
    $("distance-result").textContent = "";
    $("investigation-status").textContent = frozen ? "Read-only offline record." : "Map settings saved.";
    for (const id of ["annotation-save","motifs-all","motifs-none"]) $(id).disabled=!!frozen;
    $("investigation-reload").hidden=!!frozen;
  }
  renderInvestigation();
}
function renderInvestigation() {
  const data=state.investigation;
  if (!data || !state.result) return;
  $("map-annotations").innerHTML=data.annotations.map((a,i)=>`<p><button type="button" class="text-button" data-region-open="${i}">${escapeHTML(a.name)} · bases ${a.start+1}–${a.end}</button> · ${escapeHTML(a.orientation)} ${escapeHTML(a.source)} <button type="button" class="text-button" data-region-remove="${i}" ${frozen?'disabled':''}>Remove</button></p>`).join("");
  renderMotifChoices();
  const features=investigationFeatures();
  state.distanceFeatures=features;
  for (const id of ["distance-from","distance-to"]) {
    const node=$(id), previous=node.value;
    node.innerHTML=features.map((f,i)=>`<option value="${i}">${escapeHTML(f.name)}${f.point?'':` · ${f.start+1}–${f.end}`}</option>`).join("");
    if ([...node.options].some(o=>o.value===previous)) node.value=previous;
  }
}
function renderMotifChoices() {
  const query=$("motif-search").value.trim().toLowerCase();
  $("motif-choices").innerHTML=(state.result?.motif_scans || []).filter(m=>`${m.motif} ${m.matrix_id}`.toLowerCase().includes(query)).map(m=>`<label class="check"><input type="checkbox" data-motif-choice="${escapeHTML(m.motif)}" ${state.investigation.motifs===null || state.investigation.motifs.includes(m.motif)?'checked':''} ${frozen?'disabled':''}>${escapeHTML(m.motif)} · ${escapeHTML(m.matrix_id)} · ${m.sites} sites</label>`).join("") || '<p class="hint">No matching scanned motifs.</p>';
}
async function saveInvestigation(next) {
  if (frozen) return false;
  if (state.investigationSaving) { $("investigation-status").textContent="A save is in progress. Wait for it to finish, then retry your change."; renderMotifChoices(); return false; }
  if (state.investigationConflict) { $("investigation-status").textContent="Reload saved map settings before retrying. Your pending region fields remain available."; return false; }
  const run=state.result, id=run.app.run_id;
  state.investigationSaving=true;
  $("investigation-status").textContent="Saving map settings…";
  try {
    const saved=await api(`/api/runs/${id}/investigation`,{method:"POST",headers:{"Content-Type":"application/json"},body:JSON.stringify({revision:state.investigation.revision,annotations:next.annotations,motifs:next.motifs})});
    run.workspace={...run.workspace,investigation:saved};
    if (state.result.app.run_id===id) { state.investigation=saved; renderInvestigation(); renderMap(); $("investigation-status").textContent="Map settings saved."; }
    return true;
  } catch (failure) {
    if (state.result.app.run_id===id) { state.investigationConflict=failure.code==="investigation-conflict"; $("investigation-status").textContent=`Not saved: ${failure.message}`; renderMotifChoices(); }
    return false;
  } finally { state.investigationSaving=false; }
}
function focusRegion(region) {
  const pad=Math.max(30,region.end-region.start);
  state.start=Math.max(0,region.start-pad); state.end=Math.min(state.result.promoter.length,region.end+pad);
  renderMap(); $("promoter-landscape").scrollIntoView({block:"start"});
}
function investigationFeatures() {
  const r=state.result, features=[];
  if (!r) return features;
  // TSS-relative numbering skips zero; the point is the beginning of the +1 base.
  if (r.promoter.tss_start!==null && r.promoter.tss_start!==undefined) {
    const t=r.promoter.tss_start, offset=t<0?-t:1-t;
    if (offset>=0 && offset<=r.promoter.length) features.push({name:"Reference TSS",start:offset,end:offset,point:true});
  }
  features.push(...state.investigation.annotations.map(a=>({...a,name:`Region: ${a.name}`})));
  features.push(...(state.sequenceRecords || []).flatMap(a=>a.matches.map(m=>({name:`saRNA: ${a.name}`,start:m.offset,end:m.end}))));
  features.push(...r.candidates.map(c=>({name:c.id,start:c.offset,end:c.offset+c.target_site_dna.length})));
  features.push(...r.tracks.regulatory.map(t=>({...t,name:`Regulatory: ${t.label}`})));
  features.push(...r.tracks.motifs.filter(t=>state.investigation.motifs===null || state.investigation.motifs.includes(t.label)).map(t=>({...t,name:`Motif: ${t.label}`})));
  return features;
}
$("map-find").addEventListener("click",async()=>{
  if (!state.result) return;
  const run=state.result;
  $("map-find").disabled=true;
  try {
    const sequence=$("map-query").value.replace(/\s/g,"").toUpperCase(), orientation=$("map-query-orientation").value;
    let matches;
    if (frozen) {
      if (!/^[ACGTU]{4,10000}$/.test(sequence) || (sequence.includes("U") && sequence.includes("T"))) throw new Error("Use 4–10000 DNA or RNA bases, without mixed T and U.");
      const dna=sequence.replace(/U/g,"T"), reverse=dna.split("").reverse().map(b=>({A:"T",T:"A",C:"G",G:"C"}[b])).join("");
      matches=[];const seen=new Set();
      for (const [target,label] of [[dna,"same"],[reverse,"reverse complement"]]) {
        if ((orientation!=="both" && orientation!==label) || seen.has(target)) continue;
        seen.add(target);
        for (let start=run.promoter.sequence.indexOf(target);start>=0;start=run.promoter.sequence.indexOf(target,start+1)) {
          if(matches.length>=500)throw new Error("More than 500 matches; use a longer sequence.");
          matches.push({start,end:start+target.length,orientation:label});
        }
      }
      matches.sort((a,b)=>a.start-b.start);
    } else ({matches}=await api(`/api/runs/${run.app.run_id}/sequence-search`,{method:"POST",headers:{"Content-Type":"application/json"},body:JSON.stringify({sequence,orientation})}));
    if (state.result!==run) return;
    state.mapSearchMatches=matches;
    $("map-search-status").textContent=matches.length?`${matches.length} exact match(es). Choose a site to locate it and prepare a named annotation.`:"No exact match in the loaded promoter window. This does not establish absence elsewhere in the genome.";
    $("map-search-results").innerHTML=matches.map((m,i)=>`<button type="button" class="text-button" data-search-hit="${i}">Bases ${m.start+1}–${m.end} · ${escapeHTML(m.orientation)}</button>`).join(" ");
    renderMap();
  } catch(failure) { if(state.result===run) { state.mapSearchMatches=[]; $("map-search-results").innerHTML=""; $("map-search-status").textContent=failure.message; renderMap(); } }
  finally { $("map-find").disabled=false; }
});
$("map-search-results").addEventListener("click",event=>{
  const button=event.target.closest("[data-search-hit]");if(!button)return;
  const hit=state.mapSearchMatches[Number(button.dataset.searchHit)];
  $("annotation-start").value=hit.start+1;$("annotation-end").value=hit.end;
  $("annotation-save").closest("details").open=true;focusRegion(hit);
});
$("annotation-save").addEventListener("click",async()=>{
  const start=Number($("annotation-start").value)-1,end=Number($("annotation-end").value),name=$("annotation-name").value.trim();
  if(!name || !Number.isInteger(start) || !Number.isInteger(end) || start<0 || end<=start || end>state.result.promoter.length) { $("investigation-status").textContent="Name the region and enter valid inclusive bases within the promoter."; return; }
  const hit=(state.mapSearchMatches || []).find(m=>m.start===start && m.end===end);
  await saveInvestigation({...state.investigation,annotations:[...state.investigation.annotations,{name,start,end,source:$("annotation-source").value,orientation:hit?.orientation || "region"}]});
});
$("map-annotations").addEventListener("click",async event=>{
  const open=event.target.closest("[data-region-open]"), remove=event.target.closest("[data-region-remove]");
  if(open)focusRegion(state.investigation.annotations[Number(open.dataset.regionOpen)]);
  if(remove)await saveInvestigation({...state.investigation,annotations:state.investigation.annotations.filter((a,i)=>i!==Number(remove.dataset.regionRemove))});
});
$("investigation-reload").addEventListener("click",async()=>{
  const run=state.result;
  try { const saved=await api(`/api/runs/${run.app.run_id}/investigation`); if(state.result!==run)return; state.investigation=saved;state.investigationConflict=false;run.workspace={...run.workspace,investigation:saved};renderInvestigation();renderMap();$("investigation-status").textContent="Saved map settings reloaded. Pending region fields retained."; }
  catch(failure){$("investigation-status").textContent=failure.message;}
});
$("motif-search").addEventListener("input",renderMotifChoices);
$("motif-choices").addEventListener("change",async event=>{
  if(!event.target.matches("[data-motif-choice]"))return;
  const motifs=new Set(state.investigation.motifs ?? state.result.motif_scans.map(m=>m.motif));
  event.target.checked?motifs.add(event.target.dataset.motifChoice):motifs.delete(event.target.dataset.motifChoice);
  await saveInvestigation({...state.investigation,motifs:[...motifs]});
});
for(const [id,motifs] of [["motifs-all",null],["motifs-none",[]]]) $(id).addEventListener("click",()=>saveInvestigation({...state.investigation,motifs}));
$("distance-measure").addEventListener("click",()=>{
  const features=state.distanceFeatures || [],a=features[Number($("distance-from").value)],b=features[Number($("distance-to").value)];
  if(!a || !b)return;
  const overlap=Math.max(0,Math.min(a.end,b.end)-Math.max(a.start,b.start));
  const gap=Math.max(0,Math.max(a.start,b.start)-Math.min(a.end,b.end));
  const containsPoint=(a.point && a.start>=b.start && a.start<b.end)||(b.point && b.start>=a.start && b.start<a.end);
  $("distance-result").textContent=`${a.name} → ${b.name}: ${containsPoint?"point lies within region; ":""}${overlap?`${overlap} bp overlap`:`${gap} bp edge-to-edge gap`}; signed start displacement ${signed(b.start-a.start)} bp.`;
});

function investigationSelection(run) {
  return run.candidates.filter(c=>state.checked.size?state.checked.has(c.id):state.screenedManagedRun===run.app.run_id?false:c.rank!==null);
}
function alignedFeatures(run,candidates) {
  const data=state.investigation;
  return [
    ...candidates.map(c=>({name:c.id,kind:"saRNA target",start:c.offset,end:c.offset+c.target_site_dna.length,orientation:"same",source:`Design rank ${c.rank ?? "excluded"}`})),
    ...run.tracks.motifs.filter(t=>(data.motifs===null || data.motifs.includes(t.label))).map(t=>({...t,name:t.label,kind:"Predicted TF motif",orientation:t.strand || "unspecified",source:t.detail})),
    ...data.annotations.map(a=>({...a,kind:"Researcher region"})),
    ...(state.sequenceRecords || []).flatMap(a=>a.matches.map(m=>({name:a.name,kind:a.reference?"Reference saRNA":"Imported saRNA",start:m.offset,end:m.end,orientation:m.input_orientation,source:a.source})))
  ].sort((a,b)=>a.start-b.start || a.end-b.end);
}
function reportDocument(title,body) {
  return `<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>${escapeHTML(title)}</title><style>body{font:15px/1.5 system-ui,sans-serif;max-width:1100px;margin:32px auto;padding:0 24px;color:#18363a}table{border-collapse:collapse;width:100%;margin:16px 0}td,th{border:1px solid #ccd5d6;padding:7px;text-align:left;vertical-align:top}code,pre{font-family:ui-monospace,monospace}pre{overflow:auto;font-size:12px;line-height:1.7;background:#f4f7f7;padding:12px}td{overflow-wrap:anywhere}small{color:#496165}h2{margin-top:32px}@media print{body{margin:0;padding:0;font-size:10pt}pre{white-space:pre-wrap;overflow-wrap:anywhere}tr{break-inside:avoid}h2,h3{break-after:avoid}}</style></head><body><h1>${escapeHTML(title)}</h1>${body}</body></html>`;
}
function regionTable(run,features) {
  const p=run.promoter;
  return `<table><thead><tr><th>Feature</th><th>Type / orientation</th><th>Promoter bases (1-based, inclusive)</th><th>Genomic interval (0-based, half-open)</th><th>Source</th></tr></thead><tbody>${features.map(f=>{
    const genomic=p.start===null || p.start===undefined?"Unavailable":p.strand==="-"?`${p.contig}:${p.end-f.end}–${p.end-f.start}`:`${p.contig}:${p.start+f.start}–${p.start+f.end}`;
    return `<tr><th>${escapeHTML(f.name)}</th><td>${escapeHTML(f.kind)} · ${escapeHTML(f.orientation)}</td><td>${f.start+1}–${f.end}</td><td>${escapeHTML(genomic)}</td><td>${escapeHTML(f.source)}</td></tr>`;
  }).join("")}</tbody></table>`;
}
function alignmentDocument(run,candidates) {
  const features=alignedFeatures(run,candidates), sequence=run.promoter.sequence;
  let blocks="";
  for(let start=0;start<sequence.length;start+=80) {
    const end=Math.min(sequence.length,start+80), rows=features.filter(f=>f.end>start && f.start<end);
    const labelWidth=Math.min(64,Math.max(16,...rows.map(f=>f.name.length+2)));
    const lines=[`${"Promoter".padEnd(labelWidth)} ${sequence.slice(start,end)}`];
    for(const f of rows) {
      const lo=Math.max(start,f.start), hi=Math.min(end,f.end);
      lines.push(`${f.name.slice(0,labelWidth).padEnd(labelWidth)} ${" ".repeat(lo-start)}${sequence.slice(lo,hi)}`);
    }
    blocks+=`<h3>Promoter bases ${start+1}–${end}</h3><pre>${escapeHTML(lines.join("\n"))}</pre>`;
  }
  return reportDocument("Promoter sequence alignment",`<p>Run ${escapeHTML(run.app.run_id)} · ${escapeHTML(run.gene?.symbol || "Pasted promoter")} · ${escapeHTML(run.promoter.assembly || "assembly unavailable")}</p><p>Transcript-oriented promoter DNA with selected saRNA targets, displayed motif predictions, imported matches and named regions aligned to shared coordinates. Motif rows show promoter DNA at predicted sites; they are not TF protein sequences. All rows read in promoter orientation. Strand metadata is retained below. Includes the full promoter window.</p>${regionTable(run,features)}${blocks}`);
}
function summaryDocument(run,candidates) {
  const p=run.promoter, rules=run.configuration.rules, screens=run.screenings || [];
  const detailRows=candidates.map(c=>`<tr><th>${escapeHTML(c.id)}<br><small>Rank ${c.rank ?? "excluded"}</small></th><td><code>${escapeHTML(c.sense_strand_rna)}</code><br><code>${escapeHTML(c.antisense_strand_rna)}</code><br><small>RNA cores, 5′ → 3′; overhangs excluded</small></td><td>${escapeHTML(positionText(c))}<br>${c.gc_percent}% GC</td><td>${escapeHTML(ratingText(c))}<br>${c.rating.unavailable} unavailable preferences; ${c.rating.unchecked_gates.length} unchecked gates<br>${escapeHTML(c.rejection_summary || "No evaluated mandatory gate failed.")}</td><td>${escapeHTML(c.guide_prediction.summary)}<br>${escapeHTML(c.rating.guide_model)}</td></tr>`).join("");
  const sources=rules.flatMap(r=>r.provenance?[r.provenance]:[]);
  const candidateIds=new Set(candidates.map(c=>c.id));
  const screeningRows=(run.screening_summary || []).filter(row=>candidateIds.has(row.candidate));
  const screeningTable=`<table><tr><th>Candidate / scope</th><th>Status / reference</th><th>Exact / 18 / 17 / 16 matches</th><th>Coverage and date</th></tr>${screeningRows.map(row=>`<tr><th>${escapeHTML(row.candidate)} · ${escapeHTML(row.scope)}</th><td>${escapeHTML(row.status)}<br>${escapeHTML(row.backend)} ${escapeHTML(row.reference)}</td><td>${[19,18,17,16].map(n=>row[`matches_${n}_of_19`] ?? "unavailable").join(" / ")}</td><td>Complete enumeration: ${row.enumeration_complete ?? "unknown"}; capped details: ${row.alignments_capped ?? "unknown"}<br>${escapeHTML(row.searched_at || "not screened")}</td></tr>`).join("")}</table>`;
  const fullRules=candidates.map(c=>`<details><summary>${escapeHTML(c.id)} · individual rule results</summary><table><tr><th>Rule</th><th>Role</th><th>Outcome</th><th>Explanation</th></tr>${(c.rules || []).map(r=>`<tr><td>${escapeHTML(r.id)}</td><td>${escapeHTML(r.kind)}</td><td>${escapeHTML(r.outcome)}</td><td>${escapeHTML(r.detail)}</td></tr>`).join("")}</table></details>`).join("");
  return reportDocument("saRNA design summary",`
<p>${escapeHTML(run.gene?.symbol || "Pasted promoter")} · Run ${escapeHTML(run.app.run_id)} · ${candidates.length} selected candidates</p>
<p>Research design checklist, not a prediction of activation. Mandatory gates determine inclusion; preferences count equally, followed by the configured thermodynamic tie-break. Unavailable checks are not passes. Off-target sequence matches do not establish biological specificity.</p>
<h2>Promoter and reference</h2><p>${p.length} bases · ${escapeHTML(p.assembly || "Assembly unavailable")} · ${escapeHTML(p.contig || "Contig unavailable")} · genomic interval ${p.start ?? "?"}–${p.end ?? "?"} (0-based, half-open) · strand ${escapeHTML(p.strand || "unknown")} · first-base TSS-relative position ${p.tss_start ?? "unavailable"}.</p><p>Reference: ${escapeHTML(run.transcript?.ensembl_id || run.app.anchor?.label || "Pasted sequence")}</p><details><summary>Reference and input settings</summary><pre>${escapeHTML(JSON.stringify({transcript:run.transcript,anchor:run.app.anchor,request:run.app.request},null,2))}</pre></details>
<h2>Selected designs</h2><table><thead><tr><th>Candidate</th><th>Sense / antisense</th><th>Position / GC</th><th>Checklist</th><th>Guide prediction</th></tr></thead><tbody>${detailRows}</tbody></table>${fullRules}
<h2>Applied rules</h2><p>${escapeHTML(run.configuration.preset_name || "Saved configuration")}${run.configuration.is_custom?" · custom settings: departures from baseline are recorded below":""}. Counts can have different assessable denominators. Compare rule outcomes and missing context, and use matching configurations for comparisons.</p><table><tr><th>Rule</th><th>Role</th><th>Rationale / source</th></tr>${rules.map(r=>`<tr><td>${escapeHTML(r.name)}</td><td>${escapeHTML(r.kind)}</td><td>${escapeHTML(r.rationale || r.description || "")} ${escapeHTML(r.provenance?.display || r.provenance?.source || "")}</td></tr>`).join("")}</table><details><summary>Full rule settings and departures from baseline</summary><pre>${escapeHTML(JSON.stringify(run.configuration,null,2))}</pre></details>
<h2>Off-target evidence</h2><p>${screens.length} saved search record(s). Missing, pending, failed or capped screening evidence remains incomplete; no returned hits is not proof of safety. Sequence-only imported assessments have no assigned off-target screen.</p>${candidates.map(c=>`<p><b>${escapeHTML(c.id)}</b>: ${escapeHTML(screeningSentence(c.id))}</p>`).join("")}${screeningTable}<details><summary>Screening methods and search history</summary><pre>${escapeHTML(JSON.stringify({methods:run.screening_export_notes,searches:screens.map(s=>({id:s.id,status:s.status,scope:s.scope,backend:s.backend,error:s.error,reference:s.result?.reference,parameters:s.result?.parameters,searched_at:s.result?.searched_at}))},null,2))}</pre></details>
<h2>Regulatory context</h2><p>Includes displayed predicted TF motifs, selected targets, imported matches and researcher annotations. Proximity and motif matches do not establish TF occupancy or activation.</p>${regionTable(run,alignedFeatures(run,candidates))}
<h2>Research notes and imported assessments</h2><details><summary>Saved annotations, shortlist and imported analyses</summary><pre>${escapeHTML(JSON.stringify({workspace:run.workspace,sequence_comparisons:run.sequence_comparisons},null,2))}</pre></details>
<h2>Sources and reproducibility</h2><p>Engine ${escapeHTML(run.app.engine_version)}. Download whole-run JSON for the complete machine-readable record and retained off-target alignment details.</p><details><summary>Rule citations and software fingerprints</summary><pre>${escapeHTML(JSON.stringify({sources,software:run.reproducibility},null,2))}</pre></details>`);
}
async function exportInvestigation(kind) {
  const current=state.result;
  try {
    if(state.investigationSaving)throw new Error("Wait for map settings to finish saving before exporting.");
    if(typeof flushWorkspace==="function")await flushWorkspace(current.app.run_id);
    const run=frozen?{...current,screenings:frozenScreeningRecords(),screening_summary:frozenScreeningSummary()}:await api(`/api/runs/${current.app.run_id}/export.json`);
    if(state.result!==current)throw new Error("The active run changed. Export again from the intended run.");
    const candidates=investigationSelection(run), name=`${run.gene?.symbol || "promoter"}-${run.app.run_id.slice(0,8)}-${kind}`;
    let text,extension="html",mime="text/html";
    if(kind==="regions") {
      const p=run.promoter;
      const rows=[["name","type","promoter_start_0based","promoter_end_exclusive","orientation","assembly","contig","genomic_start_0based","genomic_end_exclusive","promoter_genomic_strand","source"],...alignedFeatures(run,candidates).map(f=>[f.name,f.kind,f.start,f.end,f.orientation,p.assembly,p.contig,p.start==null?"":p.strand==="-"?p.end-f.end:p.start+f.start,p.start==null?"":p.strand==="-"?p.end-f.start:p.start+f.end,p.strand,f.source])];
      // Prefix spreadsheet formula introducers while retaining the visible source text.
      text=rows.map(row=>row.map(value=>{let cell=String(value ?? "");if(/^[=+@\-\t\r]/.test(cell))cell="'"+cell;return '"'+cell.replace(/"/g,'""')+'"';}).join(",")).join("\r\n");extension="csv";mime="text/csv";
    } else text=kind==="summary"?summaryDocument(run,candidates):alignmentDocument(run,candidates);
    downloadBlob(new Blob([text],{type:mime+";charset=utf-8"}),`${name}.${extension}`);
  } catch(failure){error(failure.message);}
}
for(const kind of ["summary","alignment","regions"]) $("export-"+kind).addEventListener("click",()=>exportInvestigation(kind));

$("evaluate-existing").addEventListener("click",()=>{
  setWorkspaceTab("candidates",{focus:true});
  revealDisclosure($("sequence-import-controls"));
  $("sequence-import-controls").open=true;
  $("sequence-input").focus();
  $("sequence-input").scrollIntoView({block:"center"});
});
function revealDisclosure(target) {
  if (!target) return;
  for (let node=target;node;node=node.parentElement) if(node.tagName==="DETAILS")node.open=true;
}
function renderCalmSummary() {
  const run=state.result;
  if (!run) return;
  const selected=activeSelection(), checked=state.checked.size;
  const decisions=typeof screenedDecisions==="function"?screenedDecisions():[];
  const complete=decisions.filter(d=>d.status!=="pending").length;
  const unchecked=new Set(run.candidates.flatMap(c=>c.rating.unchecked_gates));
  $("candidate-summary").textContent=`${fmt(run.counts.passing)} candidates pass available gates · ${complete} with complete local screening evidence · ${checked} checked${unchecked.size?` · ${unchecked.size} gate(s) unavailable`:""}`;
  $("screening-selection-title").textContent=!selected.length?"No candidates selected":checked?`${checked} checked candidates`:`${selected.length} design-passing candidates`;
  $("screening-selection-note").textContent=selected.length?`Shared selection with Candidates and Export. Searches cover transcripts and genome.${selected.some(c=>c.rank===null)?" Candidates excluded by design gates are not screened.":""}`:"Check candidates in the candidate list to begin screening.";
  const jobs=Object.values(state.screeningJobs || {}), busy=jobs.some(j=>["running","waiting"].includes(j.status));
  const remote=screeningBackend()==="ncbi";
  const screenCount=selected.filter(c=>c.rank!==null).length;
  $("screen-selection").textContent=busy?"Screening in progress":checked?`Screen selected (${screenCount})`:`Screen all passing (${selected.length})`;
  $("screen-selection").disabled=!!frozen || !selected.some(c=>c.rank!==null) || busy || (!remote && !REFERENCE_STATUS?.ready);
  $("screening-submission-note").textContent=remote?"Remote search: clicking Screen sends these sequences to NCBI.":"";
  $("reference-compact").hidden=!!frozen || remote;
  if (run.configuration.is_custom) $("result-config").classList.add("custom-configuration");
  else $("result-config").classList.remove("custom-configuration");
}
document.addEventListener("click",event=>{
  const jump=event.target.closest("[data-open-disclosure]");
  if(jump){const target=$(jump.dataset.openDisclosure);revealDisclosure(target);target?.scrollIntoView({block:"start"});}
  const tool=event.target.closest("[data-map-tool]");
  if(tool){const target=$(tool.dataset.mapTool), open=!target.open;
    for(const panel of document.querySelectorAll(".map-tool-panel"))panel.open=panel===target && open;
    if(open)revealDisclosure(target);
  }
});
for(const panel of document.querySelectorAll(".map-tool-panel")) panel.addEventListener("toggle",()=>{
  for(const button of document.querySelectorAll("[data-map-tool]")) if(button.dataset.mapTool===panel.id){button.setAttribute("aria-expanded",String(panel.open));if(!panel.open && panel.contains(document.activeElement))button.focus({preventScroll:true});}
  if(state.result)renderMap();
});
$("screen-selection").addEventListener("click",async()=>{
  if($("screen-selection").disabled)return;
  const mode=state.checked.size?"selection":"all";
  await Promise.all([startScreening("transcript",mode),startScreening("genomic",mode)]);
});
// Opening a record from history should return to the workbench immediately.
$("saved-runs").addEventListener("change",()=>{if($("saved-runs").value)$("history-dialog").close();});
// Wrapped navigation is taller on narrow screens; anchor jumps must clear its actual height.
function syncWorkspaceScrollOffset() {
  document.documentElement.style.setProperty("--workspace-scroll-offset",`${Math.ceil($("workflow-nav").getBoundingClientRect().height)+16}px`);
}
new ResizeObserver(syncWorkspaceScrollOffset).observe($("workflow-nav"));
syncWorkspaceScrollOffset();
const groupExportButtons = $("export-batch-json").parentElement;
$("group-export-actions").append(groupExportButtons);
if (frozen) { $("open-history").hidden=true; }
setWorkspaceTab(frozen ? "candidates" : "design", {scroll:false});
if (state.result) workflowRun(state.result); else workflowStatus();
