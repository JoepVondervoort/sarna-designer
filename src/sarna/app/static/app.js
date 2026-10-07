/* Display and controls only. Scientific evaluation and ranking always run in Python. */
"use strict";
const $ = (id) => document.getElementById(id);
const escapeHTML = (s) => String(s ?? "").replace(/[&<>"']/g, c => ({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;","'":"&#39;"}[c]));
const signed = (n) => n > 0 ? `+${n}` : String(n);
const fmt = (n) => Number(n).toLocaleString();
const clone = (v) => JSON.parse(JSON.stringify(v));
const ordered = v => Array.isArray(v)?v.map(ordered):v && typeof v==="object"?Object.fromEntries(Object.keys(v).sort().map(k=>[k,ordered(v[k])])):v;
const frozen = window.SARNA_FROZEN || null;
const state = {catalogue:null, source:"gene", gene:null, draft:{}, base:{}, result:null, selected:null,
  checked:new Set(), filter:"passing", page:0, sort:"rank", descending:false, start:0, end:0, busy:false, retry:null,
  customReferences:{}, basket:[], batch:null, editingReference:null};
const roleName = {mandatory:"Gate", preference:"Preference", tiebreak:"Tiebreak", off:"Off"};
const tableNames = {"rna-xia-1998":"RNA · Xia 1998", "dna-santalucia-1998":"DNA · SantaLucia 1998", "legacy-dna-kj":"Legacy DNA · defective"};

async function api(path, options) {
  const response = await fetch(path, options);
  const data = await response.json();
  if (!response.ok) {
    const detail = Array.isArray(data.detail) ? data.detail.map(e => `${e.loc.join(".")}: ${e.msg}`).join("; ") : data.detail;
    throw Object.assign(new Error(detail || `Request failed (${response.status})`), {code:data.code});
  }
  return data;
}
function error(message, retry = null) {
  $("global-error").textContent = typeof message === "string" ? message : message.message;
  $("global-error").hidden = !message;
  state.retry = retry;
  $("error-recovery").hidden = !retry || message?.code !== "offline-cache-miss" || !!state.catalogue?.data_access?.offline_forced;
}
function syncDataAccess() {
  const forced = !!state.catalogue?.data_access?.offline_forced;
  if (forced) $("offline").checked = true;
  $("offline").disabled = forced;
  $("data-access-status").textContent = forced ?
    "Cache only — enforced by SARNA_OFFLINE. To fetch new references, unset it in the launch terminal and restart the app." : $("offline").checked ?
    "Cache only — no Ensembl requests. A gene, transcript or promoter window that has not been cached will need online retrieval." :
    "Online retrieval allowed — missing public references are fetched from Ensembl and cached. Design calculations run locally.";
  document.querySelector(".data-access").classList.toggle("cache-only", $("offline").checked);
}
function presetData() { return state.catalogue.presets.find(p => p.id === $("preset").value); }
function setPreset(id, configuration = null) {
  $("preset").value = id;
  const preset = presetData();
  state.base = {};
  for (const def of state.catalogue.rules) {
    const baseline = preset.configuration.rules.find(r => r.id === def.id);
    state.base[def.id] = {role:baseline?.kind || "off", values:Object.fromEntries((baseline || def).parameters.map(p => [p.name,p.value])), option:baseline?.selected_option?.id || ""};
  }
  state.draft = clone(state.base);
  if (configuration) {
    for (const r of configuration.rules) state.draft[r.id] = {role:r.kind, values:Object.fromEntries(r.parameters.map(p => [p.name,p.value])), option:r.selected_option?.id || ""};
    for (const r of configuration.disabled_rules) if (state.draft[r.id]) state.draft[r.id].role = "off";
  }
  $("preset-description").textContent = preset.description;
  syncQuick(); renderEditor();
}
function syncQuick() {
  const v = state.draft;
  $("gc-min").value = v["gc-content"].values["minimum-gc-fraction"] * 100;
  $("gc-max").value = v["gc-content"].values["maximum-gc-fraction"] * 100;
  $("window-start").value = v["tss-window"].values["window-start-bp"];
  $("window-end").value = v["tss-window"].values["window-end-bp"];
  for (const id of ["gc-min","gc-max"]) $(id).disabled = v["gc-content"].role === "off";
  for (const id of ["window-start","window-end"]) $(id).disabled = v["tss-window"].role === "off";
  const count = Object.keys(v).filter(id => JSON.stringify(v[id]) !== JSON.stringify(state.base[id])).length;
  $("custom-count").textContent = count ? `${count} edited ↗` : `${state.catalogue.rules.length} rules ↗`;
}
function markDraft() {
  $("draft-badge").textContent = "DRAFT";
  if (state.result) {
    $("draft-notice").hidden = false;
    $("design-button").innerHTML = 'Apply changes &amp; design <span>→</span>';
  }
  if (state.batch) $("batch-draft-notice").hidden = false;
}
function roleChoices(def) {
  if (def.id === "target-length-19") return ["mandatory"];
  if (def.default_kind === "tiebreak") return ["tiebreak","off"];
  return ["mandatory","preference","off"];
}
function renderEditor() {
  $("rule-editor").innerHTML = state.catalogue.rules.map(def => {
    const d = state.draft[def.id], b = state.base[def.id];
    const params = def.parameters.map(p => {
      const value = d.values[p.name];
      const field = p.choices ? `<select data-param="${escapeHTML(p.name)}">${p.choices.map(v => `<option value="${escapeHTML(v)}" ${v===value?"selected":""}>${escapeHTML(tableNames[v] || v)}</option>`).join("")}</select>` :
        `<input data-param="${escapeHTML(p.name)}" type="number" value="${escapeHTML(value)}" ${p.minimum!==null?`min="${p.minimum}"`:""} ${p.maximum!==null?`max="${p.maximum}"`:""} step="${Number.isInteger(p.default)?1:"any"}">`;
      return `<div><label>${escapeHTML(p.description)} ${p.unit?`(${escapeHTML(p.unit)})`:""}</label>${field}<p class="preset-baseline">Preset: ${escapeHTML(b.values[p.name])} · supported range ${escapeHTML(p.choices?.join(", ") || `${p.minimum} to ${p.maximum}`)}</p></div>`;
    }).join("");
    const selected = def.options.find(o => o.id === d.option);
    return `<article class="rule-editor-item" data-rule="${def.id}"><div class="rule-editor-heading"><div><h3>${escapeHTML(def.name)}</h3><p class="source">${escapeHTML(def.provenance.display)} <span>${escapeHTML(def.evidence_class)}</span></p><p class="preset-baseline">Role in selected preset: ${roleName[b.role]}</p></div><select data-role aria-label="Role for ${escapeHTML(def.name)}">${roleChoices(def).map(role => `<option value="${role}" ${role===d.role?"selected":""}>${roleName[role]}</option>`).join("")}</select></div><p>${escapeHTML(def.rationale)}</p><p class="consequence"><b>Changing this:</b> ${escapeHTML(def.consequence)}</p>${def.options.length?`<label class="rule-option-label">Sourced parameter alternative</label><select data-option ${d.role==="off"?"disabled":""}><option value="">Preset / manual values</option>${def.options.map(o => `<option value="${o.id}" ${o.id===d.option?"selected":""}>${escapeHTML(o.label)}</option>`).join("")}</select><p class="option-note">${selected?escapeHTML(`${selected.provenance.display}. ${selected.note}`):"Manual values are recorded as researcher choices, without assigning a paper as their source."}</p>`:""}<div class="rule-parameters">${params}</div></article>`;
  }).join("");
  for (const article of $("rule-editor").querySelectorAll("[data-rule]")) {
    if (state.draft[article.dataset.rule].role === "off") for (const field of article.querySelectorAll("[data-param]")) field.disabled = true;
  }
}

function referenceLabel(key) {
  const [kind,id] = key.split(":");
  if (kind === "custom") return state.customReferences[id]?.name || id;
  if (kind === "anchor") return state.gene?.anchors?.find(a=>a.id===id)?.label || id;
  const t = state.gene?.transcripts?.find(t=>t.id===id);
  return t ? `${t.id}${t.mane?" · MANE Select":""}` : id;
}
function populateGene(gene, selected = "") {
  state.gene = gene; $("gene-symbol").value = gene.symbol;
  const transcripts = gene.transcripts || [], anchors = gene.anchors || [];
  $("gene-info").innerHTML = `<strong>${escapeHTML(gene.symbol)} <span class="muted">${escapeHTML(gene.id)}</span></strong>${escapeHTML(gene.description)}<p class="hint">${escapeHTML(gene.assembly || "GRCh38")} · genomic ${escapeHTML(gene.strand)} strand · ${fmt(transcripts.length)} transcripts<br>${escapeHTML(gene.source)}</p>`;
  $("gene-info").hidden = false;
  const groups = gene.tss_groups || [];
  const choices = groups.length ? groups.map(g=>`<optgroup label="Annotated +1 at ${escapeHTML(g.contig)}:${fmt(g.tss+1)} · ${g.transcript_ids.length} transcript(s)">${g.transcript_ids.map(id=>{
    const t = transcripts.find(t=>t.id===id);
    return `<option value="transcript:${escapeHTML(id)}">${escapeHTML(id)}${t?.mane?" · MANE Select":t?.canonical?" · canonical":""} · ${escapeHTML(t?.biotype)}</option>`;
  }).join("")}</optgroup>`).join("") : `<optgroup label="Annotated transcripts">${transcripts.map(t=>`<option value="transcript:${escapeHTML(t.id)}">${escapeHTML(t.id)}${t.mane?" · MANE Select":t.canonical?" · canonical":""} · ${escapeHTML(t.biotype)}</option>`).join("")}</optgroup>`;
  const anchorChoices = anchors.map(a=>`<option value="anchor:${escapeHTML(a.id)}">${escapeHTML(a.label)} · ${escapeHTML(a.reference_kind)}</option>`).join("");
  const customChoices = Object.values(state.customReferences).map(a=>`<option value="custom:${escapeHTML(a.id)}">${escapeHTML(a.name)} · researcher supplied</option>`).join("");
  $("reference").innerHTML = `<option value="">Choose a promoter reference…</option>${choices}${anchorChoices?`<optgroup label="Declared literature anchors">${anchorChoices}</optgroup>`:""}${customChoices?`<optgroup label="Researcher references">${customChoices}</optgroup>`:""}`;
  $("reference").disabled = false;
  $("reference").value = [...$("reference").options].some(o=>o.value===selected) && selected ? selected : gene.default_transcript ? `transcript:${gene.default_transcript}` : "";
  $("add-promoter").disabled = !$("reference").value;
  $("save-custom-reference").disabled = false;
  $("custom-reference-locus").textContent = `GRCh38 · chr${gene.contig || "?"} · gene strand ${gene.strand}. Enter a 1-based genomic coordinate, not a TSS-relative position.`;
  renderGeneOverview(); renderBasket();
}
function renderGeneOverview() {
  const gene = state.gene;
  const groups = gene?.tss_groups || [];
  if (!gene || !Number.isFinite(gene.start) || !Number.isFinite(gene.end)) { $("gene-overview").hidden = true; return; }
  $("gene-overview").hidden = state.source !== "gene" && !frozen;
  $("gene-overview-meta").textContent = frozen ? "FROZEN REFERENCE CHOICES" : "DRAFT REFERENCE CHOICES";
  $("gene-overview-description").textContent = `${gene.symbol} · ${groups.length} annotated start coordinate(s) across ${(gene.transcripts || []).length} transcripts. Click a start-site group to inspect its transcripts below.`;
  const anchors = gene.anchors || [], custom = Object.values(state.customReferences);
  const lo = Math.min(gene.start,...anchors.map(a=>a.start),...custom.map(a=>a.tss_position-1));
  const hi = Math.max(gene.end,...anchors.map(a=>a.end),...custom.map(a=>a.tss_position));
  const span = Math.max(1,hi-lo), left=130, right=966;
  const scale = n=>left+(gene.strand === "-" ? hi-1-n : n-lo)/span*(right-left);
  const tip = value=>`data-tip="${escapeHTML(value)}" tabindex="0" aria-label="${escapeHTML(value)}"`;
  const parts=[];
  for (let i=0;i<=4;i++) {
    const coordinate = Math.round(gene.strand === "-" ? hi-1-(span-1)*i/4 : lo+(span-1)*i/4);
    parts.push(`<text x="${left+(right-left)*i/4}" y="22" text-anchor="${i===0?"start":i===4?"end":"middle"}" fill="#61767a" font-size="10">${fmt(coordinate+1)}</text>`);
  }
  parts.push(`<text x="12" y="62" fill="#61767a" font-size="10">Gene / start sites</text><line x1="${left}" x2="${right}" y1="61" y2="61" stroke="#dce5e3"/><rect x="${Math.min(scale(gene.start),scale(gene.end-1))}" y="57" width="${Math.max(2,Math.abs(scale(gene.end-1)-scale(gene.start)))}" height="8" rx="3" fill="#dce9e3" ${tip(`${gene.symbol} · chr${gene.contig}:${gene.start}–${gene.end} (0-based, half-open). Axis labels are 1-based genomic coordinates; plotted in transcript orientation.`)}/>`);
  groups.forEach((g,i)=>{
    const x=scale(g.tss), y=39+(i%3)*7;
    const help=`Annotated start at chr${g.contig}:${g.tss+1} (1-based) · ${g.transcript_ids.join(", ")}${g.mane?" · includes MANE Select":""}. Ensembl annotation does not establish activity in your cell type.`;
    parts.push(`<line x1="${x}" x2="${x}" y1="${y+4}" y2="69" stroke="${g.mane?"#087e80":"#526e68"}"/><circle cx="${x}" cy="${y}" r="${g.mane?5:4}" fill="${g.mane?"#087e80":"#526e68"}" ${tip(help)} data-tss-group="${escapeHTML(g.id)}"/>`);
  });
  parts.push('<text x="12" y="107" fill="#61767a" font-size="10">Declared anchors</text>');
  anchors.forEach((a,i)=>{
    const x1=scale(a.start),x2=scale(a.end-1), y=87+(i%3)*9;
    const key=`anchor:${a.id}`, help=`${a.label} · ${a.reference_kind}. ${a.provenance}. ${a.note} ${a.caveat || ""}`;
    parts.push(`<rect x="${Math.min(x1,x2)}" y="${y}" width="${Math.max(3,Math.abs(x2-x1))}" height="7" rx="2" fill="#ddb370" ${tip(help)} data-reference="${escapeHTML(key)}"/>`);
    if (a.reference_point !== null) parts.push(`<line x1="${scale(a.reference_point)}" x2="${scale(a.reference_point)}" y1="${y-3}" y2="${y+10}" stroke="#9a6a30" ${tip(help)} data-reference="${escapeHTML(key)}"/>`);
  });
  parts.push('<text x="12" y="145" fill="#61767a" font-size="10">Your references</text>');
  custom.forEach((a,i)=>parts.push(`<circle cx="${scale(a.tss_position-1)}" cy="${135+(i%2)*9}" r="4" fill="#8e7bb9" ${tip(`${a.name} · chr${gene.contig}:${a.tss_position} (1-based) · researcher supplied, inferred. Source: ${a.source}. ${a.note}`)} data-reference="custom:${escapeHTML(a.id)}"/>`));
  parts.push(`<text x="${right}" y="173" text-anchor="end" fill="#61767a" font-size="10">chr${escapeHTML(gene.contig)} · 1-based genomic coordinates · gene strand ${escapeHTML(gene.strand)} · 5′ → 3′</text>`);
  $("gene-overview-map").innerHTML=parts.join("");
  $("tss-group-list").innerHTML=groups.map(g=>`<details><summary>Annotated +1 · ${escapeHTML(g.contig)}:${fmt(g.tss+1)} <span class="muted">${g.transcript_ids.length} transcript(s)${g.mane?" · MANE Select":""}</span></summary><p class="hint">Shared start coordinate, not a measurement of promoter activity.</p><div>${g.transcript_ids.map(id=>`<button type="button" class="text-button" data-choose-transcript="${escapeHTML(id)}">${escapeHTML(id)}</button>`).join("")}</div></details>`).join("");
}
async function lookup() {
  const symbol = $("gene-symbol").value.trim().toUpperCase();
  if (!symbol) throw new Error("Enter a human gene symbol or Ensembl gene ID first.");
  const previousReference = state.gene?.symbol.toLowerCase()===symbol.toLowerCase() ? $("reference").value : "";
  $("lookup-gene").disabled = true; $("lookup-gene").textContent = "Finding…";
  try {
    const gene = await api(`/api/genes/${encodeURIComponent(symbol)}?offline=${$("offline").checked}`);
    if (state.gene && state.gene.symbol !== gene.symbol) { state.basket=[]; state.customReferences={}; }
    populateGene(gene,previousReference);
    markDraft(); return gene;
  } finally { $("lookup-gene").disabled = false; $("lookup-gene").textContent = "Look up"; }
}
function renderBasket() {
  $("promoter-selection").hidden = state.source !== "gene" || !state.basket.length;
  $("promoter-count").textContent = `(${state.basket.length})`;
  $("promoter-basket").innerHTML = state.basket.map((entry,i)=>`<article class="basket-item" data-basket-index="${i}"><div class="basket-heading"><button type="button" class="text-button" data-basket-action="select">${escapeHTML(referenceLabel(entry.key))}</button><button type="button" class="text-button remove-reference" data-basket-action="remove" aria-label="Remove ${escapeHTML(referenceLabel(entry.key))}">×</button></div><p class="hint">${entry.pinned?`Pinned criteria · ${escapeHTML(entry.pinned.preset_id)} · ${escapeHTML(entry.pinned.cell_line)}`:"Use current criteria at submission"}</p><div class="pair"><div><label>Upstream (bp)</label><input type="number" data-basket-span="upstream" value="${entry.upstream}" min="0" max="20000" step="1" aria-label="Upstream span for ${escapeHTML(referenceLabel(entry.key))}"></div><div><label>Downstream (bp)</label><input type="number" data-basket-span="downstream" value="${entry.downstream}" min="0" max="20000" step="1" aria-label="Downstream span for ${escapeHTML(referenceLabel(entry.key))}"></div></div><div class="basket-tools"><button type="button" class="text-button" data-basket-action="capture">Capture current criteria</button>${entry.pinned?'<button type="button" class="text-button" data-basket-action="share">Use shared criteria</button>':""}${entry.key.startsWith("custom:")?'<button type="button" class="text-button" data-basket-action="edit">Edit reference</button>':""}</div>${entry.pinned?'<details><summary class="hint">Captured settings</summary><pre class="pinned-settings">'+escapeHTML(JSON.stringify({preset:entry.pinned.preset_id,roles:entry.pinned.roles,overrides:entry.pinned.overrides,options:entry.pinned.options,cell_line:entry.pinned.cell_line,regulatory:entry.pinned.regulatory,motifs:entry.pinned.motifs,motif_shuffles:entry.pinned.motif_shuffles,chromatin:entry.pinned.chromatin},null,2))+'</pre></details>':""}</article>`).join("");
}
function addPromoter(key = $("reference").value) {
  if (!key || !state.gene) return;
  if (!state.basket.some(entry=>entry.key===key)) state.basket.push({key,upstream:Number($("upstream").value),downstream:Number($("downstream").value),pinned:null});
  renderBasket(); markDraft();
}
function saveCustomReference() {
  const name=$("custom-reference-name").value.trim(), source=$("custom-reference-source").value.trim(), note=$("custom-reference-note").value.trim();
  const coordinateText=$("custom-reference-position").value, position=Number(coordinateText);
  if (!state.gene) { error("Look up a gene before adding a genomic reference."); return; }
  if (!name || !source || !note || !coordinateText || !Number.isSafeInteger(position) || position<1) { error("Enter a name, a positive integer genomic coordinate, a source or derivation, and an uncertainty note."); return; }
  const id=state.editingReference || `researcher-${Date.now().toString(36)}-${Math.random().toString(36).slice(2,7)}`;
  state.customReferences[id]={id,name,tss_position:position,source,note};
  populateGene(state.gene,`custom:${id}`); addPromoter(`custom:${id}`);
  resetCustomEditor(); $("custom-reference-editor").open=false; error("");
}
function resetCustomEditor() {
  state.editingReference=null;
  for (const id of ["custom-reference-name","custom-reference-position","custom-reference-source","custom-reference-note"]) $(id).value="";
  $("save-custom-reference").textContent="Add researcher reference";
  $("cancel-custom-reference").hidden=true;
}
function editCustomReference(key) {
  const id=key.slice(7), value=state.customReferences[id]; if (!value) return;
  state.editingReference=id;
  for (const [field,property] of [["name","name"],["position","tss_position"],["source","source"],["note","note"]]) $("custom-reference-"+field).value=value[property];
  $("save-custom-reference").textContent="Save reference changes"; $("cancel-custom-reference").hidden=false;
  $("custom-reference-editor").open=true; $("custom-reference-name").focus();
}
function applyReference(request,key) {
  const [kind,id]=key.split(":");
  request.transcript_id=kind==="transcript"?id:null;
  request.anchor_id=kind==="anchor"?id:null;
  request.custom_reference=kind==="custom"?clone(state.customReferences[id]):null;
  return request;
}
function requestFromDraft() {
  const roles = {}, overrides = {}, options = {};
  for (const [id,d] of Object.entries(state.draft)) {
    const b = state.base[id];
    if (d.role !== b.role) roles[id] = d.role;
    if (d.role === "off") continue;
    const def = state.catalogue.rules.find(r => r.id === id);
    const option = def.options.find(o => o.id === d.option);
    if (option) options[id] = option.id;
    for (const [name,value] of Object.entries(d.values)) {
      if (value !== (option && name in option.values ? option.values[name] : b.values[name])) {
        (overrides[id] ||= {})[name] = value;
      }
    }
  }
  const [kind,reference] = $("reference").value.split(":");
  return {source:state.source, symbol:state.source==="gene"?$("gene-symbol").value.trim():"", sequence:state.source==="sequence"?$("dna-sequence").value:"",
    transcript_id:state.source==="gene" && kind==="transcript"?reference:null, anchor_id:state.source==="gene" && kind==="anchor"?reference:null,
    custom_reference:state.source==="gene" && kind==="custom"?clone(state.customReferences[reference]):null,
    upstream:Number($("upstream").value), downstream:Number($("downstream").value),
    tss_start:state.source==="sequence" && $("tss-start").value!==""?Number($("tss-start").value):null,
    cell_line:$("cell-line").value.trim() || "Not specified", preset_id:$("preset").value, roles,overrides,options,
    offline:$("offline").checked, regulatory:$("regulatory").checked, motifs:$("motifs").checked,
    motif_shuffles:Number($("motif-shuffles").value), chromatin:$("chromatin").checked};
}
async function design() {
  if (state.busy) return;
  error(""); state.busy = true; $("design-button").disabled = true; $("setup").inert = true; $("saved-runs").disabled = true;
  try {
    if (state.source === "gene") {
      if (!state.gene || state.gene.symbol.toLowerCase() !== $("gene-symbol").value.trim().toLowerCase()) await lookup();
      if (!$("reference").value) throw new Error("Select a transcript, a declared anchor, or your genomic reference. A gene without a unique MANE Select requires a deliberate choice.");
    }
    const previous = state.result;
    const queued = await api("/api/runs", {method:"POST",headers:{"Content-Type":"application/json"},body:JSON.stringify(requestFromDraft())});
    for (;;) {
      const run = await api(`/api/runs/${queued.id}`);
      $("job-status").textContent = run.stage;
      if (run.status === "error") throw Object.assign(new Error(run.error), {code:run.error_code});
      if (run.status === "complete") { state.batch=null; $("batch-comparison").hidden=true; showResult(run.result); renderDiff(previous,run.result); break; }
      await new Promise(resolve => setTimeout(resolve,700));
    }
    await refreshRuns();
  } catch (e) { error(e,"design"); $("job-status").textContent = "Run not completed. Your previous saved results are retained."; }
  finally { state.busy = false; $("design-button").disabled = false; $("setup").inert = false; $("saved-runs").disabled = false; }
}
function setBatchBusy(busy) {
  state.busy=busy; $("setup").inert=busy; $("saved-runs").disabled=busy;
  $("design-button").disabled=busy; $("design-batch").disabled=busy; $("batch-retry-online").disabled=busy;
}
function batchRequest(entry,current) {
  const request=clone(entry.pinned || current);
  Object.assign(request,{source:"gene",symbol:state.gene.symbol,sequence:"",tss_start:null,upstream:entry.upstream,downstream:entry.downstream,offline:$("offline").checked});
  return applyReference(request,entry.key);
}
async function designBatch() {
  if (state.busy || !state.basket.length || state.source!=="gene") return;
  error(""); setBatchBusy(true);
  try {
    if (!state.gene) throw new Error("Look up a gene and choose the references to compare.");
    const current=requestFromDraft(), requests=state.basket.map(entry=>batchRequest(entry,current));
    const queued=await api("/api/batches",{method:"POST",headers:{"Content-Type":"application/json"},body:JSON.stringify({label:`${state.gene.symbol} · ${requests.length} promoter references`,requests})});
    await pollBatch(queued.id); await refreshRuns();
  } catch(e) { error(e); $("job-status").textContent="Comparison not completed. Saved runs are retained."; }
  finally { setBatchBusy(false); }
}
function batchReference(run) {
  const request=run.request || {}, anchor=run.summary?.anchor || run.result?.app.anchor;
  return request.custom_reference?.name || anchor?.label || request.anchor_id || request.transcript_id || "Reference";
}
function renderBatch(batch) {
  state.batch=batch; $("batch-comparison").hidden=false; $("batch-draft-notice").hidden=true;
  $("batch-meta").textContent=`GROUP ${batch.id.slice(0,8)} · ${batch.label || batch.symbol} · ${new Date(batch.created).toLocaleString()}`;
  $("batch-status").textContent=batch.status.toUpperCase();
  const runs=batch.runs || [];
  $("batch-rows").innerHTML=runs.map(run=>{
    const req=run.request || {}, summary=run.summary || {}, counts=summary.counts || run.result?.counts;
    const first=run.result?.candidates.find(c=>c.rank===1), rating=summary.best_rating || first?.rating;
    const config=summary.configuration || (run.result?{preset_name:run.result.configuration.preset.name,is_custom:run.result.configuration.is_custom}:null);
    const label=batchReference(run), kind=req.custom_reference?"researcher supplied / inferred":summary.anchor?.reference_kind || run.result?.app.anchor?.reference_kind || "annotated transcript";
    const ratingLabel=rating && typeof rating==="object"?(rating.evaluated?`${rating.met}/${rating.evaluated}`:"Not evaluated"):rating ?? "—";
    return `<tr class="${state.result?.app.run_id===run.id?"selected":""}"><td><strong>${escapeHTML(label)}</strong><p class="hint">${escapeHTML(kind)}</p>${req.custom_reference?`<p class="hint">+1 at ${fmt(req.custom_reference.tss_position)} (1-based)</p>`:""}</td><td>${fmt(req.upstream)} upstream<br>${fmt(req.downstream)} downstream</td><td>${escapeHTML(config?.preset_name || req.preset_id)}${config?.is_custom?" · custom":""}<p class="hint">${escapeHTML(summary.cell_line || req.cell_line)}</p></td><td>${counts?`${fmt(counts.passing)} / ${fmt(counts.generated)}`:"—"}</td><td>${escapeHTML(ratingLabel)}</td><td>${run.status==="complete"?`<button type="button" class="secondary" data-open-batch-run="${escapeHTML(run.id)}">Open run</button>`:`<span class="tag">${escapeHTML(run.status.toUpperCase())}</span><p class="hint batch-error">${escapeHTML(run.error || run.stage)}</p>`}</td></tr>`;
  }).join("");
  const completed=runs.filter(r=>r.status==="complete").length, failed=runs.filter(r=>r.status==="error").length;
  $("batch-progress").textContent=`${completed} of ${runs.length} completed${failed?` · ${failed} failed; completed runs remain available`:""}${batch.retry_of?` · retry of saved group ${batch.retry_of.slice(0,8)}`:""}`;
  $("batch-retry-online").hidden=!!frozen || batch.status==="running" || !!state.catalogue?.data_access?.offline_forced || !runs.some(r=>r.status==="error" && r.error_code==="offline-cache-miss");
  $("export-batch-html").hidden=!!frozen;
  $("job-status").textContent=batch.status==="running"?"Designing the selected references separately…":`Saved promoter comparison · ${batch.id.slice(0,8)}`;
  if (!state.result) $("empty-state").hidden=true;
}
async function openBatchRun(id,restore = true) {
  const child=state.batch?.runs.find(r=>r.id===id);
  if (!child || child.status!=="complete") return;
  const result=frozen?child.result:(await api(`/api/runs/${id}`)).result;
  if (!result) throw new Error("The saved group does not contain this completed run.");
  if (!frozen && restore) restoreDraft(result);
  else if (frozen) { state.gene=result.app.gene_overview || null; renderGeneOverview(); }
  showResult(result); renderBatch(state.batch);
  if (!frozen && !restore && JSON.stringify(ordered(requestFromDraft()))!==JSON.stringify(ordered(result.app.request))) {
    markDraft(); $("batch-draft-notice").hidden=true;
  }
  if (!frozen) $("saved-runs").value=`batch:${state.batch.id}`;
}
async function pollBatch(id) {
  for (;;) {
    const batch=await api(`/api/batches/${id}`); renderBatch(batch);
    if (batch.status!=="running") {
      const first=batch.runs.find(r=>r.status==="complete");
      if (first) await openBatchRun(first.id,false);
      else { state.result=null; $("result-area").hidden=true; $("empty-state").hidden=true; }
      break;
    }
    await new Promise(resolve=>setTimeout(resolve,700));
  }
}
async function retryBatchOnline() {
  if (state.busy || !state.batch || frozen || state.catalogue?.data_access?.offline_forced) return;
  error(""); setBatchBusy(true);
  try {
    const queued=await api(`/api/batches/${state.batch.id}/retry`,{method:"POST",headers:{"Content-Type":"application/json"},body:JSON.stringify({offline:false})});
    $("offline").checked=false; syncDataAccess();
    await pollBatch(queued.id); await refreshRuns();
  } catch(e) { error(e); }
  finally { setBatchBusy(false); }
}
function restoreBatchDraft(batch) {
  const first=batch.runs.find(r=>r.status==="complete" && r.result);
  if (first) restoreDraft(first.result);
  const requests=batch.request?.requests || batch.runs.map(r=>r.request);
  state.basket=[];
  for (const request of requests) {
    if (request.custom_reference) state.customReferences[request.custom_reference.id]=clone(request.custom_reference);
    const key=request.custom_reference?`custom:${request.custom_reference.id}`:request.anchor_id?`anchor:${request.anchor_id}`:`transcript:${request.transcript_id}`;
    state.basket.push({key,upstream:request.upstream,downstream:request.downstream,pinned:clone(request)});
  }
  if (state.gene?.transcripts) populateGene(state.gene,$("reference").value);
  renderBasket();
}
async function openSavedBatch(id) {
  const batch=await api(`/api/batches/${id}`);
  renderBatch(batch);
  const first=batch.runs.find(r=>r.status==="complete");
  if (first) { await openBatchRun(first.id); restoreBatchDraft(batch); }
  else { state.result=null; $("result-area").hidden=true; $("empty-state").hidden=true; }
  $("saved-runs").value=`batch:${id}`;
}
async function exportBatch(type) {
  if (!state.batch) return;
  try {
    const name=`${state.batch.symbol || "promoters"}-${state.batch.id.slice(0,8)}.${type}`;
    if (frozen) { downloadBlob(new Blob([JSON.stringify(state.batch,null,2)],{type:"application/json"}),name); return; }
    const response=await fetch(`/api/batches/${state.batch.id}/export.${type}`);
    if (!response.ok) throw new Error(`Group export failed (${response.status})`);
    downloadBlob(await response.blob(),name);
  } catch(e) { error(e); }
}
async function refreshRuns() {
  const [runs,batches] = await Promise.all([api("/api/runs"),api("/api/batches")]);
  $("saved-runs").innerHTML = '<option value="">Open a run or group…</option>' + (batches.length?'<optgroup label="Saved promoter groups">'+batches.map(r=>`<option value="batch:${r.id}">${escapeHTML(r.label)} · ${escapeHTML(new Date(r.created).toLocaleString())} · ${escapeHTML(r.status)}</option>`).join("")+'</optgroup>':"") + '<optgroup label="Individual runs">'+runs.map(r => `<option value="${r.id}">${escapeHTML(r.label)} · ${escapeHTML(new Date(r.created).toLocaleString())} · ${escapeHTML(r.status)}</option>`).join("")+'</optgroup>';
  if (state.batch) $("saved-runs").value=`batch:${state.batch.id}`;
  else if (state.result) $("saved-runs").value = state.result.app.run_id;
}
function restoreDraft(result) {
  const request = result.app.request;
  setPreset(request.preset_id,result.configuration);
  state.source = request.source;
  setSource(request.source);
  for (const [id,key] of [["gene-symbol","symbol"],["dna-sequence","sequence"],["upstream","upstream"],["downstream","downstream"],["tss-start","tss_start"],["cell-line","cell_line"],["motif-shuffles","motif_shuffles"]]) $(id).value = request[key] ?? "";
  for (const key of ["regulatory","motifs","offline","chromatin"]) $(key).checked = request[key];
  syncDataAccess();
  if (request.source === "gene") {
    if (state.gene?.symbol !== (result.app.gene_overview?.symbol || request.symbol)) { state.basket=[]; state.customReferences={}; }
    if (request.custom_reference) state.customReferences[request.custom_reference.id]=clone(request.custom_reference);
    const ref = request.custom_reference?`custom:${request.custom_reference.id}`:request.anchor_id ? `anchor:${request.anchor_id}` : `transcript:${request.transcript_id || result.transcript?.ensembl_id}`;
    if (result.app.gene_overview) populateGene(result.app.gene_overview,ref);
    else {
      state.gene = {symbol:request.symbol};
      $("reference").innerHTML = `<option value="${escapeHTML(ref)}">${escapeHTML(result.app.anchor?.label || result.transcript?.ensembl_id)}</option>`;
      $("reference").disabled = false; $("gene-info").hidden = true; $("add-promoter").disabled=false;
      $("save-custom-reference").disabled=true; renderGeneOverview();
    }
  } else {
    state.gene=null; renderGeneOverview();
  }
  renderBasket();
}
function setSource(source) {
  state.source = source;
  $("gene-inputs").hidden = source !== "gene"; $("sequence-inputs").hidden = source !== "sequence";
  for (const b of document.querySelectorAll("[data-source]")) b.classList.toggle("active", b.dataset.source===source);
  $("gene-overview").hidden=source!=="gene" || !state.gene?.tss_groups;
  renderBasket();
}
function showResult(result) {
  state.result = result; state.selected = null; state.checked.clear(); state.page = 0;
  state.start = 0; state.end = result.promoter.length;
  $("empty-state").hidden = true; $("result-area").hidden = false; $("draft-notice").hidden = true; $("diff-notice").hidden = true;
  $("draft-badge").textContent = "SAVED";
  $("design-button").innerHTML = 'Design another run <span>→</span>';
  $("job-status").textContent = `Saved locally · ${result.app.run_id.slice(0,8)}`;
  $("result-title").textContent = result.gene ? `${result.gene.symbol} · saRNA candidates` : "Your promoter · saRNA candidates";
  $("run-meta").textContent = `RUN ${result.app.run_id.slice(0,8)} · ${new Date(result.app.created).toLocaleString()}`;
  $("result-reference").textContent = result.app.anchor ? `${result.app.anchor.label} · ${result.app.anchor.reference_kind} reference` : "Pasted promoter sequence · user-supplied reference";
  $("result-config").textContent = `${result.configuration.is_custom?"CUSTOM · ":""}${result.configuration.preset.name}`;
  $("metric-generated").textContent = fmt(result.counts.generated);
  $("metric-passing").textContent = fmt(result.counts.passing);
  const first = result.candidates.find(c => c.rank===1);
  $("metric-best").textContent = first ? `${first.rating.met}/${first.rating.evaluated}` : "—";
  $("metric-length").textContent = `${fmt(result.promoter.length)} bp`;
  renderRankingPerformance(result.ranking_performance);
  renderTopResult(result);
  const cell = result.cell_types.find(c => [c.id,c.name.toLowerCase()].includes(result.app.cell_line.toLowerCase()));
  const chromatin = cell?.chromatin;
  $("context-banner").innerHTML = `<strong>Cell context · ${escapeHTML(result.app.cell_line)}</strong><p>${chromatin?escapeHTML(chromatin.summary || "Local histone measurements loaded; inspect sources below."):"Promoter activity is unknown in this run. No cell-specific expression or histone measurement establishes that this locus is active."}</p><details><summary>Evidence gaps &amp; source notes</summary><ul>${result.app.notices.map(n => `<li>${escapeHTML(n)}</li>`).join("")}${result.app.anchor?.caveat?`<li>${escapeHTML(result.app.anchor.caveat)}</li>`:""}</ul></details>`;
  if (result.configuration.warnings.length) $("context-banner").innerHTML += `<p style="margin-top:8px"><b>Configuration notes:</b> ${result.configuration.warnings.map(escapeHTML).join(" · ")}</p>`;
  $("map-orientation").textContent = result.promoter.contig ? `chr${result.promoter.contig}:${fmt(result.promoter.start)}–${fmt(result.promoter.end)} (0-based, half-open) · DNA strand ${result.promoter.strand} · plotted 5′ → 3′ in transcript orientation` : "Plotted 5′ → 3′ as supplied · no genomic location assigned";
  $("candidate-detail").innerHTML = '<p class="hint">Select a candidate in the table or promoter map to inspect its sequences, strand prediction, and every rule verdict.</p>';
  renderMap(); renderTable();
  if (first) selectCandidate(first.id);
}
function renderDiff(previous,current) {
  if (!previous) return;
  const same = previous.promoter.sequence===current.promoter.sequence && JSON.stringify(previous.app.anchor)===JSON.stringify(current.app.anchor) && previous.promoter.tss_start===current.promoter.tss_start;
  $("diff-notice").hidden = false;
  if (!same) { $("diff-notice").textContent = "New input or reference: compare these runs as different promoter contexts, rather than attributing differences to one setting."; return; }
  const key = c => `${c.offset}:${c.target_site_dna}`;
  const old = new Map(previous.candidates.filter(c=>c.rank!==null).map(c=>[key(c),c.rank]));
  const now = new Map(current.candidates.filter(c=>c.rank!==null).map(c=>[key(c),c.rank]));
  const added = [...now.keys()].filter(k=>!old.has(k)).length, removed = [...old.keys()].filter(k=>!now.has(k)).length;
  const moved = [...now].filter(([k,r])=>old.has(k) && old.get(k)!==r).length;
  $("diff-notice").textContent = `Compared with run ${previous.app.run_id.slice(0,8)} on the same promoter: ${added} newly passing, ${removed} no longer passing, ${moved} shared candidates changed rank. Both records remain saved.`;
}
function positionText(c) { return c.tss_span ? `${signed(c.tss_span[0])}…${signed(c.tss_span[1])}` : `${c.offset+1}…${c.offset+c.target_site_dna.length}`; }
function ratingText(c) { return c.rating.evaluated ? `${c.rating.met}/${c.rating.evaluated}` : "Not evaluated"; }
function filteredCandidates() {
  const search = $("candidate-search").value.trim().toUpperCase();
  const values = state.result.candidates.filter(c => (state.filter==="all" || (state.filter==="passing" ? c.rank!==null : c.rank===null)) && (!search || c.id.toUpperCase().includes(search) || c.target_site_dna.includes(search)));
  values.sort((a,b) => {
    const value = c => state.sort==="position"?c.offset:state.sort==="gc"?c.gc_percent:state.sort==="rating"?c.rating.met:(c.rank ?? Number.MAX_SAFE_INTEGER);
    return (value(a)-value(b))*(state.descending?-1:1) || a.offset-b.offset;
  });
  return values;
}
function renderTable() {
  if (!state.result) return;
  const values = filteredCandidates(), pages = Math.max(1,Math.ceil(values.length/50));
  state.page = Math.min(state.page,pages-1);
  const page = values.slice(state.page*50,(state.page+1)*50);
  $("visible-count").textContent = `· ${fmt(values.length)}`;
  $("candidate-rows").innerHTML = page.length ? page.map(c => `<tr data-candidate="${c.id}" class="${state.selected===c.id?"selected":""}" tabindex="0" aria-label="Inspect ${c.id}"><td><input type="checkbox" data-check="${c.id}" ${state.checked.has(c.id)?"checked":""} aria-label="Select ${c.id} for export"></td><td>${c.rank || "—"}</td><td><span class="dna">${escapeHTML(c.target_site_dna)}</span><p class="candidate-id">${c.id}${c.rank===null?` · ${c.excluded_by.length} failed gates`:""}</p></td><td>${escapeHTML(positionText(c))}</td><td>${c.gc_percent}%</td><td><div class="rating-cell"><b>${ratingText(c)}</b><span class="rating-bar"><i style="width:${c.rating.percent ?? 0}%"></i></span></div>${c.rating.unavailable?`<p class="candidate-id">${c.rating.unavailable} unavailable</p>`:""}${c.rating.unchecked_gates.length?`<p class="candidate-id">${c.rating.unchecked_gates.length} unchecked gates</p>`:""}</td><td><span class="guide-badge ${escapeHTML(c.rating.guide || "unknown")}">${escapeHTML(c.rating.guide || "unknown")}</span>${!c.rating.models_agree?'<p class="disagreement">Models disagree</p>':""}</td></tr>`).join("") : '<tr><td colspan="7" style="padding:30px;text-align:center;color:var(--muted)">No candidates match this view. Inspect excluded targets or adjust your criteria and run again.</td></tr>';
  $("page-label").textContent = values.length ? `${fmt(state.page*50+1)}–${fmt(Math.min((state.page+1)*50,values.length))} of ${fmt(values.length)} · page ${state.page+1}/${pages}` : "0 candidates";
  $("page-previous").disabled = state.page===0; $("page-next").disabled = state.page>=pages-1;
  $("select-page").checked = page.length>0 && page.every(c=>state.checked.has(c.id));
  $("select-page").indeterminate = page.some(c=>state.checked.has(c.id)) && !$("select-page").checked;
  $("export-note").textContent = `${state.checked.size?`${state.checked.size} selected candidates` : "All passing candidates"}. Order CSV gives the 21-mers to synthesise, 19 nt RNA cores plus dTdT; core CSV gives the cores alone. No off-target screening has been performed.`;
}
async function selectCandidate(id) {
  const run = state.result, candidate = run.candidates.find(c=>c.id===id);
  if (!candidate) return;
  state.selected = id; renderTable(); renderMap();
  try {
    if (frozen && !candidate.rules) throw new Error("This archival record is missing the candidate’s full verdicts.");
    const full = candidate.rules ? candidate : await api(`/api/runs/${run.app.run_id}/candidates/${id}`);
    if (state.result!==run || state.selected!==id) return;
    Object.assign(candidate,full);
    const primary = full.guide_prediction.assessments.find(a=>a.table===full.rating.guide_model);
    const sequenceBox = (name,sequence) => `<div class="sequence-box"><div><span>${name} · 5′ → 3′</span><button type="button" data-copy="${escapeHTML(sequence)}">Copy</button></div><code>${escapeHTML(sequence)}</code></div>`;
    $("candidate-detail").innerHTML = `<div class="detail-heading"><div><h3>${id} <span class="muted">${full.rank?`Rank ${full.rank}`:"Excluded by gates"}</span></h3><p class="hint">${escapeHTML(positionText(full))} ${full.tss_span?"relative to reference TSS":"in supplied sequence"} · ${full.target_site_dna.length} nt core · ${full.gc_percent}% GC${full.genomic?` · genomic ${escapeHTML(full.genomic.strand)} strand`:""}</p></div><span class="tag">${ratingText(full)} PREFERENCES</span></div><div class="detail-body"><div>${sequenceBox("Promoter target DNA",full.target_site_dna)}${sequenceBox("Sense RNA core",full.sense_strand_rna)}${sequenceBox("Antisense RNA core",full.antisense_strand_rna)}<p class="hint">RNA cores shown without the dTdT overhangs required for this project’s duplex chemistry. This is a design record, with no off-target screening performed.</p></div><div><div class="guide-summary"><strong>Predicted RNA guide: ${escapeHTML(full.rating.guide || "unknown")}</strong><p>${escapeHTML(tableNames[full.rating.guide_model] || full.rating.guide_model)} · ${escapeHTML(full.guide_prediction.summary)}</p><p class="hint">This predicts duplex loading, not the identity of a promoter-associated target RNA, and not activation.</p></div><table class="model-table"><thead><tr><th>Parameter model</th><th>Guide</th><th>End difference</th></tr></thead><tbody>${full.guide_prediction.assessments.map(a=>`<tr data-tip="${escapeHTML(a.explanation)}" tabindex="0"><td>${escapeHTML(tableNames[a.table] || a.table)}</td><td>${escapeHTML(a.predicted_guide)}</td><td>${a.difference.toFixed(2)} ${escapeHTML(a.unit)}</td></tr>`).join("")}</tbody></table><details class="verdict"><summary>How the selected model reaches this prediction</summary><p>${escapeHTML(primary?.explanation || "The selected model was not assessed.")}</p></details></div></div>${full.excluded_by.length?`<p class="hint">${escapeHTML(full.rejection_summary)}</p>`:""}<div class="rule-verdicts"><h3>Every rule, explained</h3>${full.rules.map(r=>ruleVerdict(r,run.configuration)).join("")}</div>`;
  } catch (e) { error(e.message); }
}
function ruleVerdict(verdict,config) {
  const def = config.rules.find(r=>r.id===verdict.id);
  return `<details class="verdict"><summary><span class="outcome ${verdict.outcome}">${verdict.outcome==="not-applicable"?"N/A":verdict.outcome.toUpperCase()}</span><span>${escapeHTML(def?.name || verdict.id)}</span><span class="muted">${roleName[verdict.kind]}</span></summary><p>${escapeHTML(verdict.detail)}</p>${def?`<p class="source">${escapeHTML(def.provenance.display)} · ${escapeHTML(def.evidence_class)}</p><p><b>Why this rule:</b> ${escapeHTML(def.rationale)}</p><p><b>Changing it:</b> ${escapeHTML(def.consequence)}</p><p><b>Saved values:</b> ${def.parameters.map(p=>`${escapeHTML(p.name)} = ${escapeHTML(p.value)} (${escapeHTML(p.origin)})`).join("; ") || "No adjustable numeric parameters"}</p>`:""}</details>`;
}

function renderMap() {
  if (!state.result) return;
  const r = state.result, p = r.promoter, left = 133, right = 965, width = right-left;
  const scale = offset => left+(offset-state.start)/(state.end-state.start)*width;
  const inView = (start,end) => end>state.start && start<state.end;
  const tip = text => `data-tip="${escapeHTML(text)}" tabindex="0" aria-label="${escapeHTML(text)}"`;
  const parts = ['<defs><clipPath id="map-clip"><rect x="133" y="20" width="832" height="265"/></clipPath></defs>'];
  const label = (text,y,help) => parts.push(`<text class="track-label" x="12" y="${y}" ${tip(help)}>${escapeHTML(text)}</text>`);
  const rect = (start,end,y,height,color,help,attrs="") => {
    if (!inView(start,end)) return;
    const x1=scale(Math.max(start,state.start)), x2=scale(Math.min(end,state.end));
    parts.push(`<rect x="${x1}" y="${y}" width="${Math.max(2,x2-x1)}" height="${height}" rx="2" fill="${color}" ${tip(help)} ${attrs}/>`);
  };
  label("Reference / TSS",59,"An annotated transcript start or a chosen literature anchor. This is the reference for coordinates, not proof of transcription in your cell line.");
  label("Design windows",103,"Bands use the frozen settings of this run. Gates exclude failures; preferences contribute one equal point.");
  label("CpG islands",142,"DNA CpG-rich regions called using this run’s length, GC and observed/expected thresholds. Absence in a short window does not establish absence outside it.");
  label("Regulatory build",176,"Ensembl Regulatory Build pools observations across cell types. These annotations do not establish activity in the selected cell line.");
  label("TF motif matches",211,"Sequence predictions from the bundled motif matrices. A match is not measured factor occupancy. Hover for the matrix, score and strand.");
  label("Passing targets",245,"Targets that pass evaluable mandatory rules, coloured teal. Unavailable gates must still be examined in candidate details.");
  label("Excluded targets",273,"Targets with at least one failed mandatory rule. Hover or select to see which checks failed.");
  for (const y of [65,110,148,182,217,251,279]) parts.push(`<line x1="${left}" x2="${right}" y1="${y}" y2="${y}" stroke="#e8efeb"/>`);
  for (let i=0;i<=6;i++) {
    const offset=Math.round(state.start+(state.end-state.start-1)*i/6), x=scale(offset);
    // Coordinate labels only: Python already supplies scientific candidate spans.
    const continuous = p.tss_offset===null?null:offset-p.tss_offset;
    const text = continuous===null?String(offset+1):signed(continuous>=0?continuous+1:continuous);
    parts.push(`<text class="tick" text-anchor="${i===0?"start":i===6?"end":"middle"}" x="${x}" y="20">${text}</text>`);
  }
  parts.push('<g clip-path="url(#map-clip)">');
  if (p.tss_offset!==null && p.tss_offset>=state.start && p.tss_offset<state.end) {
    const x=scale(p.tss_offset);
    parts.push(`<line x1="${x}" x2="${x}" y1="29" y2="281" stroke="#52726b" stroke-dasharray="4 4" opacity=".5"/><path d="M${x-5} 37h10l-5 7z" fill="#344d55" ${tip("+1 · chosen reference. The first transcribed base in this coordinate system. Its source and uncertainty are recorded in run evidence.")}/><text x="${x+7}" y="46" fill="#344d55" font-size="10">+1 reference</text>`);
  }
  for (const t of r.tracks.alternative_tss) rect(t.start,t.end,53,12,"#8eaaa0",t.detail);
  if (p.tss_offset!==null) for (const id of ["tss-window","tss-preferred-window"]) {
    const rule=r.configuration.rules.find(q=>q.id===id); if (!rule) continue;
    const values=Object.fromEntries(rule.parameters.map(q=>[q.name,q.value]));
    const a=values[id==="tss-window"?"window-start-bp":"preferred-window-start-bp"];
    const b=values[id==="tss-window"?"window-end-bp":"preferred-window-end-bp"];
    rect(p.tss_offset+a,p.tss_offset+b+1,id==="tss-window"?84:99,11,id==="tss-window"?"#c2dfd2":"#75bfa7",`${roleName[rule.kind]}: ${rule.name}. ${signed(a)} to ${signed(b)}. Whole target containment is required. ${rule.provenance.display}`);
  }
  for (const t of r.tracks.cpg_islands) rect(t.start,t.end,134,14,"#aac7c9",t.detail);
  for (const t of r.tracks.regulatory) rect(t.start,t.end,165,14,t.kind.includes("promoter")?"#ddb370":"#efcf99",`${t.label}: ${t.detail}`);
  r.tracks.motifs.forEach((t,i)=>rect(t.start,t.end,197+(i%3)*5,4,"#a798d2",t.detail));
  for (const c of r.candidates.filter(c=>c.rank===null)) rect(c.offset,c.offset+c.target_site_dna.length,264+(c.offset%3)*3,3,"#dae2df",`${c.id}: excluded by ${c.excluded_by.join(", ")}. Position ${positionText(c)}. Click for details.`,`data-candidate="${c.id}"`);
  for (const c of r.candidates.filter(c=>c.rank!==null)) rect(c.offset,c.offset+c.target_site_dna.length,232+((c.rank-1)%3)*5,7,"#087e80",`${c.id} · rank ${c.rank} · ${ratingText(c)} evaluated preferences · ${c.rating.guide || "unknown"} guide (${c.rating.guide_model}). Position ${positionText(c)}. Click for details.`,`data-candidate="${c.id}"`);
  if (state.selected) {
    const c=r.candidates.find(c=>c.id===state.selected);
    rect(c.offset,c.offset+c.target_site_dna.length,c.rank===null?261:228,c.rank===null?19:25,"none",`${c.id} selected`,`stroke="#b97a2b" stroke-width="2" data-candidate="${c.id}"`);
  }
  parts.push('</g>'); $("promoter-map").innerHTML=parts.join("");
}
function viewport(factor,pan=0) {
  if (!state.result) return;
  const max=state.result.promoter.length, old=state.end-state.start, span=Math.min(max,Math.max(40,old*factor));
  const centre=(state.start+state.end)/2+old*pan;
  state.start=Math.max(0,Math.min(max-span,centre-span/2)); state.end=state.start+span;
  renderMap();
}
function renderEvidence() {
  const r=state.result, c=r.configuration;
  $("run-evidence-content").innerHTML = `<p><b>${escapeHTML(c.preset.name)}</b> · ${c.is_custom?"customized":"preset defaults"} · engine ${escapeHTML(r.app.engine_version)}</p><p>Run ${escapeHTML(r.app.run_id)} · ${escapeHTML(r.app.created)}</p><h3>Promoter source</h3><p>${escapeHTML(r.app.anchor?.label || "User-supplied DNA")}</p>${r.app.anchor?`<p>${escapeHTML(r.app.anchor.provenance?.display || r.app.anchor.provenance || "Source unavailable")}</p><pre style="white-space:pre-wrap;font-size:10px">${escapeHTML(JSON.stringify(r.app.anchor,null,2))}</pre>`:""}<h3>Custom departures</h3><ul>${c.customizations.map(v=>`<li>${escapeHTML(JSON.stringify(v))}</li>`).join("") || "<li>None</li>"}${c.warnings.map(v=>`<li>${escapeHTML(v)}</li>`).join("")}${c.preset.divergences.map(v=>`<li>${escapeHTML(v)}</li>`).join("")}</ul><h3>Source availability</h3><ul>${r.app.notices.map(n=>`<li>${escapeHTML(n)}</li>`).join("")}</ul>${r.regulatory_summary?`<p>${escapeHTML(r.regulatory_summary.provenance)}</p><p>${escapeHTML(r.regulatory_summary.summary)}</p>`:""}<h3>Rule effects</h3><p>Failures overlap: a candidate may fail several gates.</p><ul>${r.rejection_breakdown.map(v=>`<li>${escapeHTML(v.rule)}: ${fmt(v.rejected)} (${v.percent_of_generated}%)</li>`).join("")}</ul><h3>Frozen configuration</h3>${c.rules.map(def=>`<details class="verdict"><summary>${escapeHTML(def.name)} · ${roleName[def.kind]}</summary><p class="source">${escapeHTML(def.provenance.display)} · ${escapeHTML(def.evidence_class)}</p><p>${escapeHTML(def.rationale)}</p><p>${escapeHTML(def.consequence)}</p><p>${def.parameters.map(p=>`${escapeHTML(p.name)}: ${escapeHTML(p.value)} (${escapeHTML(p.origin)})`).join("; ")}</p></details>`).join("")}<h3>Disabled or omitted rules</h3><p>${state.catalogue.rules.filter(def=>!c.rules.some(v=>v.id===def.id)).map(def=>escapeHTML(def.name)).join(" · ")}</p>`;
  $("evidence-dialog").showModal();
}
function downloadBlob(blob,name) {
  const url=URL.createObjectURL(blob), a=document.createElement("a"); a.href=url; a.download=name; a.click(); setTimeout(()=>URL.revokeObjectURL(url),1000);
}
function csvText(result) {
  const ids=state.checked, candidates=result.candidates.filter(c=>ids.size?ids.has(c.id):c.rank!==null);
  const rows=[["run_id","candidate","rank","target_dna_5to3","sense_rna_core_5to3","antisense_rna_core_5to3","predicted_guide","guide_model","preferences_met","preferences_evaluated","core_length_nt"],...candidates.map(c=>[result.app.run_id,c.id,c.rank,c.target_site_dna,c.sense_strand_rna,c.antisense_strand_rna,c.rating.guide,c.rating.guide_model,c.rating.met,c.rating.evaluated,c.target_site_dna.length])];
  return rows.map(row=>row.map(v=>`"${String(v??"").replace(/"/g,'""')}"`).join(",")).join("\r\n");
}
function orderCsvText(result) {
  // Rows come from the Python record. Selection filters them by candidate id (column 2);
  // nothing here derives a sequence, an overhang or a verdict.
  const order=result.order, ids=state.checked;
  if (!order || !order.rows?.length) throw new Error("This record carries no order-ready duplexes.");
  const rows=ids.size?order.rows.filter(r=>ids.has(r[1])):order.rows;
  if (!rows.length) throw new Error("None of the selected candidates passed the gates, so none can be ordered.");
  return [order.columns,...rows].map(row=>row.map(v=>`"${String(v??"").replace(/"/g,'""')}"`).join(",")).join("\r\n");
}
async function exportOrder() {
  try {
    const r=state.result, name=`${r.gene?.symbol || "promoter"}-${r.app.run_id.slice(0,8)}-order.csv`;
    if (frozen) { downloadBlob(new Blob([orderCsvText(r)],{type:"text/csv"}),name); return; }
    const response=await fetch(`/api/runs/${r.app.run_id}/export.order.csv${state.checked.size?`?ids=${encodeURIComponent([...state.checked].join(","))}`:""}`);
    if (!response.ok) throw new Error(`Export failed (${response.status})`);
    downloadBlob(await response.blob(),name);
  } catch(e) { error(e.message); }
}
async function exportRun(type) {
  try {
    const r=state.result, name=`${r.gene?.symbol || "promoter"}-${r.app.run_id.slice(0,8)}.${type}`;
    if (frozen) { downloadBlob(new Blob([type==="csv"?csvText(r):JSON.stringify(r,null,2)],{type:type==="csv"?"text/csv":"application/json"}),name); return; }
    const response=await fetch(`/api/runs/${r.app.run_id}/export.${type}${type==="csv"&&state.checked.size?`?ids=${encodeURIComponent([...state.checked].join(","))}`:""}`);
    if (!response.ok) throw new Error(`Export failed (${response.status})`);
    downloadBlob(await response.blob(),name);
  } catch(e) { error(e.message); }
}
function renderRankingPerformance(measured) {
  // Shows the engine's own measured performance beside the rating it produces. Python computes
  // this from the bundled chapter 4 dataset; nothing here derives or softens the number.
  const box = $("ranking-performance");
  if (!measured) { box.hidden = true; return; }
  box.hidden = false;
  box.className = "performance-banner" + (measured.beats_chance ? "" : " at-chance");
  const dataset = measured.dataset;
  const parts = [
    `<strong>${measured.beats_chance ? "Ranking beat chance when measured" : "This ranking has not been shown to predict activity"}</strong>`,
    `<p>${escapeHTML(measured.headline)}</p>`,
    `<p class="hint">Measured against ${dataset.assayed} assayed sequences, ${dataset.activating} of which activated their gene (${escapeHTML(dataset.name)}). ${escapeHTML(measured.caveat)}</p>`,
  ];
  if (!measured.applies) parts.push(`<p class="hint warn">${escapeHTML(measured.custom_note)}</p>`);
  if (measured.circular) parts.push(`<p class="hint warn">${escapeHTML(measured.circularity_note)}</p>`);
  parts.push(`<p class="hint">Source: ${escapeHTML(measured.source)} · AUC ${measured.auc} against ${measured.chance_auc} for a coin flip, p = ${measured.auc_p_value}</p>`);
  box.innerHTML = parts.join("");
}

function referenceLabel(result) {
  // Which promoter reference this run designed against, named the way the setup panel names it.
  const anchor = result.app.anchor;
  if (anchor && anchor.id && !anchor.id.startsWith(result.transcript?.ensembl_id || "\u0000")) {
    return `${anchor.label || anchor.id}${anchor.reference_kind ? ` (${anchor.reference_kind} reference)` : ""}`;
  }
  const t = result.transcript;
  if (!t) return "pasted sequence";
  const marks = [t.is_mane_select ? "MANE Select" : null, t.biotype, t.refseq_match].filter(Boolean);
  return `${t.ensembl_id}${marks.length ? ` — ${marks.join(", ")}` : ""}`;
}

function orderRowFor(result, candidateId) {
  // The order-ready duplex Python computed for this candidate. Returns null for a candidate
  // that failed a gate: those are not in the order table, by design.
  const order = result.order;
  if (!order || !order.rows) return null;
  const row = order.rows.find(r => r[1] === candidateId);
  return row ? Object.fromEntries(order.columns.map((c, i) => [c, row[i]])) : null;
}

function renderTopResult(result) {
  // Opens on the answer a researcher came for: how many survived, what the best one is, and
  // what to order. Everything else on the page remains available below it.
  const box = $("top-result");
  const best = result.candidates.find(c => c.rank === 1);
  if (!best) {
    box.hidden = false;
    box.innerHTML = `<div class="card-heading"><div><p class="eyebrow">RESULT</p>`
      + `<h3>No candidate passed every evaluable gate</h3>`
      + `<p class="hint">${fmt(result.counts.generated)} targets were evaluated and all were excluded. `
      + `The promoter map and the Excluded view below show which rule rejected each one; the most `
      + `common reasons are in the run's rejection breakdown. This is a result, not a failure: it `
      + `usually means a gate is too strict for this promoter rather than that no saRNA exists.</p>`
      + `</div></div>`;
    return;
  }
  const order = orderRowFor(result, best.id);
  const met = best.rating.evaluated ? `${best.rating.met}/${best.rating.evaluated}` : "not evaluated";
  const duplex = order
    ? `<div class="top-duplex">`
      + `<div><div><span>Sense · 5′ → 3′</span><code>${escapeHTML(order.sense_order_5to3)}</code></div>`
      + `<button type="button" data-copy="${escapeHTML(order.sense_order_5to3)}">Copy</button></div>`
      + `<div><div><span>Antisense · 5′ → 3′</span><code>${escapeHTML(order.antisense_order_5to3)}</code></div>`
      + `<button type="button" data-copy="${escapeHTML(order.antisense_order_5to3)}">Copy</button></div>`
      + `<p class="hint">Two 21-mers: a 19 nt RNA core plus a 3′ dTdT overhang on each strand. `
      + `No off-target screening has been performed.</p></div>`
    : `<p class="hint">Order-ready sequences are not in this record.</p>`;
  box.hidden = false;
  box.innerHTML = `<div class="card-heading"><div><p class="eyebrow">START HERE</p>`
    + `<h3>${fmt(result.counts.passing)} of ${fmt(result.counts.generated)} targets passed every evaluable gate</h3>`
    + `<p class="hint">Highest-rated is <b>${escapeHTML(best.id)}</b> at ${escapeHTML(positionText(best))} `
    + `${best.tss_span ? "relative to the reference TSS" : "in the supplied sequence"}, meeting `
    + `${met} evaluated preferences, with ${escapeHTML(best.rating.guide || "an unknown")} predicted as the guide strand.</p>`
    + `</div><span class="tag">RANK 1</span></div>`
    + duplex
    + `<div class="top-actions-row">`
    + `<button type="button" id="print-duplex" class="primary">Print duplex sheet</button>`
    + `<button type="button" id="inspect-top" class="secondary">Why this one?</button>`
    + `<button type="button" id="export-order-top" class="secondary">Order CSV · all ${fmt(result.counts.passing)}</button>`
    + `</div>`;
}

function duplexSheetHTML(result, candidate) {
  // A one-page record to take to the bench. Everything comes from the run record; nothing is
  // recomputed here, so the sheet cannot disagree with the page that produced it.
  const order = orderRowFor(result, candidate.id);
  const measured = result.ranking_performance;
  const gene = result.gene ? `${result.gene.symbol} (${result.gene.ensembl_id})` : "Pasted sequence";
  const met = candidate.rating.evaluated
    ? `${candidate.rating.met} of ${candidate.rating.evaluated} evaluated preferences`
    : "no preferences could be evaluated";
  const rules = (candidate.rules || []);
  const named = r => escapeHTML(r.name || r.id);
  const gates = rules.filter(r => r.kind === "mandatory");
  const gatesPassed = gates.filter(r => r.outcome === "pass");
  const preferencesMet = rules.filter(r => r.kind === "preference" && r.outcome === "pass");
  const unavailable = rules.filter(r => r.outcome === "not-applicable");
  const row = (label, value) => `<tr><th>${label}</th><td>${value}</td></tr>`;
  return `<h1>saRNA duplex — ${escapeHTML(candidate.id)}</h1>
    <table class="sheet-table">
      ${row("Gene", escapeHTML(gene))}
      ${row("Promoter reference", escapeHTML(referenceLabel(result)))}
      ${row("Position", escapeHTML(positionText(candidate)) + (candidate.tss_span ? " relative to reference TSS" : " in supplied sequence"))}
      ${row("Genomic", candidate.genomic ? escapeHTML(`${candidate.genomic.contig}:${candidate.genomic.start}-${candidate.genomic.end} (${candidate.genomic.strand} strand, GRCh38)`) : "—")}
      ${row("Target site DNA", `<code>${escapeHTML(candidate.target_site_dna)}</code>`)}
    </table>
    <h2>Order these two oligonucleotides</h2>
    ${order ? `<table class="sheet-table">
      ${row("Sense strand", `<code>${escapeHTML(order.sense_order_5to3)}</code>`)}
      ${row("Antisense strand", `<code>${escapeHTML(order.antisense_order_5to3)}</code>`)}
      ${row("Chemistry", "19 nt RNA core + 3′ dTdT overhang on both strands; 21 nt each")}
    </table>` : "<p>This candidate did not pass the gates and has no order-ready duplex.</p>"}
    <h2>Why it was selected</h2>
    <p>Rank ${candidate.rank ?? "—"} of ${fmt(result.counts.passing)} that passed every evaluable gate,
       meeting ${escapeHTML(met)}. Predicted guide strand:
       ${escapeHTML(candidate.rating.guide || "unknown")}, by ${escapeHTML(candidate.rating.guide_model || "—")}.</p>
    <p><b>Gates passed (${gatesPassed.length} of ${gates.length}):</b> ${gatesPassed.length ? gatesPassed.map(named).join("; ") : "none"}. A gate excludes a candidate that fails it.</p>
    <p><b>Preferences met (${preferencesMet.length}):</b> ${preferencesMet.length ? preferencesMet.map(named).join("; ") : "none"}. A preference adds one point and never excludes.</p>
    ${unavailable.length ? `<p><b>Could not be evaluated (${unavailable.length}):</b> ${unavailable.map(named).join("; ")}. These counted neither for nor against it, and are not failures.</p>` : ""}
    <h2>What this does not establish</h2>
    <ul>
      <li>${measured ? escapeHTML(measured.headline) : "Ranking performance was not measured for this run."}</li>
      <li>No off-target screening has been performed. Screen these sequences independently before ordering.</li>
      <li>A rating counts criteria met. It is not a probability that this duplex activates the gene.</li>
      <li>Cell context: ${escapeHTML(result.app.cell_line || "not specified")}. A cell name alone supplies no expression evidence.</li>
    </ul>
    <p class="sheet-foot">Run ${escapeHTML(result.app.run_id)} · engine ${escapeHTML(result.app.engine_version)} ·
       configuration ${escapeHTML(result.configuration.preset.name)}${result.configuration.is_custom ? " (customised)" : ""} ·
       generated ${escapeHTML(new Date(result.app.created).toLocaleString())}</p>`;
}

async function printDuplexSheet() {
  const result = state.result;
  if (!result) return;
  const id = state.selected || result.candidates.find(c => c.rank === 1)?.id;
  if (!id) { error("There is no candidate to print."); return; }
  let candidate = result.candidates.find(c => c.id === id);
  try {
    if (!candidate.rules && !frozen) {
      candidate = await api(`/api/runs/${result.app.run_id}/candidates/${id}`);
    }
    $("duplex-sheet").innerHTML = duplexSheetHTML(result, candidate);
    window.print();
  } catch (e) { error(e.message); }
}


async function startScreening(scope) {
  // Opt-in, and the run is already readable: this adds evidence, it never gates the designs.
  const run = state.result;
  if (!run || frozen) { error("Screening needs a live run."); return; }
  for (const id of ["screen-transcript","screen-genomic"]) $(id).disabled = true;
  $("offtarget-status").textContent = "RUNNING";
  $("offtarget-progress").textContent = "Starting…";
  try {
    const started = await api(`/api/runs/${run.app.run_id}/screenings`,
      {method:"POST",headers:{"Content-Type":"application/json"},body:JSON.stringify({scope})});
    if (!started.id) throw new Error("The server did not return a screening id.");
    pollScreening(started.id, run.app.run_id);
  } catch (e) {
    error(e.message);
    $("offtarget-status").textContent = "FAILED";
    for (const id of ["screen-transcript","screen-genomic"]) $(id).disabled = false;
  }
}

async function pollScreening(screeningId, runId) {
  // The server waits NCBI's requested interval; this is only the interface catching up.
  try {
    const record = await api(`/api/screenings/${screeningId}`);
    if (state.result?.app.run_id !== runId) return;
    $("offtarget-progress").textContent = record.stage || "";
    if (record.status === "running") {
      setTimeout(() => pollScreening(screeningId, runId), 2000);
      return;
    }
    for (const id of ["screen-transcript","screen-genomic"]) $(id).disabled = false;
    if (record.status === "error") {
      $("offtarget-status").textContent = "FAILED";
      $("offtarget-progress").textContent = record.error || "The screening failed.";
      return;
    }
    $("offtarget-status").textContent = "COMPLETE";
    renderScreening(record);
  } catch (e) { error(e.message); }
}

function renderScreening(record) {
  const result = record.result, rows = result.specificity || [];
  state.screenings = state.screenings || {};
  state.screenings[record.scope] = result;
  state.screeningScope = record.scope;
  const scopeName = record.scope === "genomic"
    ? "Promoters (genomic)" : "Transcripts (RefSeq RNA)";
  const sort = state.screeningSort || "perfect_hits";
  const ordered = [...rows].sort((a, b) => {
    const value = r => r[sort] === null ? Number.POSITIVE_INFINITY : r[sort];
    return sort === "best_evalue" ? value(a) - value(b) : value(b) - value(a);
  });
  const header = (key, label, tip) =>
    `<th><button data-screensort="${key}"${sort === key ? ' class="active"' : ''} `
    + `data-tip="${escapeHTML(tip)}" tabindex="0">${label} ⇅</button></th>`;
  const genomic = record.scope === "genomic";
  $("offtarget-results").innerHTML =
    `<p class="hint"><b>${escapeHTML(scopeName)}</b> · searched `
    + `${escapeHTML(new Date(result.searched_at).toLocaleString())} · database `
    + `${escapeHTML(result.parameters.database)} · word size ${result.parameters.word_size} `
    + `· ${result.hit_count} alignments across ${rows.length} candidates. `
    + `Sequence similarity is not measured binding, and nothing here changes a rank.</p>`
    + (result.notes || []).map(n => `<p class="hint warn">${escapeHTML(n)}</p>`).join("")
    + `<div class="table-scroll"><table><thead><tr><th>Candidate</th>`
    + header("perfect_hits", "Perfect matches",
             "Full-length exact matches elsewhere. One of these is a specific unintended target.")
    + header("near_perfect_hits", "Near-perfect",
             "Alignments matching at least 17 of the 19 bases.")
    + header("longest_complementary_stretch", "Longest match",
             "The most bases matched in any single alignment.")
    + header("best_evalue", "Best E-value",
             "Lowest E-value among the alignments. Smaller means less likely by chance.")
    + header("hit_count", "Alignments", "How many alignments were returned at all.")
    + (genomic ? header("in_run_window_hits", "In your window",
         "Hits inside the same window this run designed against, near another gene's start site.")
       : "")
    + `</tr></thead><tbody>`
    + ordered.map(r => `<tr><td>${escapeHTML(r.candidate)}</td>`
        + `<td${r.perfect_hits ? ' class="flagged"' : ''}>${r.perfect_hits}</td>`
        + `<td>${r.near_perfect_hits}</td><td>${r.longest_complementary_stretch}</td>`
        + `<td>${r.best_evalue === null ? "—" : Number(r.best_evalue).toExponential(1)}</td>`
        + `<td>${r.hit_count}</td>`
        + (genomic ? `<td${r.in_run_window_hits ? ' class="flagged"' : ''}>`
                     + `${r.in_run_window_hits}</td>` : "")
        + `</tr>`).join("")
    + `</tbody></table></div>`;
}

function tooltip(target,event) {
  const text=target?.closest("[data-tip]")?.dataset.tip;
  if (!text) { $("tooltip").hidden=true; return; }
  const box=$("tooltip"); box.textContent=text; box.hidden=false;
  const rect=target.getBoundingClientRect(), x=event?.clientX ?? rect.left, y=event?.clientY ?? rect.bottom;
  box.style.left=`${Math.max(8,Math.min(x+14,innerWidth-box.offsetWidth-10))}px`;
  box.style.top=`${Math.max(8,Math.min(y+15,innerHeight-box.offsetHeight-10))}px`;
}

async function initialize() {
  try {
    state.catalogue=frozen?.catalogue || await api("/api/catalogue");
    $("engine-version").textContent=`v${state.catalogue.engine_version}`;
    $("preset").innerHTML=state.catalogue.presets.map(p=>`<option value="${p.id}">${escapeHTML(p.name)}</option>`).join("");
    $("cell-lines").innerHTML=state.catalogue.cells.map(c=>`<option value="${escapeHTML(c.name)}"></option>`).join("");
    if (!frozen) { setPreset("wang-2015"); syncDataAccess(); }
    if (frozen) {
      document.body.classList.add("frozen"); document.querySelector(".local-badge").textContent="Frozen offline record";
      document.querySelector(".intro .eyebrow").textContent=frozen.batch?"FROZEN PROMOTER COMPARISON":"FROZEN DESIGN RECORD";
      document.querySelector(".intro h1").innerHTML="A design you can inspect.<br><span>A record you can keep.</span>";
      document.querySelector(".intro > div > p:last-child").textContent="Explore saved references, sequences, settings and evidence. New designs require the local workbench.";
      if (frozen.batch) {
        for (const request of frozen.batch.request?.requests || []) if (request.custom_reference) state.customReferences[request.custom_reference.id]=clone(request.custom_reference);
        renderBatch(frozen.batch);
        const first=frozen.batch.runs.find(r=>r.status==="complete" && r.result);
        if (first) await openBatchRun(first.id,false);
      } else { state.gene=frozen.result.app.gene_overview || null; if (frozen.result.app.request.custom_reference) state.customReferences[frozen.result.app.request.custom_reference.id]=clone(frozen.result.app.request.custom_reference); renderGeneOverview(); showResult(frozen.result); }
    }
    else await refreshRuns();
  } catch(e) { error(e.message); }
}
$("design-form").addEventListener("submit",e=>{e.preventDefault();design();});
$("lookup-gene").addEventListener("click",()=>{error("");lookup().catch(e=>error(e,"lookup"));});
$("gene-symbol").addEventListener("input",()=>{state.gene=null;state.basket=[];state.customReferences={};resetCustomEditor();renderBasket();renderGeneOverview();$("reference").disabled=true;$("reference").innerHTML='<option value="">Look up this gene first</option>';$("gene-info").hidden=true;$("add-promoter").disabled=true;$("save-custom-reference").disabled=true;error("");});
$("reference").addEventListener("change",()=>{$("add-promoter").disabled=!$("reference").value;});
$("add-promoter").addEventListener("click",()=>addPromoter());
$("save-custom-reference").addEventListener("click",saveCustomReference);
$("cancel-custom-reference").addEventListener("click",resetCustomEditor);
$("design-batch").addEventListener("click",designBatch);
$("batch-retry-online").addEventListener("click",retryBatchOnline);
for (const type of ["json","html"]) $("export-batch-"+type).addEventListener("click",()=>exportBatch(type));
$("batch-rows").addEventListener("click",e=>{const button=e.target.closest("[data-open-batch-run]");if(button)openBatchRun(button.dataset.openBatchRun).catch(e=>error(e));});
$("promoter-basket").addEventListener("input",e=>{
  const field=e.target, row=field.closest("[data-basket-index]");
  if (row && field.dataset.basketSpan) { state.basket[Number(row.dataset.basketIndex)][field.dataset.basketSpan]=Number(field.value); markDraft(); }
});
$("promoter-basket").addEventListener("click",e=>{
  const button=e.target.closest("[data-basket-action]"), row=button?.closest("[data-basket-index]"); if (!row) return;
  const index=Number(row.dataset.basketIndex), entry=state.basket[index];
  if (button.dataset.basketAction==="remove") state.basket.splice(index,1);
  else if (button.dataset.basketAction==="capture") entry.pinned=clone(requestFromDraft());
  else if (button.dataset.basketAction==="share") entry.pinned=null;
  else if (button.dataset.basketAction==="edit") editCustomReference(entry.key);
  else { $("reference").value=entry.key; $("upstream").value=entry.upstream; $("downstream").value=entry.downstream; }
  renderBasket(); markDraft();
});
function chooseOverviewReference(key) {
  if (frozen || !key) return;
  $("reference").value=key; $("add-promoter").disabled=!key; markDraft();
  $("reference").focus();
}
$("tss-group-list").addEventListener("click",e=>{const button=e.target.closest("[data-choose-transcript]");if(button)chooseOverviewReference(`transcript:${button.dataset.chooseTranscript}`);});
function selectOverviewItem(target) {
  const node=target.closest("[data-reference],[data-tss-group]"); if (!node) return;
  if (node.dataset.reference) chooseOverviewReference(node.dataset.reference);
  else {
    const group=state.gene.tss_groups.find(g=>g.id===node.dataset.tssGroup);
    const details=[...$("tss-group-list").querySelectorAll("details")][state.gene.tss_groups.indexOf(group)];
    if (details) { details.open=true; details.querySelector("summary").focus(); }
  }
}
$("gene-overview-map").addEventListener("click",e=>selectOverviewItem(e.target));
$("gene-overview-map").addEventListener("keydown",e=>{if(e.key==="Enter" || e.key===" "){e.preventDefault();selectOverviewItem(e.target);}});
$("offline").addEventListener("change",()=>{syncDataAccess();error("");});
$("retry-online").addEventListener("click",()=>{
  if (state.busy || state.catalogue?.data_access?.offline_forced) return;
  const retry = state.retry;
  $("offline").checked = false; syncDataAccess(); markDraft(); error("");
  if (retry === "design") design();
  else lookup().catch(e=>error(e,"lookup"));
});
document.querySelectorAll("[data-source]").forEach(b=>b.addEventListener("click",()=>{setSource(b.dataset.source);markDraft();}));
$("design-form").addEventListener("input",markDraft);
$("design-form").addEventListener("change",markDraft);
$("preset").addEventListener("change",()=>{setPreset($("preset").value);markDraft();});
for (const [id,rule,param,divisor] of [["gc-min","gc-content","minimum-gc-fraction",100],["gc-max","gc-content","maximum-gc-fraction",100],["window-start","tss-window","window-start-bp",1],["window-end","tss-window","window-end-bp",1]]) {
  $(id).addEventListener("input",()=>{state.draft[rule].values[param]=Number($(id).value)/divisor;state.draft[rule].option="";syncQuick();renderEditor();});
}
$("open-settings").addEventListener("click",()=>$("settings-dialog").showModal());
document.querySelectorAll(".close-dialog").forEach(b=>b.addEventListener("click",()=>b.closest("dialog").close()));
$("reset-settings").addEventListener("click",()=>{setPreset($("preset").value);markDraft();});
$("rule-editor").addEventListener("change",e=>{
  const field=e.target, article=field.closest("[data-rule]"); if (!article) return;
  const id=article.dataset.rule, d=state.draft[id];
  if (field.hasAttribute("data-role")) d.role=field.value;
  else if (field.hasAttribute("data-option")) {
    d.option=field.value;
    const option=state.catalogue.rules.find(r=>r.id===id).options.find(o=>o.id===d.option);
    Object.assign(d.values,option?option.values:state.base[id].values);
  } else if (field.dataset.param) { d.values[field.dataset.param]=field.type==="number"?Number(field.value):field.value; d.option=""; }
  markDraft();syncQuick();renderEditor();
});
$("saved-runs").addEventListener("change",async()=>{
  if (!$("saved-runs").value || state.busy) return;
  try {
    const id=$("saved-runs").value;
    if (id.startsWith("batch:")) await openSavedBatch(id.slice(6));
    else { const run=await api(`/api/runs/${id}`); if (run.status!=="complete") throw new Error(run.error || run.stage); state.batch=null;$("batch-comparison").hidden=true;restoreDraft(run.result);showResult(run.result); }
    error("");
  } catch(e) { error(e); }
});
document.querySelectorAll("[data-filter]").forEach(b=>b.addEventListener("click",()=>{state.filter=b.dataset.filter;state.page=0;document.querySelectorAll("[data-filter]").forEach(v=>v.classList.toggle("active",v===b));renderTable();}));
$("candidate-search").addEventListener("input",()=>{state.page=0;renderTable();});
document.querySelectorAll("[data-sort]").forEach(b=>b.addEventListener("click",()=>{state.descending=state.sort===b.dataset.sort?!state.descending:b.dataset.sort==="rating";state.sort=b.dataset.sort;renderTable();}));
$("candidate-rows").addEventListener("click",e=>{if (!e.target.closest("input")) {const row=e.target.closest("[data-candidate]");if(row)selectCandidate(row.dataset.candidate);}});
$("candidate-rows").addEventListener("change",e=>{const id=e.target.dataset.check;if(id){e.target.checked?state.checked.add(id):state.checked.delete(id);renderTable();}});
$("candidate-rows").addEventListener("keydown",e=>{if ((e.key==="Enter" || e.key===" ") && e.target.matches("tr[data-candidate]")){e.preventDefault();selectCandidate(e.target.dataset.candidate);}});
$("promoter-map").addEventListener("click",e=>{const v=e.target.closest("[data-candidate]");if(v)selectCandidate(v.dataset.candidate);});
$("promoter-map").addEventListener("keydown",e=>{if ((e.key==="Enter" || e.key===" ") && e.target.dataset.candidate){e.preventDefault();selectCandidate(e.target.dataset.candidate);}});
$("page-previous").addEventListener("click",()=>{state.page--;renderTable();});$("page-next").addEventListener("click",()=>{state.page++;renderTable();});
$("select-page").addEventListener("change",e=>{filteredCandidates().slice(state.page*50,(state.page+1)*50).forEach(c=>e.target.checked?state.checked.add(c.id):state.checked.delete(c.id));renderTable();});
$("select-top").addEventListener("click",()=>{state.checked.clear();state.result.candidates.filter(c=>c.rank!==null).sort((a,b)=>a.rank-b.rank).slice(0,10).forEach(c=>state.checked.add(c.id));renderTable();});
$("zoom-in").addEventListener("click",()=>viewport(.5));$("zoom-out").addEventListener("click",()=>viewport(2));$("pan-left").addEventListener("click",()=>viewport(1,-.3));$("pan-right").addEventListener("click",()=>viewport(1,.3));$("zoom-reset").addEventListener("click",()=>{state.start=0;state.end=state.result.promoter.length;renderMap();});
$("candidate-detail").addEventListener("click",async e=>{const b=e.target.closest("[data-copy]");if(b){try{await navigator.clipboard.writeText(b.dataset.copy);b.textContent="Copied";}catch{error("Clipboard unavailable. Select and copy the displayed sequence.");}}});
$("run-evidence").addEventListener("click",renderEvidence);
for (const type of ["csv","json","html"]) $("export-"+type).addEventListener("click",()=>exportRun(type));
$("export-order").addEventListener("click",exportOrder);
$("screen-transcript").addEventListener("click",()=>startScreening("transcript"));
$("screen-genomic").addEventListener("click",()=>startScreening("genomic"));
$("offtarget-results").addEventListener("click",event=>{
  const button=event.target.closest("[data-screensort]");
  if (!button) return;
  state.screeningSort=button.dataset.screensort;
  const scope=state.screeningScope;
  if (scope) renderScreening({scope,result:state.screenings[scope],status:"complete"});
});
document.addEventListener("click", event => {
  const target = event.target.closest("#print-duplex,#inspect-top,#export-order-top");
  if (!target) return;
  if (target.id === "print-duplex") printDuplexSheet();
  if (target.id === "export-order-top") exportOrder();
  if (target.id === "inspect-top") {
    const best = state.result?.candidates.find(c => c.rank === 1);
    if (best) { selectCandidate(best.id); $("candidate-detail").scrollIntoView({behavior:"smooth", block:"start"}); }
  }
});
document.addEventListener("pointermove",e=>tooltip(e.target,e));document.addEventListener("focusin",e=>tooltip(e.target));document.addEventListener("focusout",()=>{$("tooltip").hidden=true;});document.addEventListener("scroll",()=>{$("tooltip").hidden=true;},true);
initialize();
