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
  customReferences:{}, basket:[], batch:null, editingReference:null, expression:null};
let sequencePending = [];
let sequenceFocus = null;
const roleName = {mandatory:"Gate", preference:"Preference", tiebreak:"Tiebreak", off:"Off"};
const tableNames = {"rna-xia-1998":"RNA · Xia 1998", "dna-santalucia-1998":"DNA · SantaLucia 1998", "legacy-dna-kj":"Legacy DNA · defective"};

// Keep activity local to the work it describes, with text for assistive technology.
function loading(id, active, message) {
  const node = $(id);
  if (message !== undefined) node.textContent = message;
  node.classList.toggle("is-loading", active);
  if (id !== "job-status") node.hidden = !active;
}

async function api(path, options) {
  const response = await fetch(path, options);
  let data;
  try { data = JSON.parse(await response.text()); }
  catch {
    throw Object.assign(new Error(`Local API returned HTTP ${response.status} `
      + `(${response.statusText || "unexpected response"}) without a readable JSON response.`),
      {status:response.status, code:"invalid-api-response"});
  }
  if (!response.ok) {
    const detail = Array.isArray(data.detail) ? data.detail.map(e => `${e.loc.join(".")}: ${e.msg}`).join("; ") : data.detail;
    throw Object.assign(new Error(detail || `Request failed (${response.status})`),
      {code:data.code, status:response.status});
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
    "Cache only, enforced by SARNA_OFFLINE. To fetch new references, unset it in the launch terminal and restart the app." : $("offline").checked ?
    "Cache only, with no Ensembl requests. A gene, transcript or promoter window that has not been cached will need online retrieval." :
    "Online retrieval allowed. Missing public references are fetched from Ensembl and cached. Design calculations run locally.";
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
  if (typeof workflowStatus === "function") workflowStatus();
}
function roleChoices(def) {
  if (def.id === "target-length-19") return ["mandatory"];
  if (def.default_kind === "tiebreak") return ["tiebreak","off"];
  return ["mandatory","preference","off"];
}
function provenanceText(p) {
  const text = escapeHTML(p?.display || "Source not recorded");
  const url = p?.url;
  return typeof url === "string" && /^https:\/\//.test(url)
    ? `<a href="${escapeHTML(url)}" target="_blank" rel="noopener noreferrer">${text}</a>` : text;
}
function evidenceLabel(def) {
  if (def.provenance?.is_warning) return "unsupported / no source";
  return ({established:"published support in source setting", empirical:"published observation or guidance", exploratory:"exploratory hypothesis"})[def.evidence_class] || def.evidence_class || "scope not recorded";
}
function ruleEvidence(def) {
  return `<p class="hint"><b>Evidence scope:</b> ${escapeHTML(def.evidence_note || "Evidence scope was not recorded in this saved run; a citation alone does not validate activation prediction.")}</p><p class="hint">${escapeHTML(def.application_note || "Gate, preference and off describe software behavior, not evidence strength.")}</p>`;
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
    return `<article class="rule-editor-item" data-rule="${def.id}"><div class="rule-editor-heading"><div><h3>${escapeHTML(def.name)}</h3><p class="source">${provenanceText(def.provenance)} <span>${escapeHTML(evidenceLabel(def))}</span></p><p class="preset-baseline">Role in selected preset: ${roleName[b.role]}</p></div><select data-role aria-label="Role for ${escapeHTML(def.name)}">${roleChoices(def).map(role => `<option value="${role}" ${role===d.role?"selected":""}>${roleName[role]}</option>`).join("")}</select></div><p>${escapeHTML(def.rationale)}</p>${ruleEvidence(def)}<p class="consequence"><b>Changing this:</b> ${escapeHTML(def.consequence)}</p>${def.options.length?`<label class="rule-option-label">Sourced parameter alternative</label><select data-option ${d.role==="off"?"disabled":""}><option value="">Preset / manual values</option>${def.options.map(o => `<option value="${o.id}" ${o.id===d.option?"selected":""}>${escapeHTML(o.label)}</option>`).join("")}</select><p class="option-note">${selected?`${provenanceText(selected.provenance)}. ${escapeHTML(selected.note)}`:"Manual values are recorded as researcher choices, without assigning a paper as their source."}</p>`:""}<div class="rule-parameters">${params}</div></article>`;
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
  // Python orders transcripts by reference_priority, so the MANE Select group and the MANE
  // transcript inside it already arrive first; the list is rendered in the order given rather
  // than re-sorted here. The marker says *why* it is first, because "on top" alone is a
  // position a reader has to guess the meaning of.
  const mark = t => t?.mane ? " · MANE Select · default"
    : t?.canonical ? " · Ensembl canonical" : "";
  const choices = groups.length ? groups.map(g=>`<optgroup label="Annotated +1 at ${escapeHTML(g.contig)}:${fmt(g.tss+1)} · ${g.transcript_ids.length} transcript(s)${g.mane?" · default reference":""}">${g.transcript_ids.map(id=>{
    const t = transcripts.find(t=>t.id===id);
    return `<option value="transcript:${escapeHTML(id)}">${escapeHTML(id)}${mark(t)} · ${escapeHTML(t?.biotype)}</option>`;
  }).join("")}</optgroup>`).join("") : `<optgroup label="Annotated transcripts">${transcripts.map(t=>`<option value="transcript:${escapeHTML(t.id)}">${escapeHTML(t.id)}${mark(t)} · ${escapeHTML(t.biotype)}</option>`).join("")}</optgroup>`;
  const anchorChoices = anchors.map(a=>`<option value="anchor:${escapeHTML(a.id)}" title="${escapeHTML(a.label)} · ${escapeHTML(a.reference_kind)}">${escapeHTML(a.id)} · ${escapeHTML(a.reference_kind)}</option>`).join("");
  const customChoices = Object.values(state.customReferences).map(a=>`<option value="custom:${escapeHTML(a.id)}">${escapeHTML(a.name)} · researcher supplied</option>`).join("");
  $("reference").innerHTML = `<option value="">Choose a promoter reference…</option>${choices}${anchorChoices?`<optgroup label="Declared literature anchors">${anchorChoices}</optgroup>`:""}${customChoices?`<optgroup label="Researcher references">${customChoices}</optgroup>`:""}`;
  $("reference").disabled = false;
  $("reference").value = [...$("reference").options].some(o=>o.value===selected) && selected ? selected : gene.default_transcript ? `transcript:${gene.default_transcript}` : "";
  $("add-promoter").disabled = !$("reference").value;
  $("save-custom-reference").disabled = false;
  $("custom-reference-locus").textContent = `GRCh38 · chr${gene.contig || "?"} · gene strand ${gene.strand}. Enter a 1-based genomic coordinate, not a TSS-relative position.`;
  renderGeneOverview(); renderBasket(); renderExpression();
}
function promoterSaRNACoverage(gene, custom, upstream, downstream) {
  if (!gene) return [];
  const known=(gene.sarna_references || []).filter(r=>r.known_location);
  const references=[...(gene.transcripts || []).map(t=>({key:`transcript:${t.id}`,label:t.id,tss:t.tss})),
    ...(gene.anchors || []).map(a=>({key:`anchor:${a.id}`,label:a.id,tss:a.reference_point})),
    ...Object.values(custom).map(a=>({key:`custom:${a.id}`,label:a.name,tss:a.tss_position-1}))];
  return references.map(reference=>{
    const start=reference.tss===null?null:gene.strand==='-'?reference.tss-downstream+1:reference.tss-upstream;
    const end=reference.tss===null?null:gene.strand==='-'?reference.tss+upstream+1:reference.tss+downstream;
    const contains=known.filter(r=>{
      const site=r.known_location;
      return start!==null && end!==null && site.assembly===gene.assembly && site.contig===gene.contig
        && site.strand===gene.strand && site.start>=start && site.end<=end;
    });
    return {...reference,start,end,contains};
  });
}
function renderKnownSaRNAPromoters() {
  const gene=state.gene, host=$("known-sarna-promoters"), note=$("known-sarna-selection");
  const known=gene?.sarna_references || [];
  if (!gene || !known.length) {host.innerHTML='';note.textContent='';return;}
  const coverage=promoterSaRNACoverage(gene,state.customReferences,Number($("upstream").value),Number($("downstream").value));
  for (const option of $("reference").options) {
    option.dataset.baseLabel ||= option.textContent;
    const row=coverage.find(r=>r.key===option.value);
    option.textContent=option.dataset.baseLabel+(row?.contains.length?` · known saRNAs: ${row.contains.map(r=>r.name).join(', ')}`:'');
  }
  const selected=coverage.find(r=>r.key===$("reference").value);
  note.textContent=selected?.contains.length?`This window contains known sites: ${selected.contains.map(r=>r.name).join(', ')}.`:'No curated saRNA sites in the selected window. See “Known saRNAs by promoter region” in the gene overview.';
  const rows=coverage.filter(r=>r.contains.length);
  host.innerHTML=`<section class="known-sarna-promoters"><h3>Known saRNAs by promoter region</h3><p class="hint">Shown before design for the current ${fmt(Number($("upstream").value))} bp upstream / ${fmt(Number($("downstream").value))} bp downstream window. Coverage uses recorded genomic sites; exact sequence matching follows when the promoter is loaded.</p><p>${known.map(r=>`<strong>${escapeHTML(r.name)}</strong>: chr${escapeHTML(r.known_location.contig)}:${fmt(r.known_location.start+1)}–${fmt(r.known_location.end)} (1-based)`).join(' · ')}</p><div class="table-scroll"><table><thead><tr><th>Promoter reference</th><th>Known saRNAs in window</th><th></th></tr></thead><tbody>${rows.map(r=>`<tr class="${r.key===$("reference").value?'selected':''}"><td>${escapeHTML(r.label)}</td><td>${r.contains.map(s=>escapeHTML(s.name)).join(', ')}</td><td>${!frozen?`<button type="button" data-known-promoter="${escapeHTML(r.key)}">${r.key===$("reference").value?'Selected':'Choose this promoter'}</button>`:''}</td></tr>`).join('') || '<tr><td colspan="3">No available reference window fully contains these known sites. Increase the window or choose a reference near the recorded sites.</td></tr>'}</tbody></table></div><p class="hint">Curated references are not a complete catalogue of assayed saRNAs. Their presence does not establish promoter activity in your cell.</p></section>`;
}
function renderGeneOverview() {
  renderKnownSaRNAPromoters();
  const gene = state.gene;
  const groups = gene?.tss_groups || [];
  if (!gene || !Number.isFinite(gene.start) || !Number.isFinite(gene.end)) { $("gene-overview").hidden = true; return; }
  $("gene-overview").hidden = state.source !== "gene" && !frozen;
  $("gene-overview-meta").textContent = frozen ? "FROZEN REFERENCE CHOICES" : "DRAFT REFERENCE CHOICES";
  $("gene-overview-description").textContent = `${gene.symbol} · ${groups.length} annotated start coordinate(s) across ${(gene.transcripts || []).length} transcripts. Click a start-site group to inspect its transcripts below.`;
  const anchors = gene.anchors || [], custom = Object.values(state.customReferences);
  const knownSaRNAs = (gene.sarna_references || state.result?.sequence_comparisons || []).filter(r=>r.reference && r.known_location && r.known_location.contig===gene.contig);
  const lo = Math.min(gene.start,...anchors.map(a=>a.start),...custom.map(a=>a.tss_position-1),...knownSaRNAs.map(r=>r.known_location.start));
  const hi = Math.max(gene.end,...anchors.map(a=>a.end),...custom.map(a=>a.tss_position),...knownSaRNAs.map(r=>r.known_location.end));
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
  parts.push('<text x="12" y="180" fill="#61767a" font-size="10">Reference saRNAs</text>');
  knownSaRNAs.forEach((record,index)=>{
    const location=record.known_location, x=scale(location.start), y=166+(index%2)*16, style=sequenceAppearance(record);
    const help=`${record.name} · chr${location.contig}:${location.start+1}–${location.end} (1-based). ${location.source}`;
    parts.push(`<circle cx="${x}" cy="${y}" r="4" fill="${style.color}" ${tip(help)} data-sarna-region="${escapeHTML(record.name)}"/><text x="${x+7}" y="${y+3}" font-size="10" fill="${style.color}" ${tip(help)} data-sarna-region="${escapeHTML(record.name)}">${escapeHTML(record.name)} · ${style.label}</text>`);
  });
  $("gene-overview-map").setAttribute('viewBox','0 0 1000 220');
  parts.push(`<text x="${right}" y="213" text-anchor="end" fill="#61767a" font-size="10">chr${escapeHTML(gene.contig)} · 1-based genomic coordinates · gene strand ${escapeHTML(gene.strand)} · 5′ → 3′</text>`);
  $("gene-overview-map").innerHTML=parts.join("");
  $("tss-group-list").innerHTML=groups.map(g=>`<details><summary>Annotated +1 · ${escapeHTML(g.contig)}:${fmt(g.tss+1)} <span class="muted">${g.transcript_ids.length} transcript(s)${g.mane?" · MANE Select":""}</span></summary><p class="hint">Shared start coordinate, not a measurement of promoter activity.</p><div>${g.transcript_ids.map(id=>`<button type="button" class="text-button" data-choose-transcript="${escapeHTML(id)}">${escapeHTML(id)}</button>`).join("")}</div></details>`).join("");
}
/* What the cell-line box actually does, shown next to the box. The control looked inert
   because it is almost inert: of the 24 (cell line x evidence axis) pairs in the registry,
   exactly one is available, AF22 chromatin, and that one needs the option switched on and the
   tracks installed. The choice is recorded with the run and never touches a gate or a rank.
   Saying so beside the field is better than leaving a researcher to infer it. */
/* The status words come from axisState, declared with the results panel below and shared on
   purpose: the setup side and the results side must lay these facts out differently, because a
   four-column table is unreadable in a 310px sidebar, but they must not be able to word them
   differently. Writing a second copy of the words here is the bug this comment prevents, since
   a later function declaration of the same name silently wins. */
function axisBrief(axis) {
  const {state, word} = axisState(axis);
  return `<li class="axis-${state}"><span class="axis-state ${state}">${word}</span> `
    + `<strong>${escapeHTML(AXIS_LABEL[axis.axis] || axis.axis)}</strong>`
    + `<span class="muted"> ${escapeHTML(axis.question || "")}</span>`
    + `<details><summary>${axis.source ? escapeHTML(axis.source) : "no source"}</summary>`
    + `<p class="hint">${escapeHTML(axis.note || "")}</p></details></li>`;
}
function cellTypeFor(typed) {
  const wanted = typed.trim().toLowerCase();
  if (!wanted) return null;
  return (state.catalogue?.cells || []).find(
    c => c.name.toLowerCase() === wanted || c.id.toLowerCase() === wanted) || null;
}
function renderCellEvidence() {
  const box = $("cell-evidence");
  if (!box) return;
  const typed = $("cell-line").value;
  const always = `<p class="hint">Recorded with the run as your experimental context. It does `
    + `not change which candidates pass, how they rank, or any sequence calculation.</p>`;
  if (!typed.trim()) {
    box.innerHTML = always + `<p class="hint">Leave it blank, pick a known line, or type your `
      + `own. Known lines carry a record of what has and has not been measured in them.</p>`;
    return;
  }
  const cell = cellTypeFor(typed);
  if (!cell) {
    box.innerHTML = always + `<p class="hint"><strong>${escapeHTML(typed.trim())}</strong> is `
      + `not in the registry, so nothing shown anywhere in this app is specific to it. No `
      + `measurement is implied for it. Add your own with the SARNA_CELL_LINES file.</p>`;
    return;
  }
  const axes = cell.data || [];
  const live = axes.filter(a => a.available);
  const chromatinLive = live.some(a => a.axis === "chromatin");
  const summary = live.length
    ? `${live.length} of ${axes.length} kinds of evidence are available for this line.`
    : `None of the ${axes.length} kinds of evidence are available for this line.`;
  box.innerHTML = always
    + `<p class="cell-name"><strong>${escapeHTML(cell.name)}</strong>`
    + (cell.rrid ? ` <span class="muted">RRID ${escapeHTML(cell.rrid)}</span>` : "")
    + (cell.experimental_context ? ` <span class="muted">used in this project's own experiments</span>` : "")
    + `</p><p class="hint">${escapeHTML(summary)}</p>`
    + `<ul class="axis-brief">${axes.map(axisBrief).join("")}</ul>`
    + (chromatinLive
        ? `<p class="hint">To use it, switch on <strong>AF22 chromatin</strong> under Data `
          + `sources and run options, and have the track files installed. Without both, no `
          + `histone measurement is loaded and the result says so rather than reporting closed `
          + `chromatin.</p>`
        : "")
    + `<details class="cell-detail"><summary>About this line</summary>`
    + `<p class="hint">${escapeHTML(cell.description)}</p>`
    + `<p class="hint">${escapeHTML(cell.provenance)}</p></details>`;
}
/* Expression evidence sits beside the reference choice and never reorders it. The picker's
   order is the annotation's (MANE Select first); a transcript being the most expressed in bulk
   tissue is a different claim, and promoting it would decide for the researcher. */
function expressionFor(id) {
  return state.expression?.transcripts?.[id] || null;
}
function annotateReferenceOptions() {
  // The option text has to fit a 310px control. Three facts compete for it: the identifier,
  // which cannot be shortened; why this one is offered first; and what GTEx measured. The
  // biotype is abbreviated and the expression reduced to one decimal so all three survive,
  // and the untruncated version goes in the option's title for anyone who wants it.
  const payload = state.expression;
  for (const option of $("reference").options) {
    if (option.dataset.baseLabel === undefined) option.dataset.baseLabel = option.textContent;
    const id = option.value.startsWith("transcript:") ? option.value.slice(11) : "";
    if (!id) { option.textContent = option.dataset.baseLabel; continue; }
    const transcript = state.gene?.transcripts?.find(t => t.id === id);
    const entry = payload?.transcripts ? payload.transcripts[id] : null;
    const mark = transcript?.mane ? " · MANE" : transcript?.canonical ? " · canonical" : "";
    const brief = compactExpression(entry);
    // protein_coding is what a reader assumes, and spelling it out is what pushes this label
    // past the width of the control. Anything else is worth the characters: the highest brain
    // value for SCN1A sits on a nonsense_mediated_decay transcript.
    const biotype = transcript?.biotype === "protein_coding"
      ? "" : ` · ${shortBiotype(transcript?.biotype)}`;
    option.textContent = id + mark + biotype + (brief ? ` · ${brief}` : "");
    option.title = option.dataset.baseLabel + (entry ? ` · ${entry.label}` : "");
  }
}
function compactExpression(entry) {
  // One decimal would print a measured 0.01 as "0.0", which in a dropdown is indistinguishable
  // from a measured zero. They are different findings, so a value too small to show at this
  // precision says so instead of rounding into the other one.
  if (!entry) return "";
  if (entry.status !== "measured") return "not in GTEx";
  const brain = entry.brain;
  if (!brain) return "no brain group";
  if (brain.median === 0) return "brain 0";
  return brain.median < 0.1 ? "brain <0.1" : `brain ${brain.median.toFixed(1)}`;
}
// Biotypes are printed in a 310px picker option alongside an identifier and a number, where
// `protein_coding_CDS_not_defined` alone overflows the control. Shortened for the option only;
// the full term stays in the comparison table and in the option's hover text.
const BIOTYPE_SHORT = {
  protein_coding: "coding",
  nonsense_mediated_decay: "NMD",
  non_stop_decay: "non-stop",
  retained_intron: "intron",
  processed_transcript: "processed",
  protein_coding_CDS_not_defined: "no CDS",
  processed_pseudogene: "pseudogene",
  unprocessed_pseudogene: "pseudogene",
  transcribed_unprocessed_pseudogene: "pseudogene",
};
const shortBiotype = (biotype) => BIOTYPE_SHORT[biotype] || (biotype || "").replace(/_/g, " ");

function renderExpression() {
  closeOrganChart();
  // Two surfaces, because one 310px column cannot hold both an answer and a comparison. The
  // sidebar answers "what about the transcript I have selected" in three lines beside the
  // control. The comparison of every transcript goes in the gene overview card, which is the
  // full width of the workspace and is where a reader is already looking at the locus.
  renderExpressionSidebar();
  renderExpressionOverview();
  annotateReferenceOptions();
}

function selectedTranscriptId() {
  const value = $("reference").value;
  return value.startsWith("transcript:") ? value.slice(11) : "";
}

function renderExpressionSidebar() {
  const box = $("expression"), payload = state.expression;
  if (state.source !== "gene" || !state.gene || !payload) { box.hidden = true; return; }
  box.hidden = false;
  if (payload.loading) {
    box.innerHTML = '<p class="hint is-loading" role="status">Loading transcript expression evidence…</p>';
    return;
  }
  if (payload.available === false) {
    box.innerHTML = `<p class="eyebrow">TRANSCRIPT EXPRESSION</p>`
      + `<p class="hint">${escapeHTML(payload.note)}</p>`
      + `<p class="hint muted">${escapeHTML(payload.detail)}</p>`
      + `<button type="button" class="text-button" id="retry-expression">Try again</button>`;
    $("retry-expression").addEventListener("click", () => loadExpression(state.gene.symbol));
    return;
  }
  const selected = selectedTranscriptId();
  const entry = selected ? expressionFor(selected) : null;
  const heading = `<p class="eyebrow">TRANSCRIPT EXPRESSION · GTEx `
    + `${escapeHTML(payload.dataset_id.replace("gtex_", ""))}</p>`;
  if (!entry) {
    box.innerHTML = heading
      + `<p class="hint">Select an annotated transcript to see what GTEx measured for it.</p>`;
    return;
  }
  const transcript = state.gene.transcripts?.find(t => t.id === selected);
  const mark = transcript?.mane ? " · MANE Select"
    : transcript?.canonical ? " · Ensembl canonical" : "";
  box.innerHTML = heading
    + `<p class="hint expression-subject">${escapeHTML(selected + mark)}</p>`
    + `<p class="expression-headline">${expressionHeadline(entry)}</p>`
    + `<p class="hint">${escapeHTML(expressionSubline(entry))}</p>`
    + (organStrip(entry, 6) ? `<p class="organ-line">${organStrip(entry, 6)}</p>` : "")
    + `<p class="hint"><a href="#gene-overview" class="text-link">`
    + `Compare all ${fmt(payload.measured_transcripts)} transcripts below</a></p>`;
}

let emojiRenders = null;
function emojiAvailable() {
  // WSL and bare Linux containers often ship no colour emoji font, and the browser then draws
  // an empty box. A box beside a measurement is worse than no symbol, so detect it once and
  // fall back to names. Measured against a noncharacter, which can never have a glyph: if the
  // two come out the same width, neither is being rendered.
  if (emojiRenders !== null) return emojiRenders;
  try {
    const context = document.createElement("canvas").getContext("2d");
    if (!context) { emojiRenders = true; return emojiRenders; }
    context.font = "16px sans-serif";
    emojiRenders = context.measureText("\u{1F9E0}").width
      !== context.measureText("\u{FFFF}").width;
  } catch { emojiRenders = true; }
  return emojiRenders;
}

function organStrip(entry, limit) {
  // Subtle by construction: symbols only, the number on hover, and nothing at all when the
  // transcript was measured at zero everywhere. An organ here means a median above zero in
  // one of its sample groups, which is why the title always carries the group and the value.
  const organs = (entry.organs || []).slice(0, limit);
  if (!organs.length) return "";
  const chips = organs.map(o => {
    const value = `${o.median.toFixed(2)} ${o.unit} in ${o.tissue_name}`;
    const donors = o.sample_count ? `, ${fmt(o.sample_count)} donors` : "";
    const caveat = o.is_body_site ? "" : ". Cultured cells, not a body site";
    const title = escapeHTML(`${o.name}: ${value}${donors}${caveat}`);
    // An organ with no honest emoji is named instead. A lookalike object would be both
    // ambiguous and undignified in a tool that reports measurements.
    return o.symbol && emojiAvailable()
      ? `<span class="organ" title="${title}" aria-label="${title}">${o.symbol}</span>`
      : `<span class="organ organ-named" title="${title}">`
        + `${escapeHTML(o.chip || o.name)}</span>`;
  }).join("");
  const hidden = (entry.other_organs || 0)
    + Math.max(0, (entry.organs || []).length - organs.length);
  const trigger = `data-expression-chart="${escapeHTML(entry.transcript_id)}" aria-haspopup="dialog" aria-controls="organ-popover"`;
  const more = hidden
    ? `<button type="button" class="organ-more" ${trigger} aria-label="Show all expression values for ${escapeHTML(entry.transcript_id)}">+${hidden} more</button>` : "";
  return `<span class="organ-strip">${chips}${more}<button type="button" class="organ-chart-button" ${trigger} aria-label="Show expression chart for ${escapeHTML(entry.transcript_id)}">View chart</button></span>`;
}

let organTrigger = null, organCloseTimer = null;
function closeOrganChart() {
  clearTimeout(organCloseTimer);
  $("organ-popover").hidden = true;
  organTrigger?.setAttribute("aria-expanded", "false");
  organTrigger = null;
}
function expressionChart(entry, payload) {
  // Old frozen records still carry the full tissue table, even without the organ summary.
  const grouped = Array.isArray(entry.all_organs);
  const rows = grouped ? entry.all_organs : (entry.all || []).filter(row=>row.median>0);
  const maximum = Math.max(...rows.map(row=>row.median), 0);
  const unit = rows[0]?.unit || entry.highest?.unit || "TPM";
  const valueText = value => value>0 && value<0.01 ? value.toPrecision(2) : value.toFixed(2);
  const bars = rows.map(row=>{
    const name = grouped ? row.name : row.tissue_name;
    const context = grouped ? row.tissue_name : "";
    return `<div class="expression-bar-row"><div class="expression-bar-label"><strong>${escapeHTML(name)}</strong><small>${escapeHTML(context)}${row.is_body_site===false?" · cultured cells":""}${row.sample_count?` · n = ${fmt(row.sample_count)} donors`:""}</small></div><div class="expression-bar-track"><span style="width:${maximum?100*row.median/maximum:0}%"></span></div><span class="expression-bar-value">${valueText(row.median)}</span></div>`;
  }).join("");
  return `<div class="expression-chart-heading"><div><p class="eyebrow">TRANSCRIPT EXPRESSION</p><h3>Expression by ${grouped?"organ":"sample group"}</h3><p class="hint">${escapeHTML(entry.transcript_id)}</p></div><button type="button" data-close-organ-chart aria-label="Close expression chart">×</button></div><p class="expression-chart-scale">Median expression · ${escapeHTML(unit)} · linear scale, 0 to ${valueText(maximum)}</p><div class="expression-bars">${bars || '<p class="hint">No sample group was measured above zero.</p>'}</div><p class="hint">${grouped?"Each bar shows the highest sample-group median within that organ, not an average across tissues.":"Each bar shows one measured sample-group median."} Only values above zero are shown. TPM means transcripts per million.</p><div class="expression-chart-source"><a href="https://gtexportal.org/home/datasets" target="_blank" rel="noopener noreferrer">Source: ${escapeHTML(payload.source || "GTEx Portal")}</a><p>${escapeHTML(payload.dataset_id)} · GENCODE ${escapeHTML(payload.annotation_release)} · transcript ${escapeHTML(entry.release_id || entry.transcript_id)}</p><p>${escapeHTML(payload.provenance)}</p><p>GTEx tissue and donor-derived cultured-cell measurements do not establish expression in AF22, HeLa, or your experimental cell line.</p></div>`;
}
function openOrganChart(trigger) {
  clearTimeout(organCloseTimer);
  const payload = state.expression;
  const entry = payload?.transcripts?.[trigger.dataset.expressionChart];
  if (!entry) return;
  const box = $("organ-popover");
  if (organTrigger!==trigger || box.hidden) {
    closeOrganChart(); organTrigger=trigger;
    box.innerHTML=expressionChart(entry,payload);
  }
  trigger.setAttribute("aria-expanded","true");
  box.hidden=false;
  const rect=trigger.getBoundingClientRect();
  box.style.left=`${Math.max(8,Math.min(rect.right-box.offsetWidth,innerWidth-box.offsetWidth-8))}px`;
  box.style.top=`${Math.max(8,Math.min(rect.bottom+8,innerHeight-box.offsetHeight-8))}px`;
}
function scheduleOrganClose() {
  clearTimeout(organCloseTimer);
  organCloseTimer=setTimeout(closeOrganChart,220);
}
document.addEventListener("pointerover",event=>{
  const trigger=event.target.closest("[data-expression-chart]");
  if (trigger) openOrganChart(trigger);
  else if (event.target.closest("#organ-popover")) clearTimeout(organCloseTimer);
});
document.addEventListener("pointerout",event=>{
  if (event.target.closest("[data-expression-chart], #organ-popover")) scheduleOrganClose();
});
document.addEventListener("focusin",event=>{
  const trigger=event.target.closest("[data-expression-chart]");
  if (trigger) openOrganChart(trigger);
  else if (!event.target.closest("#organ-popover")) closeOrganChart();
});
document.addEventListener("click",event=>{
  const trigger=event.target.closest("[data-expression-chart]");
  if (trigger) openOrganChart(trigger);
  else if (event.target.closest("[data-close-organ-chart]") || !event.target.closest("#organ-popover")) closeOrganChart();
});
document.addEventListener("keydown",event=>{if(event.key==="Escape")closeOrganChart();});
window.addEventListener("resize",closeOrganChart);
document.addEventListener("scroll",event=>{if(!event.target.closest?.("#organ-popover"))closeOrganChart();},true);

function expressionHeadline(entry) {
  // One line, and it has to answer the question the researcher came with: is this transcript
  // made in brain. A measured zero is stated as measured, never as missing.
  if (entry.status !== "measured") return `<strong>Not in GTEx ${escapeHTML("v10")}</strong>`;
  const brain = entry.brain;
  if (!brain) return `<strong>No nervous-system group reported</strong>`;
  return brain.median > 0
    ? `<strong>${brain.median.toFixed(2)} ${escapeHTML(brain.unit)}</strong> in brain`
    : `<strong>0.00 ${escapeHTML(brain.unit)}</strong> in every brain group`;
}

function expressionSubline(entry) {
  if (entry.status !== "measured") {
    return "This transcript is absent from the annotation GTEx was built on, so it was never "
      + "quantified. That is a gap, not a measurement of zero.";
  }
  const brain = entry.brain, highest = entry.highest;
  if (brain && brain.median > 0) {
    const donors = brain.sample_count ? `, ${fmt(brain.sample_count)} donors` : "";
    const elsewhere = highest && highest.tissue_id !== brain.tissue_id
      ? ` Highest anywhere: ${highest.median.toFixed(2)} in ${highest.tissue_name}.` : "";
    return `Highest of the 13 nervous-system groups: ${brain.tissue_name}${donors}.${elsewhere}`;
  }
  if (highest && highest.median > 0) {
    return `Measured at zero in all 13 nervous-system groups. Highest anywhere: `
      + `${highest.median.toFixed(2)} ${highest.unit} in ${highest.tissue_name}.`;
  }
  return `Measured in all ${fmt(entry.measured_groups)} sample groups and zero in every one. `
    + "A median of zero is a measurement, not a gap.";
}

function expressionOverviewRow(entry, selected) {
  const transcript = state.gene?.transcripts?.find(t => t.id === entry.transcript_id);
  const biotype = transcript?.biotype || "not in Ensembl's current annotation";
  const mark = transcript?.mane ? '<span class="tag-inline">MANE Select</span>'
    : transcript?.canonical ? '<span class="tag-inline muted">canonical</span>' : "";
  // A zero cell must not name a tissue. The "highest" group of an all-zero transcript is
  // whichever one sorted first, and printing "Brain - Amygdala" beside 0.00 reads as though
  // the amygdala were where it was looked for.
  const cell = (median, where, everywhere) => median > 0
    ? `<td class="num"><strong>${median.toFixed(2)}</strong>`
      + `<span class="muted"> ${escapeHTML(where)}</span></td>`
    : `<td class="num"><span class="zero">0.00</span>`
      + `<span class="muted"> in ${escapeHTML(everywhere)}</span></td>`;
  const outside = (entry.all || []).filter(m => !m.tissue_id.startsWith("Brain_"));
  const bestOutside = outside.length
    ? outside.reduce((a, b) => (b.median > a.median ? b : a)) : null;
  const cells = entry.status !== "measured"
    ? `<td colspan="2" class="not-measured">not in this GENCODE release, never quantified</td>`
    : cell(entry.brain ? entry.brain.median : 0,
           entry.brain ? entry.brain.tissue_name : "", "all 13 groups")
      + cell(bestOutside ? bestOutside.median : 0,
             bestOutside ? bestOutside.tissue_name : "", "all 41 groups");
  return `<tr class="${entry.transcript_id === selected ? "selected" : ""}">`
    + `<th scope="row"><button type="button" class="text-button" `
    + `data-choose-transcript="${escapeHTML(entry.transcript_id)}">`
    + `${escapeHTML(entry.transcript_id)}</button> ${mark}</th>`
    + `<td class="muted">${escapeHTML(biotype)}</td>${cells}`
    + `<td class="organ-cell">${organStrip(entry, 3)}</td></tr>`;
}

function renderExpressionOverview() {
  const box = $("expression-overview"), payload = state.expression;
  if (!box) return;
  if (state.source !== "gene" || !state.gene || !payload || payload.loading || payload.available === false) {
    box.hidden = true; box.innerHTML = ""; return;
  }
  box.hidden = false;
  const selected = selectedTranscriptId();
  // Picker order, not expression order. Sorting this table by expression would be the
  // recommendation ADR 0021 says the tool must not make, just moved to another surface.
  const ordered = (state.gene.transcripts || [])
    .map(t => payload.transcripts[t.id]).filter(Boolean);
  const retired = Object.values(payload.transcripts).filter(
    e => !(state.gene.transcripts || []).some(t => t.id === e.transcript_id));
  const rows = [...ordered, ...retired].map(e => expressionOverviewRow(e, selected)).join("");
  box.innerHTML = `<details class="expression-overview-detail" open>`
    + `<summary>Transcript expression · GTEx ${escapeHTML(payload.dataset_id.replace("gtex_", ""))}`
    + ` · GENCODE ${escapeHTML(payload.annotation_release)}</summary>`
    + `<p class="hint">Same order as the promoter reference list above, which is the`
    + ` annotation's order. Click an identifier to design against it.`
    + ` ${fmt(payload.measured_transcripts)} quantified,`
    + ` ${fmt(payload.silent_transcripts)} measured at zero everywhere. The two value columns`
    + ` are split so a transcript used mainly in brain can be told from one used everywhere;`
    + ` neither is a sum, and both are medians over donors.</p>`
    + `<div class="table-scroll"><table class="expression-overview-table"><thead><tr>`
    + `<th>Transcript</th><th>Biotype</th>`
    + `<th>Brain, highest of 13 groups</th><th>Highest of the other 41</th>`
    + `<th>Organs measured above zero</th>`
    + `</tr></thead><tbody>${rows}</tbody></table></div>`
    + `<p class="hint">${escapeHTML(payload.brain_note)}</p>`
    + `<p class="hint">${escapeHTML(payload.organ_note)}</p>`
    + `<details class="expression-detail"><summary>What these numbers do and do not show`
    + `</summary><ul class="caveat-list">`
    + (payload.caveat_list || payload.caveats).map(c => `<li>${escapeHTML(c)}</li>`).join("")
    + `</ul><p class="hint">${escapeHTML(payload.provenance)}</p></details></details>`;
}

async function loadExpression(symbol) {
  const asked = symbol;
  state.expression = {loading:true}; renderExpression();
  try {
    const payload = await api(`/api/genes/${encodeURIComponent(symbol)}/expression?offline=${$("offline").checked}`);
    if (state.gene?.symbol !== asked) return;   /* The user moved on while this was in flight. */
    state.expression = payload;
  } catch (failure) {
    if (state.gene?.symbol !== asked) return;
    state.expression = {available:false, transcripts:{},
      note:"No expression evidence was retrieved. This says nothing about whether these transcripts are expressed.",
      detail:failure.message};
  }
  renderExpression();
}
async function lookup() {
  const symbol = $("gene-symbol").value.trim().toUpperCase();
  if (!symbol) throw new Error("Enter a human gene symbol or Ensembl gene ID first.");
  const previousReference = state.gene?.symbol.toLowerCase()===symbol.toLowerCase() ? $("reference").value : "";
  loading("lookup-status", true, "Finding gene and promoter references…");
  $("lookup-gene").disabled = true; $("lookup-gene").textContent = "Finding…";
  try {
    const gene = await api(`/api/genes/${encodeURIComponent(symbol)}?offline=${$("offline").checked}`);
    if (state.gene && state.gene.symbol !== gene.symbol) { state.basket=[]; state.customReferences={}; }
    populateGene(gene,previousReference);
    /* Deliberately not awaited: the picker is usable before GTEx answers, and a GTEx outage
       must not stop a gene lookup from completing. */
    loadExpression(gene.symbol);
    markDraft(); return gene;
  } finally { loading("lookup-status", false); $("lookup-gene").disabled = false; $("lookup-gene").textContent = "Look up"; }
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
  $("design-ready").hidden=true;
  loading("job-status", true, "Preparing design run…");
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
  finally { loading("job-status", false); state.busy = false; $("design-button").disabled = false; $("setup").inert = false; $("saved-runs").disabled = false; }
}
function setBatchBusy(busy) {
  if (busy) $("design-ready").hidden=true;
  loading("job-status", busy, busy ? "Preparing promoter comparison…" : undefined);
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
    const ratingLabel=rating && typeof rating==="object"?(rating.evaluated?`${rating.met}/${rating.evaluated}`:"Not evaluated"):rating ?? "n/a";
    return `<tr class="${state.result?.app.run_id===run.id?"selected":""}"><td><strong>${escapeHTML(label)}</strong><p class="hint">${escapeHTML(kind)}</p>${req.custom_reference?`<p class="hint">+1 at ${fmt(req.custom_reference.tss_position)} (1-based)</p>`:""}</td><td>${fmt(req.upstream)} upstream<br>${fmt(req.downstream)} downstream</td><td>${escapeHTML(config?.preset_name || req.preset_id)}${config?.is_custom?" · custom":""}<p class="hint">${escapeHTML(summary.cell_line || req.cell_line)}</p></td><td>${counts?`${fmt(counts.passing)} / ${fmt(counts.generated)}`:"n/a"}</td><td>${escapeHTML(ratingLabel)}</td><td>${run.status==="complete"?`<button type="button" class="secondary" data-open-batch-run="${escapeHTML(run.id)}">Open run</button>`:`<span class="tag">${escapeHTML(run.status.toUpperCase())}</span><p class="hint batch-error">${escapeHTML(run.error || run.stage)}</p>`}</td></tr>`;
  }).join("");
  const completed=runs.filter(r=>r.status==="complete").length, failed=runs.filter(r=>r.status==="error").length;
  $("batch-progress").textContent=`${completed} of ${runs.length} completed${failed?` · ${failed} failed; completed runs remain available`:""}${batch.retry_of?` · retry of saved group ${batch.retry_of.slice(0,8)}`:""}`;
  $("batch-retry-online").hidden=!!frozen || batch.status==="running" || !!state.catalogue?.data_access?.offline_forced || !runs.some(r=>r.status==="error" && r.error_code==="offline-cache-miss");
  $("export-batch-html").hidden=!!frozen;
  $("job-status").textContent=batch.status==="running"?"Designing the selected references separately…":`Saved promoter comparison · ${batch.id.slice(0,8)}`;
  if (!state.result) $("empty-state").hidden=true;
  if (typeof workflowBatch === "function") workflowBatch(batch);
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
    if (typeof flushWorkspace === "function") await Promise.all(state.batch.runs.map(run=>flushWorkspace(run.id)));
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
  state.screeningEvidence = null;
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
  $("metric-best").textContent = first ? `${first.rating.met}/${first.rating.evaluated}` : "n/a";
  $("metric-length").textContent = `${fmt(result.promoter.length)} bp`;
  renderRankingPerformance(result.ranking_performance);
  renderReproducibility(result);
  renderMotifPanel(result);
  renderTopResult(result);
  const cell = result.cell_types.find(c => [c.id,c.name.toLowerCase()].includes(result.app.cell_line.toLowerCase()));
  const chromatin = cell?.chromatin;
  renderCellContext(result, cell);
  $("map-orientation").textContent = result.promoter.contig ? `chr${result.promoter.contig}:${fmt(result.promoter.start)}–${fmt(result.promoter.end)} (0-based, half-open) · DNA strand ${result.promoter.strand} · plotted 5′ → 3′ in transcript orientation` : "Plotted 5′ → 3′ as supplied · no genomic location assigned";
  $("candidate-detail").innerHTML = '<p class="hint">Select a candidate in the table or promoter map to inspect its sequences, strand prediction, and every rule verdict.</p>';
  state.sequenceRecords = result.sequence_comparisons || [];
  state.sequenceRevision = result.sequence_comparison_revision || 0;
  $("sequence-conflict-recovery").hidden=true;
  renderSequenceRecords();
  $("sequence-status").textContent="";
  if (result.candidates.some(candidate=>candidateControlMatches(candidate,result).length)) $("candidate-inspection").open=true;
  renderMap(); renderTable();
  // A screening belongs to the run, not to the browser tab that started it, so a reopened run
  // rebuilds its merged panel from the stored records rather than reading as never screened.
  // Guarded because the panel is a second script: in a frozen archive app.js runs to
  // completion before it is parsed, and a hard call here left the whole page unrendered.
  state.screenings = {}; state.screeningJobs = {}; state.screeningOpen = null;
  if (typeof renderOffTarget === "function") {
    if (frozen) for (const record of frozenScreeningRecords().filter(r=>r.status === "complete")) mergeScreening(record);
    renderOffTarget("genomic");
    renderScreeningStatus();
    if (!frozen) loadSavedScreenings(result.app.run_id);
  }
  if (typeof workflowRun === "function") workflowRun(result);
  if (first) selectCandidate(first.id, {reveal:false});
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
function ruleName(id) {
  // A rule's human name, falling back to its id. "6/7" is a number without a name; this is
  // what turns it into "missed: Position 2 is G or C".
  const rules = state.result?.configuration?.rules || [];
  return (rules.find(r => r.id === id) || {}).name || id;
}

function ratingText(c) { return c.rating.evaluated ? `${c.rating.met}/${c.rating.evaluated}` : "Not evaluated"; }
function filteredCandidates() {
  const search = $("candidate-search").value.trim().toUpperCase();
  const values = state.result.candidates.filter(c => (state.filter==="all" || (state.filter==="passing" ? c.rank!==null : c.rank===null)) && (!search || c.id.toUpperCase().includes(search) || c.target_site_dna.includes(search) || candidateControlMatches(c).some(record=>record.name.toUpperCase().includes(search))));
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
  const screeningById = new Map(typeof screenedDecisions === "function" ? screenedDecisions().map(d=>[d.candidate.id,d]) : []);
  const controlCandidates=state.result.candidates.filter(candidate=>candidateControlMatches(candidate).length);
  $("candidate-controls").innerHTML=controlCandidates.length?`<p class="hint">Known controls in this run · exact paired-core matches</p>${controlCandidates.map(candidate=>`<div class="candidate-control-summary"><strong>${escapeHTML(candidate.id)}</strong> ${candidateControlBadges(candidate)} <span>${candidate.rank===null?'Excluded by current gates':`Design rank ${candidate.rank}`}</span></div>`).join('')}` : '';
  $("visible-count").textContent = `· ${fmt(values.length)}`;
  $("candidate-rows").innerHTML = page.length ? page.map(c => `<tr data-candidate="${c.id}" class="${state.selected===c.id?"selected":""}" tabindex="0" aria-label="Inspect ${c.id}"><td><input type="checkbox" data-check="${c.id}" ${state.checked.has(c.id)?"checked":""} aria-label="Select ${c.id} for export"></td><td>${c.rank || "n/a"}</td><td><span class="dna">${escapeHTML(c.target_site_dna)}</span><p class="candidate-id">${c.id}${c.rank===null?` · ${c.excluded_by.length} failed gates`:""}</p>${twinChip(c)}${candidateControlBadges(c)}<div class="candidate-screening-note" data-candidate-screening="${c.id}">${typeof candidateScreeningSummary === "function" ? candidateScreeningSummary(screeningById.get(c.id)) : ""}</div></td><td>${escapeHTML(positionText(c))}</td><td>${c.gc_percent}%</td><td><div class="rating-cell"><b>${ratingText(c)}</b><span class="rating-bar"><i style="width:${c.rating.percent ?? 0}%"></i></span></div>${c.rating.unavailable?`<p class="candidate-id">${c.rating.unavailable} unavailable</p>`:""}${c.rating.unchecked_gates.length?`<p class="candidate-id">${c.rating.unchecked_gates.length} unchecked gates</p>`:""}${(c.rating.missed||[]).length?`<p class="missed-chip">missed: ${c.rating.missed.map(id=>escapeHTML(ruleName(id))).join(", ")}</p>`:""}</td><td><span class="guide-badge ${escapeHTML(c.rating.guide || "unknown")}">${escapeHTML(c.rating.guide || "unknown")}</span>${!c.rating.models_agree?'<p class="disagreement">Models disagree</p>':""}</td></tr>`).join("") : '<tr><td colspan="7" style="padding:30px;text-align:center;color:var(--muted)">No candidates match this view. Inspect excluded targets or adjust your criteria and run again.</td></tr>';
  $("page-label").textContent = values.length ? `${fmt(state.page*50+1)}–${fmt(Math.min((state.page+1)*50,values.length))} of ${fmt(values.length)} · page ${state.page+1}/${pages}` : "0 candidates";
  $("page-previous").disabled = state.page===0; $("page-next").disabled = state.page>=pages-1;
  $("select-page").checked = page.length>0 && page.every(c=>state.checked.has(c.id));
  $("select-page").indeterminate = page.some(c=>state.checked.has(c.id)) && !$("select-page").checked;
  renderExportNote();
  // The ticked selection is one of the batch choices, so the button's count follows the table.
  if (typeof renderBatchControls === "function") renderBatchControls();
}
function renderExportNote() {
  const selection = state.checked.size ? `${state.checked.size} checked candidates` : "all design-passing candidates (no explicit selection)";
  $("export-note").textContent = `Excel and all three CSV downloads use ${selection}. JSON and HTML preserve the whole run and its screening history. SVG follows the screening plot’s range and sort controls. Each file’s contents are listed below.`;
  if (typeof workflowStatus === "function") workflowStatus();
  if (state.workspaceTab === "export" && typeof renderExportSelection === "function") renderExportSelection();
}

/** What this run can honestly say about off-target screening, for whichever candidate.
 *
 * Three outcomes, as everywhere else: none run, run and quiet, run and found something. The
 * three places that said this each hardcoded the first one, so a screened run went on
 * reporting that nothing had been screened.
 */
function screeningSentence(candidateId = null) {
  const held = Object.values(state.screenings || {});
  if (!held.length) {
    return Object.values(state.screeningJobs || {}).some(job => ["running", "waiting"].includes(job.status))
      ? "Off-target searches are pending; no completed findings are available yet."
      : "No off-target screening has been performed.";
  }
  const scopes = held.map(h => h.scope === "genomic" ? "promoters" : "transcripts").join(" and ");
  if (!candidateId) return `Off-target screening has been run against ${scopes}; see the panel below.`;
  const rows = held.map(h => h.byCandidate[candidateId]).filter(Boolean);
  if (!rows.length) {
    return `This candidate was not included in the ${scopes} screening that was run, which is `
      + `not the same as having been screened and found clean.`;
  }
  const strong = rows.reduce((n, r) => n + (r.row.wang_homology_hits || 0), 0);
  const exact = rows.reduce((n, r) => n + (r.row.perfect_hits || 0), 0);
  return `Screened against ${scopes}: ${exact} full-length exact match`
    + `${exact === 1 ? "" : "es"} elsewhere and ${strong} alignment`
    + `${strong === 1 ? "" : "s"} at Wang 2015's 80%. Similarity is not measured binding.`;
}
async function selectCandidate(id, {reveal=true} = {}) {
  if (reveal && $("candidate-inspection")) $("candidate-inspection").open = true;
  const run = state.result, candidate = run.candidates.find(c=>c.id===id);
  if (!candidate) return;
  state.selected = id; renderTable(); renderMap();
  if (!candidate.rules && !frozen) {
    $("candidate-detail").innerHTML = '<p class="hint is-loading" role="status">Loading candidate sequences and rule explanations…</p>';
  }
  try {
    if (frozen && !candidate.rules) throw new Error("This archival record is missing the candidate’s full verdicts.");
    const full = candidate.rules ? candidate : await api(`/api/runs/${run.app.run_id}/candidates/${id}`);
    if (state.result!==run || state.selected!==id) return;
    Object.assign(candidate,full);
    renderSequenceDetail();
    const primary = full.guide_prediction.assessments.find(a=>a.table===full.rating.guide_model);
    const sequenceBox = (name,sequence) => `<div class="sequence-box"><div><span>${name} · 5′ → 3′</span><button type="button" data-copy="${escapeHTML(sequence)}">Copy</button></div><code>${escapeHTML(sequence)}</code></div>`;
    $("candidate-detail").innerHTML = `<div class="detail-heading"><div><h3>${id} <span class="muted">${full.rank?`Rank ${full.rank}`:"Excluded by gates"}</span></h3><p class="hint">${escapeHTML(positionText(full))} ${full.tss_span?"relative to reference TSS":"in supplied sequence"} · ${full.target_site_dna.length} nt core · ${full.gc_percent}% GC${full.genomic?` · genomic ${escapeHTML(full.genomic.strand)} strand`:""}</p></div><span class="tag">${ratingText(full)} PREFERENCES</span></div>${candidateControlBadges(full,run)}${duplexPanel(run,full)}<div class="detail-body"><div>${sequenceBox("Promoter target DNA",full.target_site_dna)}${sequenceBox("Sense RNA core",full.sense_strand_rna)}${sequenceBox("Antisense RNA core",full.antisense_strand_rna)}<p class="hint">The cores above are the 19 nt RNA the rules are evaluated against, without overhangs. The duplex panel at the top is what you order. ${escapeHTML(screeningSentence(id))}</p>${twinNote(full)}</div><div><div class="guide-summary"><strong>Predicted RNA guide: ${escapeHTML(full.rating.guide || "unknown")}</strong><p>${escapeHTML(tableNames[full.rating.guide_model] || full.rating.guide_model)} · ${escapeHTML(full.guide_prediction.summary)}</p><p class="hint">This predicts duplex loading, not the identity of a promoter-associated target RNA, and not activation.</p></div><table class="model-table"><thead><tr><th>Parameter model</th><th>Guide</th><th>End difference</th></tr></thead><tbody>${full.guide_prediction.assessments.map(a=>`<tr data-tip="${escapeHTML(a.explanation)}" tabindex="0"><td>${escapeHTML(tableNames[a.table] || a.table)}</td><td>${escapeHTML(a.predicted_guide)}</td><td>${a.difference.toFixed(2)} ${escapeHTML(a.unit)}</td></tr>`).join("")}</tbody></table><details class="verdict"><summary>How the selected model reaches this prediction</summary><p>${escapeHTML(primary?.explanation || "The selected model was not assessed.")}</p></details></div></div>${full.excluded_by.length?`<p class="hint">${escapeHTML(full.rejection_summary)}</p>`:""}<div class="rule-verdicts"><h3>Every rule, explained</h3>${full.rules.map(r=>ruleVerdict(r,run.configuration)).join("")}</div>`;
    if (typeof workflowCandidate === "function") workflowCandidate(candidate);
  } catch (e) {
    if (state.result!==run || state.selected!==id) return;
    $("candidate-detail").innerHTML = '<p class="hint">Candidate details could not be loaded. Select the candidate to try again.</p>';
    error(e.message);
  }
}
function ruleVerdict(verdict,config) {
  const def = config.rules.find(r=>r.id===verdict.id);
  // The group, which the engine has always carried and nothing displayed. Without it a rule
  // named "Position 20 is A or T" sits in a candidate whose core is 19 nt long and reads as a
  // contradiction; "3' flank A/U richness" says in three words that the position is outside
  // the target.
  const group = def?.group
    ? `<span class="rule-group">${escapeHTML(def.group)}</span>` : "";
  return `<details class="verdict"><summary><span class="outcome ${verdict.outcome}">${verdict.outcome==="not-applicable"?"N/A":verdict.outcome.toUpperCase()}</span><span>${escapeHTML(def?.name || verdict.id)}</span>${group}<span class="muted">${roleName[verdict.kind]}</span><span class="why">why \u203a</span></summary><p>${escapeHTML(verdict.detail)}</p>${def?`<p class="source">${provenanceText(def.provenance)} · ${escapeHTML(evidenceLabel(def))}</p><p><b>Why this rule:</b> ${escapeHTML(def.rationale)}</p>${ruleEvidence(def)}<p><b>Changing it:</b> ${escapeHTML(def.consequence)}</p><p><b>Saved values:</b> ${def.parameters.map(p=>`${escapeHTML(p.name)} = ${escapeHTML(p.value)} (${escapeHTML(p.origin)})`).join("; ") || "No adjustable numeric parameters"}</p>`:""}</details>`;
}

const REJECTION_COLOURS = ["#b4714a","#7a93c2","#6f9e7a","#b2608f","#8d8248"];
const OTHER_REJECTION = "#6e7d79";  // 3.9:1 on white; the previous grey was 2.47:1

// The track colours one mark per excluded candidate, by the first gate that excluded it, so
// the colours are shared out by that same attribution. They were shared out by
// rejection_breakdown, which counts every gate a candidate failed: on a 600 bp fixture that
// spent a colour on a rule accounting for 7 of 574 marks while the third largest actual cause,
// 27 marks, fell into the grey "other". Use rejectionCauses for anything describing the
// drawing, and rejection_breakdown only where overlapping totals are what is being reported.
function rejectionCauses(result) {
  if (result.rejection_attribution) return result.rejection_attribution;
  // A run saved before the engine recorded this has candidates but no attribution. Deriving it
  // here from the same field the track colours by keeps an archived run drawing and labelled
  // identically to a fresh one, rather than silently losing its legend and turning every mark
  // grey. Fresh runs never reach this branch.
  const excluded = (result.candidates || []).filter(c => (c.excluded_by || []).length);
  const counts = new Map();
  for (const candidate of excluded) {
    const first = candidate.excluded_by[0];
    counts.set(first, (counts.get(first) || 0) + 1);
  }
  return [...counts].sort((a, b) => b[1] - a[1]).map(([rule, rejected]) => ({
    rule,
    rejected,
    percent_of_excluded: excluded.length
      ? Math.round(1000 * rejected / excluded.length) / 10 : 0,
    percent_of_generated: (result.candidates || []).length
      ? Math.round(1000 * rejected / result.candidates.length) / 10 : 0,
  }));
}

// A candidate is a locus; the molecule follows from its target bases alone. The same 19-mer at
// several positions is therefore several rows and one product, and ticking them all buys it
// several times. The grouping is computed in Python and carried in the run, so this chip, the
// detail panel and the order CSV cannot disagree about which rows are the same molecule.
function duplexTwins(id) {
  return (state.result?.identical_duplex_groups || {})[id] || [];
}

// How many twins are named before the rest are counted.
const TWINS_NAMED = 8;
function twinNote(c) {
  const twins = duplexTwins(c.id);
  if (!twins.length) return "";
  // Named rather than counted, because the next question is always "which ones", and the
  // differing ratings are the reason they are listed instead of merged away. Capped, because
  // a tandem-repeat promoter can put the same 19-mer at dozens of positions and a wall of
  // identifiers stops being readable long before it stops being true.
  const shown = twins.slice(0, TWINS_NAMED);
  const unnamed = twins.length - shown.length;
  return `<p class="twin-note"><b>Identical duplex elsewhere in this run.</b> `
    + `${shown.map(id => `<button type="button" class="text-button" data-candidate="${
        escapeHTML(id)}">${escapeHTML(id)}</button>`).join(", ")}`
    + `${unnamed ? ` and ${unnamed} further position${unnamed === 1 ? "" : "s"}` : ""}`
    + ` order as exactly these two `
    + `21-mers. They are different promoter positions, so the flank and window rules can rate `
    + `them differently, but a vendor would synthesise one product. The order CSV repeats this `
    + `per row.</p>`;
}

function twinChip(c) {
  const twins = duplexTwins(c.id);
  if (!twins.length) return "";
  return `<p class="twin-chip" data-tip="${escapeHTML(
    `The ordered 21-mers here are identical to ${twins.join(", ")}. Those are other promoter `
    + `positions, and the flank and window rules can rate them differently, but a vendor would `
    + `make one product. Order one.`)}" tabindex="0">same duplex as ${twins.length} other`
    + `${twins.length === 1 ? "" : "s"}</p>`;
}

function rejectionPalette(result) {
  const palette = new Map();
  rejectionCauses(result).slice(0, REJECTION_COLOURS.length)
    .forEach((entry, i) => palette.set(entry.rule, REJECTION_COLOURS[i]));
  return palette;
}

function swatch(entry) {
  return `<i style="background:${entry.colour.fill};border-color:${entry.colour.edge}"></i>`;
}

function renderLegend(result) {
  // Built from TRACKS, in the order the map draws them, so each entry names the row it
  // belongs to rather than leaving the reader to match a colour by eye across seven lines.
  // Collapsed it is one line of chips; opened it says what each track means.
  const legend = $("map-legend");
  if (!legend) return;
  const rejections = rejectionCauses(result).slice(0, REJECTION_COLOURS.length);
  const palette = rejectionPalette(result);
  const motifColours = motifPalette(result);
  const motifSwatches = [...motifColours].map(([name, fill]) =>
    ({colour:{fill, edge:fill}, name}));
  const swatchesFor = track => track.key === "motifs" && motifSwatches.length
    ? [...motifSwatches, {colour:TRACK_COLOURS.motif, name:"other factor"},
       {colour:TRACK_COLOURS.motifDim, name:"dimmed, while one factor is selected"}]
    : track.swatches;
  const chips = TRACKS.map(track => {
    const marks = track.key === "excluded"
      ? rejections.map(entry =>
          `<i style="background:${palette.get(entry.rule)};border-color:${palette.get(entry.rule)}"></i>`).join("")
        + (rejections.length ? "" : `<i style="background:${OTHER_REJECTION};border-color:${OTHER_REJECTION}"></i>`)
      : swatchesFor(track).map(swatch).join("");
    const name = track.key === "motifs" && state.motifFocus
      ? `${track.name}: ${state.motifFocus} only` : track.name;
    return `<span class="legend-chip${track.key === "motifs" && state.motifFocus
      ? " focused" : ""}">${marks}${escapeHTML(name)}</span>`;
  }).join("");
  const rows = TRACKS.map(track => {
    // The excluded track's colours are per-rule and listed with their percentages in the
    // breakdown directly above the legend. Repeating all five here said the same thing twice
    // and made this the longest row in the table.
    const marks = track.key === "excluded" ? "see the breakdown above"
      : swatchesFor(track).map(entry =>
          `<span class="legend-mark">${swatch(entry)}${escapeHTML(entry.name)}</span>`).join("");
    return `<tr><th>${escapeHTML(track.name)}</th><td>${marks || "none"}</td>`
      + `<td>${escapeHTML(track.help)}</td></tr>`;
  }).join("");
  legend.innerHTML = `<summary><span class="legend-chips">${chips}</span>`
    + `<span class="legend-toggle">what am I looking at? \u203a</span></summary>`
    + `<table class="legend-table"><tbody>${rows}</tbody></table>`
    + `<p class="hint">Drag to pan \u00b7 scroll to zoom \u00b7 hover or focus any mark for its `
    + `evidence \u00b7 click a candidate to inspect it. Coordinates and annotations identify `
    + `features; they do not establish activity in your cell line.</p>`;
}

function renderRejectionLegend(result) {
  // Says which rule is actually constraining this promoter, in the same move as the colours.
  // Every number here is a share of the marks on the track, so the legend and the picture
  // cannot disagree. The overlapping count, where a candidate failing four gates is counted
  // four times, is a different measurement and is reported under "Rule effects" in the run
  // evidence, with its own caveat.
  const causes = rejectionCauses(result);
  const legend = $("rejection-legend");
  if (!legend) return;
  if (!causes.length) { legend.innerHTML = ""; return; }
  const palette = rejectionPalette(result);
  const rest = causes.slice(REJECTION_COLOURS.length);
  const restShare = Math.round(rest.reduce((sum, e) => sum + e.percent_of_excluded, 0) * 10) / 10;
  legend.innerHTML = "<span>First gate to exclude:</span>"
    + causes.slice(0, REJECTION_COLOURS.length).map(entry =>
        `<span><i style="background:${palette.get(entry.rule)}"></i>`
        + `${escapeHTML(ruleName(entry.rule))} \u00b7 ${entry.percent_of_excluded}%</span>`).join("")
    + (rest.length
        ? `<span><i style="background:${OTHER_REJECTION}"></i>`
          + `${rest.length} further ${rest.length === 1 ? "rule" : "rules"} \u00b7 ${restShare}%`
          + "</span>"
        : "")
    // Which mode the track is drawn in is a statement about what the reader is looking at, so
    // it is said rather than left to be inferred from the shape of the bars.
    + `<small>Each excluded target is counted once, against the first gate in this run's`
      + ` evaluation order that excluded it, so these shares are of the`
      + ` ${fmt(result.counts.rejected)} excluded targets and add up to 100%. Most failed`
      + ` several gates; the run evidence lists every failure per rule. `
      + `${state.excludedBinned
        ? "Excluded targets are binned by position, and each column is split by the share of "
          + "that bin each rule rejected. Every column is the same height, which carries no "
          + "meaning. Hover for counts; click a column, or zoom below " + INDIVIDUAL_MARK_SPAN
          + " bp, to see each target."
        : "Every excluded target in view is drawn individually; click one to inspect it."}`
    + "</small>";
}

// One description of the tracks, used to draw the map, to label its rows, and to build the
// legend. They were three separate lists and had drifted: the legend's swatch for the
// regulatory build was #e2b772 while the map drew #ddb370 and #efcf99, and its motif swatch
// was #a798d2 against a map that drew #6a5ca8. A legend whose colours are not the map's
// colours is worse than no legend, because it is confidently wrong.
//
// Fills stay light so bands do not dominate; each carries a darker edge, which is what gives
// the boundary the 3:1 contrast a graphic mark needs without darkening the page.
const TRACK_COLOURS = {
  tss:          {fill:"#344d55", edge:"#243840"},
  altTss:       {fill:"#8eaaa0", edge:"#5c7a70"},
  window:       {fill:"#c2dfd2", edge:"#6f9e87"},
  preferred:    {fill:"#75bfa7", edge:"#3f8a72"},
  cpg:          {fill:"#aac7c9", edge:"#5f8d90"},
  // Reddish, because a promoter is the feature this tool exists to target and it read as
  // just another shade of the same yellow.
  promoter:     {fill:"#e8b5a6", edge:"#b4573f"},
  regulatory:   {fill:"#efcf99", edge:"#b8882f"},
  // The neutral for a factor past the end of the per-factor palette. It must not *be* a
  // palette colour: while it was #6a5ca8, SOX2 fell off the end and rendered identically to
  // KLF4, and nothing on screen distinguished them. 3.6:1 on white.
  motif:        {fill:"#7f8794", edge:"#5c636e"},
  motifDim:     {fill:"#d7d1ea", edge:"#b3a9d4"},
  candidate:    {fill:"#087e80", edge:"#055e60"},
  // The highest-ranked targets. A darker shade of the same teal, because they are passing
  // candidates rather than a different kind of thing -- and colour is only one of the three
  // channels separating them: they are also full lane height and numbered.
  topCandidate: {fill:"#0b4f57", edge:"#063136"},
  referencePositive: {fill:"#7954a1", edge:"#7954a1"},
  referenceDown: {fill:"#b43f63", edge:"#b43f63"},
  imported: {fill:"#167a79", edge:"#167a79"},
};

// Rank 1 to TOP_MARKED as one hue getting paler, so the ten read as *ordered* at a glance.
// Deliberately not green-to-red: ranking is preferences met, counted equally (ADR 0006), and
// it does not beat chance on the 125 assayed sequences (AUC 0.41-0.53), so a scale that says
// "green is good, red is bad" would make a claim the data refutes. A single hue says "these
// are in order" and stops there. It is also the only option that survives colour blindness
// and does not collide with the amber regulatory track or the reddish promoter track.
const RANK_DARKEST = [11, 79, 87];      // #0b4f57, rank 1
const RANK_PALEST = [158, 196, 198];    // #9ec4c6, rank TOP_MARKED, still visible on white

function rankColour(rank) {
  // Pale fills need a darker edge or the last marks lose their outline against the page.
  const t = TOP_MARKED > 1 ? Math.min(1, Math.max(0, (rank - 1) / (TOP_MARKED - 1))) : 0;
  const mix = (a, b, amount) => Math.round(a + (b - a) * amount);
  const channel = (amount) => RANK_DARKEST.map((v, i) => mix(v, RANK_PALEST[i], amount));
  const hex = (rgb) => "#" + rgb.map(v => v.toString(16).padStart(2, "0")).join("");
  return {fill: hex(channel(t)), edge: hex(channel(t).map(v => Math.round(v * 0.55)))};
}

// How many of the highest-ranked targets the map marks out. The candidate table already
// treats the same set as a unit with its "Select top 10" button, so the two agree.
const TOP_MARKED = 10;

const TRACKS = [
  {key:"reference", y:59, name:"Reference / TSS",
   help:"An annotated transcript start or a chosen literature anchor. This is the reference for coordinates, not proof of transcription in your cell line.",
   swatches:[{colour:TRACK_COLOURS.tss, name:"chosen +1"},
             {colour:TRACK_COLOURS.altTss, name:"other annotated starts"}]},
  {key:"windows", y:103, name:"Design windows",
   help:"Bands use the frozen settings of this run. Gates exclude failures; preferences contribute one equal point.",
   swatches:[{colour:TRACK_COLOURS.window, name:"gate window"},
             {colour:TRACK_COLOURS.preferred, name:"preferred window"}]},
  {key:"cpg", y:142, name:"CpG islands",
   help:"DNA CpG-rich regions called using this run\u2019s length, GC and observed/expected thresholds. Absence in a short window does not establish absence outside it.",
   swatches:[{colour:TRACK_COLOURS.cpg, name:"island"}]},
  {key:"regulatory", y:176, name:"Regulatory build",
   help:"Ensembl Regulatory Build pools observations across cell types. These annotations do not establish activity in the selected cell line.",
   swatches:[{colour:TRACK_COLOURS.promoter, name:"promoter"},
             {colour:TRACK_COLOURS.regulatory, name:"other feature"}]},
  {key:"motifs", y:211, name:"TF motif matches",
   help:"Sequence predictions from the bundled motif matrices, one colour per factor, assigned to the most enriched first. A match is not measured factor occupancy. Hover for the matrix, score and strand.",
   swatches:[{colour:TRACK_COLOURS.motif, name:"a match"},
             {colour:TRACK_COLOURS.motifDim, name:"dimmed, while one factor is selected"}]},
  {key:"passing", y:245, name:"Passing targets",
   help:"Targets that pass every evaluable mandatory rule. The ten highest-ranked are drawn full height, numbered, and shaded from dark at rank 1 to pale at rank 10. The shading is order only. It is deliberately not a green-to-red scale, because rank is preferences met counted equally, and on 125 assayed sequences this ranking did not beat chance, so colouring rank 1 'good' would claim something the data refutes. A numbered mark is a shortlist to inspect, not a prediction of activity. Unavailable gates must still be examined in candidate details.",
   swatches:[{colour:rankColour(1), name:"rank 1"},
             {colour:rankColour(Math.round(TOP_MARKED / 2)), name:"rank 5"},
             {colour:rankColour(TOP_MARKED), name:"rank 10, palest"},
             {colour:TRACK_COLOURS.candidate, name:"passes every evaluable gate"}]},
  {key:"excluded", y:273, name:"Excluded targets",
   help:"Targets with at least one failed mandatory rule, coloured by the first rule that excluded each. Hover or select to see which checks failed.",
   swatches:[]},
  {key:"existing", y:340, name:"Existing saRNAs",
   help:"Exact sequence placements of curated and imported duplex cores. Effects describe the cited experimental setting, not new predictions of efficacy.",
   swatches:[{colour:TRACK_COLOURS.referencePositive, name:"activating reference"},
             {colour:TRACK_COLOURS.referenceDown, name:"downregulating (researcher reported)"},
             {colour:TRACK_COLOURS.imported, name:"imported, effect unassigned"}]},
];

// One colour per factor, so the motif row says *which* factor is where rather than only that
// something matched. Assigned in the panel's own order -- enriched first, by fold -- so the
// factors worth looking at get the distinct hues and any beyond the palette share the neutral
// purple. Each is at least 4.5:1 on white, well past the 3:1 a graphic mark needs, and the
// family is deliberately unlike the rejection palette so the two tracks do not read alike.
const MOTIF_COLOURS = ["#6a5ca8","#2f7d95","#9c5420","#3f7d4e","#a2436b","#6b6f24"];

function motifPalette(result) {
  const drawn = new Set((result.tracks?.motifs || []).map(t => t.label));
  const ranked = [...(result.motif_scans || [])]
    .filter(m => drawn.has(m.motif))
    .sort((a, b) => (b.exceeds_chance - a.exceeds_chance)
                 || (b.enrichment - a.enrichment) || (b.sites - a.sites));
  return new Map(ranked.slice(0, MOTIF_COLOURS.length)
    .map((m, index) => [m.motif, MOTIF_COLOURS[index]]));
}

const EXCLUDED_TOP = 256, EXCLUDED_HEIGHT = 22, EXCLUDED_BINS = 104;
// How wide a view can be before excluded targets stop being separable marks. A 19 nt target
// drawn across 832 units is distinguishable from its neighbour only while 832/span >= 3.
const INDIVIDUAL_MARK_SPAN = 280;

function drawExcluded(parts, r, rejectionColours, scale, tip, rect) {
  // A target is generated at every offset, so across the full 2,500 bp region there are about
  // three excluded candidates per pixel of the 832-unit plot. Drawn one rect each they merge
  // into a solid band that claims to show every rejection and in fact shows none of them: the
  // marks overlap, and 2,468 of them cost 14 ms a frame to redraw.
  //
  // Above the separability threshold they are therefore binned into a stacked histogram, which
  // says something the ribbon could not -- *where* along the promoter each rule does its
  // rejecting. Below it every candidate is drawn and clickable as before. Nothing is dropped
  // either way: the bin reports its own counts by rule, and clicking it zooms in far enough
  // for the individual marks to appear.
  const excluded = r.candidates.filter(c => c.rank === null);
  const cause = c => c.rating?.excluded_by_first || c.excluded_by[0] || "other";
  if (state.end - state.start <= INDIVIDUAL_MARK_SPAN) {
    for (const c of excluded) {
      rect(c.offset, c.offset + c.target_site_dna.length, 264 + (c.offset % 3) * 3, 3,
        rejectionColours.get(cause(c)) || OTHER_REJECTION,
        `${c.id}: excluded by ${c.excluded_by.join(", ")}. Position ${positionText(c)}. `
        + `Click for details.`, `data-candidate="${c.id}"`);
    }
    return false;
  }
  const width = (state.end - state.start) / EXCLUDED_BINS;
  const bins = Array.from({length: EXCLUDED_BINS}, () => new Map());
  let any = false;
  for (const c of excluded) {
    if (c.offset < state.start || c.offset >= state.end) continue;
    const index = Math.min(EXCLUDED_BINS - 1, Math.floor((c.offset - state.start) / width));
    const counts = bins[index], key = cause(c);
    counts.set(key, (counts.get(key) || 0) + 1);
    any = true;
  }
  if (!any) return false;
  // Each column is full height and split by share, not scaled by count. A target is generated
  // at every offset and on this promoter nearly all are excluded, so bin totals are constant
  // by construction: measured over SCN1A, 99 of 104 bins sat within a tenth of the largest.
  // Scaling height by count would spend the whole visual channel on that noise and draw a
  // ragged top edge implying a signal there is not one. What varies, and what the track exists
  // for, is *which* rule rejects where. Absolute totals stay in the tooltip.
  bins.forEach((counts, index) => {
    if (!counts.size) return;
    const from = state.start + index * width, to = from + width;
    const x1 = scale(from), x2 = scale(to), columnWidth = Math.max(1, x2 - x1);
    const ordered = [...counts].sort((a, b) => b[1] - a[1]);
    let total = 0;
    for (const [, n] of ordered) total += n;
    let y = EXCLUDED_TOP + EXCLUDED_HEIGHT;
    for (const [rule, n] of ordered) {
      const height = (n / total) * EXCLUDED_HEIGHT;
      y -= height;
      parts.push(`<rect x="${x1}" y="${y}" width="${columnWidth}" height="${height}" `
        + `fill="${rejectionColours.get(rule) || OTHER_REJECTION}"/>`);
    }
    // One transparent column carries the hover and the click, so the whole bar responds
    // rather than only whichever thin segment the pointer happens to land on.
    const where = r.promoter.tss_offset === null
      ? `${Math.round(from) + 1}\u2026${Math.round(to)}`
      : `${signed(Math.round(from - r.promoter.tss_offset) || 1)}\u2026`
        + `${signed(Math.round(to - r.promoter.tss_offset) || 1)}`;
    parts.push(`<rect x="${x1}" y="${EXCLUDED_TOP}" width="${columnWidth}" `
      + `height="${EXCLUDED_HEIGHT}" fill="transparent" `
      + `data-zoom-from="${Math.round(from)}" data-zoom-to="${Math.round(to)}" `
      + tip(`${where}: ${total} excluded target${total === 1 ? "" : "s"}. `
            + ordered.map(([rule, n]) =>
                `${ruleName(rule)}: ${n} (${Math.round(100 * n / total)}%)`).join(" \u00b7 ")
            + `. Click to zoom in and see each one.`) + `/>`);
  });
  return true;
}

const PASSING_TOP = 232, PASSING_STEP = 5, PASSING_HEIGHT = 7;
// 11px rather than 13: with rank also carried by the fill, a numeral that is dropped for
// crowding costs less, but five of the ten were being dropped at full extent, which is more
// than it needs to be.
const MARKED_TOP = 230, MARKED_HEIGHT = 19, MARKED_LABEL_GAP = 9;

function passingTip(c) {
  return `${c.id} · rank ${c.rank} · ${ratingText(c)} evaluated preferences · `
    + `${c.rating.guide || "unknown"} guide (${c.rating.guide_model}). `
    + `Position ${positionText(c)}. ${candidateControlMatches(c).map(r=>`Same duplex core as ${r.name} (${sequenceAppearance(r).label})`).join("; ")} Click for details.`;
}

function drawPassing(parts, r, scale, tip, rect) {
  // Every passing target used to be the same teal mark, so rank 1 and rank 900 were
  // indistinguishable and the map could not answer the question it is looked at for: where
  // are the ones I would actually order? The top ten are now separated on three channels at
  // once -- full lane height instead of a third of it, a darker fill, and the rank printed
  // above the mark -- because any one of the three alone fails somebody: colour fails a
  // colour-blind reader, height alone does not say *which* rank, and a numeral is too small
  // to find by eye at full extent.
  //
  // Rank is preferences met, counted equally (ADR 0006). The numeral is a shortlist to
  // inspect, not a measured probability, and both the legend and the banner above the map
  // say so.
  const passing = r.candidates.filter(c => c.rank !== null);
  for (const c of passing.filter(c => c.rank > TOP_MARKED)) {
    rect(c.offset, c.offset + c.target_site_dna.length,
      PASSING_TOP + ((c.rank - 1) % 3) * PASSING_STEP, PASSING_HEIGHT,
      TRACK_COLOURS.candidate, passingTip(c), `data-candidate="${c.id}"`);
  }
  const marked = passing.filter(c => c.rank <= TOP_MARKED);
  // Worst rank first, so rank 1 is drawn last and is the one on top where they overlap.
  for (const c of [...marked].sort((a, b) => b.rank - a.rank)) {
    rect(c.offset, c.offset + c.target_site_dna.length, MARKED_TOP, MARKED_HEIGHT,
      rankColour(c.rank), passingTip(c), `data-candidate="${c.id}"`);
  }
  // Numerals in rank order, skipping any that would collide with one already placed: a
  // better rank keeps its label, and the rank of an unlabelled mark is still in its tooltip.
  // Without this, two adjacent top-ten targets -- common, since neighbouring windows share
  // most of their sequence -- overprint into an unreadable smudge.
  const placed = [];
  for (const c of [...marked].sort((a, b) => a.rank - b.rank)) {
    const middle = c.offset + c.target_site_dna.length / 2;
    if (middle < state.start || middle >= state.end) continue;
    const x = scale(middle);
    if (placed.some(at => Math.abs(at - x) < MARKED_LABEL_GAP)) continue;
    placed.push(x);
    parts.push(`<text class="rank-label" x="${x.toFixed(1)}" y="${MARKED_TOP - 3}" `
      + `text-anchor="middle">${c.rank}</text>`);
  }
}

function candidateControlMatches(candidate,result=state.result) {
  if (!candidate || !result) return [];
  const target=candidate.target_site_dna;
  return (result.sequence_comparisons || []).filter(record=>{
    if (!record.reference) return false;
    const dna=record.dna;
    if (!dna || dna.length!==target.length) return false;
    const complement=dna.split('').reverse().map(base=>({A:'T',T:'A',C:'G',G:'C'})[base]).join('');
    return target===dna || target===complement;
  });
}
function candidateControlBadges(candidate,result=state.result) {
  return candidateControlMatches(candidate,result).map(record=>{
    const style=sequenceAppearance(record);
    return `<span class="sequence-effect ${style.kind}" title="Exact paired-core sequence identity; overhangs and modifications may differ. ${escapeHTML(style.description)}">Same core as ${escapeHTML(record.name)} · ${style.label}</span>`;
  }).join(' ');
}
function referenceDuplexPanel(record,match=null) {
  const analysis=match?.strand_analysis || record.strand_analysis;
  if (!analysis) return '<p class="hint">Locate this sequence or assign its strand role to display the duplex.</p>';
  const models=analysis.models, matchedCandidate=match?.candidate || record.matches?.[0]?.candidate;
  const configured=state.result?.configuration?.rules?.find(rule=>rule.id==='thermodynamic-asymmetry')
    || state.result?.configuration?.rules?.find(rule=>rule.id==='asymmetry-margin');
  const preferred=matchedCandidate?.rating.guide_model
    || configured?.parameters?.find(parameter=>parameter.name==='nearest-neighbour-table')?.value
    || 'rna-xia-1998';
  const model=models.find(model=>model.table===preferred) || models[0];
  if (!model) return '';
  const roles=strandRoles(model.predicted_guide);
  const sense=analysis.strands.find(strand=>strand.name==='sense').rna_5to3;
  const anti=analysis.strands.find(strand=>strand.name==='antisense').rna_5to3;
  const candidate={rating:{guide:model.predicted_guide,guide_model:model.table},guide_prediction:{
    end_length_nt:matchedCandidate?.guide_prediction?.end_length_nt || 4, assessments:models}};
  const tube=(name,sequence)=>`<div class="tube ${roles[name]}"><div class="tube-head"><span class="role-chip ${roles[name]}">Predicted ${roles[name]==='undecided'?'undetermined':roles[name]}</span><span>${name} · 5′ → 3′ · ${sequence.length} nt core</span><button type="button" data-copy="${escapeHTML(sequence)}">Copy</button></div><code>${escapeHTML(sequence)}</code><small>Complementary genomic DNA strand: ${escapeHTML(analysis.strands.find(strand=>strand.name===name).complementary_genomic_strand || 'unknown')}</small></div>`;
  return `<div class="duplex-panel reference-duplex"><div class="duplex-head"><h4>${escapeHTML(record.name)} · guide &amp; passenger prediction</h4><p class="hint">${escapeHTML(model.table_name)}. ${models.some(other=>other.predicted_guide!==model.predicted_guide)?'Models disagree; strand loading is uncertain.':'Prediction describes strand loading, not measured binding.'} The diagram shows paired cores only; recorded overhangs are listed separately. ${record.partner_inferred?'The partner strand is inferred.':''}</p></div>${duplexLadder({sense_core_rna_5to3:sense,antisense_core_rna_5to3:anti,core_only:true},candidate)}${endEvidence(candidate)}<div class="tubes">${tube('sense',sense)}${tube('antisense',anti)}</div></div>`;
}
function strandAnalysisPanel(record,match=null) {
  const analysis=match?.strand_analysis || record.strand_analysis;
  if (!analysis) return '<p class="hint">Assign a strand role or locate this sequence to analyse sense/antisense orientation.</p>';
  const strands=analysis.strands;
  const bases=sequence=>`<div class="indexed-bases" aria-label="${escapeHTML(sequence)}">${[...sequence].map((base,i)=>`<span class="base-${base}"><small>${i+1}</small><b>${base}</b></span>`).join('')}</div>`;
  const nucleotideCards=strands.map(strand=>`<article class="strand-analysis-card"><h4>${escapeHTML(strand.name)} RNA · 5′ → 3′</h4>${bases(strand.rna_5to3)}<p>${strand.length} nt · A ${strand.counts.A} · C ${strand.counts.C} · G ${strand.counts.G} · U ${strand.counts.U} · GC ${strand.gc_percent}%</p><p>Complementary DNA: genomic ${strand.complementary_genomic_strand?`<strong>${escapeHTML(strand.complementary_genomic_strand)}</strong> strand`:'strand unknown (no genomic orientation assigned)'}</p><div class="strand-pairing"><code>RNA 5′ ${escapeHTML(strand.rna_5to3)} 3′</code><code>DNA 3′ ${escapeHTML(strand.complementary_dna_aligned_3to5)} 5′</code></div><p class="hint">The DNA strand is shown antiparallel. Its 5′ → 3′ sequence is <code>${escapeHTML(strand.complementary_dna_5to3)}</code>.</p></article>`).join('');
  const roleSequence=role=>strands.find(strand=>strand.name===role)?.rna_5to3 || '';
  return `<section class="strand-analysis"><h4>Strand composition &amp; complementarity</h4><p class="hint">RNA cores only; overhangs are excluded. ${record.partner_inferred?'The complementary partner is inferred.':'Both cores were supplied.'}</p><div class="strand-analysis-grid">${nucleotideCards}</div><div class="table-scroll"><table><thead><tr><th>Model</th><th>Predicted guide · 5′ → 3′</th><th>Predicted passenger · 5′ → 3′</th><th>Guide complementary DNA strand</th></tr></thead><tbody>${analysis.models.map(model=>{
    const guide=strands.find(strand=>strand.name===model.predicted_guide);
    return `<tr><td title="${escapeHTML(model.explanation)}">${escapeHTML(model.table_name)}</td><td>${escapeHTML(model.predicted_guide)}<br><code>${escapeHTML(roleSequence(model.predicted_guide))}</code></td><td>${escapeHTML(model.predicted_passenger || 'undetermined')}<br><code>${escapeHTML(roleSequence(model.predicted_passenger))}</code></td><td>${guide?.complementary_genomic_strand?`Genomic ${escapeHTML(guide.complementary_genomic_strand)}`:'Undetermined'}</td></tr>`;
  }).join('')}</tbody></table></div><p class="hint">${escapeHTML(analysis.scope)} A genomic +/− label describes reference orientation, independently of guide/passenger loading. Promoter DNA association and promoter-associated RNA targeting have experimental support in different settings; this app does not determine which substrate this duplex binds. <a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC4797303/" target="_blank" rel="noopener">DNA-target evidence</a> · <a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC4783471/" target="_blank" rel="noopener">saRNA/Ago2 promoter association</a>.</p></section>`;
}
function sequenceAppearance(record) {
  if (record.reference && record.name === 'SCN1A-7') return {
    color:'#b43f63', pale:'#fff0f4', label:'↓ Downregulating', kind:'down',
    description:'Researcher-reported downregulation in HeLa and AF22 (qPCR); AF22 NES also RNA-seq.'
  };
  if (record.reference) return {color:'#7954a1',pale:'#f4effb',label:'↑ Activating',kind:'positive',
    description:'Activation reported in the cited experimental setting.'};
  return {color:'#167a79',pale:'#edf8f6',label:'Imported',kind:'imported',
    description:'Imported sequence; experimental effect not assigned.'};
}
function sequenceEffectBadge(record) {
  const style=sequenceAppearance(record);
  return `<span class="sequence-effect ${style.kind}" title="${escapeHTML(style.description)}">${style.label}</span>`;
}
function existingSaRNATrack(records,start,end,left,width,focused) {
  const scale=offset=>left+(offset-start)/(end-start)*width;
  const items=[];
  records.forEach((record,index)=>record.matches.forEach((match,j)=>{
    if (match.end<=start || match.offset>=end) return;
    const style=sequenceAppearance(record), text=`${record.name}${match.candidate?` = ${match.candidate.id}`:""} · ${style.label}`;
    const markerLeft=scale(Math.max(start,match.offset)),markerRight=scale(Math.min(end,match.end));
    const pillWidth=Math.min(width,Math.max(150,text.length*7+24));
    const x=Math.max(left,Math.min(left+width-pillWidth,(markerLeft+markerRight-pillWidth)/2));
    items.push({record,index,j,style,text,markerLeft,markerRight,pillWidth,x,
      lo:Math.min(x,markerLeft),hi:Math.max(x+pillWidth,markerRight)});
  }));
  items.sort((a,b)=>a.lo-b.lo || a.index-b.index || a.j-b.j);
  const lanes=[];
  for (const item of items) {
    let lane=lanes.findIndex(right=>right+14<=item.lo);
    if (lane<0) lane=lanes.length;
    lanes[lane]=item.hi;item.lane=lane;
  }
  const laneCount=Math.max(1,lanes.length), bottom=326+laneCount*48;
  let markup=`<rect class="map-decoration" aria-hidden="true" x="${left}" y="316" width="${width}" height="${bottom-316}" rx="8" fill="#f7f9fc" stroke="#e4e9f0"/>`;
  if (!items.length) markup+=`<text x="${left+14}" y="344" font-size="12" fill="#61767a">No existing saRNA sites in this view. See their locations in the table above.</text>`;
  for (const item of items) {
    const {record,index,j,style,text,markerLeft,markerRight,pillWidth,x}=item;
    const y=326+item.lane*48, center=(markerLeft+markerRight)/2;
    const selected=focused?.record===index && (focused.match || 0)===j;
    const help=`${record.name}: ${style.description} ${record.source}`;
    markup+=`<g class="existing-sarna-mark" role="button" tabindex="0" aria-label="${escapeHTML(help)}" data-tip="${escapeHTML(help)}" data-sequence="${index}" data-match="${j}"><rect x="${markerLeft}" y="${y}" width="${Math.max(3,markerRight-markerLeft)}" height="7" rx="2" fill="${style.color}"/><line x1="${center}" x2="${center}" y1="${y+7}" y2="${y+13}" stroke="${style.color}"/><rect class="map-decoration" aria-hidden="true" x="${x}" y="${y+13}" width="${pillWidth}" height="25" rx="7" fill="${style.pale}" stroke="${style.color}" stroke-width="${selected?2:1}"/><text x="${x+12}" y="${y+30}" font-size="12" font-weight="600" fill="${style.color}">${escapeHTML(text)}</text></g>`;
  }
  markup+=`<text x="${left+2}" y="${bottom+20}" font-size="11" fill="#7954a1">↑ Activating reference</text><text x="${left+160}" y="${bottom+20}" font-size="11" fill="#b43f63">↓ Downregulating · researcher reported</text><text x="${left+425}" y="${bottom+20}" font-size="11" fill="#167a79">Imported · effect unassigned</text>`;
  return {markup,height:bottom+34};
}
function renderMap() {
  if (!state.result) return;
  const svg = $("promoter-map");
  // One SVG unit is one CSS pixel: widening the workspace adds plot space, not larger type.
  const mapWidth = Math.max(900, Math.round(svg.parentElement.clientWidth - 24));
  MAP_PLOT = mapWidth - MAP_LEFT - 26;
  svg.setAttribute("viewBox", `0 0 ${mapWidth} 370`);
  svg.style.width = `${mapWidth}px`;
  const r = state.result, p = r.promoter, left = MAP_LEFT, right = left + MAP_PLOT;
  const width = MAP_PLOT;
  const scale = offset => left+(offset-state.start)/(state.end-state.start)*width;
  const inView = (start,end) => end>state.start && start<state.end;
  const tip = text => `data-tip="${escapeHTML(text)}" tabindex="0" aria-label="${escapeHTML(text)}"`;
  const parts = [`<defs><clipPath id="map-clip"><rect x="${left}" y="25" width="${width}" height="260"/></clipPath></defs>`];
  // Subtle lanes and a separate label gutter keep annotation colors as the data channel.
  parts.push(`<rect class="map-decoration" aria-hidden="true" x="0" y="29" width="${mapWidth}" height="256" rx="8" fill="#fbfdfc"/>`);
  for (const [top,bottom] of [[29,72],[116,153],[187,222],[254,285]]) {
    parts.push(`<rect class="map-decoration" aria-hidden="true" x="${left}" y="${top}" width="${width}" height="${bottom-top}" fill="#f2f7f5"/>`);
  }
  parts.push(`<rect class="map-decoration" aria-hidden="true" x="0" y="29" width="${left-12}" height="256" rx="8" fill="#f5f8f7"/>`);
  parts.push(`<text class="map-axis-caption" x="${left}" y="301">${p.tss_offset===null?"Position in supplied sequence (bp)":"Position relative to chosen reference (bp)"}</text>`);
  const label = (text,y,help) => parts.push(`<text class="track-label" x="12" y="${y}" `
    + `${tip(help)}>${escapeHTML(text)}</text>`);
  const emptyNote = (y, text) => parts.push(
    `<text class="track-empty" x="${left + 6}" y="${y}">${escapeHTML(text)}</text>`);
  const rect = (start,end,y,height,color,help,attrs="") => {
    if (!inView(start,end)) return;
    const x1=scale(Math.max(start,state.start)), x2=scale(Math.min(end,state.end));
    // A colour may be a palette entry carrying its edge, or a bare fill for the
    // rejection colours, which are generated rather than declared.
    const fill = typeof color === "string" ? color : color.fill;
    const edge = typeof color === "string" ? "" : ` stroke="${color.edge}" stroke-width=".75"`;
    parts.push(`<rect x="${x1}" y="${y}" width="${Math.max(2,x2-x1)}" height="${height}" rx="2" fill="${fill}"${edge} ${tip(help)} ${attrs}/>`);
  };
  for (const track of TRACKS) label(track.name, track.y, track.help);
  for (const y of [72.5,116.5,153.5,187.5,222.5,254.5,285.5]) parts.push(`<line x1="${left}" x2="${right}" y1="${y}" y2="${y}" stroke="#e8efeb"/>`);
  parts.push(`<g id="map-ticks">${tickMarkup(state.start,state.end,p)}</g>`);
  // Two nested groups on purpose: the clip stays in the untransformed coordinate system while
  // the inner group is translated during a drag. Putting both on one element leaves the clip's
  // reference space ambiguous, and browsers have disagreed about it.
  parts.push('<g clip-path="url(#map-clip)"><g id="map-pan">');
  if (p.tss_offset!==null && p.tss_offset>=state.start && p.tss_offset<state.end) {
    const x=scale(p.tss_offset);
    parts.push(`<line x1="${x}" x2="${x}" y1="29" y2="281" stroke="#52726b" stroke-dasharray="4 4" opacity=".5"/><path d="M${x-5} 37h10l-5 7z" fill="#344d55" ${tip("+1 · chosen reference. The first transcribed base in this coordinate system. Its source and uncertainty are recorded in run evidence.")}/><text x="${x+7}" y="46" fill="#344d55" font-size="11">+1 reference</text>`);
  }
  for (const t of r.tracks.alternative_tss) rect(t.start,t.end,53,12,TRACK_COLOURS.altTss,t.detail);
  if (p.tss_offset!==null) for (const id of ["tss-window","tss-preferred-window"]) {
    const rule=r.configuration.rules.find(q=>q.id===id); if (!rule) continue;
    const values=Object.fromEntries(rule.parameters.map(q=>[q.name,q.value]));
    const a=values[id==="tss-window"?"window-start-bp":"preferred-window-start-bp"];
    const b=values[id==="tss-window"?"window-end-bp":"preferred-window-end-bp"];
    rect(p.tss_offset+a,p.tss_offset+b+1,id==="tss-window"?84:99,11,id==="tss-window"?TRACK_COLOURS.window:TRACK_COLOURS.preferred,`${roleName[rule.kind]}: ${rule.name}. ${signed(a)} to ${signed(b)}. Whole target containment is required. ${rule.provenance.display}`);
  }
  for (const t of r.tracks.cpg_islands) rect(t.start,t.end,134,14,TRACK_COLOURS.cpg,t.detail);
  if (!r.tracks.cpg_islands.length) emptyNote(145,
    "No CpG island met the configured length, GC and observed/expected thresholds in this window.");
  for (const t of r.tracks.regulatory) rect(t.start,t.end,165,14,
    t.kind.includes("promoter")?TRACK_COLOURS.promoter:TRACK_COLOURS.regulatory,`${t.label}: ${t.detail}`);
  if (!r.tracks.regulatory.length) emptyNote(176, r.app.request.regulatory === false
    ? "Not retrieved. Ensembl regulatory annotations were switched off for this run."
    : "No Ensembl Regulatory Build feature overlaps this window.");
  const focused = state.motifFocus;
  const motifColours = motifPalette(r);
  r.tracks.motifs.forEach((t,i)=>{
    // Exact label match: a prefix test would make SP1 claim SP100 if one were bundled.
    const mine = !focused || t.label === focused;
    const own = motifColours.get(t.label);
    rect(t.start,t.end,197+(i%3)*5,mine?6:3,
      mine ? (own ? {fill:own, edge:own} : TRACK_COLOURS.motif) : TRACK_COLOURS.motifDim,
      t.detail,`data-motif-mark="${escapeHTML(t.label)}"`);
  });
  if (!r.tracks.motifs.length) emptyNote(212, r.app.request.motifs === false
    ? "Not scanned. Motif prediction was switched off for this run."
    : "No bundled matrix matched above its threshold in this window.");
  const rejectionColours = rejectionPalette(r);
  state.excludedBinned = drawExcluded(parts, r, rejectionColours, scale, tip, rect);
  drawPassing(parts, r, scale, tip, rect);
  if (state.selected) {
    const c=r.candidates.find(c=>c.id===state.selected);
    rect(c.offset,c.offset+c.target_site_dna.length,c.rank===null?261:228,c.rank===null?19:25,"none",`${c.id} selected`,`stroke="#b97a2b" stroke-width="2" data-candidate="${c.id}"`);
  }
  for (const candidate of r.candidates) {
    const controls=candidateControlMatches(candidate,r);
    if (!controls.length) continue;
    const style=sequenceAppearance(controls[0]);
    const y=candidate.rank===null?260:228;
    rect(candidate.offset,candidate.offset+candidate.target_site_dna.length,y,20,
      'none',`${candidate.id}: same duplex core as ${controls.map(c=>c.name+' · '+sequenceAppearance(c).label).join('; ')}`,
      `stroke="${style.color}" stroke-width="2.5" data-candidate="${candidate.id}"`);
  }
  parts.push('</g></g>');
  const existingTrack=existingSaRNATrack(state.sequenceRecords || [],state.start,state.end,left,width,sequenceFocus);
  parts.push(existingTrack.markup);
  svg.setAttribute("viewBox", `0 0 ${mapWidth} ${existingTrack.height}`);
  svg.style.height=`${existingTrack.height}px`;
  svg.innerHTML=parts.join("");
  // Only offer the grab cursor when there is somewhere to pan. At full extent the drag is a
  // no-op, and advertising it is a promise the map does not keep.
  svg.classList.toggle("pannable", state.start > 0 || state.end < p.length);
  renderRejectionLegend(r);
  renderLegend(r);
  const note = $("unevaluated-note");
  if (note) {
    const missing = r.unevaluated_rules || [];
    note.hidden = !missing.length;
    note.textContent = missing.length
      ? `${missing.length} rule${missing.length === 1 ? "" : "s"} could not be evaluated for any `
        + `candidate in this run: ${missing.map(ruleName).join(", ")}. `
        + `These counted neither for nor against anything, and they are not failures.`
      : "";
  }
}
const MAP_LEFT = 160, MAP_MIN_SPAN = 40;
let MAP_PLOT = 814;

function tickMarkup(start,end,p) {
  // Coordinate labels only: Python already supplies scientific candidate spans.
  let out = "";
  for (let i=0;i<=6;i++) {
    const offset=Math.round(start+(end-start-1)*i/6);
    const x=MAP_LEFT+(offset-start)/(end-start)*MAP_PLOT;
    const continuous = p.tss_offset===null?null:offset-p.tss_offset;
    const text = continuous===null?String(offset+1):signed(continuous>=0?continuous+1:continuous);
    out += `<text class="tick" text-anchor="${i===0?"start":i===6?"end":"middle"}" x="${x}" y="20">${text}</text>`;
  }
  return out;
}

function setViewport(start,span) {
  // One clamp, used by the buttons, the wheel and the drag, so the three cannot disagree.
  const max = state.result.promoter.length;
  const width = Math.min(max,Math.max(MAP_MIN_SPAN,span));
  state.start = Math.max(0,Math.min(max-width,start));
  state.end = state.start+width;
}

function enableMapDrag() {
  // A full redraw costs ~14.5 ms on this promoter (2,600 SVG nodes) against a 16.7 ms frame
  // budget, so redrawing on every pointermove would drop frames on any slower machine.
  // Panning is a pure horizontal translation of the plotted area, so during the drag the
  // clipped group is translated -- exact, not an approximation -- and the full redraw happens
  // once on release. Only the seven tick labels are refreshed live.
  const svg = $("promoter-map");
  let drag = null;
  const perPixel = () => svg.viewBox.baseVal.width / (svg.getBoundingClientRect().width || 1);
  new ResizeObserver(() => {
    if (state.result && Math.max(900, Math.round(svg.parentElement.clientWidth - 24)) !== svg.viewBox.baseVal.width) renderMap();
  }).observe(svg.parentElement);

  const move = event => {
    if (!drag) return;
    const dx = (event.clientX - drag.x) * drag.unit;
    drag.moved = Math.max(drag.moved, Math.abs(event.clientX - drag.x));
    // Dragging right pulls the sequence right, which reveals earlier offsets.
    const wanted = drag.start - dx * drag.span / MAP_PLOT;
    const start = Math.max(0, Math.min(drag.max - drag.span, wanted));
    drag.pending = start;
    const pan = svg.querySelector("#map-pan");
    if (pan) pan.setAttribute("transform",
      `translate(${((drag.start - start) * MAP_PLOT / drag.span).toFixed(2)},0)`);
    const ticks = svg.querySelector("#map-ticks");
    if (ticks) ticks.innerHTML = tickMarkup(start, start + drag.span, state.result.promoter);
  };
  const end = () => {
    if (!drag) return;
    const {pending, span, moved} = drag;
    drag = null;
    window.removeEventListener("pointermove", move);
    window.removeEventListener("pointerup", end);
    window.removeEventListener("pointercancel", end);
    svg.classList.remove("dragging");
    // A drag that ends over a candidate must not also select it; a click without movement must.
    state.suppressMapClick = moved > 3;
    setViewport(pending, span);
    renderMap();
  };

  svg.addEventListener("pointerdown", event => {
    if (!state.result || event.button !== 0) return;
    if (state.end - state.start >= state.result.promoter.length && state.start === 0) {
      // Nothing to pan to at full extent; leave clicks and focus alone.
      return;
    }
    event.preventDefault();
    drag = {x: event.clientX, start: state.start, span: state.end - state.start,
            max: state.result.promoter.length, unit: perPixel(), moved: 0,
            pending: state.start};
    svg.classList.add("dragging");
    window.addEventListener("pointermove", move);
    window.addEventListener("pointerup", end);
    window.addEventListener("pointercancel", end);
  });

  svg.addEventListener("wheel", event => {
    if (!state.result) return;
    event.preventDefault();
    const rect = svg.getBoundingClientRect();
    // Zoom about the cursor, so the base under the pointer stays under the pointer.
    const atX = (event.clientX - rect.left) * perPixel();
    const fraction = Math.min(1, Math.max(0, (atX - MAP_LEFT) / MAP_PLOT));
    const span = state.end - state.start;
    const focus = state.start + fraction * span;
    const wanted = Math.min(state.result.promoter.length,
      Math.max(MAP_MIN_SPAN, span * (event.deltaY > 0 ? 1.25 : 0.8)));
    setViewport(focus - fraction * wanted, wanted);
    renderMap();
  }, {passive: false});
}

function viewport(factor,pan=0) {
  if (!state.result) return;
  const old=state.end-state.start, span=Math.min(state.result.promoter.length,Math.max(MAP_MIN_SPAN,old*factor));
  setViewport((state.start+state.end)/2+old*pan-span/2, span);
  renderMap();
}
function renderEvidence() {
  const r=state.result, c=r.configuration;
  $("run-evidence-content").innerHTML = `<p><b>${escapeHTML(c.preset.name)}</b> · ${c.is_custom?"customized":"preset defaults"} · engine ${escapeHTML(r.app.engine_version)}</p><p>Run ${escapeHTML(r.app.run_id)} · ${escapeHTML(r.app.created)}</p><h3>Promoter source</h3><p>${escapeHTML(r.app.anchor?.label || "User-supplied DNA")}</p>${r.app.anchor?`<p>${escapeHTML(r.app.anchor.provenance?.display || r.app.anchor.provenance || "Source unavailable")}</p><pre style="white-space:pre-wrap;font-size:10px">${escapeHTML(JSON.stringify(r.app.anchor,null,2))}</pre>`:""}<h3>Custom departures</h3><ul>${c.customizations.map(v=>`<li>${escapeHTML(JSON.stringify(v))}</li>`).join("") || "<li>None</li>"}${c.warnings.map(v=>`<li>${escapeHTML(v)}</li>`).join("")}${c.preset.divergences.map(v=>`<li>${escapeHTML(v)}</li>`).join("")}</ul><h3>Source availability</h3><ul>${r.app.notices.map(n=>`<li>${escapeHTML(n)}</li>`).join("")}</ul>${r.regulatory_summary?`<p>${escapeHTML(r.regulatory_summary.provenance)}</p><p>${escapeHTML(r.regulatory_summary.summary)}</p>`:""}<h3>Rule effects</h3><p>Every gate each candidate failed, so failures overlap: a candidate failing four gates appears in four rows and these shares do not add up to 100%. Percentages are of all ${fmt(r.counts.generated)} generated targets. The promoter map's legend measures something else, each target counted once against the first gate that excluded it, and is listed below.</p><ul>${r.rejection_breakdown.map(v=>`<li>${escapeHTML(v.rule)}: ${fmt(v.rejected)} (${v.percent_of_generated}%)</li>`).join("")}</ul><h4>First gate to exclude, the map's attribution</h4><p>Each of the ${fmt(r.counts.rejected)} excluded targets counted once, in this run's rule evaluation order.</p><ul>${(r.rejection_attribution||[]).map(v=>`<li>${escapeHTML(v.rule)}: ${fmt(v.rejected)} (${v.percent_of_excluded}% of excluded)</li>`).join("")}</ul><h3>Frozen configuration</h3>${c.rules.map(def=>`<details class="verdict"><summary>${escapeHTML(def.name)} · ${roleName[def.kind]}</summary><p class="source">${provenanceText(def.provenance)} · ${escapeHTML(evidenceLabel(def))}</p><p>${escapeHTML(def.rationale)}</p>${ruleEvidence(def)}<p>${escapeHTML(def.consequence)}</p><p>${def.parameters.map(p=>`${escapeHTML(p.name)}: ${escapeHTML(p.value)} (${escapeHTML(p.origin)})`).join("; ")}</p></details>`).join("")}<h3>Disabled or omitted rules</h3><p>${state.catalogue.rules.filter(def=>!c.rules.some(v=>v.id===def.id)).map(def=>escapeHTML(def.name)).join(" · ")}</p>`;
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
  const order=result.screening_order || result.order, ids=state.checked;
  if (!order || !order.rows?.length) throw new Error("This record carries no order-ready duplexes.");
  const rows=ids.size?order.rows.filter(r=>ids.has(r[1])):order.rows;
  if (!rows.length) throw new Error("None of the selected candidates passed the gates, so none can be ordered.");
  return [order.columns,...rows].map(row=>row.map(v=>`"${String(v??"").replace(/"/g,'""')}"`).join(",")).join("\r\n");
}
async function exportOrder() {
  try {
    if (typeof validateScreenedExport === "function") validateScreenedExport();
    const r=state.result, name=`${r.gene?.symbol || "promoter"}-${r.app.run_id.slice(0,8)}-order.csv`;
    if (frozen) { downloadBlob(new Blob([orderCsvText(r)],{type:"text/csv"}),name); return; }
    const response=await fetch(`/api/runs/${r.app.run_id}/export.order.csv${state.checked.size?`?ids=${encodeURIComponent([...state.checked].join(","))}`:""}`);
    if (!response.ok) throw new Error(`Export failed (${response.status})`);
    downloadBlob(await response.blob(),name);
  } catch(e) { error(e.message); }
}
async function exportRun(type) {
  try {
    if ((type === "csv" || type === "xlsx") && typeof validateScreenedExport === "function") validateScreenedExport();
    if (typeof flushWorkspace === "function") await flushWorkspace(state.result.app.run_id);
    const r=state.result, name=`${r.gene?.symbol || "promoter"}-${r.app.run_id.slice(0,8)}.${type}`;
    if (frozen) { downloadBlob(new Blob([type==="csv"?csvText(r):JSON.stringify({...r,screenings:frozenScreeningRecords(),screening_summary:frozenScreeningSummary()},null,2)],{type:type==="csv"?"text/csv":"application/json"}),name); return; }
    // A selection applies to the row-wise exports; an HTML or JSON archive is the whole run.
    const selectable = type === "csv" || type === "xlsx";
    const ids = selectable && state.checked.size
      ? `?ids=${encodeURIComponent([...state.checked].join(","))}` : "";
    const response=await fetch(`/api/runs/${r.app.run_id}/export.${type}${ids}`);
    if (!response.ok) throw new Error(`Export failed (${response.status})`);
    downloadBlob(await response.blob(),name);
  } catch(e) { error(e.message); }
}
function renderReproducibility(result) {
  const box = $("reproducibility-panel"), record = result.reproducibility;
  if (!box) return;
  if (!record) {
    box.innerHTML = '<summary>Reproducibility metadata not recorded</summary><p class="hint">This historical run retains its saved results. Its exact build and reference fingerprints were not recorded and cannot be inferred from today’s installation.</p>';
    return;
  }
  const software = record.software || {};
  const revision = software.git_revision || "Git revision unknown";
  const dirty = software.git_dirty === null ? "working-tree status unknown" : software.git_dirty ? "uncommitted changes at capture" : "clean source checkout at capture";
  const sources = (record.sources || []).map(source => `<li>${escapeHTML(source.source)} · ${source.cache_hit?"cached response":"live response"} · observed ${escapeHTML(source.observed_at)} · original retrieval ${escapeHTML(source.retrieved_at || "unknown")} · annotation release ${escapeHTML(source.annotation_release || "unknown")}<br>${provenanceText({display:source.url,url:source.url})}<br>Response SHA-256: <code>${escapeHTML(source.response_sha256)}</code></li>`).join("");
  box.innerHTML = `<summary>Reproducibility: saved inputs, software &amp; sources</summary><p class="hint">${escapeHTML(record.replay_scope)} Export the whole-run JSON and run <code>sarna-designer-replay run.json</code> offline. The report compares outputs, software and environment separately.</p><p class="hint">Revision: <code>${escapeHTML(revision)}</code> · ${escapeHTML(dirty)}<br>Package SHA-256: <code>${escapeHTML(software.package_sha256)}</code><br>Promoter sequence SHA-256: <code>${escapeHTML(record.sequence_sha256)}</code><br>Python ${escapeHTML(record.environment?.python)} · ${Object.keys(record.environment?.packages || {}).length} installed package versions recorded.</p><ul class="hint">${sources || "<li>No external reference responses were recorded for this run.</li>"}</ul><ul class="hint">${record.limitations.map(note=>`<li>${escapeHTML(note)}</li>`).join("")}</ul>`;
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
  box.innerHTML = parts.slice(0,2).join("") + `<details><summary>Measurement &amp; limitations</summary>${parts.slice(2).join("")}</details>`;
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
  return `${t.ensembl_id}${marks.length ? `: ${marks.join(", ")}` : ""}`;
}

function orderRowFor(result, candidateId) {
  // The order-ready duplex Python computed for this candidate. Returns null for a candidate
  // that failed a gate: those are not in the order table, by design.
  const order = result.order;
  if (!order || !order.rows) return null;
  const row = order.rows.find(r => r[1] === candidateId);
  return row ? Object.fromEntries(order.columns.map((c, i) => [c, row[i]])) : null;
}

const RNA_PAIR = {A:"U", U:"A", G:"C", C:"G"};
const STRAND_ROLE = {
  antisense:{antisense:"guide", sense:"passenger"},
  sense:{antisense:"passenger", sense:"guide"},
};
const ROLE_NOTE = {
  guide:"Loaded into RISC. This is the strand that engages the target.",
  passenger:"Discarded after loading. It is still synthesised and still annealed.",
  undecided:"No strand is favoured: the two 5′ ends are equally stable to within the model’s tolerance.",
};

function strandRoles(guide) {
  // "ambiguous" and a missing prediction are the same display case but not the same statement;
  // the caption below the diagram distinguishes them.
  return STRAND_ROLE[guide] || {antisense:"undecided", sense:"undecided"};
}

function duplexLadder(order, candidate) {
  // Draws the molecule as it is actually annealed: antiparallel, with the 2 nt 3' dTdT
  // overhang on *each* strand, and a rung only where the bases genuinely pair. Three
  // sequences listed 5'->3' in separate boxes cannot show any of that, and the reader was
  // left to reverse-complement 19 bases in their head to check the design at all.
  const sense = order.sense_core_rna_5to3, anti = order.antisense_core_rna_5to3;
  const n = sense.length;
  if (!n || anti.length !== n) return "";
  // Antisense written 3'->5' so it lies antiparallel under the sense strand.
  const antiReversed = [...anti].reverse();
  const paired = [...sense].map((b, i) => RNA_PAIR[b] === antiReversed[i]);
  const roles = strandRoles(candidate.rating.guide);
  const endLen = Math.min(candidate.guide_prediction.end_length_nt || 4, n);
  const cell = (text, cls="") => `<span class="${cls}">${escapeHTML(text)}</span>`;
  // The two ends the prediction compares, marked on the drawing rather than described beside it.
  const endClass = i => i < endLen ? " end-five" : i >= n - endLen ? " end-three" : "";
  const senseCells = [...sense].map((b, i) => cell(b, "b" + endClass(i)));
  const antiCells = antiReversed.map((b, i) => cell(b, "b" + endClass(i)));
  const rungs = paired.map((ok, i) =>
    cell(ok ? "│" : "×", (ok ? "rung" : "rung mismatch") + endClass(i)));
  const dT = () => order.core_only ? [] : [cell("dT", "ov"), cell("dT", "ov")];
  const blank = () => order.core_only ? [] : [cell(" ", "pad"), cell(" ", "pad")];
  const row = (label, role, terminus5, cells, terminus3) =>
    `<div class="ladder-row ${role}">`
    + `<span class="ladder-label">${escapeHTML(label)}<i>${escapeHTML(role)}</i></span>`
    + `<span class="term">${terminus5}</span><span class="bases">${cells.join("")}</span>`
    + `<span class="term">${terminus3}</span></div>`;
  return `<div class="ladder" role="img" aria-label="${escapeHTML(
      `Duplex: sense ${sense} paired with antisense ${anti}${order.core_only ? ", paired RNA cores; overhangs omitted" : ", each with a 3 prime dTdT overhang"}`)}">`
    // Sense reads 5'->3' left to right; its dTdT therefore sits on the right.
    + row("Sense", roles.sense, "5′", [...blank(), ...senseCells, ...dT()], "3′")
    + `<div class="ladder-row rungs"><span class="ladder-label"></span><span class="term"></span>`
    + `<span class="bases">${[...blank(), ...rungs, ...blank()].join("")}</span>`
    + `<span class="term"></span></div>`
    // Antisense reads 3'->5' left to right, so its dTdT sits on the left and its 5' end,
    // the one the loading rule is about, is at the right-hand end of the drawing.
    + row("Antisense", roles.antisense, "3′", [...dT(), ...antiCells, ...blank()], "5′")
    + `</div>`;
}

function endEvidence(candidate) {
  // The measurement behind the guide call, as numbers. Until the ends were serialised this
  // existed only inside one English sentence, so the interface could state the verdict and
  // not show the subtraction it rests on (ADR 0009).
  const pred = candidate.guide_prediction;
  const model = candidate.rating.guide_model;
  const a = pred.assessments.find(x => x.table === model) || pred.assessments[0];
  if (!a || !a.five_prime) return "";
  const box = (side, end, cls) =>
    `<div class="end-box ${cls}"><span>Target ${escapeHTML(end.end)} end · ${escapeHTML(side)}</span>`
    + `<b>${end.energy.toFixed(2)} ${escapeHTML(end.unit)}</b>`
    + `<small>${end.stacks.map(([s, e]) => `${escapeHTML(s)} ${e.toFixed(2)}`).join(" + ")}</small></div>`;
  const weaker = a.three_prime.energy > a.five_prime.energy ? "3′" : "5′";
  return `<div class="end-evidence">`
    + box("left of the diagram", a.five_prime, "left")
    + `<div class="end-verdict"><b>${a.difference >= 0 ? "+" : ""}${a.difference.toFixed(2)} ${escapeHTML(a.unit)}</b>`
    + `<small>E(3′) − E(5′)</small></div>`
    + box("right of the diagram", a.three_prime, "right")
    + `<p class="hint">Less negative is less stable. The ${weaker} end is the weaker of the two. `
    + `RISC loads whichever strand has the more weakly paired 5′ end (Schwarz et al. 2003; `
    + `Khvorova et al. 2003), and the antisense strand’s 5′ end pairs with the target’s `
    + `3′ end, the right-hand end above. ${escapeHTML(
        a.predicted_guide === "ambiguous"
          ? "Here the two ends are equal to within the model’s tolerance, so no strand is favoured."
          : `Hence ${a.predicted_guide} is predicted to be loaded.`)}</p>`
    + `</div>`;
}

function duplexPanel(result, candidate) {
  // One panel answering the two questions a researcher has at the bench: what do I order,
  // and which tube is the guide.
  const order = orderRowFor(result, candidate.id);
  if (!order) {
    return `<div class="duplex-panel"><p class="hint">This candidate failed a gate, so it has no `
      + `order-ready duplex. The rule verdicts below say which gate and why.</p></div>`;
  }
  const roles = strandRoles(candidate.rating.guide);
  const tube = (name, role, sequence) =>
    `<div class="tube ${role}"><div class="tube-head"><span class="role-chip ${role}">`
    + `${escapeHTML(role === "undecided" ? "not predicted" : role)}</span>`
    + `<span>${escapeHTML(name)} · 5′ → 3′ · ${order.ordered_length_nt} nt</span>`
    + `<button type="button" data-copy="${escapeHTML(sequence)}">Copy</button></div>`
    + `<code>${escapeHTML(sequence)}</code>`
    + `<small>${escapeHTML(ROLE_NOTE[role])}</small></div>`;
  const warning = candidate.rating.guide === "sense"
    ? `<p class="duplex-warning">The <b>sense</b> strand is predicted to be loaded. The duplex is `
      + `designed for the antisense strand to engage the promoter, so this prediction runs against `
      + `the design intent. The candidate is not withheld; the prediction is reported.</p>`
    : "";
  return `<div class="duplex-panel">`
    + `<div class="duplex-head"><h4>The duplex you would order</h4>`
    + `<p class="hint">Two ${order.ordered_length_nt}-mers, annealed 1:1. Each is a `
    + `${order.core_length_nt} nt RNA core plus a 3′ ${escapeHTML(order.overhang)} overhang, `
    + `so the duplex carries a 2 nt 3′ overhang at both ends. Order as RNA with a 3′ dTdT; `
    + `the overhang is deoxythymidine, not a ribo-UU.</p></div>`
    + duplexLadder(order, candidate)
    + endEvidence(candidate)
    + warning
    + `<div class="tubes">`
    + tube("Sense", roles.sense, order.sense_order_5to3)
    + tube("Antisense", roles.antisense, order.antisense_order_5to3)
    + `</div></div>`;
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
  const controlBadges = candidateControlBadges(best,result);
  const met = best.rating.evaluated ? `${best.rating.met}/${best.rating.evaluated}` : "not evaluated";
  const roles = strandRoles(best.rating.guide);
  const line = (name, role, sequence) =>
    `<div><div><span class="role-chip ${role}">${escapeHTML(
        role === "undecided" ? "not predicted" : role)}</span>`
    + `<span>${escapeHTML(name)} · 5′ → 3′</span><code>${escapeHTML(sequence)}</code></div>`
    + `<button type="button" data-copy="${escapeHTML(sequence)}">Copy</button></div>`;
  const duplex = order
    ? `${controlBadges}<div class="top-duplex">`
      + line("Sense", roles.sense, order.sense_order_5to3)
      + line("Antisense", roles.antisense, order.antisense_order_5to3)
      + `<p class="hint">Two ${order.ordered_length_nt}-mers: a ${order.core_length_nt} nt RNA core `
      + `plus a 3′ dTdT overhang on each strand, annealed 1:1. `
      + `${roles.antisense === "guide"
          ? "The antisense strand is predicted to be loaded, which is the design intent."
          : roles.sense === "guide"
          ? "The sense strand is predicted to be loaded, against the design intent. See the candidate below."
          : "Neither strand is favoured for loading by the selected model."} `
      + `${escapeHTML(screeningSentence())} `
      + `<button type="button" id="inspect-top-duplex" class="text-button">See how they pair ›</button></p></div>`
    : `<p class="hint">Order-ready sequences are not in this record.</p>`;
  box.hidden = false;
  box.innerHTML = `<div class="card-heading"><div><p class="eyebrow">START HERE</p>`
    + `<h3>${fmt(result.counts.passing)} of ${fmt(result.counts.generated)} targets passed every evaluable gate</h3>`
    + `<p class="hint">Highest-rated is <b>${escapeHTML(best.id)}</b> at ${escapeHTML(positionText(best))} `
    + `${best.tss_span ? "relative to the reference TSS" : "in the supplied sequence"}, meeting `
    + `${met} evaluated preferences, with ${escapeHTML(best.rating.guide || "an unknown")} predicted as the guide strand.</p>`
    + `</div><span class="tag">RANK 1</span></div>`
    + duplex
    + `<div id="top-offtarget-summary">${typeof topOffTargetSummary === "function" ? topOffTargetSummary() : ""}</div>`
    + `<div class="top-actions-row">`
    + `<button type="button" id="print-duplex" class="primary">Print duplex sheet</button>`
    + `<button type="button" id="inspect-top" class="secondary">Why this one?</button>`
    + `<button type="button" id="export-order-top" class="secondary">Review in Export →</button>`
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
  // Which tube is which, on the rows that name the tubes. "Predicted guide strand: antisense"
  // further down is the same fact, three paragraphs from where somebody is labelling them.
  const sheetRoles = strandRoles(candidate.rating.guide);
  // Spelled out rather than upper-cased in CSS: this page is printed, and a role that depends
  // on a stylesheet is a role that can go missing on paper.
  const sheetLabel = role =>
    role === "undecided" ? "STRAND NOT PREDICTED" : role.toUpperCase();
  // A screening is a separate dated record, so the sheet can only speak for what this page
  // has. Saying "none performed" after one has run is the error the order CSV just stopped
  // making, and it would be worse here: this is the page that goes to the bench.
  const screened = Object.entries(state.screenings || {})
    .map(([scope, held]) => {
      // Each candidate carries the date of the search that covered it, because batches mean
      // one scope can hold several searches made on different days.
      const mine = held.byCandidate?.[candidate.id];
      if (!mine) return null;
      return `${scope} (${mine.search?.parameters?.database || "?"}, `
        + `${String(mine.search?.searched_at || "").slice(0, 10)}): `
        + `${mine.row.perfect_hits} full-length exact, `
        + `${mine.row.wang_homology_hits} at Wang 2015's 80%, `
        + `of ${mine.row.hit_count} alignments`;
    }).filter(Boolean);
  return `<h1>saRNA duplex ${escapeHTML(candidate.id)}</h1>
    <table class="sheet-table">
      ${row("Gene", escapeHTML(gene))}
      ${row("Promoter reference", escapeHTML(referenceLabel(result)))}
      ${row("Position", escapeHTML(positionText(candidate)) + (candidate.tss_span ? " relative to reference TSS" : " in supplied sequence"))}
      ${row("Genomic", candidate.genomic ? escapeHTML(`${candidate.genomic.contig}:${candidate.genomic.start}-${candidate.genomic.end} (${candidate.genomic.strand} strand, GRCh38)`) : "not placed on the genome")}
      ${row("Target site DNA", `<code>${escapeHTML(candidate.target_site_dna)}</code>`)}
    </table>
    <h2>Order these two oligonucleotides</h2>
    ${order ? `<table class="sheet-table">
      ${row(`Sense strand <i>${escapeHTML(sheetLabel(sheetRoles.sense))}</i>`, `<code>${escapeHTML(order.sense_order_5to3)}</code>`)}
      ${row(`Antisense strand <i>${escapeHTML(sheetLabel(sheetRoles.antisense))}</i>`, `<code>${escapeHTML(order.antisense_order_5to3)}</code>`)}
      ${row("Chemistry", `${order.core_length_nt} nt RNA core + 3′ dTdT overhang on both strands; ${order.ordered_length_nt} nt each. Anneal 1:1.`)}
    </table>` : "<p>This candidate did not pass the gates and has no order-ready duplex.</p>"}
    <h2>Why it was selected</h2>
    <p>Rank ${candidate.rank ?? "n/a"} of ${fmt(result.counts.passing)} that passed every evaluable gate,
       meeting ${escapeHTML(met)}. Predicted guide strand:
       ${escapeHTML(candidate.rating.guide || "unknown")}, by ${escapeHTML(candidate.rating.guide_model || "an unnamed model")}.</p>
    <p><b>Gates passed (${gatesPassed.length} of ${gates.length}):</b> ${gatesPassed.length ? gatesPassed.map(named).join("; ") : "none"}. A gate excludes a candidate that fails it.</p>
    <p><b>Preferences met (${preferencesMet.length}):</b> ${preferencesMet.length ? preferencesMet.map(named).join("; ") : "none"}. A preference adds one point and never excludes.</p>
    ${unavailable.length ? `<p><b>Could not be evaluated (${unavailable.length}):</b> ${unavailable.map(named).join("; ")}. These counted neither for nor against it, and are not failures.</p>` : ""}
    <h2>What this does not establish</h2>
    <ul>
      <li>${measured ? escapeHTML(measured.headline) : "Ranking performance was not measured for this run."}</li>
      <li>${screened.length
        ? `Off-target screening: ${screened.map(escapeHTML).join("; ")}. `
          + `Sequence similarity is not measured binding, and one search of one database on `
          + `one day does not establish specificity.`
        : "No off-target screening has been performed. Screen these sequences independently "
          + "before ordering."}</li>
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




function foldText(scan) {
  // A matrix with no sites has no fold to report: 0 found against 0 expected is undefined, not
  // zero, and printing "null" or "0x" would each claim more than the scan knows.
  if (!scan.sites && !scan.expected) return "no sites found";
  if (scan.enrichment === null || scan.enrichment === undefined) return "n/a";
  return `${scan.enrichment}\u00d7${scan.exceeds_chance ? "" : " (within chance)"}`;
}

// Keyed by the engine's own axis strings, which are hyphenated. An underscore here left
// "tf-occupancy" rendering raw beside two properly named rows.
const AXIS_LABEL = {
  expression:"Expression", chromatin:"Chromatin", "tf-occupancy":"Factor occupancy",
};

function axisState(axis) {
  // Three outcomes here too: loaded, not loaded but reachable, and nothing and no route.
  // Rendering the last two alike would turn "nobody has looked" into "there is nothing there".
  //
  // "Available" rather than "Loaded": this describes what exists for the cell line, not what
  // this run read. A run with the chromatin option off still faces an available axis, and
  // reading the two as the same thing would credit the run with a measurement it did not make.
  //
  // Shared by the setup sidebar and the results banner so the two can lay the same facts out
  // differently without being able to word them differently.
  return axis.available ? {state: "loaded", word: "Available"}
    : axis.reachable ? {state: "reachable", word: "Not wired in"}
    : {state: "none", word: "No route"};
}
function cellAxisRow(axis) {
  const {state, word} = axisState(axis);
  return `<tr class="axis-${state}"><th>${escapeHTML(AXIS_LABEL[axis.axis] || axis.axis)}`
    + `<span>${escapeHTML(axis.question || "")}</span></th>`
    + `<td><span class="axis-state ${state}">${word}</span></td>`
    + `<td>${axis.source ? escapeHTML(axis.source) : "no source"}</td>`
    + `<td>${escapeHTML(axis.note || "")}</td></tr>`;
}

function renderCellContext(result, cell) {
  // The registry carries each line's Cellosaurus accession, its citation and what evidence
  // exists on three axes. All of it reached the browser and none of it was drawn: the banner
  // printed the name the researcher typed and nothing else. A cell line named without its
  // accession cannot be checked against anybody else's record, which is the whole point of
  // recording one.
  const chromatin = cell?.chromatin;
  const typed = escapeHTML(result.app.cell_line);
  if (!cell) {
    $("context-banner").innerHTML = `<strong>Cell context \u00b7 ${typed}</strong>`
      + `<p>No record is held for this name, so nothing here is specific to it. A cell name on `
      + `its own supplies no expression, chromatin or occupancy evidence, and the ranking above `
      + `does not use one. Known lines are offered in the cell-line box; your own can be added `
      + `with SARNA_CELL_LINES.</p>`
      + configurationNotes(result);
    return;
  }
  const accession = cell.rrid
    ? `<a class="accession" href="https://www.cellosaurus.org/${encodeURIComponent(cell.rrid)}" `
      + `target="_blank" rel="noopener noreferrer" `
      + `data-tip="Cellosaurus record for this line: its origin, its known problems and its `
      + `cross-references. Opens in a new tab so the identity can be checked independently.">`
      + `${escapeHTML(cell.rrid)} \u2197</a>`
    : `<span class="accession none" data-tip="No Cellosaurus accession is recorded for this `
      + `line, so its identity cannot be checked against an external record.">no accession</span>`;
  const axes = (cell.data_axes || []).map(cellAxisRow).join("");
  $("context-banner").innerHTML =
    `<strong>Cell context \u00b7 ${escapeHTML(cell.name)}</strong> ${accession}`
    + `<p>${chromatin
        ? escapeHTML(chromatin.summary || "Local histone measurements loaded; sources below.")
        : "Promoter activity is unknown in this run. No cell-specific expression or histone "
          + "measurement establishes that this locus is active."}</p>`
    + `<p class="cell-source">${escapeHTML(cell.source || "")}`
    + `${cell.experimental_context
        ? " \u00b7 this project ran experiments in this line"
        : " \u00b7 declared for reference; nothing in this project was measured in this line"}</p>`
    + `<details><summary>What is known about ${escapeHTML(cell.name)}, and what is not</summary>`
    + `<p class="hint">${escapeHTML(cell.description || "")}</p>`
    + (axes ? `<table class="axis-table"><thead><tr><th>Axis</th><th>Status</th><th>Source</th>`
        + `<th>Detail</th></tr></thead><tbody>${axes}</tbody></table>`
        + `<p class="hint">This describes the cell line, not this run: an available axis is one `
        + `the tool can read for these cells, and whether this run read it is in the run's own `
        + `evidence. None of it changes a rank.</p>` : "")
    + `</details>`
    + `<details><summary>Evidence gaps &amp; source notes</summary><ul>`
    + result.app.notices.map(n => `<li>${escapeHTML(n)}</li>`).join("")
    + (result.app.anchor?.caveat ? `<li>${escapeHTML(result.app.anchor.caveat)}</li>` : "")
    + `</ul></details>`
    + configurationNotes(result);
}

function configurationNotes(result) {
  return result.configuration.warnings.length
    ? `<p style="margin-top:8px"><b>Configuration notes:</b> `
      + `${result.configuration.warnings.map(escapeHTML).join(" \u00b7 ")}</p>`
    : "";
}

function renderMotifPanel(result) {
  // The map draws 144 identical marks, which honestly says "motifs are everywhere" and is
  // useless for finding one. The signal the scan already computed -- how many sites a factor
  // has against how many a shuffled background predicts -- was in the record and on no screen.
  const panel = $("motif-panel");
  if (!panel) return;
  const scans = result.motif_scans || [];
  if (!scans.length) { panel.hidden = true; panel.innerHTML = ""; return; }
  const ranked = [...scans].sort((a, b) =>
    (b.exceeds_chance - a.exceeds_chance) || (b.enrichment - a.enrichment) || (b.sites - a.sites));
  const focused = state.motifFocus;
  const colours = motifPalette(result);
  panel.hidden = false;
  panel.innerHTML =
    `<h3>Transcription factor motifs <span class="muted">· ${scans.length} matrices scanned`
    + `</span></h3>`
    + `<p class="hint">Sites found against sites expected from ${result.app.request.motif_shuffles} `
    + `shuffled backgrounds of the same composition. A match is a sequence match: most are not `
    + `bound in any given cell type, and none of this is a measurement in your cell line. `
    + `Select a row to show only that factor on the map.</p>`
    + `<div class="table-scroll"><table><thead><tr><th>Factor</th><th>Matrix</th><th>Sites</th>`
    + `<th>Expected</th><th>Fold</th><th>Assay</th><th>Why this matrix is bundled</th></tr></thead>`
    + `<tbody>` + ranked.map(m => {
        const on = focused === m.motif;
        return `<tr data-motif="${escapeHTML(m.motif)}" class="${on ? "selected" : ""}" `
          + `tabindex="0" aria-label="Show only ${escapeHTML(m.motif)} on the map">`
          + `<td><b>${colours.has(m.motif)
              ? `<i class="motif-swatch" style="background:${colours.get(m.motif)}"></i>` : ""}`
          + `${escapeHTML(m.motif)}</b></td>`
          + `<td class="chance">${escapeHTML(m.matrix_id)}</td>`
          + `<td>${m.sites}</td><td class="chance">${m.expected}</td>`
          + `<td class="${m.exceeds_chance ? "enriched" : "chance"}">`
          + `${foldText(m)}</td>`
          + `<td class="chance">${escapeHTML(m.data_type || "n/a")}`
          + `${m.in_vivo ? "" : " (in vitro)"}</td>`
          + `<td class="relevance">${escapeHTML(m.relevance || "")}</td></tr>`;
      }).join("")
    + `</tbody></table></div>`
    + (focused ? `<p class="hint"><b>Showing ${escapeHTML(focused)} only.</b> `
        + `<button type="button" class="text-button" id="clear-motif-focus">Show all factors`
        + `</button></p>` : "");
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
  loading("saved-status", true, "Loading workbench…");
  try {
    state.catalogue=frozen?.catalogue || await api("/api/catalogue");
    $("engine-version").textContent=`v${state.catalogue.engine_version}`;
    $("preset").innerHTML=state.catalogue.presets.map(p=>`<option value="${p.id}">${escapeHTML(p.name)}</option>`).join("");
    $("cell-lines").innerHTML=state.catalogue.cells.map(c=>`<option value="${escapeHTML(c.name)}"></option>`).join("");
    renderCellEvidence();
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
  finally { loading("saved-status", false); }
}
$("design-form").addEventListener("submit",e=>{e.preventDefault();design();});
$("lookup-gene").addEventListener("click",()=>{error("");lookup().catch(e=>error(e,"lookup"));});
$("gene-symbol").addEventListener("input",()=>{state.gene=null;state.expression=null;$("expression").hidden=true;state.basket=[];state.customReferences={};resetCustomEditor();renderBasket();renderGeneOverview();$("reference").disabled=true;$("reference").innerHTML='<option value="">Look up this gene first</option>';$("gene-info").hidden=true;$("add-promoter").disabled=true;$("save-custom-reference").disabled=true;error("");});
$("reference").addEventListener("change",()=>{$("add-promoter").disabled=!$("reference").value;renderExpression();renderKnownSaRNAPromoters();});
$("cell-line").addEventListener("input",renderCellEvidence);
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
  renderExpression();   /* Setting .value in script fires no change event. */
  renderKnownSaRNAPromoters();
  $("reference").focus();
}
$("tss-group-list").addEventListener("click",e=>{const button=e.target.closest("[data-choose-transcript]");if(button)chooseOverviewReference(`transcript:${button.dataset.chooseTranscript}`);});
$("expression").addEventListener("click",e=>{const button=e.target.closest("[data-choose-transcript]");if(button)chooseOverviewReference(`transcript:${button.dataset.chooseTranscript}`);});
$("expression-overview").addEventListener("click",e=>{const button=e.target.closest("[data-choose-transcript]");if(button)chooseOverviewReference(`transcript:${button.dataset.chooseTranscript}`);});
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
  loading("saved-status", true, "Opening saved results…");
  $("saved-runs").disabled=true;
  try {
    const id=$("saved-runs").value;
    if (id.startsWith("batch:")) await openSavedBatch(id.slice(6));
    else { const run=await api(`/api/runs/${id}`); if (run.status!=="complete") throw new Error(run.error || run.stage); state.batch=null;$("batch-comparison").hidden=true;restoreDraft(run.result);showResult(run.result); }
    error("");
  } catch(e) { error(e); }
  finally { loading("saved-status", false); $("saved-runs").disabled=state.busy; }
});
document.querySelectorAll("[data-filter]").forEach(b=>b.addEventListener("click",()=>{state.filter=b.dataset.filter;state.page=0;document.querySelectorAll("[data-filter]").forEach(v=>v.classList.toggle("active",v===b));renderTable();}));
$("candidate-search").addEventListener("input",()=>{state.page=0;renderTable();});
document.querySelectorAll("[data-sort]").forEach(b=>b.addEventListener("click",()=>{state.descending=state.sort===b.dataset.sort?!state.descending:b.dataset.sort==="rating";state.sort=b.dataset.sort;renderTable();}));
$("candidate-rows").addEventListener("click",e=>{if (!e.target.closest("input,[data-screening-candidate]")) {const row=e.target.closest("[data-candidate]");if(row)selectCandidate(row.dataset.candidate);}});
$("candidate-rows").addEventListener("change",e=>{const id=e.target.dataset.check;if(id){e.target.checked?state.checked.add(id):state.checked.delete(id);renderTable();}});
$("candidate-rows").addEventListener("keydown",e=>{if ((e.key==="Enter" || e.key===" ") && e.target.matches("tr[data-candidate]")){e.preventDefault();selectCandidate(e.target.dataset.candidate);}});
// The detail panel names the candidates ordering as this same duplex; they are buttons so
// the twin can be inspected, which is the point of naming it rather than counting it.
$("candidate-detail").addEventListener("click",e=>{const v=e.target.closest("[data-candidate]");if(v)selectCandidate(v.dataset.candidate);});
$("promoter-map").addEventListener("click",e=>{
  if(state.suppressMapClick){state.suppressMapClick=false;return;}
  const v=e.target.closest("[data-candidate]");
  if(v){selectCandidate(v.dataset.candidate);return;}
  // A column of the excluded-target histogram zooms to the targets it is counting, which is
  // how the individual rejections stay reachable at full-region zoom.
  const bin=e.target.closest("[data-zoom-from]");
  if(bin){
    const from=+bin.dataset.zoomFrom, to=+bin.dataset.zoomTo;
    setViewport((from+to)/2-INDIVIDUAL_MARK_SPAN/2, INDIVIDUAL_MARK_SPAN);
    renderMap();
  }
});
enableMapDrag();
$("promoter-map").addEventListener("keydown",e=>{if ((e.key==="Enter" || e.key===" ") && e.target.dataset.candidate){e.preventDefault();selectCandidate(e.target.dataset.candidate);}});
$("page-previous").addEventListener("click",()=>{state.page--;renderTable();});$("page-next").addEventListener("click",()=>{state.page++;renderTable();});
$("select-page").addEventListener("change",e=>{filteredCandidates().slice(state.page*50,(state.page+1)*50).forEach(c=>e.target.checked?state.checked.add(c.id):state.checked.delete(c.id));renderTable();});
$("select-top").addEventListener("click",()=>{state.checked.clear();state.result.candidates.filter(c=>c.rank!==null).sort((a,b)=>a.rank-b.rank).slice(0,10).forEach(c=>state.checked.add(c.id));renderTable();});
$("zoom-in").addEventListener("click",()=>viewport(.5));$("zoom-out").addEventListener("click",()=>viewport(2));$("pan-left").addEventListener("click",()=>viewport(1,-.3));$("pan-right").addEventListener("click",()=>viewport(1,.3));$("zoom-reset").addEventListener("click",()=>{state.start=0;state.end=state.result.promoter.length;renderMap();});
$("candidate-detail").addEventListener("click",async e=>{const b=e.target.closest("[data-copy]");if(b){try{await navigator.clipboard.writeText(b.dataset.copy);b.textContent="Copied";}catch{error("Clipboard unavailable. Select and copy the displayed sequence.");}}});
$("run-evidence").addEventListener("click",renderEvidence);
for (const type of ["csv","json","html","xlsx"]) $("export-"+type).addEventListener("click",()=>exportRun(type));
$("export-order").addEventListener("click",exportOrder);
$("motif-panel").addEventListener("click",event=>{
  if (event.target.closest("#clear-motif-focus")) {
    state.motifFocus=null; renderMotifPanel(state.result); renderMap(); return;
  }
  const row=event.target.closest("[data-motif]");
  if (!row) return;
  state.motifFocus = state.motifFocus===row.dataset.motif ? null : row.dataset.motif;
  renderMotifPanel(state.result); renderMap();
});
$("motif-panel").addEventListener("keydown",event=>{
  const row=event.target.closest("[data-motif]");
  if (!row || (event.key!=="Enter" && event.key!==" ")) return;
  event.preventDefault();
  state.motifFocus = state.motifFocus===row.dataset.motif ? null : row.dataset.motif;
  renderMotifPanel(state.result); renderMap();
});
document.addEventListener("click", event => {
  const target = event.target.closest("#print-duplex,#inspect-top,#inspect-top-duplex,#export-order-top");
  if (!target) return;
  if (target.id === "print-duplex") printDuplexSheet();
  if (target.id === "export-order-top") {
    if (typeof setWorkspaceTab === "function") setWorkspaceTab("export", {focus:true});
    else exportOrder();
  }
  if (target.id === "inspect-top" || target.id === "inspect-top-duplex") {
    const best = state.result?.candidates.find(c => c.rank === 1);
    if (best) { selectCandidate(best.id); $("candidate-detail").scrollIntoView({behavior:"smooth", block:"start"}); }
  }
});
document.addEventListener("pointermove",e=>tooltip(e.target,e));document.addEventListener("focusin",e=>tooltip(e.target));document.addEventListener("focusout",()=>{$("tooltip").hidden=true;});document.addEventListener("scroll",()=>{$("tooltip").hidden=true;},true);
initialize();

// Existing duplexes retain their supplied identity; Python owns strand matching and predictions.
function renderSequenceRecords() {
  const records = state.sequenceRecords || [];
  const opened=new Set([...$("sequence-records").querySelectorAll('details[data-record-key][open]')].map(node=>node.dataset.recordKey));
  $("sequence-analysis-notice").textContent=state.result?.sequence_comparison_notice || '';
  $("sequence-recalculate").hidden=!!frozen;
  $("sequence-import-controls").hidden = !!frozen;
  $("sequence-records").innerHTML = records.length ? `<div class="table-scroll"><table><thead><tr><th>Reference / import</th><th>Placement</th><th>Source / observation</th><th></th></tr></thead><tbody>${records.map((r,i)=>`<tr><td><button type="button" data-sequence="${i}">${escapeHTML(r.name)}</button>${r.reference?' · reference':''}<br>${sequenceEffectBadge(r)}</td><td>${r.matches.length?r.matches.map((m,j)=>`<button type="button" data-sequence="${i}" data-match="${j}">Bases ${m.offset+1}–${m.end}</button>`).join(' '):escapeHTML(r.status)}${!r.matches.length&&r.known_location?`<br><small>Known site: chr${escapeHTML(r.known_location.contig)}:${fmt(r.known_location.start+1)}–${fmt(r.known_location.end)} · ${escapeHTML(r.known_location.assembly)} · ${escapeHTML(r.known_location.strand)} strand (1-based)</small>${!frozen?`<br><button type="button" data-sarna-region="${escapeHTML(r.name)}">Choose this promoter</button>`:''}`:''}${r.matches.length>1?' · multiple matches':''}</td><td>${escapeHTML(r.source)}</td><td>${!r.reference&&!frozen?`<button type="button" data-sequence-remove="${i}">Remove</button>`:''}</td></tr><tr class="reference-analysis-row"><td colspan="4">${referenceDuplexPanel(r)}<details class="reference-strand-disclosure" data-record-key="${escapeHTML(JSON.stringify([r.name,r.sequence,r.role]))}"><summary>Analyse both strands · ${escapeHTML(r.name)}</summary>${strandAnalysisPanel(r)}</details></td></tr>`).join('')}</tbody></table></div>` : '<p class="hint">No existing saRNAs added for this run.</p>';
  for (const node of $("sequence-records").querySelectorAll('details[data-record-key]')) node.open=opened.has(node.dataset.recordKey);
  renderSequenceDetail();
}
function renderSequenceDetail() {
  const record = (state.sequenceRecords || [])[sequenceFocus?.record];
  if (!record) { $("sequence-detail").innerHTML=''; return; }
  const match = record.matches[sequenceFocus.match || 0];
  const selected = state.result?.candidates.find(c=>c.id===state.selected);
  const duplex = (name,sense,antisense)=>`<h4>${escapeHTML(name)}</h4><p>Sense 5′ <code>${escapeHTML(sense)}</code> 3′<br>Antisense 5′ <code>${escapeHTML(antisense)}</code> 3′</p>`;
  $("sequence-detail").innerHTML = `<h3>${escapeHTML(record.name)} ${sequenceEffectBadge(record)}</h3><p class="hint">Original: ${escapeHTML(record.sequence)} · ${escapeHTML(record.role)} · ${escapeHTML(record.direction)}. ${record.partner_inferred?'Complementary partner inferred.':'Both cores supplied and complementarity checked.'} Overhang: ${escapeHTML(record.overhang || 'unknown')}; partner overhang: ${escapeHTML(record.partner_overhang || 'unknown')}.</p>` + (match ? duplex('Located duplex',match.sense,match.antisense) + `<p class="hint">Input orientation: ${escapeHTML(match.input_orientation)}. ${match.candidate?`Current run: ${escapeHTML(match.candidate.id)} · rank ${match.candidate.rank ?? 'excluded'} · preferences ${match.candidate.rating.met}/${match.candidate.rating.evaluated} · excluded by ${escapeHTML(match.candidate.excluded_by?.join(', ') || 'none')}`:'No candidate of this length in the current run; rule rating unavailable.'}</p><div class="table-scroll"><table><thead><tr><th>Model</th><th>Predicted guide</th><th>Predicted passenger</th><th>End difference</th></tr></thead><tbody>${match.thermodynamics.map(a=>`<tr title="${escapeHTML(a.explanation)}"><td>${escapeHTML(a.table_name)}</td><td>${escapeHTML(a.predicted_guide)}</td><td>${a.predicted_guide==='antisense'?'sense':a.predicted_guide==='sense'?'antisense':'undetermined'}</td><td>${a.difference.toFixed(2)} ${escapeHTML(a.unit)}</td></tr>`).join('')}</tbody></table></div><p class="hint">Predictions use paired cores; matched candidates retain the run’s model settings. Other lengths use four-base ends. Different model predictions remain separate. They do not establish activation or measured strand loading.</p>` + (selected ? duplex('Selected new design: '+selected.id, selected.sense_strand_rna,selected.antisense_strand_rna) + `<p>Predicted guide: ${escapeHTML(selected.rating.guide)} · preferences ${selected.rating.met}/${selected.rating.evaluated}. Target overlap: ${Math.max(0,Math.min(match.end,selected.offset+selected.target_site_dna.length)-Math.max(match.offset,selected.offset))} bases. GC: existing ${(100*(match.target.match(/[GC]/g)||[]).length/match.target.length).toFixed(1)}%, selected ${selected.gc_percent}%. Aligned core differences: ${match.target.length===selected.target_site_dna.length?[...match.target].filter((base,i)=>base!==selected.target_site_dna[i]).length:"different lengths"}.</p>` : '<p class="hint">Select a new candidate on the map or in the table to compare its duplex here.</p>') + '<button type="button" id="sequence-nearby">Find overlapping / nearby designs</button>' : `<p>${escapeHTML(record.status)}. ${record.known_location?`Known site: chr${escapeHTML(record.known_location.contig)}:${fmt(record.known_location.start+1)}–${fmt(record.known_location.end)} (${escapeHTML(record.known_location.assembly)}, 1-based). ${escapeHTML(record.known_location.source)}. The current map covers a different region; use “Choose this promoter” to update your draft, then submit the design when ready.`:'The sequence may lie outside this region or differ from the reference. No location is assumed.'}</p>`);
  $("sequence-detail").insertAdjacentHTML("beforeend",referenceDuplexPanel(record,match)+strandAnalysisPanel(record,match));
  $("sequence-nearby")?.addEventListener('click',()=>{
    const nearby = state.result.candidates.filter(c=>Math.abs(c.offset-match.offset)<=100).sort((a,b)=>Math.abs(a.offset-match.offset)-Math.abs(b.offset-match.offset)).slice(0,10);
    $("sequence-detail").insertAdjacentHTML('beforeend',`<p>Within 100 bases: ${nearby.map(c=>`<button type="button" data-candidate="${escapeHTML(c.id)}">${escapeHTML(c.id)}</button>`).join(' ') || 'No generated candidates'}</p>`);
  });
}
function previewSequences(records,errors=[]) {
  sequencePending=records;
  $("sequence-status").textContent=errors.join(' ');
  $("sequence-preview").innerHTML=records.length?`<p>${records.length} sequences ready to add:</p><ul>${records.map(r=>`<li>${escapeHTML(r.name)}: <code>${escapeHTML(r.sequence)}</code> · ${escapeHTML(r.role)}</li>`).join('')}</ul><button type="button" id="sequence-confirm">Add these sequences to this run</button>`:'';
  $("sequence-confirm")?.addEventListener('click', async()=>{
    try { await saveSequenceRecords([...(state.sequenceRecords||[]).filter(r=>!r.reference),...sequencePending]);sequencePending=[];$("sequence-preview").innerHTML='';$("sequence-input").value=''; }
    catch(e) { $("sequence-status").textContent=e.message; }
  });
}
function sequenceInputs(records) {
  return records.map(record=>Object.fromEntries(['name','sequence','antisense','role','direction','source']
    .filter(key=>record[key]!==undefined).map(key=>[key,record[key]])));
}
function applySavedSequences(data) {
  state.sequenceRecords=[...data.references,...data.records];
  state.sequenceRevision=data.revision;
  state.result.sequence_comparisons=state.sequenceRecords;
  state.result.sequence_comparison_revision=data.revision;
  state.result.sequence_comparison_notice=data.notice || '';
  sequenceFocus=null;
  renderSequenceRecords();renderMap();renderTable();
}
async function saveSequenceRecords(records) {
  if (state.sequenceSaving) throw new Error('A sequence save is already in progress. Retry when it finishes.');
  const id=state.result.app.run_id;
  state.sequenceSaving=true;
  try {
    const data=await api(`/api/runs/${id}/sequence-comparisons`,{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({revision:state.sequenceRevision,records:sequenceInputs(records)})});
    if (state.result.app.run_id!==id) throw new Error('The run changed while saving. Sequences were saved to the original run; pending input is retained.');
    applySavedSequences(data);
    $("sequence-conflict-recovery").hidden=true;
    $("sequence-status").textContent='Comparison saved with this run.';
  } catch (failure) {
    if (failure.code==='sequence-comparison-conflict' && state.result.app.run_id===id) $("sequence-conflict-recovery").hidden=false;
    throw failure;
  } finally {state.sequenceSaving=false;}
}
$("sequence-reload").addEventListener('click',async()=>{
  const id=state.result.app.run_id;
  try {
    const data=await api(`/api/runs/${id}/sequence-comparisons`);
    if (state.result.app.run_id!==id) return;
    applySavedSequences(data);
    $("sequence-conflict-recovery").hidden=true;
    $("sequence-status").textContent='Saved sequences reloaded. Pending input is retained; add it again when ready.';
  } catch (failure) {$("sequence-status").textContent=failure.message;}
});
$("sequence-recalculate").addEventListener('click',async()=>{
  const id=state.result.app.run_id;
  if (state.sequenceSaving) return;
  state.sequenceSaving=true;
  $("sequence-recalculate").disabled=true;
  try {
    const data=await api(`/api/runs/${id}/sequence-comparisons/recalculate`,{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({revision:state.sequenceRevision})});
    if (state.result.app.run_id!==id) return;
    applySavedSequences(data);
    $("sequence-status").textContent='Comparisons recalculated explicitly with current software. The prior revision remains recorded.';
  } catch (failure) {
    if (failure.code==='sequence-comparison-conflict') $("sequence-conflict-recovery").hidden=false;
    $("sequence-status").textContent=failure.message;
  } finally {state.sequenceSaving=false;$("sequence-recalculate").disabled=false;}
});
$("sequence-add").addEventListener('click',async()=>{
  try {
    const lines=$("sequence-input").value.trim().split(/\n/).filter(Boolean);
    const role=$("sequence-role").value, direction=$("sequence-direction").value;
    const text='name,sequence,role,direction\n'+lines.map((line,i)=>{const fields=line.split(',');return `${fields.length>1?fields[0]:'Sequence '+(i+1)},${fields.at(-1).trim()},${role},${direction}`;}).join('\n');
    const bytes=new TextEncoder().encode(text);
    const data=await api('/api/sequences/import',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({filename:'manual.csv',content:btoa(String.fromCharCode(...bytes))})});
    previewSequences(data.records,data.errors);
  } catch(e) {$("sequence-status").textContent=e.message;}
});
$("sequence-file").addEventListener('change',async(event)=>{
  try {
    const file=event.target.files[0];if(!file)return;
    if(file.size>2000000)throw new Error('Upload at most 2 MB.');
    const bytes=new Uint8Array(await file.arrayBuffer());let binary='';for(const byte of bytes)binary+=String.fromCharCode(byte);
    const data=await api('/api/sequences/import',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({filename:file.name,content:btoa(binary)})});
    previewSequences(data.records,data.errors);
  } catch(e) {$("sequence-status").textContent=e.message;} finally {event.target.value='';}
});
$("sequence-template").addEventListener('click',()=>{
  const link=document.createElement('a');link.href=URL.createObjectURL(new Blob(['name,sequence,role,direction,sense,antisense,source\nExample,CUUAACUCGCCUUCCAAAA,sense,5to3,,,Lab notebook\n'],{type:'text/csv'}));link.download='sarna-import-template.csv';link.click();URL.revokeObjectURL(link.href);
});
document.addEventListener('click',async(event)=>{
  const remove=event.target.closest('[data-sequence-remove]');
  if(remove){try{await saveSequenceRecords(state.sequenceRecords.filter((r,i)=>!r.reference&&i!==Number(remove.dataset.sequenceRemove)));}catch(e){$("sequence-status").textContent=e.message;}return;}
  const button=event.target.closest('[data-sequence]');
  if(button){sequenceFocus={record:Number(button.dataset.sequence),match:Number(button.dataset.match||0)};const match=state.sequenceRecords[sequenceFocus.record].matches[sequenceFocus.match];if(match){setViewport(match.offset-60,Math.max(160,match.end-match.offset+120));renderMap();}renderSequenceDetail();}
  const candidateButton=event.target.closest('#sequence-detail [data-candidate]');
  if(candidateButton) await selectCandidate(candidateButton.dataset.candidate);
  if(event.target.closest('[data-candidate]'))setTimeout(renderSequenceDetail,0);
});

$("promoter-map").addEventListener('keydown',event=>{
  if ((event.key==='Enter'||event.key===' ')&&event.target.hasAttribute('data-sequence')) {event.preventDefault();event.target.dispatchEvent(new MouseEvent('click',{bubbles:true}));}
});

// A curated locus outside the loaded sequence needs a new promoter run, not local panning.
async function designAroundSaRNA(name) {
  if (frozen || state.busy) return;
  const records=state.gene?.sarna_references || state.sequenceRecords || [];
  const record=records.find(r=>r.name===name && r.known_location);
  if (!record) return;
  try {
    if (!state.gene?.sarna_references) await lookup();
    const location=record.known_location;
    let key=location.anchor_id?`anchor:${location.anchor_id}`:`transcript:${location.nearest_transcript}`;
    if (![...$("reference").options].some(option=>option.value===key)) await lookup();
    if (![...$("reference").options].some(option=>option.value===key)) {
      const available=promoterSaRNACoverage(state.gene,state.customReferences,Number($("upstream").value),Number($("downstream").value))
        .find(row=>row.contains.some(reference=>reference.name===name)
          && [...$("reference").options].some(option=>option.value===row.key));
      if (available) key=available.key;
    }
    if (![...$("reference").options].some(option=>option.value===key)) throw new Error('The reference promoter is unavailable in the current annotation. Select a reference near the displayed genomic site.');
    $("reference").value=key;
    chooseOverviewReference(key);
    if (typeof setWorkspaceTab==='function') setWorkspaceTab('design',{scroll:false});
    $("sequence-status").textContent=`Reference for ${name} selected. Current criteria and window sizes are retained; submit the design when ready.`;
    $("job-status").hidden=false;
    $("job-status").textContent=`Reference for ${name} selected. Review your settings, then click Design saRNAs.`;
  } catch (failure) {$("sequence-status").textContent=failure.message;error(failure.message);}
}
document.addEventListener('click',event=>{
  const button=event.target.closest('[data-sarna-region]');
  if (button) designAroundSaRNA(button.dataset.sarnaRegion);
});
$("gene-overview-map").addEventListener('keydown',event=>{
  if ((event.key==='Enter'||event.key===' ')&&event.target.hasAttribute('data-sarna-region')) {event.preventDefault();designAroundSaRNA(event.target.dataset.sarnaRegion);}
});

for (const id of ['upstream','downstream']) $(id).addEventListener('input',renderKnownSaRNAPromoters);
$("known-sarna-promoters").addEventListener('click',event=>{
  const button=event.target.closest('[data-known-promoter]');
  if (button) chooseOverviewReference(button.dataset.knownPromoter);
});
