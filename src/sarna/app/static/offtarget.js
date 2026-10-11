/* Off-target screening: the panel, its batches, and where each alignment lands.
 *
 * Display only. Every threshold shown here is computed in Python and arrives in the record;
 * nothing in this file decides what counts as a strong match (ADR 0019, ADR 0023).
 *
 * The panel merges several screening records per scope, because screening now runs in batches
 * and a run can hold a first fifteen, a next fifteen and a re-screen months later. Each row
 * therefore carries the date of the search it came from rather than the panel claiming one.
 */
"use strict";

const SCOPE_BUTTON = {transcript:"screen-transcript", genomic:"screen-genomic"};
const SCOPE_NAME = {transcript:"Transcripts (RefSeq RNA)", genomic:"Promoters (genomic)"};
const BATCH_SIZES = [1, 5, 10, 15, 25];

// One hue, darker is closer, per ADR 0022. Green to red would read as good to bad, and an
// alignment is neither: it is a measurement the researcher judges.
const BANDS = [
  {key:"exact", label:"exact", fill:"#1d5c52"},
  {key:"strong", label:"at Wang's 80%", fill:"#3f8d7d"},
  {key:"partial", label:"partial", fill:"#9cc6ba"},
  {key:"weak", label:"weak", fill:"#dde9e4"},
];

const otTip = (text) =>
  `data-tip="${escapeHTML(text)}" tabindex="0" aria-label="${escapeHTML(text)}"`;

function screenings() {
  state.screenings = state.screenings || {};
  return state.screenings;
}

/* ---------------------------------------------------------------- merging records */

/** Fold one record into the merged view for its scope, newest search winning per candidate. */
function mergeScreening(record) {
  const result = record.result || {}, scope = record.scope;
  const held = screenings()[scope] || {scope, searches:[], byCandidate:{}};
  if (held.searches.some(s => s.id === record.id)) return held;
  const search = {
    id: record.id,
    searched_at: result.searched_at,
    parameters: result.parameters || {},
    notes: result.notes || [],
    request_id: result.request_id,
    window: result.window || null,
  };
  held.searches = [...held.searches, search];
  const hits = {}, placements = {};
  for (const hit of result.hits || []) (hits[hit.candidate] = hits[hit.candidate] || []).push(hit);
  for (const p of result.proximities || []) {
    (placements[p.hit.candidate] = placements[p.hit.candidate] || []).push(p);
  }
  for (const row of result.specificity || []) {
    held.byCandidate = {...held.byCandidate, [row.candidate]: {
      row, search, hits: hits[row.candidate] || [], placements: placements[row.candidate] || [],
    }};
  }
  screenings()[scope] = held;
  return held;
}

/** Every candidate id this run has screened in one scope, for "the next fifteen". */
function alreadyScreened(scope) {
  return new Set(Object.keys(screenings()[scope]?.byCandidate || {}));
}

/* ---------------------------------------------------------------- starting a search */

/** Passing candidates in rank order: the pool every batch is drawn from. */
function screenablePool() {
  return (state.result?.candidates || [])
    .filter(c => c.rank !== null && c.rank !== undefined)
    .sort((a, b) => a.rank - b.rank)
    .map(c => c.id);
}

function batchChoice(scope) {
  state.screenBatch = state.screenBatch || {};
  return state.screenBatch[scope] ?? 15;
}

/** Resolve a chooser into the explicit ids the server validates. One path, not one per button. */
function resolveBatch(scope, mode) {
  const pool = screenablePool(), done = alreadyScreened(scope);
  const remaining = pool.filter(id => !done.has(id));
  if (mode === "selection") {
    const ticked = pool.filter(id => state.checked.has(id));
    return {ids: ticked, why: "the candidates you ticked"};
  }
  if (mode === "rest") return {ids: remaining, why: "every candidate not yet screened"};
  if (mode === "all") return {ids: pool, why: "every candidate that passed the gates"};
  const size = batchChoice(scope);
  const take = done.size ? remaining : pool;
  return {ids: take.slice(0, size), why: done.size ? `the next ${size}` : `the top ${size}`};
}

async function startScreening(scope, mode = "batch") {
  // Opt-in, and the run is already readable: this adds evidence, it never gates the designs.
  const run = state.result;
  if (!run || frozen) { error("Screening needs a live run."); return; }
  state.screeningJobs = state.screeningJobs || {};
  if (state.screeningJobs[scope]?.status === "running") return;
  const {ids, why} = resolveBatch(scope, mode);
  if (!ids.length) {
    error(mode === "selection"
      ? "No candidates are ticked. Tick some in the table, or choose a batch size."
      : "Every candidate that passed the gates has already been screened in this scope.");
    return;
  }
  // Several jobs rather than one oversized submission: NCBI gets one at a time, and each
  // batch becomes its own dated record (ADR 0019).
  const chunks = [];
  for (let i = 0; i < ids.length; i += 25) chunks.push(ids.slice(i, i + 25));
  state.screeningJobs[scope] = {status:"running", stage:"Starting…", done:0, total:chunks.length};
  renderScreeningStatus();
  for (const [index, chunk] of chunks.entries()) {
    try {
      const started = await api(`/api/runs/${run.app.run_id}/screenings`,
        {method:"POST", headers:{"Content-Type":"application/json"},
         body:JSON.stringify({scope, candidates:chunk})});
      if (!started.id) throw new Error("The server did not return a screening id.");
      const label = chunks.length > 1 ? `batch ${index + 1} of ${chunks.length}` : why;
      await pollScreening(started.id, run.app.run_id, scope, label, index, chunks.length);
    } catch (e) {
      error(e.message);
      state.screeningJobs[scope] = {status:"error", stage:e.message};
      renderScreeningStatus();
      return;
    }
    if (state.result?.app.run_id !== run.app.run_id) return;
  }
  state.screeningJobs[scope] = {status:"complete", stage:"Finished"};
  renderScreeningStatus();
}

/** Resolve once the screening finishes; the server waits NCBI's own interval, this catches up. */
function pollScreening(screeningId, runId, scope, label, index, total) {
  return new Promise((resolve, reject) => {
    const tick = async () => {
      try {
        const record = await api(`/api/screenings/${screeningId}`);
        if (state.result?.app.run_id !== runId) { resolve(); return; }
        if (record.status === "running") {
          state.screeningJobs[scope] = {
            status:"running", stage:`${label}: ${record.stage || "working"}`,
            done:index, total,
          };
          renderScreeningStatus();
          setTimeout(tick, 2000);
          return;
        }
        if (record.status === "error") {
          reject(new Error(record.error || "The screening failed."));
          return;
        }
        mergeScreening(record);
        renderOffTarget(scope);
        // The three other places that state what has been screened are rendered elsewhere, so
        // refresh them here rather than let them go on saying none was performed.
        if (typeof renderTopResult === "function") renderTopResult(state.result);
        if (state.selected && typeof selectCandidate === "function") selectCandidate(state.selected);
        resolve();
      } catch (e) { reject(e); }
    };
    tick();
  });
}

/** Rebuild the merged panel for a saved run, so screenings survive a reload. */
async function loadSavedScreenings(runId) {
  try {
    const listed = await api(`/api/runs/${runId}/screenings`);
    let any = null;
    for (const entry of listed.filter(s => s.status === "complete").reverse()) {
      mergeScreening(await api(`/api/screenings/${entry.id}`));
      any = entry.scope;
    }
    if (any) { renderOffTarget(any); renderScreeningStatus(); }
  } catch { /* A run with no screenings is the ordinary case, not an error worth showing. */ }
}

/* ---------------------------------------------------------------- the controls */

function renderScreeningStatus() {
  const jobs = state.screeningJobs || {};
  const states = Object.values(jobs).map(j => j.status);
  const done = Object.keys(screenings()).length;
  $("offtarget-status").textContent = states.includes("running") ? "RUNNING"
    : states.includes("error") ? (done ? "PARTIAL" : "FAILED")
    : done ? "COMPLETE" : "NOT RUN";
  $("offtarget-progress").textContent = states.includes("running")
    ? "Searches run at NCBI and report here when they finish. You can keep working meanwhile."
    : "";
  for (const scope of Object.keys(SCOPE_BUTTON)) {
    const lane = $(`screen-state-${scope}`), job = jobs[scope];
    if (lane) {
      lane.textContent = job ? (job.stage || job.status) : "";
      lane.className = "screen-state " + (job ? job.status : "")
        + (job?.status === "running" ? " is-loading" : "");
    }
    const button = $(SCOPE_BUTTON[scope]);
    if (button) button.disabled = job?.status === "running";
  }
  renderBatchControls();
}

function renderBatchControls() {
  const host = $("offtarget-batch");
  if (!host) return;
  const pool = screenablePool();
  if (!pool.length || frozen) { host.innerHTML = ""; return; }
  const scope = state.screeningScope || "transcript";
  const done = alreadyScreened(scope), remaining = pool.length - done.size;
  const size = batchChoice(scope);
  const ticked = pool.filter(id => state.checked.has(id)).length;
  host.innerHTML =
    `<div class="batch-row"><span class="batch-label">How many at a time</span>`
    + `<span class="segmented">`
    + BATCH_SIZES.map(n => `<button type="button" data-batch="${n}"`
        + `${n === size ? ' class="active"' : ""}>${n}</button>`).join("")
    + `</span>`
    + `<button type="button" data-batch-mode="all" class="secondary"`
    + ` ${otTip("Screen every candidate that passed the gates. Submitted in sequential jobs of "
        + "25, each its own dated record, so NCBI is never asked for all of them at once.")}>`
    + `All ${pool.length}</button>`
    + (done.size && remaining > 0
        ? `<button type="button" data-batch-mode="rest" class="secondary">`
          + `The rest (${remaining})</button>` : "")
    + (ticked
        ? `<button type="button" data-batch-mode="selection" class="secondary">`
          + `Ticked (${ticked})</button>` : "")
    + `</div>`
    + `<p class="hint">${!done.size
        ? `Nothing screened yet in this scope. Fewer candidates finish sooner: a genomic `
          + `search of fifty can take longer than NCBI will hold it.`
        : remaining > 0
        ? `${done.size} of ${pool.length} screened in this scope. The next search covers `
          + `${Math.min(size, remaining)} more unless you choose otherwise.`
        : `All ${pool.length} screened in this scope. Screening again searches the same `
          + `sequences against whatever the database holds today, as a second dated record.`
      }</p>`;
}

/* ---------------------------------------------------------------- the results */

/** Which band an alignment falls in, from thresholds the record carries. */
function band(hit, row, core) {
  if (hit.identities >= core && hit.gaps === 0 && hit.alignment_length === core) return "exact";
  if (hit.identities >= (row.strong_threshold || core)) return "strong";
  if (hit.identities >= Math.ceil(core * 0.6)) return "partial";
  return "weak";
}

function bandCounts(entry, core) {
  const counts = {exact:0, strong:0, partial:0, weak:0};
  for (const hit of entry.hits) counts[band(hit, entry.row, core)] += 1;
  return counts;
}

/** A compact strip of one candidate's alignments, binned by how much of the core matched. */
function identityStrip(counts, total) {
  if (!total) return `<span class="strip empty" ${otTip("No alignments were returned.")}></span>`;
  const widest = Math.max(total, 1);
  return `<span class="strip">` + BANDS.map(b => counts[b.key]
      ? `<i style="width:${(counts[b.key] / widest) * 100}%;background:${b.fill}" `
        + otTip(`${counts[b.key]} ${b.label}`) + `></i>`
      : "").join("") + `</span>`;
}

/** The strongest alignment, named by gene rather than by accession. */
function worstHit(entry) {
  let best = null;
  for (const hit of entry.hits) if (!best || hit.identities > best.identities) best = hit;
  return best;
}

/** A readable gene name out of a BLAST description, falling back to the accession. */
function hitName(hit, placement) {
  if (placement?.nearest_gene) return placement.nearest_gene;
  const title = (hit.title || "").replace(/^Homo sapiens\s+/i, "");
  const symbol = title.match(/\(([A-Z0-9orf-]{2,12})\),/);
  if (symbol) return symbol[1];
  return title.split(",")[0].slice(0, 48) || hit.accession;
}

function renderOffTarget(scope) {
  state.screeningScope = scope || state.screeningScope;
  const held = screenings()[state.screeningScope];
  const host = $("offtarget-results");
  if (!host) return;
  if (!held) { host.innerHTML = ""; renderBatchControls(); return; }
  const core = state.result?.candidates?.[0]?.target_site_dna?.length || 19;
  const entries = Object.values(held.byCandidate);
  const genomic = state.screeningScope === "genomic";
  host.innerHTML =
    scopeTabs()
    + summaryLine(held, entries, core, genomic)
    + notesFor(held)
    + resultsTable(held, entries, core, genomic)
    + (genomic ? tssPlot(entries) : "");
  renderBatchControls();
}

function scopeTabs() {
  const present = Object.keys(screenings());
  if (present.length < 2) return "";
  return `<div class="segmented scope-tabs">` + present.map(s =>
    `<button type="button" data-scope-tab="${s}"`
    + `${s === state.screeningScope ? ' class="active"' : ""}>`
    + `${escapeHTML(SCOPE_NAME[s] || s)}</button>`).join("") + `</div>`;
}

function summaryLine(held, entries, core, genomic) {
  const strong = entries.filter(e => e.row.wang_homology_hits > 0).length;
  const exact = entries.filter(e => e.row.perfect_hits > 0).length;
  const quiet = entries.filter(e => e.row.hit_count === 0).length;
  const promiscuous = entries.filter(e => e.row.promiscuous).length;
  const threshold = entries[0]?.row.strong_threshold || Math.ceil(core * 0.8);
  const inWindow = genomic
    ? entries.filter(e => (e.row.in_run_window_hits || 0) > 0).length : 0;
  const dates = held.searches.map(s => new Date(s.searched_at));
  const when = dates.length
    ? new Date(Math.max(...dates)).toLocaleDateString() : "an unknown date";
  const databases = [...new Set(held.searches.map(s => s.parameters.database))].join(", ");
  const findings = [
    exact ? `<b>${exact}</b> with a full-length exact match elsewhere` : "",
    strong ? `<b>${strong}</b> matching at least ${threshold} of ${core} bases, Wang 2015's 80%` : "",
    promiscuous ? `<b>${promiscuous}</b> occurring in many places` : "",
    inWindow ? `<b>${inWindow}</b> with a hit near another gene's start site, inside this run's own window` : "",
    quiet ? `<b>${quiet}</b> with no alignments at all` : "",
  ].filter(Boolean);
  return `<div class="ot-summary"><p><b>${escapeHTML(SCOPE_NAME[state.screeningScope])}</b> · `
    + `${entries.length} candidate${entries.length === 1 ? "" : "s"} screened across `
    + `${held.searches.length} search${held.searches.length === 1 ? "" : "es"}, latest `
    + `${escapeHTML(when)} · ${escapeHTML(databases)}</p>`
    + `<p>${findings.length ? findings.join(". ") + "." : "Nothing found."}</p>`
    + `<p class="hint">Sequence similarity is not measured binding, and nothing here changes a `
    + `rank. Wang 2015 applies 80% to control duplexes rather than to candidate off-targets, so `
    + `that column is this project's reading of the paper's only stated threshold. The paper's `
    + `own caveat stands: carefully selected saRNAs cannot avoid partial homology with other `
    + `sequences, and putative off-target genes should be screened experimentally.</p></div>`;
}

function notesFor(held) {
  const notes = [...new Set(held.searches.flatMap(s => s.notes))];
  return notes.map(n => `<p class="hint warn">${escapeHTML(n)}</p>`).join("");
}

function resultsTable(held, entries, core, genomic) {
  const sort = state.screeningSort || "wang_homology_hits";
  const value = (e) => {
    const v = e.row[sort];
    return v === null || v === undefined ? (sort === "best_evalue" ? Infinity : -1) : v;
  };
  const ordered = [...entries].sort((a, b) =>
    sort === "best_evalue" ? value(a) - value(b) : value(b) - value(a));
  const quiet = ordered.filter(e => e.row.hit_count === 0);
  const found = ordered.filter(e => e.row.hit_count > 0);
  const header = (key, label, tip) =>
    `<th><button data-screensort="${key}"${sort === key ? ' class="active"' : ""} `
    + `data-tip="${escapeHTML(tip)}" tabindex="0">${label} &#8645;</button></th>`;
  const threshold = entries[0]?.row.strong_threshold || Math.ceil(core * 0.8);
  return `<div class="table-scroll"><table class="ot-table"><thead><tr>`
    + `<th>Candidate</th><th>Strongest match</th><th>Alignments by identity</th>`
    + header("perfect_hits", "Exact",
             `Full-length exact matches elsewhere, all ${core} bases with no gaps. One of these `
             + `is a specific unintended target.`)
    + header("wang_homology_hits", `Wang 80%`,
             `Alignments matching at least ${threshold} of ${core} bases. 80% is the one `
             + `homology threshold Wang 2015 states, where it rejects a control duplex.`)
    + header("best_evalue", "Best E-value",
             "Lowest E-value among the alignments. Smaller means less likely by chance.")
    + header("hit_count", "Alignments", "How many alignments were returned at all.")
    + (genomic ? header("in_run_window_hits", "In your window",
         "Hits inside the same window this run designed against, near another gene's start "
         + "site. A statement about this run's configuration, not a definition of a promoter.")
       : "")
    + `</tr></thead><tbody>`
    + found.map(e => candidateRow(e, core, genomic)).join("")
    + (quiet.length
        ? `<tr class="ot-quiet"><td colspan="${genomic ? 8 : 7}">`
          + `<details><summary>${quiet.length} candidate${quiet.length === 1 ? "" : "s"} `
          + `returned no alignments. That is not proof of specificity: it is one search, of one `
          + `database, on one day.</summary><p>`
          + quiet.map(e => escapeHTML(e.row.candidate)).join(", ")
          + `</p></details></td></tr>`
        : "")
    + `</tbody></table></div>`;
}

function candidateRow(entry, core, genomic) {
  const row = entry.row, counts = bandCounts(entry, core);
  const best = worstHit(entry);
  const placement = best
    ? entry.placements.find(p => p.hit.accession === best.accession) : null;
  const open = state.screeningOpen === row.candidate;
  const flag = (n) => n ? ' class="flagged"' : "";
  const main = `<tr data-screenrow="${escapeHTML(row.candidate)}"`
    + `${open ? ' class="open"' : ""}>`
    + `<td>${escapeHTML(row.candidate)}`
    + (row.alignments_capped ? `<p class="candidate-id" ${otTip(
        "More alignments were returned than were kept, so every count in this row is a floor.")
      }>list capped</p>` : "")
    + (row.on_target_hits ? `<p class="candidate-id" ${otTip(
        `${row.on_target_hits} alignment at this candidate's own target locus was removed from `
        + `these findings. A sequence matching the promoter it was designed against is the `
        + `intended target, not an off-target one.`)}>own locus removed</p>` : "")
    + `</td>`
    + `<td>${best
        ? `<b>${escapeHTML(hitName(best, placement))}</b> <span class="muted">`
          + `${best.identities}/${core}</span>`
          + (placement && placement.distance_to_tss !== null
             ? `<p class="candidate-id">${signed(placement.distance_to_tss)} bp from its TSS</p>`
             : "")
        : `<span class="muted">none</span>`}</td>`
    + `<td>${identityStrip(counts, row.hit_count)}`
    + (row.promiscuous ? `<p class="candidate-id">occurs in many places</p>` : "")
    + `</td>`
    + `<td${flag(row.perfect_hits)}>${row.perfect_hits}</td>`
    + `<td${flag(row.wang_homology_hits)}>${row.wang_homology_hits}</td>`
    + `<td>${row.best_evalue === null || row.best_evalue === undefined
        ? "n/a" : Number(row.best_evalue).toExponential(1)}</td>`
    + `<td>${row.hit_count}</td>`
    + (genomic ? `<td${flag(row.in_run_window_hits)}>${row.in_run_window_hits || 0}</td>` : "")
    + `</tr>`;
  return main + (open ? alignmentDetail(entry, core, genomic) : "");
}

/* ---------------------------------------------------------------- one candidate's hits */

/** Group a candidate's alignments by what they landed in, not one row per HSP. */
function groupHits(entry) {
  const groups = new Map();
  for (const hit of entry.hits) {
    const placement = entry.placements.find(
      p => p.hit.accession === hit.accession && p.hit.subject_start === hit.subject_start);
    const key = hit.accession;
    const held = groups.get(key) || {accession:hit.accession, title:hit.title, hits:[], placement};
    held.hits = [...held.hits, hit];
    held.placement = held.placement || placement;
    groups.set(key, held);
  }
  return [...groups.values()].sort((a, b) =>
    Math.max(...b.hits.map(h => h.identities)) - Math.max(...a.hits.map(h => h.identities)));
}

/** Draw which part of the core matched, and where the mismatches fall.
 *
 * No seed region is marked. Wang 2015 names none, and positions 2 to 8 are an siRNA
 * convention; drawing it here would put a borrowed mechanism on screen with this tool's
 * authority behind it (ADR 0023).
 */
function footprint(hit, core) {
  const width = 150, unit = width / core;
  const from = Math.min(hit.query_start, hit.query_end) - 1;
  const to = Math.max(hit.query_start, hit.query_end);
  const mismatched = new Set(hit.mismatch_positions || []);
  const counted = Math.max(0, (to - from) - hit.identities);
  const cells = [];
  for (let i = 0; i < core; i += 1) {
    const inside = i >= from && i < to;
    // Three states, not two. Outside the alignment, matched, and mismatched: drawing the
    // last two alike makes 16 of 19 and 19 of 19 identical pictures, which is what the
    // first version of this panel did.
    const fill = !inside ? "#e8eeeb" : mismatched.has(i + 1) ? "#c08a2e" : "#3f8d7d";
    cells.push(`<rect x="${(i * unit).toFixed(2)}" y="0" width="${(unit - 0.6).toFixed(2)}" `
      + `height="11" rx="1" fill="${fill}"/>`);
  }
  const where = hit.has_alignment_detail
    ? (mismatched.size
        ? `, mismatched at base${mismatched.size === 1 ? "" : "s"} `
          + [...mismatched].join(", ") + " of the core"
        : " with no mismatch in the aligned span")
    : (counted
        ? `, ${counted} mismatched somewhere in the aligned span; this report did not carry `
          + `the aligned sequences, so which bases is not known`
        : ", and this report did not carry the aligned sequences");
  const label = `Query bases ${from + 1} to ${to} of ${core} aligned, `
    + `${hit.identities} matching${where}. The ${hit.subject_strand} strand of the subject. `
    + `No seed region is marked: Wang 2015 names none.`;
  return `<svg class="footprint" viewBox="0 0 ${width} 11" width="${width}" height="11" `
    + `role="img" ${otTip(label)}>${cells.join("")}</svg>`;
}

function alignmentDetail(entry, core, genomic) {
  const groups = groupHits(entry);
  const shown = entry.row.promiscuous && !state.screeningExpandAll
    ? groups.slice(0, 5) : groups;
  const body = shown.map(group => {
    const strongest = group.hits.reduce((a, b) => a.identities >= b.identities ? a : b);
    const p = group.placement;
    return `<tr><td><b>${escapeHTML(hitName(strongest, p))}</b>`
      + `<p class="candidate-id">${escapeHTML(group.accession)}`
      + `${group.hits.length > 1 ? ` · ${group.hits.length} alignments` : ""}</p></td>`
      + `<td class="ot-title">${escapeHTML((group.title || "").slice(0, 90))}</td>`
      + `<td>${strongest.identities}/${core}</td>`
      + `<td>${footprint(strongest, core)}</td>`
      + `<td>${escapeHTML(strongest.subject_strand)}</td>`
      + `<td>${Number(strongest.evalue).toExponential(1)}</td>`
      + (genomic ? `<td>${p
          ? `${escapeHTML(p.contig)}${p.on_chromosome === false ? "" : `:${fmt(p.position)}`}`
            + (p.on_chromosome === false
               ? `<p class="candidate-id">not an assembled chromosome, so it could not be `
                 + `placed</p>`
               : p.distance_to_tss === null
               ? `<p class="candidate-id">no annotated transcript within 10 kb</p>`
               : `<p class="candidate-id">${signed(p.distance_to_tss)} bp from `
                 + `${escapeHTML(p.nearest_gene || p.nearest_transcript || "an unnamed transcript")}`
                 + `${p.inside_run_window ? " · inside your window" : ""}</p>`)
          : `<span class="muted">unplaced</span>`}</td>` : "")
      + `</tr>`;
  }).join("");
  return `<tr class="ot-detail"><td colspan="${genomic ? 8 : 7}">`
    + `<p class="hint">${escapeHTML(entry.row.summary || "")}</p>`
    + `<div class="table-scroll"><table><thead><tr>`
    + `<th>Where</th><th>What it is</th><th>Matched</th>`
    + `<th ${otTip("Which part of the core aligned, base by base. Teal matched, gold did not, "
        + "pale was outside the alignment. No seed region is marked: Wang 2015 names none, and "
        + "positions 2 to 8 are an siRNA convention this tool does not import.")}>Footprint `
      + `<span class="fp-key"><i style="background:#3f8d7d"></i>`
      + `<i style="background:#c08a2e"></i><i style="background:#e8eeeb"></i></span></th>`
    + `<th>Strand</th><th>E-value</th>${genomic ? "<th>Position</th>" : ""}`
    + `</tr></thead><tbody>${body}</tbody></table></div>`
    + (shown.length < groups.length
        ? `<p class="hint">Showing the 5 strongest of ${groups.length} places this sequence `
          + `aligns. <button type="button" class="link" data-expand-all="1">Show all</button></p>`
        : "")
    + `</td></tr>`;
}

/* ---------------------------------------------------------------- the genomic picture */

/** Every genomic hit against its signed distance to the nearest annotated start site.
 *
 * No threshold line is drawn and nothing is classified as "a promoter hit". The only band is
 * the window this run designed against, which is a statement about the run (ADR 0019).
 */
function tssPlot(entries) {
  const placed = entries.flatMap(e => e.placements.filter(p => p.distance_to_tss !== null)
    .map(p => ({...p, candidate: e.row.candidate})));
  if (!placed.length) {
    return `<p class="hint">No genomic hit could be placed against an annotated start site, so `
      + `there is nothing to plot. Distance is reported only where a transcript was found `
      + `within 10 kb.</p>`;
  }
  const width = 960, height = 148, left = 52, right = width - 52;
  const span = Math.max(1000, ...placed.map(p => Math.abs(p.distance_to_tss)));
  // Symmetric log: distances run from tens of bases to tens of kilobases, and a linear axis
  // would pile every near hit onto one pixel at the centre.
  const scale = (d) => {
    const mid = (left + right) / 2, half = (right - left) / 2;
    const t = Math.sign(d) * (Math.log10(Math.abs(d) + 1) / Math.log10(span + 1));
    return mid + t * half;
  };
  // From the screening record, which states the window it measured against, rather than from
  // the run: the flag beside each hit is meaningless without the window that produced it.
  const window = (held => held?.searches.map(s => s.window).find(Boolean) || null)(
    screenings()[state.screeningScope]);
  const parts = [];
  if (window) {
    const x1 = scale(window[0]), x2 = scale(window[1]);
    parts.push(`<rect x="${Math.min(x1, x2)}" y="22" width="${Math.abs(x2 - x1)}" height="74" `
      + `fill="#eef5f1" ${otTip(`The window this run designed against: ${signed(window[0])} to `
        + `${signed(window[1])} relative to the TSS. A hit inside it is a sequence that would `
        + `have been a legal target for that other gene too.`)}/>`);
  }
  parts.push(`<line x1="${left}" x2="${right}" y1="96" y2="96" stroke="#dce5e3"/>`);
  parts.push(`<line x1="${scale(0)}" x2="${scale(0)}" y1="18" y2="100" stroke="#9fb3ae" `
    + `stroke-dasharray="3 3"/>`);
  parts.push(`<text x="${scale(0)}" y="113" text-anchor="middle" fill="#61767a" `
    + `font-size="10">TSS</text>`);
  for (const d of [-span, -span / 10, span / 10, span]) {
    const x = scale(d);
    if (x < left || x > right) continue;
    parts.push(`<text x="${x}" y="113" text-anchor="middle" fill="#61767a" font-size="10">`
      + `${fmt(Math.round(d))}</text>`);
  }
  parts.push(`<text x="${(left + right) / 2}" y="131" text-anchor="middle" fill="#61767a" `
    + `font-size="10">bp from the nearest annotated start site, upstream left</text>`);
  placed.forEach((p, i) => {
    // Vertical offset is spacing so overlapping marks stay countable. It means nothing, and
    // the caption says so rather than leaving a reader to find a pattern in it.
    const y = 34 + (i % 5) * 12;
    const strong = p.hit.identities >= Math.ceil(p.hit.alignment_length * 0.8);
    parts.push(`<circle cx="${scale(p.distance_to_tss).toFixed(1)}" cy="${y}" `
      + `r="${p.inside_run_window ? 5 : 4}" fill="${strong ? "#1d5c52" : "#9cc6ba"}" `
      + `${p.inside_run_window ? 'stroke="#c08a2e" stroke-width="1.5"' : ""} `
      + otTip(`${p.candidate} aligns ${p.hit.identities} of ${p.hit.alignment_length} bases on `
        + `${p.contig}:${fmt(p.position)}, ${Math.abs(p.distance_to_tss).toLocaleString()} bp `
        + `${p.distance_to_tss < 0 ? "upstream of" : "downstream of"} the start site of `
        + `${p.nearest_gene || p.nearest_transcript || "an unnamed transcript"}`
        + `${p.inside_run_window ? ", inside the window this run designed against" : ""}. `
        + `Proximity to an annotated start site does not establish that the promoter is active `
        + `or that the duplex binds it.`) + `/>`);
  });
  return `<div class="ot-plot"><h4>Where the genomic hits sit, relative to the nearest `
    + `annotated start site</h4>`
    + `<svg viewBox="0 0 ${width} ${height}" role="img" `
    + `aria-label="Genomic off-target hits by signed distance to the nearest annotated `
    + `transcription start site">${parts.join("")}</svg>`
    + `<p class="hint">Distance is on a log scale either side of the start site, so near hits `
    + `stay separable from far ones. Vertical position is spacing only and carries no meaning. `
    + `Darker marks reach Wang's 80%${window ? ", and the shaded band is the window this run "
        + "designed against" : ""}. Nothing here is classified as a promoter hit: this project `
    + `has not defended a definition of promoter proximity, and inventing one would turn a `
    + `convenience into a scientific claim.</p></div>`;
}

/* ---------------------------------------------------------------- events */

document.addEventListener("click", (event) => {
  const batch = event.target.closest("[data-batch]");
  if (batch) {
    state.screenBatch = state.screenBatch || {};
    state.screenBatch[state.screeningScope || "transcript"] = Number(batch.dataset.batch);
    renderBatchControls();
    return;
  }
  const mode = event.target.closest("[data-batch-mode]");
  if (mode) { startScreening(state.screeningScope || "transcript", mode.dataset.batchMode); return; }
  const tab = event.target.closest("[data-scope-tab]");
  if (tab) { renderOffTarget(tab.dataset.scopeTab); return; }
  const sort = event.target.closest("[data-screensort]");
  if (sort) { state.screeningSort = sort.dataset.screensort; renderOffTarget(); return; }
  if (event.target.closest("[data-expand-all]")) {
    state.screeningExpandAll = true; renderOffTarget(); return;
  }
  const row = event.target.closest("[data-screenrow]");
  if (row) {
    const id = row.dataset.screenrow;
    state.screeningOpen = state.screeningOpen === id ? null : id;
    state.screeningExpandAll = false;
    renderOffTarget();
  }
});

$("screen-transcript").addEventListener("click", () => startScreening("transcript"));
$("screen-genomic").addEventListener("click", () => startScreening("genomic"));

// This file is a second script, so a frozen archive has already run app.js and shown its run
// by the time this executes. Catch up rather than leave the panel blank: an archive of a
// screened run that showed no screening would be a record quietly missing its own evidence.
if (state.result) {
  for (const record of (frozen?.screenings || [])) mergeScreening(record);
  renderOffTarget(frozen?.screenings?.[0]?.scope || "transcript");
  renderScreeningStatus();
}
