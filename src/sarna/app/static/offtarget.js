/* Off-target screening: the panel, its batches, and where each alignment lands.
 *
 * Display only. Every threshold shown here is computed in Python and arrives in the record;
 * nothing in this file decides what counts as a strong match. Selection assessments come
 * from the separate, explicitly unvalidated Python project policy.
 *
 * The panel merges several screening records per scope, because screening now runs in batches
 * and a run can hold a first fifteen, a next fifteen and a re-screen months later. Each row
 * therefore carries the date of the search it came from rather than the panel claiming one.
 */
"use strict";

const SCOPE_BUTTON = {transcript:"screen-transcript", genomic:"screen-genomic"};
const SCOPE_NAME = {transcript:"Transcript matches (RefSeq RNA)",
                    genomic:"Genomic matches (GRCh38)"};
const BATCH_SIZES = [1, 5, 10, 15, 25];
const MATCH_COLOURS = {"19":"#81572b", "18":"#327968", "17":"#97bcb1"};


// One hue, darker is closer, per ADR 0022. Green to red would read as good to bad, and an
// alignment is neither: it is a measurement the researcher judges.
const BANDS = [
  {key:"exact", label:"exact", fill:"#1d5c52"},
  {key:"strong", label:"at the reporting threshold", fill:"#3f8d7d"},
  {key:"partial", label:"partial", fill:"#9cc6ba"},
  {key:"weak", label:"weak", fill:"#dde9e4"},
];

const otTip = (text) =>
  `data-tip="${escapeHTML(text)}" tabindex="0" aria-label="${escapeHTML(text)}"`;

function screenings() {
  state.screenings = state.screenings || {};
  return state.screenings;
}

const SCREENING_HTML = new WeakMap();
const SCREENING_GROUPS = new WeakMap();
const SCREENING_ANNOTATIONS = new Map();
const SCREENING_POLLS = new Map();

function screeningHTML(node, html) {
  if (!node || SCREENING_HTML.get(node) === html) return;
  node.innerHTML = html;
  SCREENING_HTML.set(node, html);
}

function localScreeningDelay(record) {
  const native = record.progress?.backend === "local" || record.result?.backend === "local";
  const requested = Number(record.poll_after_seconds || (native ? 1 : 5)) * 1000;
  if (native) return Math.max(1000, Number.isFinite(requested) ? requested : 1000, document.hidden ? 15000 : 0);
  // These requests only read the local store. NCBI's separate scheduler is unchanged.
  return Math.max(5000, Number.isFinite(requested) ? requested : 5000,
    record.progress?.phase === "waiting" ? 15000 : 0,
    document.hidden ? 30000 : 0);
}

function screeningLimit(entry) {
  if (!entry.row.alignments_capped) return "";
  const parameters = entry.search.parameters || {};
  const sequences = new Set(entry.hits.map(h => h.accession.split(".")[0])).size;
  if (parameters.max_hits && sequences >= parameters.max_hits) {
    return `NCBI returned the maximum ${parameters.max_hits} database sequences for this candidate. Additional matches may exist.`;
  }
  if (parameters.max_alignments_per_candidate && entry.row.hit_count >= parameters.max_alignments_per_candidate) {
    return `The app retained at most ${parameters.max_alignments_per_candidate} alignments for this candidate. The list may be incomplete.`;
  }
  return "This saved result reached a reporting limit. Its returned matches may be incomplete.";
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
    annotations: result.gene_annotations || {},
    backend: result.backend || "ncbi",
    reference: result.reference || null,
    search_policy: result.search_policy || null,
  };
  held.searches = [...held.searches, search];
  const hits = {}, placements = {};
  for (const hit of result.hits || []) (hits[hit.candidate] = hits[hit.candidate] || []).push(hit);
  for (const p of result.proximities || []) {
    (placements[p.hit.candidate] = placements[p.hit.candidate] || []).push(p);
  }
  held.byEvidence = held.byEvidence || {};
  const evidenceKey = search.backend === "local" ? `local:${search.reference?.release || ""}` : "ncbi";
  const evidence = held.byEvidence[evidenceKey] || {};
  for (const row of result.specificity || []) {
    const entry = {row, search, hits:hits[row.candidate] || [], placements:placements[row.candidate] || []};
    const earlier = evidence[row.candidate];
    if (!earlier || !(new Date(earlier.search.searched_at) > new Date(search.searched_at)))
      evidence[row.candidate] = entry;
    const previous = held.byCandidate[row.candidate];
    if (previous && new Date(previous.search.searched_at) > new Date(search.searched_at)) continue;
    held.byCandidate = {...held.byCandidate, [row.candidate]:entry};
  }
  held.byEvidence[evidenceKey] = evidence;
  screenings()[scope] = held;
  return held;
}

/** Every candidate id this run has screened in one scope, for "the next fifteen". */
function selectedEvidenceKey() {
  if (screeningBackend() !== "local") return "ncbi";
  return `local:${REFERENCE_STATUS?.installed?.release || REFERENCE_STATUS?.catalog_release || ""}`;
}

function currentScreeningJob(scope) {
  const job = state.screeningJobs?.[scope];
  return job && (frozen || (job.backend || job.progress?.backend || "ncbi") === screeningBackend()) ? job : null;
}

function alreadyScreened(scope) {
  const held = screenings()[scope];
  if (frozen) return new Set(Object.keys(held?.byCandidate || {}));
  const key = selectedEvidenceKey();
  if (held?.byEvidence) return new Set(Object.keys(held.byEvidence[key] || {}));
  return new Set(Object.entries(held?.byCandidate || {}).filter(([,entry]) => {
    const source = entry.search.backend === "local" ? `local:${entry.search.reference?.release || ""}` : "ncbi";
    return source === key;
  }).map(([id])=>id));
}

/* A paused worker is not a paused NCBI job. Check the saved RID while this run is open,
 * without creating another search or letting stale timers follow the user to a new run. */
const AUTO_SCREEN_CHECKS = new Map();
const AUTO_SCREEN_CHECK_MS = 60000;

function scheduleScreeningCheck(scope) {
  const job = currentScreeningJob(scope), runId = state.result?.app.run_id;
  if (frozen || !runId || !job?.id || job.status !== "waiting" || job.retrievalIssue) return;
  const key = `${runId}:${scope}:${job.id}`;
  if (AUTO_SCREEN_CHECKS.has(key)) return;
  const due = Date.now() + AUTO_SCREEN_CHECK_MS;
  const timer = setTimeout(() => {
    AUTO_SCREEN_CHECKS.delete(key);
    const current = state.screeningJobs?.[scope];
    if (state.result?.app.run_id !== runId || current?.id !== job.id
        || current.status !== "waiting" || current.retrievalIssue) return;
    resumeScreening(scope);
  }, AUTO_SCREEN_CHECK_MS);
  // Node-based state tests must not wait a real minute for a browser-only timer.
  if (typeof timer?.unref === "function") timer.unref();
  AUTO_SCREEN_CHECKS.set(key, {timer, due});
}

function cancelScreeningCheck(scope) {
  const job = state.screeningJobs?.[scope], runId = state.result?.app.run_id;
  const key = `${runId}:${scope}:${job?.id}`;
  const pending = AUTO_SCREEN_CHECKS.get(key);
  if (pending) clearTimeout(pending.timer);
  AUTO_SCREEN_CHECKS.delete(key);
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
  const backend = screeningBackend();
  // Opt-in, and the run is already readable: this adds evidence, it never gates the designs.
  const run = state.result;
  if (!run || frozen) { error("Screening needs a live run."); return; }
  if (screeningBackend() === "local" && !REFERENCE_STATUS?.ready) {
    $("local-reference-panel").open = true;
    if (typeof setWorkspaceTab === "function") setWorkspaceTab("screening");
    $("local-reference-panel").scrollIntoView({behavior:"smooth", block:"center"});
    error("Prepare the local human references before searching, or explicitly choose NCBI BLAST.");
    return;
  }
  state.screeningJobs = state.screeningJobs || {};
  if (["running", "waiting"].includes(currentScreeningJob(scope)?.status)) return;
  cancelScreeningCheck(scope);
  state.screeningScope = scope;
  renderOffTarget(scope);
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
  const chunkSize = backend === "local" ? 50 : 25;
  for (let i = 0; i < ids.length; i += chunkSize) chunks.push(ids.slice(i, i + chunkSize));
  state.screeningJobs[scope] = {backend, status:"running", stage:"Starting…", done:0, total:chunks.length, progress:backend === "local" ? {backend:"local", phase:"reading"} : null};
  renderScreeningStatus();
  for (const [index, chunk] of chunks.entries()) {
    try {
      const started = await api(`/api/runs/${run.app.run_id}/screenings`,
        {method:"POST", headers:{"Content-Type":"application/json"},
         body:JSON.stringify({scope, candidates:chunk, backend})});
      if (screeningBackend() !== backend) return;
      if (!started.id) throw new Error("The server did not return a screening id.");
      const label = chunks.length > 1 ? `batch ${index + 1} of ${chunks.length}` : why;
      const outcome = await pollScreening(started.id, run.app.run_id, scope, label, index, chunks.length, backend);
      if (outcome !== "complete") return;
    } catch (e) {
      if (state.result?.app.run_id !== run.app.run_id || screeningBackend() !== backend) return;
      error(e.message);
      state.screeningJobs[scope] = {backend, status:"error", stage:e.message};
      renderScreeningStatus();
      return;
    }
    if (state.result?.app.run_id !== run.app.run_id || screeningBackend() !== backend) return;
  }
  state.screeningJobs[scope] = {backend, status:"complete", stage:"Finished"};
  renderScreeningStatus();
}

/** Resolve once the screening finishes; the server waits NCBI's own interval, this catches up. */
function pollScreening(screeningId, runId, scope, label, index, total, backend = screeningBackend()) {
  const key = `${runId}:${scope}:${screeningId}`;
  if (SCREENING_POLLS.has(key)) return SCREENING_POLLS.get(key);
  const pending = new Promise((resolve, reject) => {
    const tick = async () => {
      if (state.result?.app.run_id !== runId) { resolve("changed-run"); return; }
      if (screeningBackend() !== backend) { resolve("changed-backend"); return; }
      try {
        const record = await api(`/api/screenings/${screeningId}`);
        if (state.result?.app.run_id !== runId) { resolve("changed-run"); return; }
      if (screeningBackend() !== backend) { resolve("changed-backend"); return; }
        if (record.status === "running") {
          state.screeningJobs[scope] = {backend,
            status:"running", stage:`${label}: ${record.stage || "working"}`,
            done:index, total, progress:record.progress || (record.result?.backend === "local" ? {backend:"local", phase:"reading"} : null), id:screeningId,
          };
          renderScreeningStatus();
          setTimeout(tick, localScreeningDelay(record));
          return;
        }
        if (record.status === "waiting") {
          // Not finished and not failed. NCBI still holds the search, so the interface offers
          // to go back for it rather than reporting a failure and losing the job.
          state.screeningJobs[scope] = {backend,
            status:"waiting", stage:record.error || "Still running at NCBI.",
            done:index, total, id:screeningId, progress:record.progress || null,
            retrievalIssue:record.stage === "Report retrieval needs retry",
          };
          renderScreeningStatus();
          resolve("waiting");
          return;
        }
        if (record.status === "error") {
          reject(new Error(record.error || "The screening failed."));
          return;
        }
        if (record.status !== "complete") throw new Error("Unexpected screening status; results were not interpreted.");
        mergeScreening(record);
        state.screeningJobs[scope] = {backend, status:"complete", stage:"Search finished"};
        renderOffTarget(scope);
        renderScreeningStatus();
        // The three other places that state what has been screened are rendered elsewhere, so
        // refresh them here rather than let them go on saying none was performed.
        if (typeof renderTopResult === "function") renderTopResult(state.result);
        if (state.selected && typeof selectCandidate === "function") selectCandidate(state.selected);
        resolve("complete");
      } catch (e) { reject(e); }
    };
    tick();
  });
  const watched = pending.finally(() => {
    if (SCREENING_POLLS.get(key) === watched) SCREENING_POLLS.delete(key);
  });
  SCREENING_POLLS.set(key, watched);
  return watched;
}

/** Go back to NCBI for a search that outran its wait, by the request id it already has. */
async function resumeScreening(scope) {
  const job = currentScreeningJob(scope);
  const runId = state.result?.app.run_id;
  if (!job?.id || job.status !== "waiting") return;
  cancelScreeningCheck(scope);
  state.screeningJobs[scope] = {...job, status:"running", stage:"Checking NCBI again"};
  renderScreeningStatus();
  try {
    await api(`/api/screenings/${job.id}/resume`, {
      method:"POST", headers:{"Content-Type":"application/json"}, body:JSON.stringify({}),
    });
    if (state.result?.app.run_id !== runId || screeningBackend() !== "ncbi") return;
    await pollScreening(job.id, runId, scope, "collecting", 0, 1);
  } catch (e) {
    if (state.result?.app.run_id !== runId || screeningBackend() !== "ncbi") return;
    if (["screening-already-collecting", "screening-already-complete"].includes(e.code)) {
      // Another browser tab claimed the same RID. Follow its worker instead of submitting
      // or falsely labelling the search failed.
      try {
        await pollScreening(job.id, runId, scope, "collecting", 0, 1);
      } catch (pollError) {
        if (state.result?.app.run_id !== runId || screeningBackend() !== "ncbi") return;
        state.screeningJobs[scope] = {...job, status:"error", stage:pollError.message};
      }
    } else {
      error(e.message);
      state.screeningJobs[scope] = {...job, status:"error", stage:e.message};
    }
  }
  renderScreeningStatus();
}

/** Rebuild the merged panel for a saved run, so screenings survive a reload. */
async function loadSavedScreenings(runId) {
  const engine = screeningBackend();
  const current = () => state.result?.app.run_id === runId && screeningBackend() === engine;
  try {
    const listed = await api(`/api/runs/${runId}/screenings`);
    if (!current()) return;
    let any = null;
    for (const entry of listed.filter(s => s.status === "complete").reverse()) {
      if (screenings()[entry.scope]?.searches.some(search=>search.id === entry.id)) { any = entry.scope; continue; }
      const record = await api(`/api/screenings/${entry.id}`);
      if (!current()) return;
      mergeScreening(record);
      any = entry.scope;
    }
    state.screeningJobs = state.screeningJobs || {};
    const active = listed.filter(s => ["running", "waiting"].includes(s.status));
    for (const entry of [...active, ...listed]) {
      if (!current()) return;
      if (currentScreeningJob(entry.scope)) continue;
      const saved = entry.status !== "complete"
        ? await api(`/api/screenings/${entry.id}`) : null;
      if (!current()) return;
      const backend = saved?.result?.backend || saved?.progress?.backend
        || screenings()[entry.scope]?.searches.find(search=>search.id === entry.id)?.backend || "ncbi";
      if (backend !== engine) continue;
      state.screeningJobs[entry.scope] = {backend,
        status:entry.status, id:entry.id, stage:entry.error || entry.stage,
        retrievalIssue:entry.stage === "Report retrieval needs retry",
      };
      if (entry.status === "running") {
        pollScreening(entry.id, runId, entry.scope, "restored search", 0, 1).catch(e => {
          if (!current()) return;
          state.screeningJobs[entry.scope] = {backend, status:"error", id:entry.id, stage:e.message};
          renderScreeningStatus();
        });
      }
    }
    renderOffTarget(any || state.screeningScope || "transcript");
    renderScreeningStatus();
    if (typeof renderTopResult === "function") renderTopResult(state.result);
  } catch (e) {
    if (!current()) return;
    $("offtarget-progress").textContent = `Saved screening history could not be loaded: ${e.message}`;
    $("offtarget-status").textContent = "HISTORY UNAVAILABLE";
  }
}

/* ---------------------------------------------------------------- the controls */

/* The bar. NCBI gives no completion percentage and never has: a search is WAITING or READY,
 * and the only forecast is RTOE, its own estimate at submission. So the bar shows a true
 * fraction where one exists (placing hits, and batches finished) and otherwise shows elapsed
 * against that estimate, labelled as an estimate and never allowed past nine tenths. Once the
 * estimate is passed it stops pretending and says so. A bar that sat at 99% for twenty
 * minutes would be worse than no bar. */
const PHASES = [
  {key:"submitting", label:"Sending to NCBI"},
  {key:"waiting", label:"NCBI is searching"},
  {key:"reading", label:"Reading alignments"},
  {key:"placing", label:"Placing hits"},
];

function phaseFraction(progress) {
  if (!progress) return {fraction:null, note:""};
  if (progress.phase === "placing" && progress.to_place) {
    return {fraction: progress.placed / progress.to_place,
            note: `${progress.placed} of ${progress.to_place}`};
  }
  if (progress.phase === "waiting") {
    const elapsed = Number(progress.elapsed_seconds || 0);
    const estimate = Number(progress.estimate_seconds || 0);
    return {fraction: null,
            note: estimate > 0
              ? `${elapsed}s, longer than the ${estimate}s NCBI estimated`
              : `${elapsed}s elapsed, with no estimate from NCBI`};
  }
  return {fraction: null, note: ""};
}

function progressBar(job) {
  if (!job || job.status !== "running") return "";
  const progress = job.progress || null;
  const phase = progress?.phase || "submitting";
  if (progress?.backend === "local" || job.backend === "local") return `<p class="hint">${escapeHTML(job.stage || "Searching the local reference…")} · sequences remain on this computer</p><span class="bar indeterminate"><i></i></span>`;
  const reached = Math.max(0, PHASES.findIndex(p => p.key === phase));
  const {fraction, note} = phaseFraction(progress);
  const batch = job.total > 1 ? `Batch ${(job.done || 0) + 1} of ${job.total} · ` : "";
  const label = PHASES[reached]?.label || "Working";
  const bar = fraction === null
    ? `<span class="bar indeterminate"><i></i></span>`
    : `<span class="bar"><i style="width:${Math.round(fraction * 100)}%"></i></span>`;
  return `<div class="screen-progress" role="group" aria-label="Screening progress">`
    + `<div class="phase-row">`
    + PHASES.map((p, index) => `<span class="phase${index < reached ? " done" : ""}`
        + `${index === reached ? " now" : ""}">${escapeHTML(p.label)}</span>`).join("")
    + `</div>${bar}`
    + `<p class="hint">${escapeHTML(batch + label)}`
    + `${note ? ` · ${escapeHTML(note)}` : ""}`
    + `${fraction === null && phase === "waiting"
        ? ". NCBI reports only started or finished, so there is no percentage to show."
        : ""}</p></div>`;
}

function renderScreeningStatus() {
  renderScreenedShortlist();
  const jobs = Object.fromEntries(Object.keys(SCOPE_BUTTON).map(scope=>[scope,currentScreeningJob(scope)]).filter(([,job])=>job));
  const states = Object.values(jobs).map(j => j.status);
  const done = Object.keys(screenings()).length;
  const pool = screenablePool();
  const coverage = Object.keys(SCOPE_BUTTON).map(scope => {
    const count = alreadyScreened(scope).size;
    const job = jobs[scope];
    const status = job?.status === "running" ? "search running"
      : job?.status === "waiting" ? "waiting at NCBI"
      : job?.status === "error" ? "search failed" : count ? "results available" : "not started";
    return `<div class="coverage-item"><span>${scope === "transcript" ? "Transcript" : "Genomic"} search · ${frozen ? "saved" : screeningBackend() === "local" ? "Local" : "NCBI"}</span>`
      + `<strong aria-label="${count} of ${pool.length} candidates have completed results">${count} / ${pool.length}</strong>`
      + `<span>Completed results · ${escapeHTML(status)}</span></div>`;
  });
  const coverageHost = $("screening-coverage");
  screeningHTML(coverageHost, coverage.join(""));
  $("offtarget-status").textContent = states.includes("running") ? "RUNNING"
    : states.includes("waiting") ? "AT NCBI"
    : states.includes("error") ? (done ? "PARTIAL" : "FAILED")
    : done ? "RESULTS AVAILABLE" : "NOT RUN";
  $("offtarget-progress").textContent = states.includes("running")
    ? "Searches run in the background and report here when they finish. You can keep working meanwhile."
    : "";
  for (const scope of Object.keys(SCOPE_BUTTON)) {
    const lane = $(`screen-state-${scope}`), job = jobs[scope];
    if (lane) {
      lane.textContent = job?.status === "waiting"
        ? (job.retrievalIssue ? "Report retrieval needs retry" : "Waiting for NCBI results")
        : job ? (job.stage || job.status) : "";
      lane.className = "screen-state " + (job ? job.status : "");
    }
    const bar = $(`screen-bar-${scope}`);
    if (bar) {
      if (job?.status === "waiting") scheduleScreeningCheck(scope);
      screeningHTML(bar, progressBar(job)
        + (job?.status === "waiting"
          ? `<div class="screen-wait-note" role="status"><p>${job.retrievalIssue
              ? "The downloaded report could not be read. Retry collection of the saved search."
              : "NCBI has not finished this search. The app will check again automatically in about a minute while this run is open."}</p>`
            + `<p class="hint">${job.retrievalIssue
              ? "Your saved request ID is retained; this retries retrieval rather than submitting a new search."
              : "You can keep working. If you close this run, checking resumes when you reopen it."} `
            + `<button type="button" class="link" data-resume="${scope}">Check NCBI again</button></p>`
            + `<details><summary>Search details</summary><p>${escapeHTML(job.stage || "")}</p>`
            + `<p>Browser status requests stay on this computer (every 5–15 seconds while visible). NCBI status checks are separately limited to once a minute per search.</p>`
            + `<p>NCBI keeps request IDs for about a day; an expired search must be submitted again.</p></details></div>`
          : ""));
    }
    const button = $(SCOPE_BUTTON[scope]);
    if (button) {
      const remaining = screenablePool().filter(id => !alreadyScreened(scope).has(id)).length;
      const amount = Math.min(batchChoice(scope), remaining);
      const name = scope === "transcript" ? "transcript" : "genomic";
      if (button.firstChild) button.firstChild.textContent = ["running", "waiting"].includes(job?.status)
        ? `${scope === "transcript" ? "Transcript" : "Genomic"} search pending `
        : amount ? `Search next ${amount} ${name} candidate${amount === 1 ? "" : "s"} `
        : `All ${name} candidates searched `;
      button.disabled = ["running", "waiting"].includes(job?.status) || !remaining
        || (screeningBackend() === "local" && !REFERENCE_STATUS?.ready);
    }
  }
  renderBatchControls();
  if (state.result && typeof renderExportNote === "function") renderExportNote();
}

const SCOPE_WARNING = {
  genomic: "Genomic searches can take substantially longer than transcript searches. NCBI controls "
    + "the queue; its estimate is not a completion time. You can return to this saved run later.",
};

function renderBatchControls() {
  const host = $("offtarget-batch");
  if (!host) return;
  const pool = screenablePool();
  if (!pool.length || frozen) { screeningHTML(host, ""); return; }
  const scope = state.screeningScope || "transcript";
  const done = alreadyScreened(scope), remaining = pool.length - done.size;
  const size = batchChoice(scope);
  const busy = ["running", "waiting"].includes(currentScreeningJob(scope)?.status);
  const disabled = busy || (screeningBackend() === "local" && !REFERENCE_STATUS?.ready) ? " disabled" : "";
  const ticked = pool.filter(id => state.checked.has(id)).length;
  screeningHTML(host,
    `<div class="batch-row"><span class="batch-label">${scope === "genomic" ? "Genomic" : "Transcript"} batch size</span>`
    + `<span class="segmented">`
    + BATCH_SIZES.map(n => `<button type="button" data-batch="${n}"`
        + `${n === size ? ' class="active"' : ""} aria-pressed="${n === size}">${n}</button>`).join("")
    + `</span>`
    + `<button type="button" data-batch-mode="all" class="secondary"${disabled}`
    + ` ${otTip("Screen every candidate that passed the gates. Local jobs cover up to 50 candidates; "
        + "remote NCBI jobs cover up to 25. Each job retains a dated record.")}>`
    + `${done.size ? "Re-screen all" : "Search all"} ${pool.length}</button>`
    + (done.size && remaining > 0
        ? `<button type="button" data-batch-mode="rest" class="secondary"${disabled}>`
          + `Search remaining ${remaining}</button>` : "")
    + (ticked
        ? `<button type="button" data-batch-mode="selection" class="secondary"${disabled}>`
          + `Search selected (${ticked})</button>` : "")
    + `</div>`
    + ((screeningBackend() === "ncbi" && SCOPE_WARNING[scope]) ? `<p class="hint">${escapeHTML(SCOPE_WARNING[scope])}</p>` : "")
    + `<p class="hint">${busy
        ? "A submitted search is pending. No new search is needed; completed results will appear here."
        : !done.size
        ? `Nothing screened with the selected engine/reference yet in this scope. Start with the candidates you plan to review. `
          + (screeningBackend() === "local" ? "Local searches reuse the prepared index and completed alignments." : "NCBI queue time varies; a smaller batch does not guarantee a faster result.")
        : remaining > 0
        ? `${done.size} of ${pool.length} screened in this scope. The next search covers `
          + `${Math.min(size, remaining)} more unless you choose otherwise.`
        : `All ${pool.length} screened in this scope. Screening again searches the same `
          + (screeningBackend() === "local" ? "sequences against the installed, pinned reference as a second dated record; cached alignments can be reused." : "sequences against whatever the database holds today, as a second dated record.")
      }</p>`);
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
  if (entry.row.enumeration_complete) return {exact:entry.row.perfect_hits, strong:entry.row.wang_homology_hits-entry.row.perfect_hits, partial:entry.row.hit_count-entry.row.wang_homology_hits, weak:0};
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
  state.screeningScope = scope || state.screeningScope || "genomic";
  const genomic=state.screeningScope === "genomic";
  const choices=screeningSortChoices(19,genomic);
  if(!choices.some(([key])=>key === (state.screeningSort || "design_rank"))) {
    state.screeningSort="design_rank";state.screeningSortDescending=false;
  }
  const host=$("offtarget-results");
  if(!host) return;
  const decisions=screeningComparisonDecisions();
  renderScreeningComparison(decisions);
  const held=screenings()[state.screeningScope];
  const core=state.result?.candidates?.[0]?.target_site_dna?.length || 19;
  const entries=decisions.map(d=>comparisonEntry(state.screeningScope,d.candidate.id)
    || {row:{candidate:d.candidate.id,enumeration_complete:false},hits:[],placements:[],missing:true});
  const actual=entries.filter(e=>!e.missing);
  const shownHeld=held ? {...held,byCandidate:Object.fromEntries(actual.map(e=>[e.row.candidate,e])),
    searches:[...new Map(actual.map(e=>[e.search.id,e.search])).values()]} : null;
  const tableOpen=host.querySelector?.(".screening-detail-table")?.open;
  const plotOpen=host.querySelector?.("details.ot-plot")?.open;
  const methodsOpen=host.querySelector?.(".ot-methods")?.open;
  const mismatchOpen=host.querySelector?.('[data-comparison-plot="mismatches"]')?.open;
  const notesOpen=host.querySelector?.(".ot-record-notes")?.open;
  host.innerHTML=`<details class="screening-detail-table" ${tableOpen ? "open" : ""}><summary>Detailed sequence matches · ${genomic ? "genome" : "transcripts"} · ${actual.length} screened / ${entries.length} displayed candidates</summary>`
    + scopeTabs()
    + (held && actual.length ? summaryLine(shownHeld,actual,core,genomic) : '<p class="hint">No completed results in this scope for the displayed candidates.</p>')
    + resultsTable(held,entries,core,genomic)+`</details>`
    + (shownHeld?.searches.length ? notesFor(shownHeld) : "")
    + mismatchPositionPlot(decisions)
    + (genomic && actual.length ? tssPlot(actual) : "");
  if(plotOpen && host.querySelector("details.ot-plot")) host.querySelector("details.ot-plot").open=true;
  if(methodsOpen && host.querySelector(".ot-methods")) host.querySelector(".ot-methods").open=true;
  if(notesOpen && host.querySelector(".ot-record-notes")) host.querySelector(".ot-record-notes").open=true;
  if(mismatchOpen && host.querySelector('[data-comparison-plot="mismatches"]')) host.querySelector('[data-comparison-plot="mismatches"]').open=true;
  renderBatchControls();
}

function scopeTabs() {
  const present = Object.keys(SCOPE_BUTTON);
  return `<div class="segmented scope-tabs">` + present.map(s =>
    `<button type="button" data-scope-tab="${s}"`
    + `${s === state.screeningScope ? ' class="active"' : ""}>`
    + `${escapeHTML(SCOPE_NAME[s] || s)}</button>`).join("") + `</div>`;
}

function summaryLine(held, entries, core, genomic) {
  const sources = [...new Map(entries.map(e => [e.search.id, e.search])).values()];
  const localSources = sources.filter(s => s.backend === "local");
  const remoteSources = sources.filter(s => s.backend !== "local");
  const hasLocal = localSources.length > 0, hasRemote = remoteSources.length > 0;
  const exact = entries.filter(e => e.row.perfect_hits > 0).length;
  const quiet = entries.filter(e => e.row.hit_count === 0).length;
  const threshold = entries[0]?.row.strong_threshold || Math.ceil(core * 0.8);
  const capped = entries.filter(e => e.row.alignments_capped).length;
  const dates = sources.map(s => new Date(s.searched_at)).filter(d => !isNaN(d));
  const when = dates.length ? new Date(Math.max(...dates)).toLocaleDateString() : "date unavailable";
  const databases = [...new Set(sources.map(s => s.parameters.database))].join(", ");
  const local = sources.every(s => s.backend === "local");
  const byDatabase = sources.every(s => s.parameters.species_restriction === "database");
  const returned = entries.reduce((total,e) => total + e.row.hit_count, 0);
  const matching = entries.reduce((total,e) => total + (e.row.wang_homology_hits || 0), 0);
  const pool = screenablePool().length;
  const limit = capped === 1 ? screeningLimit(entries.find(e=>e.row.alignments_capped))
    : capped ? `${capped} candidate result lists may be incomplete. Each candidate shows its search limit.` : "";
  return `<div class="ot-summary"><p><b>${escapeHTML(SCOPE_NAME[state.screeningScope])}</b> · `
    + `latest saved ${escapeHTML(when)}</p>`
    + (!frozen && entries.some(e => (e.search.backend === "local" ? `local:${e.search.reference?.release || ""}` : "ncbi") !== selectedEvidenceKey()) ? `<p class="hint">Saved findings include another engine or reference. Coverage counters refer to the selected engine/reference; each row names the method used.</p>` : "")
    + (hasLocal && hasRemote ? `<p class="hint">Mixed search methods: some candidates have local results and others have NCBI BLAST results. Each candidate names its method; their search coverage differs.</p>` : "")
    + `<div class="ot-result-flow" aria-label="Search coverage and returned sequence matches">`
    + `<div><strong>${entries.length}<small> / ${pool}</small></strong><span>saRNA candidates searched</span>`
    + `<small>${Math.max(0,pool-entries.length)} not yet searched in this scope</small></div>`
    + `<div><strong>${returned}</strong><span>alignments returned</span><small>${local ? "Full-length alignments with 0–3 mismatches" : "Includes weaker sequence matches"}</small></div>`
    + `<div><strong>${matching}</strong><span>alignments with ≥${threshold}/${core} matching bases</span>`
    + `<small>A reporting threshold, not a biological risk score</small></div></div>`
    + `<p class="hint">${exact} candidate${exact === 1 ? "" : "s"} with full-length exact matches · `
    + `${quiet} with no reported alignments. Binding, gene-expression effects and specificity remain untested.</p>`
    + (limit ? `<p class="ot-limit"><b>Result list may be incomplete.</b> ${escapeHTML(limit)}</p>` : "")
    + `<details class="ot-methods"><summary>Methods, evidence and search limits</summary>`
    + (hasLocal ? `<p><b>Local Bowtie search.</b> Full-length 19-base alignments, 0–3 mismatches, both strands, no gaps. All eligible alignments are counted; up to 200 details per candidate are stored. Shorter/seed-only matches, bulges, individual variants and biological activity are not assessed. Alternate and patch contigs can represent overlapping loci.</p><p>Reference: ${escapeHTML(databases)}. Tool: ${escapeHTML(localSources.map(s=>"Bowtie " + s.reference?.bowtie_version).join(", "))}. Nearest eligible transcript start within 10 kb uses the pinned RefSeq annotation. These TSS proxies are annotated 5′ boundaries of non-partial mRNA, lncRNA, transcript or primary-transcript models, not measured cell-specific starts. Mature small-RNA feature starts are excluded. Sources, checksums and parameters are included in the JSON record. <a href="https://bowtie-bio.sourceforge.net/manual.shtml" target="_blank" rel="noopener">Bowtie manual</a> · <a href="https://www.ncbi.nlm.nih.gov/refseq/" target="_blank" rel="noopener">NCBI RefSeq</a>.</p>` : "")
    + `<p>Database: ${escapeHTML(databases)}. ${byDatabase
      ? (local ? "The reference contains human sequence." : "The genomic database contains human GRCh38 sequence.")
      : "Human filtering was requested using txid9606. The Entrez filter is not listed in the current supported URL API; species restriction is not independently verified."}</p>`
    + `<p>The ≥${threshold}/${core} column counts matching bases${local ? " across the full core" : ", including partial or gapped alignments in remote BLAST results"}. `
    + `It is a reporting threshold, not a validated saRNA risk cutoff. Applying 80% to candidate `
    + `off-targets is a project interpretation; Wang 2015 discusses it for control duplexes. `
    + `<a href="https://doi.org/10.14440/jbm.2015.39" target="_blank" rel="noopener">Wang et al., 2015</a>.</p>`
    + (hasRemote ? `<p>Remote BLAST counts describe retained alignments, not affected genes. The server returns at most `
    + `${remoteSources[0]?.parameters.max_hits || 100} database sequences per query, and the client retains `
    + `at most ${remoteSources[0]?.parameters.max_alignments_per_candidate || 200} alignments. `
    + `A limit reached means counts may be lower bounds. A negative result is not proof of specificity.</p>` : "")
    + `<p>This is a sequence-similarity search, not a comprehensive seed-risk screen. Seed-mediated `
    + `off-target silencing has experimental evidence in siRNA; that evidence does not validate `
    + `prediction of saRNA activation. <a href="https://doi.org/10.1261/rna.25706" target="_blank" `
    + `rel="noopener">Jackson et al., 2006</a>. No biological risk score is assigned.</p>`
    + (hasLocal ? `<p>Local transcript chromosome context comes from accession mappings in the pinned RefSeq GFF. Missing mappings are marked unavailable; no live annotation requests are needed.</p>` : "")
    + (hasRemote ? `<p>Remote transcript chromosomes are optional Ensembl gene-symbol annotations, not mappings of transcript alignment coordinates to the genome.</p>` : "")
    + `<p>Descriptions are source annotations, including predicted functions, rather than functional validation.</p></details></div>`;
}

function notesFor(held) {
  const current = [...new Map(Object.values(held.byCandidate).map(e=>[e.search.id,e.search])).values()];
  const notes = [...new Set(current.flatMap(s => s.notes))];
  const important = notes.filter(n => /unavailable|not placed|weaker ones/i.test(n) && !/reporting limit|floors/i.test(n));
  const other = notes.filter(n => !important.includes(n));
  return important.map(n => `<p class="hint warn">${escapeHTML(n)}</p>`).join("")
    + (other.length ? `<details class="ot-record-notes"><summary>Search record notes</summary>`
      + other.map(n => `<p class="hint">${escapeHTML(n)}</p>`).join("") + `</details>` : "");
}

function screeningTableOrder(entries, sort="design_rank", descending=false) {
  const candidates = new Map((state.result?.candidates || []).map(c=>[c.id,c]));
  const rank = e => candidates.get(e.row.candidate)?.rank ?? Infinity;
  const value = e => {
    const v = sort === "design_rank" ? rank(e)
      : sort === "distinct_loci" ? (e.row.distinct_loci ?? distinctLoci(e.hits))
      : sort === "distinct_sequences" ? (e.row.distinct_sequences ?? new Set(e.hits.map(h=>h.accession.split(".")[0])).size)
      : e.row[sort];
    return Number.isFinite(v) && v >= 0 ? v : null;
  };
  return [...entries].sort((a,b)=>{
    // Incomplete or missing counts must not appear to be the fewest-match findings.
    if (sort !== "design_rank") {
      const incomplete = e => !!e.row.alignments_capped || e.row.enumeration_complete === false;
      if (incomplete(a) !== incomplete(b)) return incomplete(a) ? 1 : -1;
    }
    const av=value(a), bv=value(b);
    if ((av===null) !== (bv===null)) return av===null ? 1 : -1;
    if (av!==bv) return (descending ? -1 : 1)*(av-bv);
    if (rank(a)!==rank(b)) return rank(a)-rank(b);
    return a.row.candidate.localeCompare(b.row.candidate,undefined,{numeric:true});
  });
}
function resultsTable(held, entries, core, genomic) {
  const threshold = entries[0]?.row.strong_threshold || Math.ceil(core * 0.8);
  const sort = state.screeningSort || "design_rank", descending = !!state.screeningSortDescending;
  const ordered = screeningTableOrder(entries,sort,descending);
  const header = (key, label, tip) =>
    `<th${sort === key ? ` aria-sort="${descending ? "descending" : "ascending"}"` : ""}><button data-screensort="${key}"${sort === key ? ' class="active"' : ""} `
    + `data-tip="${escapeHTML(tip)}" tabindex="0">${label} ${sort === key ? descending ? "↓" : "↑" : "↕"}</button></th>`;
  return `<div class="table-scroll"><table class="ot-table"><thead><tr>`
    + `<th>Candidate</th>`
    + header("design_rank", "Design rank", "The original design rank shown in Candidates. Screening does not change it.")
    + `<th>Design score</th><th>Strongest match</th><th>Alignments by identity</th>`
    + header("perfect_hits", "Exact alignments",
             `Full-length exact matches elsewhere, all ${core} bases with no gaps. One of these `
             + `is a reported sequence match; biological targeting is untested.`)
    + header("wang_homology_hits", `Alignments ≥${threshold}/${core}`,
             `Retained alignments matching at least ${threshold} of ${core} bases. A reporting `
             + `threshold, not a validated biological risk cutoff.`)
    + header(genomic ? "distinct_loci" : "distinct_sequences", genomic ? "Distinct loci" : "Transcripts",
             "Distinct reported transcript accessions; genomic rows show distinct alignment loci. Not a count of affected genes.")
    + header("hit_count", "Returned alignments", "How many alignments were returned at all.")
    + (genomic ? header("in_run_window_hits", "In your window",
         "Hits inside the same window this run designed against, near another gene's start "
         + "site. A statement about this run's configuration, not a definition of a promoter.")
       : "")
    + `</tr></thead><tbody>`
    + ordered.map(e => e.missing ? missingScreeningRow(e.row.candidate,genomic) : candidateRow(e, core, genomic)).join("")
    + `</tbody></table></div>`;
}

function missingScreeningRow(id,genomic) {
  const c=state.result?.candidates.find(c=>c.id===id);
  return `<tr data-screenrow="${escapeHTML(id)}" data-missing-evidence="true"><td>${escapeHTML(id)}</td><td>${c?.rank ?? "Not available"}</td><td>${c?.rating?.met ?? "?"} / ${c?.rating?.evaluated ?? "?"}</td><td colspan="${genomic ? 7 : 6}">Not screened against this reference in this scope. Missing evidence is not zero.</td></tr>`;
}

function candidateRow(entry, core, genomic) {
  const row = entry.row, counts = bandCounts(entry, core);
  const candidate = state.result?.candidates?.find(c=>c.id===row.candidate);
  const rating = candidate?.rating;
  const score = Number.isInteger(rating?.met) && Number.isInteger(rating?.evaluated) && rating.evaluated > 0
    ? `${rating.met} / ${rating.evaluated}` : "Not available";
  const best = worstHit(entry);
  const placement = best
    ? entry.placements.find(p => p.hit.accession === best.accession && p.hit.subject_start === best.subject_start
      && p.hit.subject_end === best.subject_end && p.hit.subject_strand === best.subject_strand) : null;
  const open = state.screeningOpen === row.candidate;
  const flag = (n) => n ? ' class="match-count"' : "";
  const main = `<tr data-screenrow="${escapeHTML(row.candidate)}"`
    + ` tabindex="0" role="button" aria-expanded="${open}"${open ? ' class="open"' : ""}>`
    + `<td>${escapeHTML(row.candidate)}`
    + `<p class="candidate-id" ${otTip(entry.search.parameters.database || "Reference unavailable")}>${entry.search.backend === "local" ? "Local · Bowtie" : "Remote · NCBI BLAST"}</p>`
    + (row.alignments_capped ? `<p class="candidate-id" ${otTip(
        "A server or client reporting limit was reached. Counts may be lower bounds.")
      }>list may be incomplete</p>` : "")
    + (row.on_target_hits ? `<p class="candidate-id" ${otTip(
        `${row.on_target_hits} alignment at this candidate's own target locus was removed from `
        + `these findings. A sequence matching the promoter it was designed against is the `
        + `intended target, not an off-target one.`)}>own locus removed</p>` : "")
    + (row.hit_count === 0 ? '<p class="candidate-id">No reported alignments; specificity untested</p>' : "")
    + `</td><td data-label="Design rank">${candidate?.rank ?? "Not available"}</td>`
    + `<td data-label="Design score"><strong>${score}</strong><p class="candidate-id">preferences met / evaluated</p></td>`
    + `<td class="ot-strongest">${best
        ? `<b>${escapeHTML(hitName(best, placement))}</b> <span class="muted">`
          + `${best.identities}/${core}</span>`
          + (placement && placement.distance_to_tss !== null
             ? `<p class="candidate-id">${signed(placement.distance_to_tss)} bp from its TSS</p>`
             : "")
        : `<span class="muted">none</span>`}</td>`
    + `<td class="ot-distribution">${identityStrip(counts, row.hit_count)}`
    + (row.promiscuous ? `<p class="candidate-id">multiple reported alignments</p>` : "")
    + `</td>`
    + `<td data-label="Exact matches"${flag(row.perfect_hits)}>${row.perfect_hits}</td>`
    + `<td data-label="Alignments at reporting threshold"${flag(row.wang_homology_hits)}>${row.wang_homology_hits}</td>`
    + `<td data-label="Distinct ${genomic ? "loci" : "transcripts"}">${genomic ? (row.distinct_loci ?? distinctLoci(entry.hits))
        : (row.distinct_sequences ?? new Set(entry.hits.map(h => h.accession.split(".")[0])).size)}${!genomic && row.named_genes !== undefined ? `<p class="candidate-id">${row.named_genes} named genes${row.unnamed_alignments ? " + unnamed matches" : ""}</p>` : ""}</td>`
    + `<td data-label="Returned alignments">${row.hit_count}</td>`
    + (genomic ? `<td data-label="In your window"${flag(row.in_run_window_hits)}>${row.in_run_window_hits === null ? "Unknown" : row.in_run_window_hits ?? "Unknown"}</td>` : "")
    + `</tr>`;
  return main + (open ? alignmentDetail(entry, core, genomic) : "");
}

/* ---------------------------------------------------------------- one candidate's hits */

/** Group a candidate's alignments by what they landed in, not one row per HSP. */
function groupHits(entry) {
  const genomic = state.screeningScope === "genomic";
  const cached = SCREENING_GROUPS.get(entry);
  if (cached?.genomic === genomic) return cached.groups;
  const placementKey = h => `${h.accession}:${h.subject_start}:${h.subject_end}:${h.subject_strand}`;
  const placements = new Map(entry.placements.map(p => [placementKey(p.hit), p]));
  const groups = new Map();
  for (const hit of entry.hits) {
    const placement = placements.get(placementKey(hit));
    const key = state.screeningScope === "genomic"
      ? `${hit.accession}:${hit.subject_start}:${hit.subject_end}:${hit.subject_strand}` : hit.accession;
    const held = groups.get(key) || {accession:hit.accession, title:hit.title, hits:[], placement};
    held.hits = [...held.hits, hit];
    held.placement = held.placement || placement;
    groups.set(key, held);
  }
  const ordered = [...groups.values()].sort((a, b) =>
    Math.max(...b.hits.map(h => h.identities)) - Math.max(...a.hits.map(h => h.identities)));
  SCREENING_GROUPS.set(entry, {genomic, groups:ordered});
  return ordered;
}

function distinctLoci(hits) {
  return new Set(hits.map(h => `${h.accession}:${Math.min(h.subject_start,h.subject_end)}:`
    + `${Math.max(h.subject_start,h.subject_end)}:${h.subject_strand}`)).size;
}

/** Show the actual aligned DNA, including insertions; never invent bases in old reports. */
function footprint(hit, core) {
  const candidate = state.result?.candidates?.find(c => c.id === hit.candidate);
  const full = candidate?.target_site_dna || "";
  const from = Math.min(hit.query_start, hit.query_end) - 1;
  const to = Math.max(hit.query_start, hit.query_end);
  const rawQuery = hit.query_sequence || "", rawSubject = hit.subject_sequence || "";
  const exactBases = rawQuery && rawSubject && rawQuery.length === rawSubject.length;
  const columns = [];
  for (let i=0; i<Math.max(0,from); i++) columns.push({q:full[i] || "?", s:"·", kind:"outside", pos:i+1});
  let pos = from;
  if (exactBases) {
    for (let i=0; i<rawQuery.length; i++) {
      const q = rawQuery[i].toUpperCase(), subject = rawSubject[i].toUpperCase();
      if (q !== "-") pos++;
      columns.push({q, s:subject, pos:q === "-" ? "+" : pos,
        kind:q === "-" || subject === "-" ? "gap" : q === subject ? "match" : "mismatch"});
    }
  } else {
    const mismatches = new Set(hit.mismatch_positions || []);
    for(let i=Math.max(0,from); i<to; i++) columns.push({q:full[i] || "?", s:"?", pos:i+1,
      kind:hit.has_alignment_detail ? (mismatches.has(i+1) ? "mismatch" : "match") : "unknown"});
  }
  for(let i=Math.max(to,0); i<core; i++) columns.push({q:full[i] || "?", s:"·", kind:"outside", pos:i+1});
  const cell = (c, value) => `<span class="fp-base ${c.kind}" title="${escapeHTML(`Position ${c.pos}: ${c.q} / ${c.s}; ${c.kind}`)}">${escapeHTML(String(value))}</span>`;
  return `<div class="footprint-wrap"><div class="footprint" role="group" aria-label="DNA sequence alignment, query and subject in BLAST orientation">`
    + `<div class="fp-row fp-position"><span>Position</span>${columns.map(c => cell(c,c.pos)).join("")}</div>`
    + `<div class="fp-row"><span>Query 5′</span>${columns.map(c => cell(c,c.q)).join("")}<small>3′</small></div>`
    + `<div class="fp-row"><span>Subject</span>${columns.map(c => cell(c,c.s)).join("")}</div></div>`
    + `<p class="fp-note">DNA alphabet: T corresponds to U in RNA. Subject shown in query alignment orientation `
    + `(${escapeHTML(hit.subject_strand)} strand). ${exactBases ? "− gap · outside alignment." : "? aligned bases unavailable in this saved report; rerun to retrieve them."}</p></div>`;
}

async function loadHitAnnotations(entry, accessions) {
  if (frozen) return;
  const runId = state.result?.app.run_id;
  const wanted = [...new Set(accessions)].filter(a => !entry.search.annotations[a]);
  const key = accession => `${runId}:${entry.search.id}:${accession}`;
  for (let index=0; index<wanted.length; index+=5) {
    if (state.result?.app.run_id !== runId || state.screeningOpen !== entry.row.candidate
        || state.screeningScope !== "transcript") return;
    const visible = new Set(visibleHitGroups(entry)
      .slice(0,state.screeningExpandAll ? undefined : 5).map(g=>g.accession));
    const batch = wanted.slice(index,index+5).filter(a => visible.has(a)
      && !entry.search.annotations[a] && !SCREENING_ANNOTATIONS.has(key(a)));
    if (!batch.length) continue;
    batch.forEach(a=>SCREENING_ANNOTATIONS.set(key(a), true));
    try {
      const annotations = await api(`/api/screenings/${entry.search.id}/annotations`, {
        method:"POST", headers:{"Content-Type":"application/json"},
        body:JSON.stringify({accessions:batch}),
      });
      if (state.result?.app.run_id !== runId) return;
      Object.assign(entry.search.annotations, annotations);
    } catch (e) {
      if (state.result?.app.run_id !== runId) return;
      for (const accession of batch) entry.search.annotations[accession] = {
        status:"unavailable", reason:"Chromosome annotation could not be loaded.",
      };
    } finally {
      batch.forEach(a=>SCREENING_ANNOTATIONS.delete(key(a)));
    }
    if (state.screeningOpen === entry.row.candidate && state.screeningScope === "transcript") renderOffTarget();
  }
}

function chromosomeLabel(group, entry, genomic) {
  if (genomic) {
    const p = group.placement;
    const accession = group.accession.split(".")[0];
    const number = /^NC_0000(\d{2})$/.exec(accession)?.[1];
    const chromosome = number && Number(number) <= 24
      ? (Number(number) === 23 ? "X" : Number(number) === 24 ? "Y" : String(Number(number))) : null;
    return chromosome ? `Chr ${chromosome}` : (p?.on_chromosome ? `Chr ${escapeHTML(p.contig)}` : entry.search.backend === "local" ? `Contig ${escapeHTML(p?.contig || group.accession)}` : "Unmapped contig");
  }
  const annotation = entry.search.annotations[group.accession];
  if (!annotation) return frozen ? "Not recorded" : "Loading…";
  return annotation.status === "available"
    ? `<span ${otTip(`${annotation.source}; gene context, not transcript-alignment coordinates. ${annotation.retrieved_at || ""}`)}>Chr ${escapeHTML(annotation.chromosome)}<small>${entry.search.backend === "local" ? "RefSeq transcript mapping" : "gene annotation"}</small></span>`
    : `<span class="muted" ${otTip(annotation.reason || "Annotation unavailable")}>Unavailable</span>`;
}

function visibleHitGroups(entry, core = state.result?.candidates?.find(c=>c.id === entry.row.candidate)?.target_site_dna?.length || 19) {
  const groups = groupHits(entry);
  if (state.screeningMatchFilter !== "threshold") return groups;
  return groups.filter(g=>g.hits.some(h=>h.identities >= (entry.row.strong_threshold || Math.ceil(core*0.8))));
}

function alignmentDetail(entry, core, genomic) {
  const groups = visibleHitGroups(entry, core);
  const shown = !state.screeningExpandAll ? groups.slice(0, 5) : groups;
  const body = shown.map(group => {
    const strongest = group.hits.reduce((a, b) => a.identities >= b.identities ? a : b);
    const p = group.placement;
    const context = !genomic ? "" : !p ? "TSS context not annotated"
      : p.annotation_status === "unavailable" ? "TSS annotation unavailable"
      : p.on_chromosome === false && entry.search.backend !== "local" ? "Contig not mapped to an assembled chromosome"
      : p.distance_to_tss === null ? (entry.search.backend === "local" ? "No eligible annotated transcript start within 10 kb" : "No annotated TSS within 10 kb")
      : `${signed(p.distance_to_tss)} bp from ${p.nearest_gene || p.nearest_transcript || "nearest TSS"}${entry.search.backend === "local" ? " · annotated TSS proxy" : ""}`
        + (p.inside_run_window ? " · inside this run’s window" : "");
    return `<article class="ot-hit"><div class="ot-hit-info"><div class="ot-hit-heading">`
      + `<b>${escapeHTML(hitName(strongest,p))}</b><span>${chromosomeLabel(group,entry,genomic)}</span></div>`
      + `<p class="candidate-id">${escapeHTML(group.accession)}`
      + `${group.hits.length > 1 ? ` · ${group.hits.length} alignments` : ""}</p>`
      + `<p class="ot-description">${escapeHTML(group.title || "Description unavailable")}</p>`
      + `<p class="ot-hit-metrics"><b>${strongest.identities}/${core}</b> bases matched · `
      + `${escapeHTML(strongest.subject_strand)} subject strand · E-value ${strongest.evalue == null ? "not applicable (local alignment)" : Number(strongest.evalue).toExponential(1)}</p>`
      + (genomic ? `<p class="hint">${p?.on_chromosome ? `Chr ${escapeHTML(p.contig)}:${fmt(p.position)} · ` : ""}${escapeHTML(context)}</p>` : "")
      + `</div><div class="ot-hit-alignment">${footprint(strongest,core)}</div></article>`;
  }).join("");
  const missingAnnotations = !genomic && Object.values(entry.search.annotations).some(a=>a.status === "unavailable");
  const threshold = entry.row.strong_threshold || Math.ceil(core*0.8);
  const matching = entry.row.wang_homology_hits ?? 0;
  const weak = Math.max(0,entry.row.hit_count-matching);
  return `<tr class="ot-detail"><td colspan="${genomic ? 10 : 9}">`
    + `<p class="ot-candidate-explanation"><b>1 saRNA candidate searched → ${entry.row.hit_count} alignments returned → `
    + `${matching} with ≥${threshold}/${core} matching bases.</b> `
    + `${weak ? `${weak} weaker alignment${weak === 1 ? " is" : "s are"} also included in the returned list.`
      : "All returned alignments meet this reporting threshold."}</p>`
    + `<p class="hint">These are sequence findings, not confirmed biological targets.</p>`
    + (entry.row.details_retained < entry.row.hit_count ? `<p class="hint">All ${entry.row.hit_count} alignments were counted. The strongest ${entry.row.details_retained} alignment details are retained for this candidate; the detail list is not the full search output.</p>` : "")
    + (entry.row.alignments_capped ? `<p class="ot-limit">${escapeHTML(screeningLimit(entry))}</p>` : "")
    + `<div class="ot-match-filter" role="group" aria-label="Filter displayed alignments">`
    + `<button type="button" data-match-filter="all" aria-pressed="${state.screeningMatchFilter !== "threshold"}">All returned (${entry.row.hit_count})</button>`
    + `<button type="button" data-match-filter="threshold" aria-pressed="${state.screeningMatchFilter === "threshold"}">≥${threshold}/${core} matches (${matching})</button>`
    + `<span>Display only; counts and candidate ranks stay unchanged.</span></div>`
    + (!groups.length ? `<p class="hint">No returned alignments meet this display filter. Select “All returned” to inspect weaker matches; a negative result does not establish specificity.</p>` : "")
    + `<p class="fp-legend"><span class="match">Match</span><span class="mismatch">Mismatch</span>`
    + `<span class="gap">Gap</span><span class="outside">Outside alignment</span></p>`
    + `<div class="ot-hit-list">${body}</div>`
    + (shown.length < groups.length
      ? `<p class="hint">Showing the 5 strongest of ${groups.length} reported ${genomic ? "loci" : "sequences"}. `
        + `<button type="button" class="link" data-expand-all="1">Show all</button></p>` : "")
    + (missingAnnotations && !frozen ? `<p class="hint"><button type="button" class="link" data-retry-annotations="1">Retry unavailable chromosome annotations</button></p>` : "")
    + `</td></tr>`;
}

/* ---------------------------------------------------------------- the genomic picture */

/** Every genomic hit against its signed distance to the nearest annotated start site.
 *
 * No threshold line is drawn and nothing is classified as "a promoter hit". The only band is
 * the window this run designed against, which is a statement about the run (ADR 0019).
 */
function tssPlot(entries) {
  const subset = entries.some(e => e.row.details_retained !== undefined && e.row.details_retained < e.row.hit_count);
  const placed = entries.flatMap(e => e.placements.filter(p => p.distance_to_tss !== null)
    .map(p => ({...p, candidate: e.row.candidate})));
  const total = entries.reduce((n,e)=>n + e.row.hit_count, 0);
  const retained = entries.reduce((n,e)=>n + (e.row.details_retained ?? e.hits?.length ?? e.placements.length), 0);
  const percent = total ? (100 * placed.length / total).toLocaleString(undefined, {maximumFractionDigits:1}) : "0";
  const scopeNote = `<p class="hint"><b>${fmt(placed.length)} of ${fmt(total)} reported genomic alignments (${percent}%) are shown</b>: displayed alignments within 10 kb of the nearest eligible annotated start site anywhere in the reference, which may belong to another gene.</p>`
    + `<p class="hint">${fmt(retained)} alignment details retained${subset ? " (up to 200 per candidate, prioritising the closest sequence matches)" : ""}. ${subset ? "Other alignments may also lie near a start site; their locations are not shown. This percentage describes plot coverage, not the fraction of all alignments near a TSS." : "This is annotation proximity, not evidence of binding or gene regulation."}</p>`;
  if (!placed.length) {
    return scopeNote + `<p class="hint">No displayed genomic alignment has an eligible annotated start nearby, so `
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
  const window = entries.map(e=>e.search?.window).find(Boolean) || null;
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
  const full19 = placed.every(p=>p.hit.alignment_length === 19);
  placed.forEach((p, i) => {
    // Vertical offset is spacing so overlapping marks stay countable. It means nothing, and
    // the caption says so rather than leaving a reader to find a pattern in it.
    const y = 34 + (i % 5) * 12;
    const colour = !full19 ? (p.hit.identities >= Math.ceil(p.hit.alignment_length * 0.8) ? "#1d5c52" : "#9cc6ba") : p.hit.identities >= 19 ? "#81572b" : p.hit.identities === 18 ? "#327968" : p.hit.identities === 17 ? "#699d91" : "#b4d0c8";
    parts.push(`<circle cx="${scale(p.distance_to_tss).toFixed(1)}" cy="${y}" `
      + `r="${p.inside_run_window ? 5 : 4}" fill="${colour}" `
      + `${p.inside_run_window ? 'stroke="#c08a2e" stroke-width="1.5"' : ""} `
      + otTip(`${p.candidate} aligns ${p.hit.identities} of ${p.hit.alignment_length} bases on `
        + `${p.contig}:${fmt(p.position)}, ${Math.abs(p.distance_to_tss).toLocaleString()} bp `
        + `${p.distance_to_tss < 0 ? "upstream of" : "downstream of"} the start site of `
        + `${p.nearest_gene || p.nearest_transcript || "an unnamed transcript"}`
        + `${p.inside_run_window ? ", inside the window this run designed against" : ""}. `
        + `Proximity to an annotated start site does not establish that the promoter is active `
        + `or that the duplex binds it.`) + `/>`);
  });
  return `<details class="ot-plot"><summary>Genomic alignment locations · ${fmt(placed.length)} / ${fmt(total)} shown (${percent}%)</summary><h4>Where the genomic hits sit, relative to the nearest `
    + `annotated start site</h4>`
    + scopeNote
    + `<svg viewBox="0 0 ${width} ${height}" role="img" `
    + `aria-label="Genomic off-target hits by signed distance to the nearest annotated `
    + `transcription start site">${parts.join("")}</svg>`
    + `<p class="hint">Distance is on a log scale either side of the start site, so near hits `
    + `stay separable from far ones. Vertical position is spacing only and carries no meaning. `
    + `${full19 ? "Colour indicates matching bases: 19 (brown), 18 (dark green), 17 (green), 16 (pale green)" : "Darker marks reach the reporting threshold"}${window ? ", and the shaded band is the window this run "
        + "designed against" : ""}. Annotation proximity does not establish promoter activity, binding or an off-target effect.</p></details>`;
}

/* ---------------------------------------------------------------- events */

document.addEventListener("click", (event) => {
  const batch = event.target.closest("[data-batch]");
  if (batch) {
    state.screenBatch = state.screenBatch || {};
    state.screenBatch[state.screeningScope || "transcript"] = Number(batch.dataset.batch);
    renderScreeningStatus();
    return;
  }
  const mode = event.target.closest("[data-batch-mode]");
  if (mode) { startScreening(state.screeningScope || "transcript", mode.dataset.batchMode); return; }
  const tab = event.target.closest("[data-scope-tab]");
  if (tab) { renderOffTarget(tab.dataset.scopeTab); return; }
  const matchFilter = event.target.closest("[data-match-filter]");
  if (matchFilter) {
    state.screeningMatchFilter = matchFilter.dataset.matchFilter;
    state.screeningExpandAll = false;
    renderOffTarget();
    const entry = comparisonEntry(state.screeningScope,state.screeningOpen);
    if (entry && state.screeningScope === "transcript") {
      loadHitAnnotations(entry, visibleHitGroups(entry).slice(0,5).map(g=>g.accession));
    }
    return;
  }
  const sort = event.target.closest("[data-screensort]");
  if (sort) {
    const previous = state.screeningSort || "design_rank";
    state.screeningSortDescending = previous === sort.dataset.screensort ? !state.screeningSortDescending : false;
    state.screeningSort = sort.dataset.screensort;
    renderOffTarget(); return;
  }
  if (event.target.closest("#download-screening-svg")) { downloadScreeningSvg(); return; }
  if (event.target.closest("#screening-order-direction")) {
    state.screeningSortDescending = !state.screeningSortDescending;
    renderOffTarget(); return;
  }
  if (event.target.closest("[data-expand-all]")) {
    state.screeningExpandAll = true; renderOffTarget();
    const entry = comparisonEntry(state.screeningScope,state.screeningOpen);
    if (entry && state.screeningScope === "transcript") loadHitAnnotations(entry, visibleHitGroups(entry).map(g => g.accession));
    return;
  }
  if (event.target.closest("[data-retry-annotations]")) {
    const entry = comparisonEntry(state.screeningScope,state.screeningOpen);
    if (entry) {
      for (const [key,value] of Object.entries(entry.search.annotations)) {
        if (value.status === "unavailable") delete entry.search.annotations[key];
      }
      renderOffTarget();
      loadHitAnnotations(entry, visibleHitGroups(entry).slice(0,state.screeningExpandAll ? undefined : 5).map(g=>g.accession));
    }
    return;
  }
  const resume = event.target.closest("[data-resume]");
  if (resume) { resumeScreening(resume.dataset.resume); return; }
  const row = event.target.closest("[data-screenrow]");
  if (row) {
    const id = row.dataset.screenrow;
    state.screeningOpen = state.screeningOpen === id ? null : id;
    state.screeningExpandAll = false;
    renderOffTarget();
    if (state.screeningOpen && state.screeningScope === "transcript") {
      const entry = comparisonEntry("transcript",id);
      if (entry) {
        const groups = visibleHitGroups(entry);
        loadHitAnnotations(entry, groups.slice(0,5).map(g => g.accession));
      }
    }
  }
});

document.addEventListener("keydown", event => {
  if ((event.key === "Enter" || event.key === " ") && event.target.matches("[data-screenrow]")) {
    event.preventDefault(); event.target.click();
  }
});

let REFERENCE_STATUS = null;
let REFERENCE_POLL = null;
function screeningBackend() { return $("screening-backend")?.value || "ncbi"; }
function renderReferenceStatus() {
  const info = REFERENCE_STATUS;
  if (!info || !$("reference-status")) return;
  const busy = info.status === "running";
  const plan = $("reference-download-plan");
  if (plan) {
    const names = {genomic:"Human genome", transcript:"RefSeq transcripts", annotation:"Gene annotations", assembly:"Chromosome map", tools:"Search tools"};
    const size = bytes => bytes >= 1e9 ? `${(bytes/1e9).toFixed(2)} GB` : bytes >= 1e6 ? `${(bytes/1e6).toFixed(1)} MB` : `${(bytes/1e3).toFixed(0)} kB`;
    plan.innerHTML = `<div class="reference-plan-heading"><strong>${info.ready ? "Installed · no download needed" : "One-time reference download"}</strong><span>${size(info.download_bytes || 0)} total</span></div>`
      + `<div class="reference-assets">${(info.download_plan || []).map(asset=>`<div><span>${escapeHTML(names[asset.key] || asset.key)}</span><strong>${size(asset.bytes)}</strong></div>`).join("")}</div>`
      + `<p class="hint">Allow ${(info.required_free_bytes/1024**3).toFixed(0)} GiB free for preparation. Download time depends on your connection; building indexes happens afterwards and can take an hour or more. Keep the app running. This is a one-time setup per reference release.</p>`;
  }
  if (!info.ready && !state.referenceSetupShown) {
    $("local-reference-panel").open = true;
    state.referenceSetupShown = true;
  }
  $("reference-summary").textContent = "Local off-target screening · " + (busy ? "preparing references…" : info.ready ? "references ready" : info.status === "error" ? "setup needs retry" : "reference setup");
  const label = info.installed ? `${info.installed.release} · Bowtie ${info.installed.bowtie_version}` : "No local reference installed";
  $("reference-status").textContent = `${label}. ${info.message || ""}` + (info.file ? ` · ${info.file}` : "")
    + (info.total ? ` · ${(info.downloaded / 1e6).toFixed(0)} / ${(info.total / 1e6).toFixed(0)} MB` : "")
    + (busy && info.elapsed_seconds >= 60 ? ` · ${Math.floor(info.elapsed_seconds/60)} min in this step` : "");
  const progress = $("reference-progress");
  progress.hidden = !busy;
  if (info.total) { progress.max = info.total; progress.value = info.downloaded || 0; }
  else progress.removeAttribute("value");
  $("reference-install").disabled = busy || !info.supported || (info.ready && !info.update_available);
  $("reference-install").textContent = info.ready ? (info.update_available ? "Download reference update" : "References ready · reused automatically")
    : ["paused", "error"].includes(info.status) ? "Resume / retry reference setup" : "Download & prepare local references";
  $("reference-pause").hidden = !busy;
  $("open-reference-setup").hidden = info.ready || screeningBackend() !== "local";
  $("local-screening-note").textContent = screeningBackend() === "ncbi"
    ? "Remote cross-check: candidate sequences are sent to NCBI. Queue time and reporting limits apply."
    : info.ready ? `Ready: ${info.installed.release}. Full-length 19-base matches with up to 3 mismatches; completed searches are cached.`
    : info.supported ? "Set up local references above before searching. Downloads happen only when you click the setup button."
    : "Automatic local setup currently supports Windows x64 and Linux x64. You may explicitly choose the remote NCBI cross-check.";
  $("screening-privacy").textContent = screeningBackend() === "ncbi" ? "Submitting a remote search sends candidate sequences to NCBI." : "Local searches keep candidate sequences on this computer.";
  if (state.result) renderScreeningStatus();
}
async function loadReferenceStatus() {
  if (frozen) return;
  clearTimeout(REFERENCE_POLL);
  try { REFERENCE_STATUS = await api("/api/references"); renderReferenceStatus(); }
  catch (e) { $("reference-status").textContent = e.message; }
  if (REFERENCE_STATUS?.status === "running") REFERENCE_POLL = setTimeout(loadReferenceStatus, 2000);
}
$("reference-install")?.addEventListener("click", async () => {
  $("reference-install").disabled = true;
  try { await api("/api/references/setup", {method:"POST", headers:{"Content-Type":"application/json"}, body:"{}"}); await loadReferenceStatus(); }
  catch (e) { error(e.message); $("reference-install").disabled = false; }
});
$("reference-pause")?.addEventListener("click", async () => {
  try { await api("/api/references/pause", {method:"POST", headers:{"Content-Type":"application/json"}, body:"{}"}); await loadReferenceStatus(); }
  catch (e) { error(e.message); }
});
$("reference-check")?.addEventListener("click", loadReferenceStatus);
$("open-reference-setup")?.addEventListener("click", () => {
  if (typeof setWorkspaceTab === "function") setWorkspaceTab("screening");
  $("local-reference-panel").open = true;
  $("local-reference-panel").scrollIntoView({behavior:"smooth", block:"center"});
});
$("screening-backend")?.addEventListener("change", () => {
  state.screeningEvidence=null;
  for (const scope of Object.keys(SCOPE_BUTTON)) cancelScreeningCheck(scope);
  state.screeningJobs = {};
  renderReferenceStatus();
  if (state.result) {
    renderScreeningStatus(); renderOffTarget(state.screeningScope || "transcript");
    loadSavedScreenings(state.result.app.run_id);
  }
});
if (frozen) {
  $("local-reference-panel").hidden = true;
  $("screening-backend").closest("label").hidden = true;
  $("open-reference-setup").hidden = true;
  $("local-screening-note").textContent = "Saved observations. Each candidate names its search method; reference versions and dates remain recorded.";
  $("screening-privacy").textContent = "This archive displays saved findings and does not submit screening requests.";
}
else loadReferenceStatus();

$("screen-transcript").addEventListener("click", () => startScreening("transcript"));
$("screen-genomic").addEventListener("click", () => startScreening("genomic"));

// This file is a second script, so a frozen archive has already run app.js and shown its run
// by the time this executes. Catch up rather than leave the panel blank: an archive of a
// screened run that showed no screening would be a record quietly missing its own evidence.
if (state.result) {
  for (const record of frozenScreeningRecords().filter(r=>r.status === "complete")) mergeScreening(record);
  renderOffTarget("genomic");
  renderScreeningStatus();
}

/* Off-target findings accompany the original design rank. Historical project-rule
 * assessments remain review notes; they never remove a design or block export.
 */
function screenedDecisions() {
  const candidates = (state.result?.candidates || []).filter(c=>c.rank !== null);
  let key = `local:${REFERENCE_STATUS?.installed?.release || REFERENCE_STATUS?.catalog_release || ""}`;
  if (frozen) {
    const shared = Object.keys(screenings().genomic?.byEvidence || {}).filter(k=>k.startsWith("local:") && screenings().transcript?.byEvidence?.[k]);
    key = shared.at(-1) || "local:";
  }
  return candidates.map(candidate => {
    const evidence = Object.keys(SCOPE_BUTTON).map(scope=>({scope, entry:screenings()[scope]?.byEvidence?.[key]?.[candidate.id]}));
    const reasons = [], pending = [];
    for (const {scope,entry} of evidence) {
      const assessment = entry?.row.selection_assessment;
      if (!entry) pending.push(`${scope}: not screened with this local reference`);
      else if (!assessment || assessment.policy?.id !== "conservative-full-length-v1" || !entry.row.enumeration_complete || entry.row.alignments_capped) pending.push(`${scope}: refresh screening to obtain complete selection evidence; local cached searches are reused`);
      else if (assessment.status === "excluded") reasons.push(...assessment.reasons.map(reason=>`${scope}: ${reason}`));
      else if (assessment.status !== "eligible") pending.push(...assessment.reasons.map(reason=>`${scope}: ${reason}`));
    }
    return {candidate, evidence, reasons:[...reasons,...pending], status:pending.length ? "pending" : reasons.length ? "review" : "complete"};
  }).sort((a,b)=>a.candidate.rank-b.candidate.rank);
}
function offTargetCounts(entry) {
  const c = entry?.row.identity_counts;
  if (!c || ["19","18","17"].some(k=>!Number.isInteger(c[k]) || c[k] < 0)) return '<span class="hint">Not available</span>';
  const complete = entry.row.enumeration_complete && !entry.row.alignments_capped;
  return `<span class="ot-counts">${fmt(c["19"])} / ${fmt(c["18"])} / ${fmt(c["17"])}</span>`
    + (complete ? "" : '<small>Incomplete evidence</small>');
}
function candidateScreeningSummary(decision) {
  if (!decision) return "";
  return decision.evidence.map(e=>`<span>${e.scope === "genomic" ? "Genome" : "RNA"} 19/18/17: ${offTargetCounts(e.entry)}</span>`).join("")
    + (decision.status === "pending" ? '<span>Screening evidence incomplete</span>' : decision.status === "review" ? '<span class="screening-evidence review">Review sequence matches</span>' : "")
    + `<button type="button" class="text-button" data-screening-candidate="${escapeHTML(decision.candidate.id)}">View screening details</button>`;
}
function topOffTargetSummary() {
  const best = screenedDecisions().find(d=>d.candidate.rank === 1);
  if (!best) return "";
  return `<details class="top-screening-note"><summary>Off-target sequence matches · ${best.status === "pending" ? "evidence incomplete" : best.status === "review" ? "review findings" : "both scopes recorded"}</summary>`
    + `<p>Full-length reference alignments at 19 / 18 / 17 matching bases: `
    + best.evidence.map(e=>`${e.scope === "genomic" ? "Genome" : "Transcripts"}: ${offTargetCounts(e.entry)}`).join(" · ")
    + `.</p><p class="hint">Sequence similarity is a review signal. It does not establish binding or a biological effect. Details are in Off-target screening.</p><button type="button" class="text-button" data-screening-candidate="${escapeHTML(best.candidate.id)}">View screening details →</button></details>`;
}
/* Comparisons use saved measurements only; no chart changes rank or selection. */
function screeningEvidenceKey() {
  const keys = [...new Set(Object.values(screenings()).flatMap(h=>Object.keys(h.byEvidence || {})))];
  if (state.screeningEvidence && keys.includes(state.screeningEvidence)) return state.screeningEvidence;
  if (!frozen) return selectedEvidenceKey();
  const shared = keys.filter(k=>k.startsWith("local:") && screenings().genomic?.byEvidence?.[k] && screenings().transcript?.byEvidence?.[k]);
  return shared.at(-1) || keys.filter(k=>k.startsWith("local:")).at(-1) || keys.at(-1) || "local:";
}
function comparisonEntry(scope,id) {
  const held = screenings()[scope];
  return held?.byEvidence ? held.byEvidence[screeningEvidenceKey()]?.[id] : held?.byCandidate?.[id];
}
function screeningComparisonDecisions() {
  let decisions = (state.result?.candidates || []).filter(c=>c.rank !== null)
    .sort((a,b)=>a.rank-b.rank).map(candidate=>({candidate,
      evidence:["transcript","genomic"].map(scope=>({scope,entry:comparisonEntry(scope,candidate.id)}))}));
  if (state.screeningPlotLimit === 10 || state.screeningPlotLimit === 15) decisions = decisions.slice(0,state.screeningPlotLimit);
  else decisions = decisions.filter(d=>d.evidence.some(e=>e.entry));
  const entries = decisions.map(d=>comparisonEntry(state.screeningScope || "genomic",d.candidate.id)
    || {row:{candidate:d.candidate.id,enumeration_complete:false},hits:[],placements:[],missing:true});
  const order = screeningTableOrder(entries,state.screeningSort || "design_rank",!!state.screeningSortDescending);
  const byId = new Map(decisions.map(d=>[d.candidate.id,d]));
  return order.map(e=>byId.get(e.row.candidate));
}
function screeningSortChoices(core=19,genomic=true) {
  return [["design_rank","Design rank · same as Candidates"],
    ["hit_count","Fewest reported sequence matches first"],
    ["perfect_hits","Fewest exact alignments first"],
    ["wang_homology_hits",`Fewest alignments ≥${Math.ceil(core*.8)}/${core} first`],
    [genomic ? "distinct_loci" : "distinct_sequences",genomic ? "Fewest distinct reported loci first" : "Fewest transcript records first"],
    ...(genomic ? [["in_run_window_hits","Fewest matches in the design window first"]] : [])];
}
function screeningComparisonControls() {
  const genomic = (state.screeningScope || "genomic") === "genomic", choices = screeningSortChoices(19,genomic);
  const sort = state.screeningSort || "design_rank", descending = !!state.screeningSortDescending;
  const evidenceKeys=[...new Set([screeningEvidenceKey(),...Object.values(screenings()).flatMap(h=>Object.keys(h.byEvidence || {}))])];
  const evidenceSelect=`<label>Evidence set <select id="screening-evidence-set">${evidenceKeys.map(key=>`<option value="${escapeHTML(key)}" ${key===screeningEvidenceKey() ? "selected" : ""}>${escapeHTML(key.startsWith("local:") ? `Local · ${key.slice(6) || "current reference"}` : "NCBI BLAST · reported results")}</option>`).join("")}</select></label>`;
  return `<div class="screening-order">${evidenceSelect}<label>Order using <select id="screening-sort-scope"><option value="genomic" ${genomic ? "selected" : ""}>Genome counts</option><option value="transcript" ${!genomic ? "selected" : ""}>Transcript counts</option></select></label><label>Order by <select id="screening-order">${choices.map(([key,label])=>`<option value="${key}" ${sort===key ? "selected" : ""}>${escapeHTML(descending ? label.replace(/^Fewest /,"Most ") : label)}</option>`).join("")}</select></label><button type="button" class="secondary" id="screening-order-direction">${sort === "design_rank" ? descending ? "Reverse rank order ↓" : "Rank 1 first ↑" : descending ? "Most matches first ↓" : "Fewest matches first ↑"}</button><label>Candidate range <select id="screening-plot-limit"><option value="all" ${!state.screeningPlotLimit ? "selected" : ""}>All screened designs</option><option value="10" ${state.screeningPlotLimit===10 ? "selected" : ""}>Design ranks 1–10</option><option value="15" ${state.screeningPlotLimit===15 ? "selected" : ""}>Design ranks 1–15</option></select></label></div><p class="hint">The plot and detailed table share this order and range. ${sort === "design_rank" ? descending ? "Reverse design-rank order." : "Same default order as Candidates." : "Count ordering uses the chosen scope; incomplete or missing evidence is placed last."} Both bars stay together. Fewer sequence matches do not establish biological specificity.</p>`;
}
function matchCountData(entry) {
  if (!entry || !Number.isInteger(entry.row.hit_count) || entry.row.hit_count<0) return null;
  const total=entry.row.hit_count, c=entry.row.identity_counts;
  const valid=c && ["19","18","17","16"].every(k=>Number.isInteger(c[k]) && c[k]>=0)
    && ["19","18","17","16"].reduce((n,k)=>n+c[k],0)<=total;
  const counts=Object.fromEntries(["19","18","17","16"].map(k=>[k,valid ? c[k] : 0]));
  counts.other=total-Object.values(counts).reduce((a,b)=>a+b,0);
  return {counts,total,complete:!!entry.row.enumeration_complete && !entry.row.alignments_capped};
}
function matchCountPlot(decisions) {
  if (!decisions.length) return `<section class="paired-match-plot" data-comparison-plot="counts"><h4>Sequence matches per saRNA</h4><p class="hint">Screen candidates to compare their transcript and genome matches. Missing evidence is not a zero result.</p></section>`;
  const max=Math.max(1,...decisions.flatMap(d=>d.evidence.map(e=>matchCountData(e.entry)?.total || 0)));
  const width=Math.max(740,70+decisions.length*98), height=394, left=58, base=282, plotHeight=222;
  const step=(width-left-20)/decisions.length, barWidth=Math.min(28,step*.28);
  const colours={...MATCH_COLOURS,"16":"#d6e5df",other:"#9aa6ad"};
  const legend=["19","18","17","16","other"].map(key=>[key,colours[key]]).map(([key,colour])=>`<span><i style="background:${colour}"></i>${key==="other" ? "Other / unclassified" : `${key}/19 bases`}</span>`).join("");
  const parts=[`<title>Reported sequence matches by candidate, transcripts and genome</title><desc>Two stacked bars per candidate on one linear scale. Colours represent 19, 18, 17 and 16 matching bases across a 19-base core. Grey represents other or unclassified reported alignments. Missing evidence is marked with a question mark, not zero.</desc><defs><pattern id="missing-screening-pattern" width="6" height="6" patternUnits="userSpaceOnUse"><path d="M0 6L6 0" stroke="#9aa6ad" stroke-width="1"/></pattern></defs>`];
  for(let tick=0;tick<=4;tick++) {
    const y=base-plotHeight*tick/4, count=Math.round(max*tick/4);
    parts.push(`<line x1="${left}" x2="${width-12}" y1="${y}" y2="${y}" stroke="#e1e9e5"/><text x="${left-8}" y="${y+4}" text-anchor="end" font-size="11" fill="#526b75">${fmt(count)}</text>`);
  }
  parts.push(`<text x="${left}" y="25" font-size="12" fill="#254c43">Reported reference alignments · linear scale 0–${fmt(max)}</text>`);
  decisions.forEach((d,i)=>{
    const centre=left+step*(i+.5), score=d.candidate.rating;
    const detail=d.evidence.map(e=>{
      const data=matchCountData(e.entry);
      return `${e.scope === "genomic" ? "Genome" : "Transcripts"}: ${data ? `${data.total} reported alignments; ${["19","18","17","16"].map(k=>`${data.counts[k]} at ${k}/19`).join(", ")}; ${data.counts.other} other/unclassified; ${data.complete ? "complete enumeration" : "reported counts only / incomplete"}; ${e.entry.search?.reference?.release || e.entry.search?.parameters?.database || "reference not recorded"}; ${e.entry.search?.searched_at || "date not recorded"}` : "not screened / counts unavailable"}`;
    }).join(". ");
    const tip=`${d.candidate.id}, design rank ${d.candidate.rank}, score ${score?.met ?? "?"}/${score?.evaluated ?? "?"} preferences. Target DNA 5′–3′: ${d.candidate.target_site_dna || "not available"}. ${detail}. Click to open this candidate’s screening details.`;
    parts.push(`<g class="paired-candidate" data-plot-candidate="${escapeHTML(d.candidate.id)}" data-screening-candidate="${escapeHTML(d.candidate.id)}" role="button" ${otTip(tip)}><title>${escapeHTML(tip)}</title>`);
    for(const [index,scope] of ["transcript","genomic"].entries()) {
      const e=d.evidence.find(e=>e.scope===scope), data=matchCountData(e?.entry);
      const x=centre+(index===0 ? -barWidth-3 : 3);
      parts.push(`<g data-bar-scope="${scope}" data-screening-candidate="${escapeHTML(d.candidate.id)}" data-screening-scope="${scope}" data-total="${data?.total ?? "unavailable"}">`);
      if (!data) parts.push(`<rect x="${x}" y="${base-28}" width="${barWidth}" height="28" fill="url(#missing-screening-pattern)" stroke="#9aa6ad"/><text x="${x+barWidth/2}" y="${base-34}" text-anchor="middle" font-size="12" fill="#526b75">?</text>`);
      else {
        let y=base;
        for(const key of ["19","18","17","16","other"]) {
          const h=plotHeight*data.counts[key]/max;
          if(h>0) { y-=h;parts.push(`<rect data-match-bases="${key}" data-count="${data.counts[key]}" x="${x}" y="${y}" width="${barWidth}" height="${h}" fill="${colours[key]}"><title>${escapeHTML(`${scope}: ${data.counts[key]} ${key==="other" ? "other/unclassified alignments" : `alignments at ${key}/19`}`)}</title></rect>`); }
        }
        if(!data.total) parts.push(`<line x1="${x}" x2="${x+barWidth}" y1="${base}" y2="${base}" stroke="#327968" stroke-width="2"/>`);
        if(!data.complete) parts.push(`<rect x="${x}" y="${base-Math.max(3,plotHeight*data.total/max)}" width="${barWidth}" height="${Math.max(3,plotHeight*data.total/max)}" fill="none" stroke="#81572b" stroke-dasharray="3 2"/>`);
        parts.push(`<text x="${x+barWidth/2}" y="${Math.max(42,y-7)}" text-anchor="middle" font-size="10" fill="#254c43">${fmt(data.total)}${data.complete ? "" : "*"}</text>`);
      }
      parts.push(`<text x="${x+barWidth/2}" y="${base+17}" text-anchor="middle" font-size="10" fill="#526b75">${index===0 ? "RNA" : "DNA"}</text></g>`);
    }
    parts.push(`<text x="${centre}" y="${base+39}" text-anchor="middle" font-size="11" font-weight="600" fill="#087e80">${escapeHTML(d.candidate.id)}</text><text x="${centre}" y="${base+56}" text-anchor="middle" font-size="10" fill="#526b75">Rank ${d.candidate.rank} · ${score?.met ?? "?"}/${score?.evaluated ?? "?"}</text></g>`);
  });
  parts.push(`<text x="${left}" y="${height-15}" font-size="11" fill="#526b75">RNA = transcript reference · DNA = genome reference · ? unavailable · * counts not established complete</text>`);
  return `<section class="paired-match-plot" data-comparison-plot="counts"><div class="plot-heading"><div><h4>Sequence matches per saRNA</h4><p class="hint">Transcripts and genome side by side. Hover or focus a candidate for counts, rank, score and search details; click to inspect.</p></div><button type="button" class="secondary" id="download-screening-svg">Download plot · SVG</button></div><div class="plot-legend">${legend}</div><div class="paired-plot-scroll"><svg id="screening-match-svg" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${width} ${height}" width="${width}" height="${height}" font-family="system-ui,sans-serif" role="group" aria-label="Paired stacked transcript and genome match counts">${parts.join("")}</svg></div><p class="hint">${decisions.length} candidates shown. Counts are reference alignments, including isoforms and alternate records; these are not independent genes or measured biological off-targets. Stacks include 16/19 matches. Grey covers reported alignments without compatible full-length identity counts. Missing evidence is not zero.</p></section>`;
}
function renderScreeningComparison(decisions=screeningComparisonDecisions()) {
  const host=$("screening-comparisons");
  if(host) screeningHTML(host,screeningComparisonControls()+matchCountPlot(decisions));
}
function downloadScreeningSvg() {
  const svg=$("screening-match-svg");
  if(!svg) return;
  const copy=svg.cloneNode(true);
  const metadata=document.createElementNS("http://www.w3.org/2000/svg","metadata");
  metadata.textContent=JSON.stringify({run:state.result?.app.run_id,reference:screeningEvidenceKey(),sort:state.screeningSort || "design_rank",sort_scope:state.screeningScope,descending:!!state.screeningSortDescending,range:state.screeningPlotLimit || "all screened",interpretation:"Reference alignments; not measured biological off-targets."});
  copy.prepend(metadata);
  const dimensions=copy.getAttribute("viewBox").split(/\s+/).map(Number), originalHeight=dimensions[3];
  const colours={...MATCH_COLOURS,"16":"#d6e5df",other:"#9aa6ad"};
  const legend=document.createElementNS("http://www.w3.org/2000/svg","g");
  ["19","18","17","16","other"].forEach((key,i)=>{
    const rect=document.createElementNS("http://www.w3.org/2000/svg","rect");
    Object.entries({x:58+i*126,y:originalHeight+8,width:12,height:12,fill:colours[key]}).forEach(([k,v])=>rect.setAttribute(k,v));
    const text=document.createElementNS("http://www.w3.org/2000/svg","text");
    Object.entries({x:76+i*126,y:originalHeight+18,"font-size":11,fill:"#254c43"}).forEach(([k,v])=>text.setAttribute(k,v));
    text.textContent=key === "other" ? "Other / unclassified" : `${key}/19 bases`;
    legend.append(rect,text);
  });
  const note=document.createElementNS("http://www.w3.org/2000/svg","text");
  Object.entries({x:58,y:originalHeight+43,"font-size":11,fill:"#526b75"}).forEach(([k,v])=>note.setAttribute(k,v));
  note.textContent="Reference alignments, not independent genes or measured biological effects. Retained details can be a subset.";
  legend.append(note);copy.append(legend);dimensions[3]+=60;
  copy.setAttribute("viewBox",dimensions.join(" "));copy.setAttribute("height",dimensions[3]);
  downloadBlob(new Blob([new XMLSerializer().serializeToString(copy)],{type:"image/svg+xml"}),`${state.result?.gene?.symbol || "promoter"}-screening.svg`);
}

function preferenceHeatmap(decisions) {
  const configured = state.result?.configuration?.rules || [];
  const ids = [...new Set([...configured.filter(r=>r.kind === "preference").map(r=>r.id),...decisions.flatMap(d=>(d.candidate.rules || d.candidate.preference_verdicts || []).filter(r=>r.kind === "preference").map(r=>r.id))])];
  const names = new Map(configured.map(r=>[r.id,r.name || r.id]));
  if (!ids.length) return `<details class="candidate-plot" data-comparison-plot="criteria"><summary>Design-criteria heatmap</summary><p class="hint">No preference criteria are available in this saved run.</p></details>`;
  const cells = decisions.map(d=>`<tr><th scope="row"><button class="text-button" type="button" data-screened-inspect="${escapeHTML(d.candidate.id)}">#${d.candidate.rank} ${escapeHTML(d.candidate.id)}</button>${typeof candidateControlBadges === "function" ? candidateControlBadges(d.candidate) : ""}</th>${ids.map(id=>{
    const r=(d.candidate.rules || d.candidate.preference_verdicts || []).find(r=>r.id===id && r.kind === "preference");
    const outcome = r?.outcome === "pass" ? "pass" : r?.outcome === "fail" ? "fail" : "unknown";
    const label = outcome === "pass" ? "Met" : outcome === "fail" ? "Missed" : "Not evaluated";
    const detail = `${d.candidate.id}, ${names.get(id)||id}: ${label}. ${r?.detail || "No saved verdict is available."}`;
    return `<td class="criterion-cell ${outcome}"><button type="button" class="criterion-mark" ${otTip(detail)}>${outcome === "pass" ? "✓" : outcome === "fail" ? "×" : "?"}</button></td>`;
  }).join("")}</tr>`).join("");
  return `<details class="candidate-plot" data-comparison-plot="criteria"><summary>Design-criteria heatmap · ${ids.length} preference criteria</summary><p class="hint">Saved preference verdicts for the displayed top designs: ✓ met, × missed, ? not evaluated or unavailable. Hover or focus a cell for its explanation. Mandatory design gates and tie-breakers are separate from these preferences; no verdict is inferred from the overall score.</p><div class="plot-scroll"><table class="criteria-heatmap"><caption>Design preferences by candidate</caption><thead><tr><th scope="col">Candidate</th>${ids.map(id=>`<th scope="col">${escapeHTML(names.get(id)||id)}</th>`).join("")}</tr></thead><tbody>${cells}</tbody></table></div></details>`;
}
function mismatchMapData(entry, core) {
  const retained = entry?.hits || [];
  const usable = retained.filter(h=>h.alignment_length===core && h.gaps===0
    && Math.min(h.query_start,h.query_end)===1 && Math.max(h.query_start,h.query_end)===core
    && h.has_alignment_detail && Array.isArray(h.mismatch_positions)
    && h.mismatch_positions.every(p=>Number.isInteger(p) && p>=1 && p<=core)
    && new Set(h.mismatch_positions).size === core-h.identities);
  return {retained:retained.length, usable:usable.sort((a,b)=>b.identities-a.identities || String(a.accession).localeCompare(String(b.accession)) || a.subject_start-b.subject_start)};
}
function mismatchPositionPlot(decisions) {
  const chosen = decisions.find(d=>d.candidate.id===state.mismatchPlotCandidate) || decisions[0];
  if (!chosen) return "";
  const scope = state.mismatchPlotScope === "transcript" ? "transcript" : "genomic";
  const entry = chosen.evidence.find(e=>e.scope===scope)?.entry;
  const core = chosen.candidate.target_site_dna?.length || 19;
  const data = mismatchMapData(entry,core), shown = data.usable.slice(0,40);
  const candidateSelect = `<label>Candidate <select id="mismatch-plot-candidate">${decisions.map(d=>`<option value="${escapeHTML(d.candidate.id)}" ${d===chosen ? "selected" : ""}>#${d.candidate.rank} ${escapeHTML(d.candidate.id)}</option>`).join("")}</select></label>`;
  const scopeSelect = `<label>Reference <select id="mismatch-plot-scope"><option value="genomic" ${scope==="genomic" ? "selected" : ""}>Genome</option><option value="transcript" ${scope==="transcript" ? "selected" : ""}>Transcripts</option></select></label>`;
  const total = entry?.row.hit_count;
  const note = !entry ? "No screening evidence is available for this candidate in this scope."
    : `${fmt(shown.length)} alignments shown out of ${fmt(data.retained)} retained details and ${fmt(total)} reported alignments. ${fmt(data.usable.length)} retained details have complete, ungapped, full-length mismatch positions.${entry.row.alignments_capped || !entry.row.enumeration_complete ? " Search evidence is incomplete." : ""} Only the strongest 40 usable retained alignments are drawn; this is not a complete map of all matched loci.`;
  const header = Array.from({length:core},(_,i)=>`<th scope="col">${i+1}<small>${escapeHTML(chosen.candidate.target_site_dna?.[i] || "?")}</small></th>`).join("");
  const body = shown.map(h=>{
    const mismatches = new Set(h.mismatch_positions);
    const locus = `${h.accession}:${Math.min(h.subject_start,h.subject_end)}–${Math.max(h.subject_start,h.subject_end)} (${h.subject_strand})`;
    return `<tr><th scope="row" title="${escapeHTML(locus)}">${escapeHTML(locus)}<small>${h.identities}/${core} bases</small></th>${Array.from({length:core},(_,i)=>`<td class="mismatch-map-cell ${mismatches.has(i+1) ? "mismatch" : "match"}" title="${escapeHTML(`${locus}; target position ${i+1}: ${mismatches.has(i+1) ? "mismatch" : "match"}`)}">${mismatches.has(i+1) ? "×" : "·"}</td>`).join("")}</tr>`;
  }).join("");
  return `<details class="candidate-plot" data-comparison-plot="mismatches"><summary>Mismatch-position map · inspect a candidate</summary><div class="plot-controls">${candidateSelect}${scopeSelect}</div><p class="hint">${escapeHTML(note)}</p><p class="hint">Columns follow the promoter target DNA, 5′ → 3′ (T corresponds to U in RNA), not the antisense guide orientation. × mismatch; · match. Gapped, partial or missing-position alignments are omitted. This does not estimate seed effects or biological specificity.</p>${shown.length ? `<div class="plot-scroll"><table class="mismatch-heatmap"><caption>Mismatch positions for ${escapeHTML(chosen.candidate.id)} · ${scope === "genomic" ? "genome" : "transcripts"}</caption><thead><tr><th scope="col">Reference alignment</th>${header}</tr></thead><tbody>${body}</tbody></table></div>` : `<p class="hint">No usable mismatch-position details to draw. Missing details do not mean zero matches.</p>`}</details>`;
}
function renderScreenedShortlist() {
  const host = $("screened-shortlist");
  if (!host || !state.result) return;
  const plotHosts = [host, $("design-comparison"), $("screening-comparisons"), $("offtarget-results")].filter(Boolean);
  const openPlots = plotHosts.flatMap(h=>[...(h.querySelectorAll?.("[data-comparison-plot][open]") || [])]).map(node=>node.dataset.comparisonPlot);
  const decisions = screenedDecisions(), limit = state.screenedShortlistLimit || 10;
  const chosen = decisions.slice(0,limit), pending = decisions.filter(d=>d.status === "pending");
  const byId = new Map(decisions.map(d=>[d.candidate.id,d]));
  for (const note of document.querySelectorAll?.("[data-candidate-screening]") || []) {
    screeningHTML(note, candidateScreeningSummary(byId.get(note.dataset.candidateScreening)));
  }
  const topNote = $("top-offtarget-summary");
  if (topNote) screeningHTML(topNote, topOffTargetSummary());
  screeningHTML(host, `<div class="screened-heading"><div><p class="screened-eyebrow">DESIGN RANK &amp; SCREENING</p><h3>Top saRNAs with off-target findings</h3><p class="hint">Ordered by original design rank. Off-target findings stay visible for review; they do not remove candidates.</p></div><label>Show top <select id="screened-shortlist-limit" aria-label="Top saRNA list size"><option value="10" ${limit===10?"selected":""}>10</option><option value="15" ${limit===15?"selected":""}>15</option></select></label></div>`
    + `<p class="screened-context hint">${chosen.length} shown · ${decisions.length - pending.length} with both scopes complete · ${pending.length} with incomplete evidence. Counts below are reference alignments at <b>19 / 18 / 17 matching bases</b>, including isoforms and alternate records; they are not counts of independent genes.</p>`
    + (chosen.length ? `<div class="screened-table-wrap"><table class="screened-table"><thead><tr><th>Design rank</th><th>Candidate</th><th>Design score</th><th>Transcripts<br><small>19 / 18 / 17 bases</small></th><th>Genome<br><small>19 / 18 / 17 bases</small></th><th>Screening evidence</th><th></th></tr></thead><tbody>${chosen.map(d=>`<tr><td><span class="screened-number">${d.candidate.rank}</span></td><td><strong>${escapeHTML(d.candidate.id)}</strong>${typeof candidateControlBadges === "function" ? candidateControlBadges(d.candidate) : ""}</td><td>${d.candidate.rating?.met ?? "n/a"} / ${d.candidate.rating?.evaluated ?? "n/a"}</td><td>${offTargetCounts(d.evidence.find(e=>e.scope==="transcript").entry)}</td><td>${offTargetCounts(d.evidence.find(e=>e.scope==="genomic").entry)}</td><td><span class="screening-evidence ${d.status === "review" ? "review" : ""}">${d.status === "pending" ? "Incomplete" : d.status === "review" ? "Review findings" : "Recorded"}</span>${d.reasons.length ? `<details><summary>Why?</summary><p class="hint">${d.reasons.map(escapeHTML).join(" · ")}</p></details>` : ""}</td><td><button class="secondary" type="button" data-screened-inspect="${escapeHTML(d.candidate.id)}">Inspect</button></td></tr>`).join("")}</tbody></table></div>` : `<div class="screened-empty"><p>No candidates passed the design gates. Review the design results or adjust your criteria.</p></div>`)
    + `<details class="screened-policy"><summary>How to interpret screening</summary><p>The former project rule flagged any 18/19 or 19/19 alignment, or five or more 17/19 alignments in either scope. These unvalidated cutoffs now provide review notes only. Missing or incomplete screening is labelled explicitly. Design rank and export selection are independent of these findings.</p><p>The intended genomic locus is removed where coordinates are known. Transcript matches, including those to the intended gene, remain visible because potential RNA silencing matters. Sequence counts do not measure binding, gene activation, silencing or safety. Seed effects, gaps, variants and cell-specific activity are not assessed.</p><p>Background: <a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC4379447/" target="_blank" rel="noopener">Wang et al., 2015</a>; <a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC4797303/" target="_blank" rel="noopener">Meng et al., 2016</a>. These papers do not validate the former numeric cutoffs.</p></details>`
    + `<div class="screened-actions">${pending.length && !frozen ? `<button id="screened-screen-all" class="secondary" type="button" ${!REFERENCE_STATUS?.ready || screeningBackend() !== "local" || Object.keys(SCOPE_BUTTON).some(scope=>["running","waiting"].includes(currentScreeningJob(scope)?.status)) ? "disabled" : ""}>Screen all locally · both scopes</button>` : ""}<button id="screened-use" class="secondary" type="button" ${chosen.length ? "" : "disabled"}>Select these ${chosen.length} for export</button><button id="screened-audit" class="secondary" type="button">Download screening record</button><span class="hint">Research selection; experimental validation remains necessary.</span></div>`);
  if ($("design-comparison")) screeningHTML($("design-comparison"), chosen.length ? `<div class="candidate-comparisons card"><h4>Compare design preferences · top ${chosen.length}</h4>${preferenceHeatmap(chosen)}</div>` : "");
  renderScreeningComparison();
  if (typeof workflowStatus === "function") workflowStatus();
  if (state.workspaceTab === "export" && typeof renderExportSelection === "function") renderExportSelection();
  for (const key of openPlots) {
    for (const target of plotHosts) {
      const plot = target.querySelector?.(`[data-comparison-plot="${key}"]`);
      if (plot) plot.open = true;
    }
  }
}
function handleComparisonChange(event) {
  if (event.target.id === "mismatch-plot-candidate" || event.target.id === "mismatch-plot-scope") {
    if (event.target.id === "mismatch-plot-candidate") state.mismatchPlotCandidate = event.target.value;
    else state.mismatchPlotScope = event.target.value;
    renderOffTarget();
    $(event.target.id)?.focus();
    return;
  }
  if (event.target.id !== "screened-shortlist-limit") return;
  state.screenedShortlistLimit = Number(event.target.value) === 15 ? 15 : 10;
  renderScreenedShortlist();
}
for (const id of ["screened-shortlist","design-comparison","screening-comparisons","offtarget-results"]) {
  $(id)?.addEventListener("change", handleComparisonChange);
}
function handleComparisonClick(event) {
  if (event.target.id === "screened-screen-all") {
    if (typeof setWorkspaceTab === "function") setWorkspaceTab("screening");
    if (screeningBackend() === "local") Promise.all([startScreening("transcript","all"),startScreening("genomic","all")]);
    return;
  }
  const inspect = event.target.closest("[data-screened-inspect]");
  if (inspect) { if (typeof setWorkspaceTab === "function") setWorkspaceTab("candidates"); selectCandidate(inspect.dataset.screenedInspect); $("candidate-layout").scrollIntoView({block:"start"}); }
  if (event.target.id === "screened-use") {
    state.checked = new Set(screenedDecisions().slice(0,state.screenedShortlistLimit || 10).map(d=>d.candidate.id));
    state.screenedManagedSelection = new Set(state.checked);
    state.screenedManagedRun = state.result.app.run_id;
    renderTable();
    event.target.textContent = `${state.checked.size} selected for export`;
  }
  if (event.target.id === "screened-audit") {
    const decisions = screenedDecisions().map(d=>({candidate:d.candidate.id, original_rank:d.candidate.rank, design_score:d.candidate.rating, status:d.status, reasons:d.reasons, evidence:d.evidence.map(e=>({scope:e.scope, search_id:e.entry?.search.id, searched_at:e.entry?.search.searched_at, reference:e.entry?.search.reference?.release, counts:e.entry?.row.identity_counts, assessment:e.entry?.row.selection_assessment}))}));
    downloadBlob(new Blob([JSON.stringify({run_id:state.result.app.run_id, generated_at:new Date().toISOString(), limit:state.screenedShortlistLimit || 10, order:"Original design rank; off-target findings are advisory", decisions},null,2)],{type:"application/json"}),"screened-selection.json");
  }
}
for (const id of ["screened-shortlist","design-comparison","screening-comparisons","offtarget-results"]) {
  $(id)?.addEventListener("click", handleComparisonClick);
}

function validateScreenedExport() {
  if (state.screenedManagedRun !== state.result?.app.run_id) return;
  if (!state.checked.size) {
    throw new Error("Select at least one candidate for export.");
  }
}

function handleScreeningOrderChange(event) {
  const id=event.target.id;
  if(id === "screening-evidence-set") { state.screeningEvidence=event.target.value;renderOffTarget(); }
  else if(id === "screening-sort-scope") renderOffTarget(event.target.value);
  else if(id === "screening-order") {
    state.screeningSort=event.target.value;state.screeningSortDescending=false;renderOffTarget();
  } else if(id === "screening-plot-limit") {
    state.screeningPlotLimit=event.target.value === "all" ? null : Number(event.target.value);renderOffTarget();
  } else return;
  $(id)?.focus();
}
for(const id of ["offtarget-results","screening-comparisons"]) $(id)?.addEventListener("change",handleScreeningOrderChange);
document.addEventListener("keydown",event=>{
  const candidate=event.target.closest?.("svg [data-screening-candidate]");
  if(candidate && ["Enter"," "].includes(event.key)) {
    event.preventDefault();
    if(typeof viewScreeningCandidate === "function") viewScreeningCandidate(candidate.dataset.screeningCandidate,candidate.dataset.screeningScope);
  }
});

function frozenScreeningRecords() {
  return state.result?.screenings || frozen?.screening_records || frozen?.screenings || [];
}
function frozenScreeningSummary() {
  return state.result?.screening_summary || frozen?.screening_summary || [];
}
async function exportScreeningCsv() {
  try {
    validateScreenedExport();
    const r=state.result, name=`${r.gene?.symbol || "promoter"}-${r.app.run_id.slice(0,8)}-screening.csv`;
    if(frozen) {
      const ids=new Set((state.result.candidates || []).filter(c=>state.checked.size ? state.checked.has(c.id) : c.rank!==null).map(c=>c.id));
      const rows=frozenScreeningSummary().filter(row=>ids.has(row.candidate));
      if(!rows.length) throw new Error("This older offline report has no screening-summary export. Open the run in the live app and download a new report.");
      const columns=Object.keys(rows[0]), quote=v=>`"${(v===true ? "True" : v===false ? "False" : String(v??"")).replace(/"/g,'""')}"`;
      const csv=[columns,...rows.map(row=>columns.map(key=>row[key]))].map(row=>row.map(quote).join(",")).join("\r\n");
      downloadBlob(new Blob([csv],{type:"text/csv"}),name);return;
    }
    const response=await fetch(`/api/runs/${r.app.run_id}/export.screening.csv${state.checked.size ? `?ids=${encodeURIComponent([...state.checked].join(","))}` : ""}`);
    if(!response.ok) throw new Error(`Screening export failed (${response.status})`);
    downloadBlob(await response.blob(),name);
  } catch(e) { error(e.message); }
}
$("export-analysis")?.addEventListener("click",exportScreeningCsv);
$("export-plot")?.addEventListener("click",()=>{
  if(!$("screening-match-svg")) {
    if(typeof setWorkspaceTab === "function") setWorkspaceTab("screening",{focus:true});
    error("Screen candidates before downloading a sequence-match plot.");return;
  }
  downloadScreeningSvg();
});
