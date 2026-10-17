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
const SCOPE_NAME = {transcript:"Transcript matches (RefSeq RNA)",
                    genomic:"Genomic matches (GRCh38)"};
const BATCH_SIZES = [1, 5, 10, 15, 25];

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
  const requested = Number(record.poll_after_seconds || 5) * 1000;
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
  };
  held.searches = [...held.searches, search];
  const hits = {}, placements = {};
  for (const hit of result.hits || []) (hits[hit.candidate] = hits[hit.candidate] || []).push(hit);
  for (const p of result.proximities || []) {
    (placements[p.hit.candidate] = placements[p.hit.candidate] || []).push(p);
  }
  for (const row of result.specificity || []) {
    const previous = held.byCandidate[row.candidate];
    if (previous && new Date(previous.search.searched_at) > new Date(search.searched_at)) continue;
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

/* A paused worker is not a paused NCBI job. Check the saved RID while this run is open,
 * without creating another search or letting stale timers follow the user to a new run. */
const AUTO_SCREEN_CHECKS = new Map();
const AUTO_SCREEN_CHECK_MS = 60000;

function scheduleScreeningCheck(scope) {
  const job = state.screeningJobs?.[scope], runId = state.result?.app.run_id;
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
  // Opt-in, and the run is already readable: this adds evidence, it never gates the designs.
  const run = state.result;
  if (!run || frozen) { error("Screening needs a live run."); return; }
  state.screeningJobs = state.screeningJobs || {};
  if (["running", "waiting"].includes(state.screeningJobs[scope]?.status)) return;
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
      const outcome = await pollScreening(started.id, run.app.run_id, scope, label, index, chunks.length);
      if (outcome !== "complete") return;
    } catch (e) {
      if (state.result?.app.run_id !== run.app.run_id) return;
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
  const key = `${runId}:${scope}:${screeningId}`;
  if (SCREENING_POLLS.has(key)) return SCREENING_POLLS.get(key);
  const pending = new Promise((resolve, reject) => {
    const tick = async () => {
      if (state.result?.app.run_id !== runId) { resolve("changed-run"); return; }
      try {
        const record = await api(`/api/screenings/${screeningId}`);
        if (state.result?.app.run_id !== runId) { resolve("changed-run"); return; }
        if (record.status === "running") {
          state.screeningJobs[scope] = {
            status:"running", stage:`${label}: ${record.stage || "working"}`,
            done:index, total, progress:record.progress || null, id:screeningId,
          };
          renderScreeningStatus();
          setTimeout(tick, localScreeningDelay(record));
          return;
        }
        if (record.status === "waiting") {
          // Not finished and not failed. NCBI still holds the search, so the interface offers
          // to go back for it rather than reporting a failure and losing the job.
          state.screeningJobs[scope] = {
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
        state.screeningJobs[scope] = {status:"complete", stage:"Search finished"};
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
  const job = state.screeningJobs?.[scope];
  const runId = state.result?.app.run_id;
  if (!job?.id || job.status !== "waiting") return;
  cancelScreeningCheck(scope);
  state.screeningJobs[scope] = {...job, status:"running", stage:"Checking NCBI again"};
  renderScreeningStatus();
  try {
    await api(`/api/screenings/${job.id}/resume`, {
      method:"POST", headers:{"Content-Type":"application/json"}, body:JSON.stringify({}),
    });
    if (state.result?.app.run_id !== runId) return;
    await pollScreening(job.id, runId, scope, "collecting", 0, 1);
  } catch (e) {
    if (state.result?.app.run_id !== runId) return;
    if (["screening-already-collecting", "screening-already-complete"].includes(e.code)) {
      // Another browser tab claimed the same RID. Follow its worker instead of submitting
      // or falsely labelling the search failed.
      try {
        await pollScreening(job.id, runId, scope, "collecting", 0, 1);
      } catch (pollError) {
        if (state.result?.app.run_id !== runId) return;
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
  const current = () => state.result?.app.run_id === runId;
  try {
    const listed = await api(`/api/runs/${runId}/screenings`);
    if (!current()) return;
    let any = null;
    for (const entry of listed.filter(s => s.status === "complete").reverse()) {
      const record = await api(`/api/screenings/${entry.id}`);
      if (!current()) return;
      mergeScreening(record);
      any = entry.scope;
    }
    state.screeningJobs = state.screeningJobs || {};
    const active = listed.filter(s => ["running", "waiting"].includes(s.status));
    for (const entry of [...active, ...listed]) {
      if (!current()) return;
      if (state.screeningJobs[entry.scope]) continue;
      state.screeningJobs[entry.scope] = {
        status:entry.status, id:entry.id, stage:entry.error || entry.stage,
        retrievalIssue:entry.stage === "Report retrieval needs retry",
      };
      if (entry.status === "running") {
        pollScreening(entry.id, runId, entry.scope, "restored search", 0, 1).catch(e => {
          if (!current()) return;
          state.screeningJobs[entry.scope] = {status:"error", id:entry.id, stage:e.message};
          renderScreeningStatus();
        });
      }
    }
    renderOffTarget(any || state.screeningScope || "transcript");
    renderScreeningStatus();
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
  const jobs = state.screeningJobs || {};
  const states = Object.values(jobs).map(j => j.status);
  const done = Object.keys(screenings()).length;
  const pool = screenablePool();
  const coverage = Object.keys(SCOPE_BUTTON).map(scope => {
    const count = alreadyScreened(scope).size;
    const job = jobs[scope];
    const status = job?.status === "running" ? "search running"
      : job?.status === "waiting" ? "waiting at NCBI"
      : job?.status === "error" ? "search failed" : count ? "results available" : "not started";
    return `<div class="coverage-item"><span>${scope === "transcript" ? "Transcript" : "Genomic"} search</span>`
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
    ? "Searches run at NCBI and report here when they finish. You can keep working meanwhile."
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
      button.disabled = ["running", "waiting"].includes(job?.status) || !remaining;
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
  const busy = ["running", "waiting"].includes(state.screeningJobs?.[scope]?.status);
  const disabled = busy ? " disabled" : "";
  const ticked = pool.filter(id => state.checked.has(id)).length;
  screeningHTML(host,
    `<div class="batch-row"><span class="batch-label">${scope === "genomic" ? "Genomic" : "Transcript"} batch size</span>`
    + `<span class="segmented">`
    + BATCH_SIZES.map(n => `<button type="button" data-batch="${n}"`
        + `${n === size ? ' class="active"' : ""} aria-pressed="${n === size}">${n}</button>`).join("")
    + `</span>`
    + `<button type="button" data-batch-mode="all" class="secondary"${disabled}`
    + ` ${otTip("Screen every candidate that passed the gates. Submitted in sequential jobs of "
        + "25, each its own dated record, so NCBI is never asked for all of them at once.")}>`
    + `${done.size ? "Re-screen all" : "Search all"} ${pool.length}</button>`
    + (done.size && remaining > 0
        ? `<button type="button" data-batch-mode="rest" class="secondary"${disabled}>`
          + `Search remaining ${remaining}</button>` : "")
    + (ticked
        ? `<button type="button" data-batch-mode="selection" class="secondary"${disabled}>`
          + `Search selected (${ticked})</button>` : "")
    + `</div>`
    + (SCOPE_WARNING[scope] ? `<p class="hint">${escapeHTML(SCOPE_WARNING[scope])}</p>` : "")
    + `<p class="hint">${busy
        ? "A submitted search is pending. No new search is needed; completed results will appear here."
        : !done.size
        ? `Nothing screened yet in this scope. Start with the candidates you plan to review. `
          + `NCBI queue time varies; a smaller batch does not guarantee a faster result.`
        : remaining > 0
        ? `${done.size} of ${pool.length} screened in this scope. The next search covers `
          + `${Math.min(size, remaining)} more unless you choose otherwise.`
        : `All ${pool.length} screened in this scope. Screening again searches the same `
          + `sequences against whatever the database holds today, as a second dated record.`
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
  if (!held) {
    const pending = ["running", "waiting"].includes(state.screeningJobs?.[state.screeningScope]?.status);
    host.innerHTML = scopeTabs() + `<p class="ot-empty hint">${pending
      ? "Search submitted. Results will appear here when collection finishes; there are no findings to interpret yet."
      : "No results recorded in this scope yet."}</p>`;
    renderBatchControls(); return;
  }
  const methodsOpen = host.querySelector?.(".ot-methods")?.open;
  const notesOpen = host.querySelector?.(".ot-record-notes")?.open;
  const core = state.result?.candidates?.[0]?.target_site_dna?.length || 19;
  const entries = Object.values(held.byCandidate);
  const genomic = state.screeningScope === "genomic";
  host.innerHTML =
    scopeTabs()
    + summaryLine(held, entries, core, genomic)
    + notesFor(held)
    + resultsTable(held, entries, core, genomic)
    + (genomic ? tssPlot(entries) : "");
  if (methodsOpen) host.querySelector(".ot-methods").open = true;
  if (notesOpen && host.querySelector(".ot-record-notes")) host.querySelector(".ot-record-notes").open = true;
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
  const exact = entries.filter(e => e.row.perfect_hits > 0).length;
  const quiet = entries.filter(e => e.row.hit_count === 0).length;
  const threshold = entries[0]?.row.strong_threshold || Math.ceil(core * 0.8);
  const capped = entries.filter(e => e.row.alignments_capped).length;
  const dates = held.searches.map(s => new Date(s.searched_at)).filter(d => !isNaN(d));
  const when = dates.length ? new Date(Math.max(...dates)).toLocaleDateString() : "date unavailable";
  const databases = [...new Set(held.searches.map(s => s.parameters.database))].join(", ");
  const byDatabase = held.searches.every(s => s.parameters.species_restriction === "database");
  const returned = entries.reduce((total,e) => total + e.row.hit_count, 0);
  const matching = entries.reduce((total,e) => total + (e.row.wang_homology_hits || 0), 0);
  const pool = screenablePool().length;
  const limit = capped === 1 ? screeningLimit(entries.find(e=>e.row.alignments_capped))
    : capped ? `${capped} candidate result lists may be incomplete. Each candidate shows its search limit.` : "";
  return `<div class="ot-summary"><p><b>${escapeHTML(SCOPE_NAME[state.screeningScope])}</b> · `
    + `latest ${escapeHTML(when)}</p>`
    + `<div class="ot-result-flow" aria-label="Search coverage and returned sequence matches">`
    + `<div><strong>${entries.length}<small> / ${pool}</small></strong><span>saRNA candidates searched</span>`
    + `<small>${Math.max(0,pool-entries.length)} not yet searched in this scope</small></div>`
    + `<div><strong>${returned}</strong><span>alignments returned</span><small>Includes weaker sequence matches</small></div>`
    + `<div><strong>${matching}</strong><span>alignments with ≥${threshold}/${core} matching bases</span>`
    + `<small>A reporting threshold, not a biological risk score</small></div></div>`
    + `<p class="hint">${exact} candidate${exact === 1 ? "" : "s"} with full-length exact matches · `
    + `${quiet} with no reported alignments. Binding, gene-expression effects and specificity remain untested.</p>`
    + (limit ? `<p class="ot-limit"><b>Result list may be incomplete.</b> ${escapeHTML(limit)}</p>` : "")
    + `<details class="ot-methods"><summary>Methods, evidence and search limits</summary>`
    + `<p>Database: ${escapeHTML(databases)}. ${byDatabase
      ? "The genomic database contains human GRCh38 sequence."
      : "Human filtering was requested using txid9606. The Entrez filter is not listed in the current supported URL API; species restriction is not independently verified."}</p>`
    + `<p>The ≥${threshold}/${core} column counts matching bases, including partial or gapped alignments. `
    + `It is a reporting threshold, not a validated saRNA risk cutoff. Applying 80% to candidate `
    + `off-targets is a project interpretation; Wang 2015 discusses it for control duplexes. `
    + `<a href="https://doi.org/10.14440/jbm.2015.39" target="_blank" rel="noopener">Wang et al., 2015</a>.</p>`
    + `<p>Counts describe retained alignments, not affected genes. The server returns at most `
    + `${held.searches[0]?.parameters.max_hits || 100} database sequences per query, and the client retains `
    + `at most ${held.searches[0]?.parameters.max_alignments_per_candidate || 200} alignments. `
    + `A limit reached means counts may be lower bounds. A negative result is not proof of specificity.</p>`
    + `<p>This is a sequence-similarity search, not a comprehensive seed-risk screen. Seed-mediated `
    + `off-target silencing has experimental evidence in siRNA; that evidence does not validate `
    + `prediction of saRNA activation. <a href="https://doi.org/10.1261/rna.25706" target="_blank" `
    + `rel="noopener">Jackson et al., 2006</a>. No biological risk score is assigned.</p>`
    + `<p>Transcript chromosomes are optional Ensembl gene-symbol annotations, not mappings of `
    + `transcript alignment coordinates to the genome. Descriptions are source annotations, `
    + `including predicted functions, rather than functional validation.</p></details></div>`;
}

function notesFor(held) {
  const notes = [...new Set(held.searches.flatMap(s => s.notes))];
  const important = notes.filter(n => /unavailable|not placed|weaker ones/i.test(n) && !/reporting limit|floors/i.test(n));
  const other = notes.filter(n => !important.includes(n));
  return important.map(n => `<p class="hint warn">${escapeHTML(n)}</p>`).join("")
    + (other.length ? `<details class="ot-record-notes"><summary>Search record notes</summary>`
      + other.map(n => `<p class="hint">${escapeHTML(n)}</p>`).join("") + `</details>` : "");
}

function resultsTable(held, entries, core, genomic) {
  const sort = state.screeningSort || "wang_homology_hits";
  const value = (e) => {
    const v = sort === "distinct_loci" ? (e.row.distinct_loci ?? distinctLoci(e.hits))
      : sort === "distinct_sequences" ? (e.row.distinct_sequences ?? new Set(e.hits.map(h=>h.accession.split(".")[0])).size)
      : e.row[sort];
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
    ? entry.placements.find(p => p.hit.accession === best.accession && p.hit.subject_start === best.subject_start
      && p.hit.subject_end === best.subject_end && p.hit.subject_strand === best.subject_strand) : null;
  const open = state.screeningOpen === row.candidate;
  const flag = (n) => n ? ' class="match-count"' : "";
  const main = `<tr data-screenrow="${escapeHTML(row.candidate)}"`
    + ` tabindex="0" role="button" aria-expanded="${open}"${open ? ' class="open"' : ""}>`
    + `<td>${escapeHTML(row.candidate)}`
    + (row.alignments_capped ? `<p class="candidate-id" ${otTip(
        "A server or client reporting limit was reached. Counts may be lower bounds.")
      }>list may be incomplete</p>` : "")
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
    + `<p class="fp-note">DNA alphabet: T corresponds to U in RNA. Subject shown in BLAST alignment orientation `
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
    return chromosome ? `Chr ${chromosome}` : (p?.on_chromosome ? `Chr ${escapeHTML(p.contig)}` : "Unmapped contig");
  }
  const annotation = entry.search.annotations[group.accession];
  if (!annotation) return frozen ? "Not recorded" : "Loading…";
  return annotation.status === "available"
    ? `<span ${otTip(`${annotation.source}; gene context, not transcript-alignment coordinates. ${annotation.retrieved_at || ""}`)}>Chr ${escapeHTML(annotation.chromosome)}<small>gene annotation</small></span>`
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
      : p.on_chromosome === false ? "Contig not mapped to an assembled chromosome"
      : p.distance_to_tss === null ? "No annotated TSS within 10 kb"
      : `${signed(p.distance_to_tss)} bp from ${p.nearest_gene || p.nearest_transcript || "nearest TSS"}`
        + (p.inside_run_window ? " · inside this run’s window" : "");
    return `<article class="ot-hit"><div class="ot-hit-info"><div class="ot-hit-heading">`
      + `<b>${escapeHTML(hitName(strongest,p))}</b><span>${chromosomeLabel(group,entry,genomic)}</span></div>`
      + `<p class="candidate-id">${escapeHTML(group.accession)}`
      + `${group.hits.length > 1 ? ` · ${group.hits.length} alignments` : ""}</p>`
      + `<p class="ot-description">${escapeHTML(group.title || "Description unavailable")}</p>`
      + `<p class="ot-hit-metrics"><b>${strongest.identities}/${core}</b> bases matched · `
      + `${escapeHTML(strongest.subject_strand)} subject strand · E-value ${Number(strongest.evalue).toExponential(1)}</p>`
      + (genomic ? `<p class="hint">${p?.on_chromosome ? `Chr ${escapeHTML(p.contig)}:${fmt(p.position)} · ` : ""}${escapeHTML(context)}</p>` : "")
      + `</div><div class="ot-hit-alignment">${footprint(strongest,core)}</div></article>`;
  }).join("");
  const missingAnnotations = !genomic && Object.values(entry.search.annotations).some(a=>a.status === "unavailable");
  const threshold = entry.row.strong_threshold || Math.ceil(core*0.8);
  const matching = entry.row.wang_homology_hits ?? 0;
  const weak = Math.max(0,entry.row.hit_count-matching);
  return `<tr class="ot-detail"><td colspan="${genomic ? 8 : 7}">`
    + `<p class="ot-candidate-explanation"><b>1 saRNA candidate searched → ${entry.row.hit_count} alignments returned → `
    + `${matching} with ≥${threshold}/${core} matching bases.</b> `
    + `${weak ? `${weak} weaker alignment${weak === 1 ? " is" : "s are"} also included in the returned list.`
      : "All returned alignments meet this reporting threshold."}</p>`
    + `<p class="hint">These are sequence findings, not confirmed biological targets.</p>`
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
    + `Darker marks reach the reporting threshold${window ? ", and the shaded band is the window this run "
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
    const entry = screenings()[state.screeningScope]?.byCandidate[state.screeningOpen];
    if (entry && state.screeningScope === "transcript") {
      loadHitAnnotations(entry, visibleHitGroups(entry).slice(0,5).map(g=>g.accession));
    }
    return;
  }
  const sort = event.target.closest("[data-screensort]");
  if (sort) { state.screeningSort = sort.dataset.screensort; renderOffTarget(); return; }
  if (event.target.closest("[data-expand-all]")) {
    state.screeningExpandAll = true; renderOffTarget();
    const entry = screenings()[state.screeningScope]?.byCandidate[state.screeningOpen];
    if (entry && state.screeningScope === "transcript") loadHitAnnotations(entry, visibleHitGroups(entry).map(g => g.accession));
    return;
  }
  if (event.target.closest("[data-retry-annotations]")) {
    const entry = screenings()[state.screeningScope]?.byCandidate[state.screeningOpen];
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
      const entry = screenings().transcript?.byCandidate[id];
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
