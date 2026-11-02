# Visual and functional UI audit

Audited on 2 November 2026 with Chromium, cached SCN1A reference data and isolated test runs.
The fixes are included in version 0.3.2.

## Findings and fixes

- **Focused field overlaps its label:** nested labels did not separate their text from the
  input's focus outline. Nested controls now start on a new block with explicit spacing;
  label line height and focus-outline offset are defined.
- **Crowded rank numbers:** the map previously estimated spacing by centre distance, which
  could let a two-digit rank collide with its neighbour. The map now measures rendered text
  bounds and keeps the higher-priority label when space is insufficient. Every target still
  retains its rank in its tooltip and candidate record.
- **Stray existing-sequence track label:** the label was drawn even when its track was hidden,
  allowing it to share a row with the literature annotation. Label visibility now follows
  track visibility.
- **Long imported names:** reference-pill text could exceed its background. Names now fit
  the measured available width with an ellipsis; full names remain in accessible descriptions.
- **Clipped import-field focus outlines:** a nested card clipped outlines at its edges.
  The import panel now allows the outlines to remain visible.
- **History quick-open label spacing:** the saved-run dropdown needed additional separation
  from its label. This is now explicit.
- **Mobile navigation covers a heading:** a fixed scroll offset was too small when tabs
  wrapped. Scroll targets and the desktop inspector now use the navigation's measured height.

## Repeatable visual checks

```bash
.venv/bin/python scripts/verify_workbench_browser.py \
  --chrome /path/to/chrome --gene SCN1A --visual-audit \
  --screenshots out/ui-audit/visual
```

The audit checks 105 scenes across widths of 1920, 1440, 1280, 1024, 768, 390 and 320 CSS pixels.
It covers Design and expanded reference evidence; the rule editor; Candidates; each map tool;
references and literature; expanded candidate inspection; Screening; Export; and History.
It includes dense motif/regulatory tracks and a deliberately long imported sequence name.

The completed audit recorded **270 label/focus checks and 1,309 SVG text checks**, with no
reported overlaps, clipped field focus rings, out-of-bounds map text, page-wide horizontal
overflow, or headings covered by sticky navigation in those scenes. Scientific tracks and
candidate tables can scroll within their own containers on small screens.

`layout-audit.json` and screenshots are written to the supplied output directory. These are
UI verification artifacts; they do not establish biological efficacy or specificity.

The separate standard browser workflow verifies reference selection, design changes,
sequence import, selection, screening fixtures, history, exports, offline reports and mobile
layout. Run it without `--visual-audit` to exercise that complete workflow.

## Final verification after restart

- Regression suite: 2,076 passed, 11 skipped, 23 deselected.
- Full cached-SCN1A browser workflow: passed, including offline search, exports and mobile layout.
- Visual geometry audit: all 105 scenes passed.
- Ruff, strict type checks and whitespace checks: passed.

Results are retained under `out/ui-audit/`: `tests.log`, `workflow.log`, `visual.log`,
`visual/layout-audit.json` and the corresponding screenshots. The offline browser test waits
for the interaction handlers to be ready before testing search and declares its search
orientation explicitly.
