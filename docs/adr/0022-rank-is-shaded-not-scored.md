# ADR 0022: Rank is shaded in one hue, never green to red

**Status:** accepted
**Date:** 2026-10-09

## Context

The top ten saRNAs need to be distinguishable in the promoter map: full lane height, a rank
numeral and a fill that separates them from other passing targets. Identical dark marks became
crowded, prompting consideration of a green-to-red colour scale.

The request is reasonable and the underlying problem is real. The difficulty is what a
green-to-red ramp asserts. Rank here is **preferences met, counted equally** (ADR 0006). On the
125 assayed sequences in chapter 4 it does not beat chance, AUC 0.463 against 0.5 for a coin
flip, p = 0.709. A scale whose convention is "green is good, red is bad" would tell a reader
that rank 1 is a better bet than rank 10, which is precisely the claim the only measurement
this project has refutes. The banner above the map already says so in words; the colours would
contradict the banner.

Two practical problems compound it. The map already uses amber for the Ensembl Regulatory Build
track and a reddish fill for promoters, so mid-rank and low-rank marks would read as regulatory
features. And red against green is the least separable pairing for the roughly 8% of men with a
colour vision deficiency, in a picture whose whole purpose is telling marks apart.

Options put to the user:

- **One hue, dark to pale.** Says "ordered", says nothing about quality.
- **Green to red**, as proposed.
- **Viridis.** Ordered, perceptually uniform, colourblind-safe, but less immediately readable
  and its bright yellow end competes with the rest of the picture.

## Decision

**Rank 1 to 10 is shaded along a single teal hue, `#0b4f57` at rank 1 to `#9ec4c6` at rank 10.**
Chosen by the user from the three options above.

- The shading encodes **order only**. The legend and the track help both say so, and the help
  says explicitly why it is not a green-to-red scale.
- Rank stays encoded on three channels at once: lane height, the numeral, and now the shade.
  Any one alone fails somebody, and the shade is the one that survives at full extent where
  numerals have to be dropped for crowding.
- The endpoints are chosen so the palest mark still reads against white. Measured: contrast
  ratios run 9.25:1 at rank 1 to 1.88:1 at rank 10, monotone, and every mark carries a darker
  edge stroke so the pale end keeps its outline.
- The ramp is **computed** from `TOP_MARKED` by one function, `rankColour`. The legend and the
  browser check both call it rather than listing shades, because a hand-copied palette is how a
  legend and a map drift apart.

## Consequences

- The ten are distinguishable from each other and from the rest at a glance, which is what was
  asked for, without the picture claiming an efficacy ordering the data does not support.
- A reader who expects green-to-red has to read the legend to learn that darker means
  better-ranked. That cost is accepted: the alternative is a familiar scale that misleads.
- Adjacent ranks are hard to separate by eye at full extent, since ten steps across one hue is
  a fine gradient on a six-pixel mark. The numeral and the tooltip carry the exact rank, and
  zooming widens the marks.
- `scripts/verify_workbench_browser.py` asserts that the marks drawn are exactly the top ten,
  that each carries the colour its own rank maps to, that all ten shades are distinct, and that
  relative luminance increases monotonically with rank. A ramp that stopped being ordered, or
  that drifted from the legend, fails the check rather than shipping.
- Nothing about the ranking itself changed. This is presentation, and ADR 0006 still governs
  what a rank is.
