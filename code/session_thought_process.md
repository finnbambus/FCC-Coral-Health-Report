# Analysis Session — Thought Process Timeline

**Project:** FCC Coral Health Report · **Site:** Farallon (S2P2)

---

## 1 — Starting point: what do we have?

Opened the project with photoquadrat survey data: species, colony area, plot, site, date. The first question was simply *what analysis makes sense here?* The data is longitudinal (repeated measures on the same colonies), multi-species, and has a natural hierarchy (colonies nested in plots). That immediately pointed toward **Linear Mixed Models** as the backbone.

---

## 2 — Unique colony IDs

The dataset uses `genet_id` locally per plot — the same number in two plots means two different physical colonies. We created `colony_id = site + plot_number + genet_id` to make every colony globally unique. This was a prerequisite for correct random-effect grouping.

---

## 3 — Split and merge detection

Some colonies split into fragments or merged between surveys, causing multiple rows per `colony_id × date`. We detected these as `n_entries` changes over time, flagged them as "split"/"merged" events, and summed areas within a survey date so each colony has one row per time point.

---

## 4 — Observation gap removal

A colony that disappears for one survey and reappears later has an ambiguous identity — we can't confirm it is the same organism. These colonies are removed by comparing each colony's observed survey dates (within its first–last active span) against all survey dates recorded for the same plot. Any colony missing an expected date is excluded.

---

## 5 — Response variable: area vs log(area)

Initial models used raw area as the response. Assumption checks showed a strong funnel pattern in residuals-vs-fitted and heavy-tailed QQ residuals — classic signs of right-skewed, heteroscedastic data. Switching to `log(area)` resolved both issues and is ecologically appropriate (growth is multiplicative).

---

## 6 — Size predictor: from categories to a continuous z-score

**First approach:** split each species' initial log-area range into thirds and label colonies small/medium/large. This worked mechanically but throws away information by discretising a continuous variable.

**Revised approach:** `size_z` = z-score of log(initial area), computed within each species. This preserves the full continuous gradient, avoids arbitrary bin boundaries, and keeps the interpretation species-relative (a +1 on size_z means 1 SD above the species mean, regardless of species). `size_class` labels are retained for descriptive plots only.

The ecological rationale: coral growth is fundamentally continuous — a single coefficient for the size_z slope is more powerful and more interpretable than two category contrasts.

*Edge case:* species with a single colony (or all-equal areas) are assigned `size_z = 0` (species mean).

---

## 7 — Model building strategy

Built the base LMM with `time_days × class + time_days × size_z`. The interaction terms allow each species and size class to have its own growth slope. Ran assumption checks first, confirmed log-transform made them pass, then interpreted results.

Because `size_z` was significant, the crossed model (`time_days × class × size_z`) was also run conditionally. It revealed that SSID has the most extreme size stratification (interaction +1.72), and that PAST and PSTR also have significantly steeper size gradients than CNAT.

---

## 8 — Post-hoc all-vs-all species comparisons

The dummy-coded model only compares each species to CNAT (the reference). To understand *all* pairwise species differences we added `emmeans` + `multcomp::cld()` (Compact Letter Display). This also handles the Tukey multiple-comparisons correction automatically. With `size_z` continuous, no post-hoc is needed for size — the slope coefficient captures the full effect.

**Finding that prompted the next step:** significant size_z effects raised the question of whether smaller colonies are also more likely to *lose tissue* or *die*.

---

## 9 — Partial mortality

We filled the `partial_mortality` column from the data itself: if a colony's area decreased between surveys, partial mortality = % area lost. If area increased or was the first observation, partial mortality = NA.

This created two new response variables:
- **Occurrence** (binary: did partial mortality happen at all?) — GLMM, binomial
- **Severity** (how much area was lost, given it occurred?) — LMM on log(% area lost)

---

## 10 — Complete colony mortality

Partial mortality doesn't capture colonies that disappear entirely. We detected complete mortality by comparing each colony's last survey date against its plot's last survey date — if the colony dropped out early, it died. Added these as terminal rows with 100% mortality for the severity model.

This produced a third response variable: **colony fate** (binary: died vs alive) — GLMM, binomial.

---

## 11 — Power-based species sample size check

With multiple species and limited field data, some species have very few colonies. We introduced a two-tier check:

| Tier | Threshold | Action |
|------|-----------|--------|
| Hard floor | n < 5 colonies | Species removed — model parameters cannot be estimated |
| Power flag | n < 26 colonies (d=0.8, power=0.80, two-sample t-test lower bound) | Species kept, results flagged ⚠ in output |

The two-sample t-test power analysis is a *conservative* lower bound — LMMs with repeated measures gain power from having multiple observations per colony, so the actual required n is lower. Using it as a flag rather than a strict filter allows the analysis to proceed while keeping readers informed about which results have limited replication.

After applying the hard floor, the final species set was: **CNAT, PAST, PSTR, SSID, SINT**.

**Reference species switch (CNAT → PAST):** PAST was chosen as the model reference level because it has the highest colony count (lowest SE on all contrasts) and is ecologically the most common species at this site. Using PAST as reference makes the intercept more stable and gives more interpretable contrasts — other species are compared against the most abundant and resilient species rather than the most mortality-prone one (CNAT).

---

## 12 — Results summary dataframe

All significant results are collected into `results_sig` via two helper functions:
- `extract_fixed()`: pulls fixed-effect coefficients, SE, t/z statistics and p-values from any lmer/glmer model
- `extract_pairs()`: pulls Tukey pairwise contrasts from emmeans objects

A `power_note` column flags any term whose name contains a low-power species.

---

## Key analytical decisions in one view

| Decision | Chosen approach | Alternative considered |
|----------|----------------|------------------------|
| Response variable | `log(area)` | Raw area — failed assumption checks |
| Size predictor | Continuous `size_z` (within-species z-score of log area) | Categorical thirds — discards information |
| Size model structure | General size_z effect first; crossed (species × size_z) if significant | Species-nested — premature complexity |
| Species comparison | emmeans + Tukey CLD for all pairwise | Dummy-coded contrasts vs reference only |
| Complete mortality | Disappearance before plot's last date | Relying on `partial_mortality` column only |
| Observation gaps | Remove colonies with absent-then-present pattern | Keep and impute — identity unreliable |
| Underpowered species | Keep with ⚠ flag (hard floor = 5 colonies) | Remove entirely — loses too many species |
| Random effects | `(1\|plot_id) + (1\|colony_id)` | `(1\|colony_id)` only — misses plot-level clustering |

---

## Current findings (Farallon S2P2, post power-filter species)

- **Size is the dominant driver** across all outcomes: larger initial size → larger colonies, slower relative growth, lower partial mortality severity, lower complete mortality risk.
- **SSID is the most size-variable species** in relative terms (size_z interaction +1.72 in crossed model).
- **PAST is the most resilient** — significantly lower mortality probability than CNAT despite similar size.
- **All post-hoc contrasts and several species effects are low-power** — the patterns are ecologically consistent and effect sizes are substantial, but wider replication is needed for confirmatory inference.
