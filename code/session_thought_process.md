# Analysis Session — Thought Process Timeline

**Project:** FCC Coral Health Report · **Sites:** S1P1 Juanillo · S2P2 Farallon

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

## 5 — Response variable and model family choice

**Initial approach:** `lmer(log(area) ~ ...)` — log-transforming before fitting. This is statistically sound in terms of normality but creates a systematic visualisation problem: `predict()` on the log scale gives the geometric mean (median), not the arithmetic mean. Back-transforming with `exp()` systematically underestimates observed arithmetic means, especially when random effect variances are large.

**The Jensen's inequality problem:** For a log-normal distribution, `E[Y] = exp(μ + σ²/2)`, not `exp(μ)`. We attempted to apply the correction `exp(predict + 0.5 × σ²_residual)` (Jensen's inequality), then expanded it to `exp(predict + 0.5 × (σ²_residual + σ²_colony + σ²_plot))` to include random effect variances — but this still underestimated by ~80% because `size_z` (coefficient ≈ 1.09) adds another massive variance term: `exp(β_sz² / 2) ≈ 1.8`.

**Final solution — two-level model separation:**
1. **Inference models:** Switched from `lmer(log(y))` to `glmer(y, family = Gamma(link="log"))`. Gamma GLMM with log link is the correct model for positive continuous data — it models the arithmetic mean directly on the log scale without the Jensen's inequality issue in the residuals.
2. **Visualisation models:** Separate `glm(Gamma)` fit on per-plot arithmetic means. Since the plotted points ARE per-plot arithmetic means, a GLM on those means produces trend lines that pass through the data by construction, regardless of size_z distribution or random effect variance.

This cleanly separates inferential models (individual colonies, full mixed structure, for p-values and coefficients) from visualisation models (per-plot means, for accurate trend lines in plots).

---

## 6 — Numerical stability: rescaling time

Gamma GLMMs with `time_days` (range 0–700+) interacting with binary zone predictors (0/1) produced convergence failures — the optimizer reported "Rescale variables? large eigenvalue ratio." Solution: `time_sc = time_days / 100`, used in all Gamma GLMMs. Coefficients scale accordingly (multiply by 3.65 for log-unit annual rate). The lmer models for size effects retain `time_days` since they converge without issues.

---

## 7 — Size predictor: from categories to a continuous z-score

**First approach:** split each species' initial log-area range into thirds and label colonies small/medium/large. This worked mechanically but throws away information by discretising a continuous variable.

**Revised approach:** `size_z` = z-score of log(initial area), computed within each species. This preserves the full continuous gradient, avoids arbitrary bin boundaries, and keeps the interpretation species-relative (+1 on size_z = 1 SD above species mean).

*Edge case:* species with a single colony (or all-equal areas) are assigned `size_z = 0`.

---

## 8 — Model building strategy

Built the growth GLMM with `time_sc × class + time_sc × size_z + site`. Interactions allow each species and colony size to have its own growth slope. A separate size-only LMM (`log(area) ~ time_days × size_z`) is used purely for the size-effect visualisation (02b) where the `time_days` coefficient scale is more interpretable.

---

## 9 — Partial mortality

We filled `partial_mortality` from the data itself: if a colony's area decreased between surveys, partial mortality = % area lost. If area increased or was the first observation, partial mortality = NA.

Two response variables:
- **Severity** (how much area was lost, given it occurred?) — LMM on log(% area lost)
- **Trajectory over time** — Gamma GLM on per-plot means for visualisation

---

## 10 — Complete colony mortality and trajectory filtering

Partial mortality doesn't capture colonies that disappear entirely. We detected complete mortality by comparing each colony's last survey date against its plot's last survey date — if the colony dropped out early, it died. Terminal rows with 100% mortality are added for the severity model.

**Trajectory plots (02d, 03d):** Complete mortality events (100%) were initially included in the temporal trajectory plots. However, they cluster at the start of the time series (early colony loss) and severely skew the trajectory line upward early on, misrepresenting the ongoing partial tissue-loss dynamic. These events are excluded from `pm_partial` which feeds only the trajectory plots and their visualisation models. The full `pm_events` (including complete deaths) is retained for the severity model and fate model.

---

## 11 — Power-based species sample size check

With multiple species and limited field data, some species have very few colonies. We introduced a two-tier check:

| Tier | Threshold | Action |
|------|-----------|--------|
| Hard floor | n < 5 colonies | Species removed — model parameters cannot be estimated |
| Power flag | n < 26 colonies (d=0.8, power=0.80, two-sample t-test lower bound) | Species kept, results flagged ⚠ in output |

After applying the hard floor, the final species set was: **PAST, SINT, PPOR, PSTR, SSID, OFAV, DLAB, CNAT**.

**Reference species:** PAST chosen as model reference level — highest colony count, lowest SE, ecologically most abundant and resilient.

---

## 12 — Results summary dataframe

All significant results collected into `results_sig` via:
- `extract_fixed()`: pulls fixed-effect coefficients, SE, z/t statistics and p-values from any lmer/glmer model (column detection handles both t values and z values automatically)
- A `power_note` column flags terms involving low-power species

---

## Key analytical decisions in one view

| Decision | Chosen approach | Alternative considered |
|----------|----------------|------------------------|
| Inference model family | `glmer(Gamma, log link)` | `lmer(log(y))` — back-transformation bias |
| Visualisation trend lines | `glm(Gamma)` on per-plot means | GLMM predictions with variance corrections — still biased due to size_z |
| Size predictor | Continuous `size_z` (within-species z-score of log area) | Categorical thirds — discards information |
| Time scaling | `time_sc = time_days / 100` for GLMMs | Raw `time_days` — causes eigenvalue instability |
| Complete deaths in trajectories | Excluded from 02d/03d (pm_partial) | Included — skews early trajectory by clustering at t=0 |
| Complete mortality detection | Disappearance before plot's last date | Relying on `partial_mortality` column only |
| Observation gaps | Remove colonies with absent-then-present pattern | Keep and impute — identity unreliable |
| Underpowered species | Keep with ⚠ flag (hard floor = 5 colonies) | Remove entirely — loses too many species |
| Trajectory model structure | Additive `time_sc + class` (02d) | `time_sc * class` — degenerate Hessian with sparse mortality data |

---

## Current findings (both sites)

- **S2P2 Farallon outperforms S1P1 Juanillo in colony survival:** 78% lower odds of complete mortality (OR = 0.22).
- **Cover is increasing at both sites:** Juanillo shows a steeper rate of increase, converging toward Farallon (zone × time p = 0.008).
- **Colony size is the dominant cross-model predictor:** larger initial size → larger area (×3 per SD), slower relative growth (ontogenetic convergence), lower mortality severity, lower mortality probability.
- **SSID uniquely accelerates over time** (time_sc:SSID z = 2.40, p = 0.016).
- **DLAB grows fastest at baseline but decelerates** (time_sc:DLAB z = −2.45, p = 0.014).
- **PPOR and CNAT carry the highest mortality risk** vs PAST (both p ≤ 0.040).
- **Species composition and initial size differ between sites** (χ² p < 0.001; Wilcoxon p = 0.035).
