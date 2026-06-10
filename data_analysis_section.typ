// ============================================================
// Section 2.4 — Data Analysis
// Paste this block into your main .typ document after section 2.3
// ============================================================

== Data Analysis

#block(
  fill: rgb("#f5f5ea"),
  inset: 10pt,
  radius: 4pt,
  width: 100%,
)[
  *Data Availability* \
  The R code (including detailed documentation) and raw data used in this project,
  as well as the Typst scripts that compile this PDF document, can be found in the
  project's GitHub repository under: \
  #align(center)[#link("https://github.com/finnbambus/FCC-Coral-Health-Report")]
]

All analysis was performed in R using the packages `lme4`, `lmerTest`, `DHARMa`,
`dplyr`, `tidyr`, `ggplot2`, `patchwork`, `viridisLite`, `scales`, `pwr`, and `here`.

=== Data Preparation

Exported `.csv` files from both sites were combined and cleaned prior to analysis.
Because `genet_id` values are assigned locally per plot, globally unique colony
identifiers were constructed by concatenating site, plot number, and `genet_id`
(e.g. `S1P1_1_G42`). Where a colony split into fragments or merged with another
between surveys, multiple rows were present for the same colony and date; these events
were detected as changes in per-colony entry count over time, flagged in a notes field,
and resolved by summing fragment areas so that each colony had a single area value per
survey date.

Colonies that disappeared for one or more surveys within their active monitoring span
and reappeared at a later date were excluded, as identity across such gaps cannot be
confirmed. Expected survey dates for each colony were derived from all dates recorded
for its plot; any colony with fewer observed dates than expected within its
first-to-last span was removed.

Partial mortality was calculated directly from the area time series: if a colony's area
decreased between consecutive surveys, partial mortality was recorded as the percentage
area lost; otherwise the value was set to missing. Colony fate was determined by
comparing each colony's last observed date against the final survey date of its plot —
colonies that disappeared before the plot's last survey were classified as having died.

=== Species Filtering

Species with fewer than five colonies across both sites were excluded because model
parameters cannot be reliably estimated at that sample size. For the remaining species,
a power threshold was calculated using a two-sample _t_-test power analysis (`pwr`
package; _d_ = 0.8, _α_ = 0.05, 1 − _β_ = 0.80), yielding a minimum of 26 colonies
per species for adequate power. Species below this threshold were retained but their
results are flagged as low-power throughout. _Porites astreoides_ (PAST) was set as the
model reference level, as it had the highest colony count and the lowest standard error
across all survey dates.

=== Statistical Models

All longitudinal models account for the nested structure of the data — colonies
repeatedly measured within plots — using random effects. Time was expressed in days
from the first survey date and, for Generalized Linear Mixed Models (GLMMs), rescaled
by dividing by 100 (`time_sc`) to avoid numerical instability from large predictor
ratios. Initial colony size was expressed as a within-species z-score of
log-transformed initial area (`size_z`), which places each colony relative to its
species mean and avoids arbitrary size-category boundaries.

*Species overview.* Mean percentage cover per species was smoothed using
`smooth.spline` (5 degrees of freedom) fitted to per-survey arithmetic means and
predicted on a weekly grid for visualization. Colony abundance trends were modelled
with a Poisson GLMM (log link):

#block(inset: (left: 1.5em))[
  `n_colonies ~ time_sc × class + (1 | plot_id)`
]

where per-plot colony counts are the response and plot identity is a random intercept.
Counts were completed to a full plot × survey-date × species grid before fitting, so
that a species absent from a surveyed plot on a given date enters as a structural zero
rather than a missing row; the trend therefore reflects abundance across all plots, not
abundance conditional on presence. A count family with a log link also constrains fitted
values to be non-negative.

*Colony growth.* The primary inference model for colony area was a Gamma GLMM with log
link:

#block(inset: (left: 1.5em))[
  `area ~ time_sc × class + time_sc × size_z + site + (1 | plot_id) + (1 | colony_id)`
]

A Gamma family was chosen over log-transforming the response because the log link
models the arithmetic mean directly, avoiding the systematic back-transformation bias
introduced by Jensen's inequality when random-effect variances are large. Fitting used
the `bobyqa` optimizer (maximum 200,000 iterations). This model, and the size–growth
model below, were fitted to follow-up observations only — each colony's baseline survey
was excluded from the response. Because `size_z` is derived from a colony's first
recorded area, retaining the baseline would let a colony partly predict its own initial
value and mechanically inflate the size terms (a regression-to-the-mean artifact);
fitting on subsequent surveys separates the predictor from the response. For trajectory
visualization only, a separate Gamma GLM was fitted to per-plot arithmetic means
(`mean_area ~ time_sc × class`); because the plotted points are themselves per-plot
means, this GLM produces trend lines that pass through the data by construction and are
not subject to random-effect variance inflation. The size–growth relationship was
additionally quantified with an LMM on log-transformed area
(`log(area) ~ time_days × size_z + (1 | plot_id) + (1 | colony_id)`), fitted on the same
follow-up observations, to obtain annual growth rate coefficients on an interpretable
scale.

*Partial mortality.* Partial mortality severity (percentage area lost, given that a
mortality event occurred) was modelled as:

#block(inset: (left: 1.5em))[
  `log(partial_mortality) ~ size_z + class + (1 | colony_id) + (1 | plot_id)` (LMM)
]

This model was fitted to partial-mortality events only. Terminal complete-mortality
events, which are coded as 100% loss, were excluded because a whole-colony death is a
distinct process and a ceiling value that would bias the severity slope; they are shown
separately in the corresponding figure. Temporal trends in partial mortality were not
formally tested but were summarised for visualization with a Gamma GLM fitted to
per-plot arithmetic means (`mean_pm ~ time_sc + class`). Complete-mortality events were
likewise excluded from these trajectory plots because they cluster early in the time
series and would bias the trend line away from the ongoing tissue-loss dynamic.

*Site comparison.* When data from both sites were present, the following additional
models were fitted: total coral cover by site
(`pct_cover ~ time_sc × zone + (1 | plot_id)`, Gamma GLMM); colony abundance by site
(`n_colonies ~ time_sc × zone + (1 | plot_id)`, Poisson GLMM); and colony growth by site
(`area ~ time_sc × zone + (1 | colony_id)`, Gamma GLMM). Temporal trends in partial
mortality by site were summarised for visualization with a Gamma GLM on per-plot means
(`mean_pm ~ time_sc × zone`), consistent with the species-level trajectory plots.
Differences in species composition between sites were tested with a chi-square test of
independence on colony counts per species per site; because several species have small
expected counts, a Monte-Carlo simulated _p_-value (10⁵ replicates) was used in place of
the asymptotic approximation. Differences in initial colony size were tested with a
Wilcoxon rank-sum test on log-transformed initial area. Colony fate (the probability
that a colony died completely) was modelled with a logistic GLM:
`P(died) ~ size_z + class + site`.

=== Model Validation

Assumptions of the primary Gamma GLMM were evaluated using simulation-based scaled
residuals generated by the `DHARMa` package (500 simulations). Residuals were inspected
visually via a residuals-versus-fitted plot, a uniform quantile–quantile plot, and a
Pearson residuals plot. Random-effect normality was assessed with a normal Q-Q plot of
colony-level random intercepts. Visualization GLMs fitted to per-plot means are used
solely for plotting and were not subject to formal inference checks.
