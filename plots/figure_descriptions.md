# Figure Descriptions

---

**Figure 00 — LMM model assumption checks.**
Four diagnostic plots for the linear mixed-effects model of colony growth (log area ~ time × species + size_z + random effects for plot and colony). Top left: residuals versus fitted values, with a LOESS smoother, assessing linearity and homoscedasticity — random scatter around zero indicates no systematic pattern. Top right: normal Q-Q plot of model residuals against theoretical quantiles — points following the reference line indicate approximate normality. Bottom left: scale-location plot (√|residuals| vs. fitted) assessing homogeneity of variance across the fitted range. Bottom right: normal Q-Q plot of colony-level random intercepts — deviations from the line indicate non-normality in the random-effects distribution.

---

**Figure 01a — Species composition of total coral cover over time.**
Stacked area chart showing mean percentage cover (% of 25 m² plot area) contributed by each coral species across survey dates. Survey dates span from the first to the most recent monitoring event across both sites. Each colour represents one species (codes: SSID = *Siderastrea siderea*, DLAB = *Diploria labyrinthiformis*, PAST = *Porites astreoides*, etc.). Stack height reflects total coral cover; the relative areas show compositional shifts over time.

---

**Figure 01b — Total coral cover trajectory per plot.**
Line plot showing total coral cover (% of 25 m²) over time for each monitored plot, coloured and labelled by site and plot number (e.g. Juanillo - P1). The accent-coloured line shows the population-level trend from a linear mixed-effects model (LMM: log(% cover) ~ time, random intercept per plot); the shaded band is the 95% confidence interval derived from fixed-effect variance. Individual plot trajectories illustrate within- and between-site variability around the overall trend.

---

**Figure 02a — Species growth trajectories.**
Faceted line plot (one panel per species, 2-row layout) showing observed mean colony area (± SE, points) and LMM-predicted trajectory (solid line) over time. Predictions are made at mean initial size (size_z = 0) for a reference site, pooling across plots. Free y-axes allow comparison of trajectory shape within each species independent of absolute size differences.

---

**Figure 02b — Per-species annual growth rates.**
Dot-and-error-bar plot of estimated annual growth rates (% change in colony area per year) derived from LMM slope coefficients (exp(β × 365) − 1). Error bars represent propagated standard error, accounting for covariance between the reference species (PAST) slope and species interaction terms. Significance stars (Wald test: * p < 0.05, ** p < 0.01, *** p < 0.001) are shown above bars. Species are ordered by growth rate.

---

**Figure 03 — Growth vs. tissue loss balance per species.**
Horizontal diverging bar chart showing mean annual gained tissue (positive, right) and lost tissue (negative, left) as a percentage of colony area per year. Each species' bar pair is derived by averaging per-colony interval-level rates, weighted equally per colony. The filled dot with error bars shows the net annual change ± SE across colonies. A vertical accent line marks zero net change. Species are ordered by net change.

---

**Figure 04a — Colony initial size vs. growth rate.**
Scatter plot of each colony's empirical annual growth rate (% per year, log-linear estimate between first and last observation) against its initial size z-score (within-species standardised log area). Points are coloured by species. The accent-coloured line shows the model-estimated size effect on growth rate based on the time_days:size_z interaction term from the LMM (reference species PAST), illustrating the overall negative relationship between initial colony size and subsequent growth.

---

**Figure 04b — Initial size class distribution per species.**
Horizontal stacked bar chart showing the proportion of observations in each size class (small / medium / large, defined as equal-width thirds of the within-species log-area range) for each species. Bars sum to 100%. Size classes are colour-coded from dark (small) to light (large).

---

**Figure 05a — Colony mortality risk by species.**
Dot-and-error-bar plot of estimated probability of complete colony death (P(died)) per species, derived from a binomial GLM (fate_binary ~ size_z + class). Estimates and 95% confidence intervals are on the response (probability) scale via emmeans. Compact letter display (Tukey-adjusted) indicates species not significantly different in mortality risk share a letter.

---

**Figure 05b — Mortality severity by initial colony size.**
Scatter plot of percentage area lost (partial or complete mortality events) against initial size z-score. Partial mortality events (area decrease between surveys) and complete deaths (colony disappears before final survey date, coded as 100% loss) are shown as jittered points. The accent-coloured line is the LMM-predicted severity for the reference species (PAST): log(% area lost) ~ size_z, with random effects for plot and colony.

---

**Figure 06a — Colony area by species and size class.**
Bubble matrix with species on the y-axis and initial size class (small / medium / large) on the x-axis. Each bubble is sized proportional to mean colony area (cm²); the number inside the bubble is the colony count; the label to the left (→) shows the mean area in cm². A dashed line separates individual species from the "All species" aggregate row at the top. Significance stars next to species codes indicate a significant pairwise difference in mean size vs. the reference species PAST (Tukey-adjusted emmeans).

---

**Figure 06b — Pairwise species size comparison (Tukey).**
Dot-and-error-bar plot of estimated geometric mean colony area (cm²) per species from the LMM, with 95% confidence intervals on the log-back-transformed scale. Species are ordered by mean area. Compact letter display (Tukey-adjusted) indicates species sharing a letter are not significantly different in mean size.

---

**Figure 07a — Total coral cover over time by site.**
Line plot of mean total coral cover (% of 25 m²) per site over time. The shaded ribbon shows ± 1 SE across plots within each site. Sites are distinguished by colour (viridis palette). Allows direct visual comparison of cover trajectories between zones.

---

**Figure 07b — Species composition by site.**
Horizontal stacked bar chart showing the proportion of unique colonies belonging to each species at each site. Proportions sum to 100% per site. Colour codes match the species palette used throughout. Allows comparison of community composition between zones.

---

**Figure 07c — Initial colony size distribution by site.**
Violin plot with overlaid box plot showing the distribution of log-transformed initial colony area (cm², at each colony's first observation) per site. The violin shows the full density; the box shows median and interquartile range (outliers suppressed). Allows comparison of the size structure of the monitored community between zones.
