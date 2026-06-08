// Figure 00 — standalone
#figure(
  image("figures/plots/00_lmm_assumption_checks.jpeg", width: 100%),
  caption: [
    Four diagnostic plots for the linear mixed-effects model of colony growth
    (log area \~ time × species + size\_z + random effects for plot and colony).
    Top left: residuals versus fitted values, with a LOESS smoother, assessing
    linearity and homoscedasticity — random scatter around zero indicates no
    systematic pattern. Top right: normal Q-Q plot of model residuals against
    theoretical quantiles — points following the reference line indicate
    approximate normality. Bottom left: scale-location plot (√|residuals|
    vs.\ fitted) assessing homogeneity of variance across the fitted range.
    Bottom right: normal Q-Q plot of colony-level random intercepts — deviations
    from the line indicate non-normality in the random-effects distribution.
  ],
) <fig-00_lmm_assumption_checks>

// Figures 01a & 01b
#v(0.5em)
#grid(
  columns: (1fr, 1fr),
  gutter: 0.5em,
  [#figure(
    image("figures/plots/01a_community_cover_stacked.jpeg", width: 100%),
    caption: [
      Stacked area chart showing mean percentage cover (% of 25 m² plot area)
      contributed by each coral species across survey dates. Survey dates span
      from the first to the most recent monitoring event across both sites. Each
      colour represents one species (codes: SSID = _Siderastrea siderea_, DLAB =
      _Diploria labyrinthiformis_, PAST = _Porites astreoides_, etc.). Stack
      height reflects total coral cover; the relative areas show compositional
      shifts over time.
    ],
  ) <fig-01a_community_cover_stacked>],
  [#figure(
    image("figures/plots/01b_community_cover_trajectory.jpeg", width: 100%),
    caption: [
      Line plot showing total coral cover (% of 25 m²) over time for each
      monitored plot, coloured and labelled by site and plot number (e.g.\
      Juanillo -- P1). The accent-coloured line shows the population-level trend
      from a linear mixed-effects model (LMM: log(% cover) \~ time, random
      intercept per plot); the shaded band is the 95% confidence interval derived
      from fixed-effect variance. Individual plot trajectories illustrate within-
      and between-site variability around the overall trend.
    ],
  ) <fig-01b_community_cover_trajectory>],
)

// Figures 02a, 02b & 03 — full width then pair
#v(0.5em)
#[#figure(
  image("figures/plots/02a_growth_trajectories.jpeg", width: 100%),
  caption: [
    Faceted line plot (one panel per species, 2-row layout) showing observed
    mean colony area (± SE, points) and LMM-predicted trajectory (solid line)
    over time. Predictions are made at mean initial size (size\_z = 0) for a
    reference site, pooling across plots. Free y-axes allow comparison of
    trajectory shape within each species independent of absolute size
    differences.
  ],
) <fig-02a_growth_trajectories>]
#v(0.5em)
#grid(
  columns: (1fr, 1fr),
  gutter: 0.5em,
  [#figure(
    image("figures/plots/02b_growth_rates.jpeg", width: 100%),
    caption: [
      Dot-and-error-bar plot of estimated annual growth rates (% change in colony
      area per year) derived from LMM slope coefficients (exp(β × 365) − 1).
      Error bars represent propagated standard error, accounting for covariance
      between the reference species (PAST) slope and species interaction terms.
      Significance stars (Wald test: \* p < 0.05, \*\* p < 0.01,
      \*\*\* p < 0.001) are shown above bars. Species are ordered by growth rate.
    ],
  ) <fig-02b_growth_rates>],
  [#figure(
    image("figures/plots/03_turnover_balance.jpeg", width: 100%),
    caption: [
      Horizontal diverging bar chart showing mean annual gained tissue (positive,
      right) and lost tissue (negative, left) as a percentage of colony area per
      year. Each species' bar pair is derived by averaging per-colony
      interval-level rates, weighted equally per colony. The filled dot with error
      bars shows the net annual change ± SE across colonies. A vertical accent
      line marks zero net change. Species are ordered by net change.
    ],
  ) <fig-03_turnover_balance>],
)

// Figures 04a & 04b
#v(0.5em)
#grid(
  columns: (1fr, 1fr),
  gutter: 0.5em,
  [#figure(
    image("figures/plots/04a_size_vs_growth.jpeg", width: 100%),
    caption: [
      Scatter plot of each colony's empirical annual growth rate (% per year,
      log-linear estimate between first and last observation) against its initial
      size z-score (within-species standardised log area). Points are coloured by
      species. The accent-coloured line shows the model-estimated size effect on
      growth rate based on the time\_days:size\_z interaction term from the LMM
      (reference species PAST), illustrating the overall negative relationship
      between initial colony size and subsequent growth.
    ],
  ) <fig-04a_size_vs_growth>],
  [#figure(
    image("figures/plots/04b_size_distribution.jpeg", width: 100%),
    caption: [
      Horizontal stacked bar chart showing the proportion of observations in each
      size class (small / medium / large, defined as equal-width thirds of the
      within-species log-area range) for each species. Bars sum to 100%. Size
      classes are colour-coded from dark (small) to light (large).
    ],
  ) <fig-04b_size_distribution>],
)

// Figures 05a & 05b
#v(0.5em)
#grid(
  columns: (1fr, 1fr),
  gutter: 0.5em,
  [#figure(
    image("figures/plots/05a_mortality_risk_by_species.jpeg", width: 100%),
    caption: [
      Dot-and-error-bar plot of estimated probability of complete colony death
      (P(died)) per species, derived from a binomial GLM
      (fate\_binary \~ size\_z + class). Estimates and 95% confidence intervals
      are on the response (probability) scale via emmeans. Compact letter display
      (Tukey-adjusted) indicates species not significantly different in mortality
      risk share a letter.
    ],
  ) <fig-05a_mortality_risk_by_species>],
  [#figure(
    image("figures/plots/05b_mortality_severity_by_size.jpeg", width: 100%),
    caption: [
      Scatter plot of percentage area lost (partial or complete mortality events)
      against initial size z-score. Partial mortality events (area decrease
      between surveys) and complete deaths (colony disappears before final survey
      date, coded as 100% loss) are shown as jittered points. The accent-coloured
      line is the LMM-predicted severity for the reference species (PAST):
      log(% area lost) \~ size\_z, with random effects for plot and colony.
    ],
  ) <fig-05b_mortality_severity_by_size>],
)

// Figures 06a & 06b
#v(0.5em)
#grid(
  columns: (1fr, 1fr),
  gutter: 0.5em,
  [#figure(
    image("figures/plots/06a_species_size_bubble.jpeg", width: 100%),
    caption: [
      Bubble matrix with species on the y-axis and initial size class
      (small / medium / large) on the x-axis. Each bubble is sized proportional
      to mean colony area (cm²); the number inside the bubble is the colony count;
      the label to the left (→) shows the mean area in cm². A dashed line
      separates individual species from the "All species" aggregate row at the
      top. Significance stars next to species codes indicate a significant
      pairwise difference in mean size vs.\ the reference species PAST
      (Tukey-adjusted emmeans).
    ],
  ) <fig-06a_species_size_bubble>],
  [#figure(
    image("figures/plots/06b_species_size_tukey.jpeg", width: 100%),
    caption: [
      Dot-and-error-bar plot of estimated geometric mean colony area (cm²) per
      species from the LMM, with 95% confidence intervals on the
      log-back-transformed scale. Species are ordered by mean area. Compact
      letter display (Tukey-adjusted) indicates species sharing a letter are not
      significantly different in mean size.
    ],
  ) <fig-06b_species_size_tukey>],
)

// Figures 07a, 07b & 07c — full width then pair
#v(0.5em)
#[#figure(
  image("figures/plots/07a_site_cover_trajectory.jpeg", width: 100%),
  caption: [
    Line plot of mean total coral cover (% of 25 m²) per site over time. The
    shaded ribbon shows ± 1 SE across plots within each site. Sites are
    distinguished by colour (viridis palette). Allows direct visual comparison
    of cover trajectories between zones.
  ],
) <fig-07a_site_cover_trajectory>]
#v(0.5em)
#grid(
  columns: (1fr, 1fr),
  gutter: 0.5em,
  [#figure(
    image("figures/plots/07b_site_species_composition.jpeg", width: 100%),
    caption: [
      Horizontal stacked bar chart showing the proportion of unique colonies
      belonging to each species at each site. Proportions sum to 100% per site.
      Colour codes match the species palette used throughout. Allows comparison
      of community composition between zones.
    ],
  ) <fig-07b_site_species_composition>],
  [#figure(
    image("figures/plots/07c_site_initial_size.jpeg", width: 100%),
    caption: [
      Violin plot with overlaid box plot showing the distribution of
      log-transformed initial colony area (cm², at each colony's first
      observation) per site. The violin shows the full density; the box shows
      median and interquartile range (outliers suppressed). Allows comparison of
      the size structure of the monitored community between zones.
    ],
  ) <fig-07c_site_initial_size>],
)
