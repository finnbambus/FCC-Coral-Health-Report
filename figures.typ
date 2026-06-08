// Figure 01 — full width (landscape)
#figure(
  image("figures/plots/01_species_overview.jpeg", width: 100%),
  caption: [
    Two-panel figure with a single shared legend centered at the top. Left panel:
    smoothed stacked area chart of mean % cover per species over time — stack
    height shows total coral cover, band widths show compositional change. Right
    panel: per-plot colony counts per survey date (background scatter) with
    LMM-predicted trends per species. Both panels use the same species colour
    palette.
  ],
) <fig-01_species_overview>

// Figure 02 — full width (landscape, faceted)
#v(0.5em)
#figure(
  image("figures/plots/02_growth_trajectories.jpeg", width: 100%),
  caption: [
    Faceted plot (one panel per species, 2-row layout) showing per-plot mean
    colony area ± SE (faded background pointranges, coloured by site) and a
    Gamma GLM trend line (fit on per-plot means). Free y-axes allow within-species
    comparison of trajectory shape independent of absolute size. Inference uses
    the full Gamma GLMM (model\_growth); the GLM trend lines are for
    visualisation only.
  ],
) <fig-02_growth_trajectories>

// Figure 03 — full width (landscape, faceted)
#v(0.5em)
#figure(
  image("figures/plots/03_mortality_trajectories.jpeg", width: 100%),
  caption: [
    Faceted plot (one panel per species, 2-row layout) showing per-plot mean
    partial mortality (% area lost, complete deaths excluded) ± SE and a Gamma
    GLM trend line. Complete colony deaths are excluded to reflect ongoing
    tissue-loss dynamics only. Y-axes scale to the observed data range.
  ],
) <fig-03_mortality_trajectories>

// Figures 05 & 06 — side by side (square format)
#v(0.5em)
#grid(
  columns: (1fr, 1fr),
  gutter: 0.5em,
  [#figure(
    image("figures/plots/04_growth_by_size.jpeg", width: 100%),
    caption: [
      Scatter plot of empirical annual growth rate against initial size z-score,
      coloured by species. The accent-coloured line is the LMM-estimated
      size-growth relationship (PAST reference), illustrating ontogenetic growth
      slowdown: larger colonies grow proportionally slower.
    ],
  ) <fig-04_growth_by_size>],
  [#figure(
    image("figures/plots/05_mortality_by_size.jpeg", width: 100%),
    caption: [
      Scatter plot of percentage area lost against initial size z-score for all
      mortality events. Points are coloured by event type (partial mortality vs
      complete death). The accent-coloured line is the LMM-predicted severity for
      PAST. Larger initial colonies lose significantly less proportional tissue.
    ],
  ) <fig-05_mortality_by_size>],
)

// Figure 06 — full width (landscape, multi-site)
#v(0.5em)
#figure(
  image("figures/plots/06_site_cover_abundance.jpeg", width: 100%),
  caption: [
    Two-panel figure with a shared site legend centered at the bottom. Left
    panel: per-plot total coral cover (%) over time with Gamma GLMM trend lines
    per site — Juanillo starts lower but shows a steeper rate of increase
    (zone × time p = 0.008). Right panel: per-plot colony counts over time with
    LMM trend lines per site.
  ],
) <fig-06_site_cover_abundance>

// Figure 07 — full width (landscape, multi-site)
#v(0.5em)
#figure(
  image("figures/plots/07_site_trajectories.jpeg", width: 100%),
  caption: [
    Two-panel figure with a shared site legend centered at the bottom. Left
    panel: mean colony area over time with both zones plotted on the same axes —
    Gamma GLM trend lines plus per-plot mean ± SE pointranges coloured by zone.
    Right panel: same structure for partial mortality (% area lost, complete
    deaths excluded).
  ],
) <fig-07_site_trajectories>

// Figure 08 — full width (multi-site)
#v(0.5em)
#figure(
  image("figures/plots/08_species_by_site.jpeg", width: 100%),
  caption: [
    Horizontal stacked bar chart showing the proportion of unique colonies per
    species at each site. Chi-square test confirms significant compositional
    differences between sites (χ² = 33.44, p < 0.001).
  ],
) <fig-08_species_by_site>
