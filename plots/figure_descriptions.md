# Figure Descriptions

Figures are numbered by their order of export in `code/analysis.r`. Figures 04–08 are only produced when `multi_site = TRUE`.

---

**Figure 01 — Species cover and abundance overview.**
Two-panel landscape figure with a single shared legend centered at the top. Left panel: smoothed stacked area chart of mean % cover per species over time — stack height shows total coral cover, band widths show compositional change. Right panel: per-plot colony counts per survey date (background scatter) with LMM-predicted trends per species. Both panels use the same species colour palette; the colour legend from the right panel is collected and centered at the top of the combined figure. Landscape format (1890 × 945 px).

---

**Figure 02 — Colony growth trajectories per species.**
Faceted plot (one panel per species, 2-row layout) showing per-plot mean colony area ± SE (faded background pointranges, coloured by site) and a Gamma GLM trend line (fit on per-plot means). Free y-axes allow within-species comparison of trajectory shape independent of absolute size. Inference uses the full Gamma GLMM (model_growth); the GLM trend lines are for visualisation only. Landscape format (1890 × 945 px).

---

**Figure 03 — Partial mortality trajectories per species.**
Faceted plot (one panel per species, 2-row layout) showing per-plot mean partial mortality (% area lost, complete deaths excluded) ± SE and a Gamma GLM trend line. Complete colony deaths excluded to reflect ongoing tissue-loss dynamics only. Y-axes scale to the observed data range. Landscape format (1890 × 945 px).

---

**Figure 04 — Site cover and colony abundance by site.**
Two-panel landscape figure with a shared site legend centered at the bottom. Left panel: per-plot total coral cover (%) over time with Gamma GLMM trend lines per site — Juanillo starts lower but shows a steeper rate of increase (zone × time p = 0.008). Right panel: per-plot colony counts over time with LMM trend lines per site. Landscape format (1890 × 945 px).

---

**Figure 05 — Colony growth rate by initial size.**
Scatter plot of empirical annual growth rate against initial size z-score, coloured by species. The accent-coloured line is the LMM-estimated size-growth relationship (PAST reference), illustrating ontogenetic growth slowdown: larger colonies grow proportionally slower. Square format (945 × 945 px).

---

**Figure 06 — Mortality severity by initial colony size.**
Scatter plot of percentage area lost against initial size z-score for all mortality events. Points are coloured by event type (partial mortality vs complete death). The accent-coloured line is the LMM-predicted severity for PAST. Larger initial colonies lose significantly less proportional tissue. Square format (945 × 945 px).

---

**Figure 07 — Colony growth and partial mortality trajectories by site.**
Two-panel landscape figure with a shared site legend centered at the bottom. Left panel: mean colony area over time with both zones plotted on the same axes — Gamma GLM trend lines plus per-plot mean ± SE pointranges coloured by zone. Right panel: same structure for partial mortality (% area lost, complete deaths excluded). Landscape format (1890 × 945 px).

---

**Figure 08 — Species composition by site.**
Horizontal stacked bar chart showing the proportion of unique colonies per species at each site. Chi-square test confirms significant compositional differences (χ² = 33.44, p < 0.001). Square format (945 × 945 px).
