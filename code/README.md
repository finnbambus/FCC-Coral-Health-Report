# FCC Coral Health Report — Statistical Results Summary

**Sites:** S1P1 Juanillo · S2P2 Farallon · **Data:** `combined_data.csv`
**Species in analysis (n colonies):** PAST 229 · SINT 17 · PPOR 16 · PSTR 16 · SSID 13 · OFAV 7 · DLAB 6 · CNAT 5
**Excluded:** species with < 5 colonies (PCLV, SRAD removed).
**Note:** Species with n < 26 flagged ⚠ (below d = 0.8, 1−β = 0.80 power threshold) — interpret with caution.
**Reference level:** *Porites astreoides* (PAST) — largest colony count, lowest SE.

---

## Model structure

| Model | Type | Response | Fixed effects | Random effects |
|-------|------|----------|---------------|----------------|
| GLMM: colony growth | `glmer` (Gamma) | `area` [log link] | `time_sc × class + time_sc × size_z + site` | `(1\|plot_id) + (1\|colony_id)` |
| LMM: growth by size | `lmer` | `log(area)` | `time_days × size_z` | same |
| LMM: mortality severity | `lmer` | `log(% area lost)` | `size_z + class` | `(1\|colony_id) + (1\|plot_id)` |
| GLM: colony fate | `glm` (binomial) | P(colony died) | `size_z + class + site` | — |
| GLMM: site cover | `glmer` (Gamma) | `pct_cover` [log link] | `time_sc × zone` | `(1\|plot_id)` |
| GLMM: growth by site | `glmer` (Gamma) | `area` [log link] | `time_sc × zone` | `(1\|colony_id)` |
| Chi-square: species composition | `chisq.test` | colony counts | by site × species | — |
| Wilcoxon: initial size by site | `wilcox.test` | `log(area)` | by site | — |

`time_sc` = `time_days / 100` — rescaled for numerical stability in Gamma GLMMs. Coefficients for `time_sc` describe change per 100 days; multiply by 3.65 for annual rate on log scale.
`size_z` = within-species z-score of log(initial area). Reference level: **PAST** (species), `size_z = 0` (average initial size).
Colonies with observation gaps (present → absent → present) were excluded prior to modelling.

**Trajectory visualisation models** (02a, 02d, 03c, 03d) use separate `glm(Gamma)` fitted on per-plot arithmetic means — these match the plotted points and are not used for inference.

---

## Results

### Colony growth — GLMM fixed effects (Gamma, log link)

| Term | Estimate | SE | z | p | sig | note |
|------|----------|----|---|---|-----|------|
| (Intercept) | 3.589 | 0.036 | 99.34 | <0.001 | *** | PAST baseline |
| time_sc | 0.0683 | 0.0062 | 10.95 | <0.001 | *** | ~+28%/100 days (~+25%/yr) |
| size_z | 1.085 | 0.024 | 45.36 | <0.001 | *** | ~3× area per +1 SD |
| time_sc:size_z | −0.0461 | 0.0055 | −8.41 | <0.001 | *** | growth convergence |
| siteS2P2 | +0.094 | 0.039 | 2.42 | 0.016 | * | S2P2 colonies ~10% larger |
| classColpophyllia natans | −0.754 | 0.176 | −4.28 | <0.001 | *** | ⚠ very low power |
| classOrbicella faveolata | +1.099 | 0.170 | 6.48 | <0.001 | *** | ⚠ very low power |
| classPorites porites | −1.023 | 0.114 | −8.97 | <0.001 | *** | ⚠ low power |
| classPseudodiploria strigosa | −0.867 | 0.115 | −7.53 | <0.001 | *** | ⚠ low power |
| classDiploria labyrinthiformis | +0.647 | 0.170 | 3.81 | <0.001 | *** | ⚠ very low power |
| classSiderastrea siderea | −0.290 | 0.125 | −2.32 | 0.020 | * | ⚠ low power |
| time_sc:classDiploria labyrinthiformis | −0.0952 | 0.0388 | −2.45 | 0.014 | * | ⚠ very low power |
| time_sc:classSiderastrea siderea | +0.0701 | 0.0292 | 2.40 | 0.016 | * | ⚠ low power |

> **Read:** PAST colonies grow ~25%/yr (time_sc β = 0.068 → ×3.65 → 0.249 log-units/yr → exp(0.249)−1 ≈ +28%). Colony size is the dominant predictor: each +1 SD in initial size ≈ ×3 colony area. Larger initial colonies grow proportionally slower (ontogenetic convergence). DLAB shows faster baseline growth but decelerates more strongly over time. SSID uniquely accelerates. S2P2 colonies start ~10% larger on average.

---

### Growth by size (LMM) — used for size-effect visualisation

| Term | Estimate | SE | t | p | sig |
|------|----------|----|---|---|-----|
| (Intercept) | 3.537 | 0.058 | 60.89 | <0.001 | *** |
| time_days | 0.0006 | 0.0001 | 10.49 | <0.001 | *** |
| size_z | 1.096 | 0.035 | 31.58 | <0.001 | *** |
| time_days:size_z | −0.0005 | 0.0001 | −9.09 | <0.001 | *** |

> **Read:** Consistent with the Gamma GLMM. Used to draw the size-effect prediction line in 02b only.

---

### Mortality severity — LMM fixed effects

| Term | Estimate | SE | t | p | sig |
|------|----------|----|---|---|-----|
| (Intercept) | 3.061 | 0.202 | 15.13 | <0.001 | *** |
| size_z | −0.223 | 0.060 | −3.69 | <0.001 | *** |
| classPorites porites | +0.560 | 0.282 | 1.99 | 0.047 | * | ⚠ low power |

> **Read:** Larger initial colonies lose significantly less proportional tissue when mortality occurs (~20% less per +1 SD). PPOR shows significantly higher severity than PAST. No other species differ significantly.

---

### Colony fate — GLM fixed effects (binomial)

| Term | Estimate | SE | z | p | sig | note |
|------|----------|----|---|---|-----|------|
| (Intercept) | −1.108 | 0.235 | −4.71 | <0.001 | *** | PAST, S1P1 baseline |
| size_z | −1.007 | 0.201 | −5.02 | <0.001 | *** | OR = 0.37/SD |
| siteS2P2 | −1.517 | 0.360 | −4.21 | <0.001 | *** | OR = 0.22 |
| classPorites porites | +2.008 | 0.621 | +3.24 | 0.001 | ** | ⚠ low power |
| classColpophyllia natans | +2.128 | 1.038 | +2.05 | 0.040 | * | ⚠ very low power |

> **Read:** Larger initial size strongly predicts survival — each +1 SD reduces odds of death by 63% (OR = 0.37). Colonies at S2P2 Farallon have 78% lower mortality odds than S1P1 Juanillo (OR = 0.22). PPOR and CNAT have significantly higher complete mortality probability than PAST.

---

### Site cover — GLMM (Gamma, log link)

| Term | Estimate | SE | z | p | sig |
|------|----------|----|---|---|-----|
| (Intercept) | 3.657 | 0.108 | 33.92 | <0.001 | *** |
| time_sc | 0.0644 | 0.0063 | 10.21 | <0.001 | *** |
| zoneJuanillo | −0.449 | 0.152 | −2.96 | 0.003 | ** |
| time_sc:zoneJuanillo | +0.037 | 0.014 | 2.64 | 0.008 | ** |

> **Read:** Significant positive cover trend at both sites (p < 0.001). Juanillo starts with lower cover (~−45% vs Farallon baseline) but shows a significantly steeper rate of increase (interaction p = 0.008), converging toward Farallon over time.

---

### Growth by site — GLMM (Gamma, log link)

| Term | Estimate | SE | z | p | sig |
|------|----------|----|---|---|-----|
| (Intercept) | 3.656 | 0.107 | 34.12 | <0.001 | *** |
| time_sc | 0.0371 | 0.014 | 2.64 | 0.008 | ** |
| zoneJuanillo | −0.449 | 0.152 | −2.96 | 0.003 | ** |

> **Read:** S2P2 Farallon colonies are significantly larger at baseline. Growth rate is positive at Farallon; the zone × time interaction is not significant, indicating similar growth trajectories at both sites after accounting for baseline differences.

---

### Species composition by site — Chi-square

| Statistic | Value |
|-----------|-------|
| χ² | 33.44 |
| df | 7 |
| p | < 0.001 *** |

> **Read:** Species composition differs significantly between sites. Interpret with caution due to small cell counts for rare species.

---

### Initial colony size by site — Wilcoxon rank-sum

| Statistic | Value |
|-----------|-------|
| W | 13051.5 |
| p | 0.035 * |

> **Read:** Colonies at the two sites have significantly different initial size distributions. S2P2 Farallon hosts larger colonies on average at the start of monitoring.

---

## Summary narrative

- **S2P2 Farallon outperforms S1P1 Juanillo in colony survival:** 78% lower odds of complete mortality (OR = 0.22), the strongest effect in the fate model.
- **Total cover is increasing at both sites:** Farallon starts higher; Juanillo shows a steeper rate of increase (interaction p = 0.008), suggesting recovery.
- **Colony size is the dominant cross-model predictor:** larger initial size → larger area (×3 per SD), slower relative growth, lower partial mortality severity (~20% less tissue lost per event), and lower complete mortality odds (OR = 0.37 per SD).
- **SSID uniquely accelerates over time:** positive `time_sc:SSID` interaction (z = 2.40, p = 0.016).
- **DLAB grows fastest at baseline but decelerates:** negative `time_sc:DLAB` interaction (z = −2.45, p = 0.014).
- **OFAV and DLAB start significantly larger than PAST; PPOR and PSTR start smaller.**
- **PPOR and CNAT carry the highest mortality risk** (both p ≤ 0.040 vs PAST in fate model).
- **Species composition and initial size differ between sites** (χ² p < 0.001; Wilcoxon p = 0.035).
- **Power caveat:** species with n < 26 flagged ⚠. Effect sizes are ecologically substantial and consistent across models; wider replication needed for confirmatory inference.

---

## Variable definitions

| Variable | Definition |
|----------|------------|
| `area` | Colony planar area (cm²) |
| `size_z` | Z-score of log(initial area) within species — 0 = species mean size, +1 = 1 SD above mean |
| `time_sc` | `time_days / 100` — used in Gamma GLMMs for numerical stability |
| `partial_mortality` | % area lost relative to previous survey (only when area decreased); complete death coded as 100% |
| `colony_fate` | `"died"` if colony absent before plot's last survey date |
| `time_days` | Days since first survey |
| `colony_id` | Unique ID: `site_plotnumber_genet_id` |
| `zone` | Site identifier (S1P1 = Juanillo, S2P2 = Farallon) |
| `PLOT_AREA_CM2` | 250,000 cm² (5 × 5 m plot) — used to convert area to % cover |
