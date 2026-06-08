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
| LMM: total cover (single-site) | `lmer` | `log(% cover)` | `time_days` | `(1\|plot_id)` |
| LMM: colony growth | `lmer` | `log(area)` | `time_days × class + time_days × size_z` | `(1\|plot_id) + (1\|colony_id)` |
| LMM: growth (species × size) | `lmer` | `log(area)` | `time_days × class × size_z` | same |
| GLMM: mortality occurrence | `glmer` (binomial) | P(any mortality) | `size_z + class` | same |
| LMM: mortality severity | `lmer` | `log(% area lost)` | `size_z + class` | same |
| GLMM: colony fate | `glmer` (binomial) | P(colony died) | `size_z + class + site` | `(1\|plot_id)` |
| LMM: site cover comparison | `lmer` | `log(% cover)` | `time_days × zone` | `(1\|plot_id)` |
| Chi-square: species composition | `chisq.test` | colony counts | by site × species | — |
| Wilcoxon: initial size by site | `wilcox.test` | `log(area)` | by site | — |
| t-test: net turnover vs zero | one-sample `t.test` | net change %/yr | per species | — |

`size_z` = within-species z-score of log(initial area). Reference level: **PAST** (species), `size_z = 0` (average initial size).
Colonies with observation gaps (present → absent → present) were excluded prior to modelling.

---

## Results

### Total cover trend — LMM (single-site, pooled)

| Term | Estimate | SE | t | p | sig |
|------|----------|----|---|---|-----|
| time_days | 0.0002 | — | 1.620 | 0.120 | ns |

> **Read:** No significant trend in total coral cover when pooling across both sites. The site-level comparison (below) reveals a significant positive trend when site is modelled explicitly.

---

### Colony growth — LMM fixed effects

| Term | Estimate | SE | t | p | sig | note |
|------|----------|----|---|---|-----|------|
| (Intercept) | 3.592 | 0.061 | — | <0.001 | *** | PAST baseline |
| time_days | 0.0006 | 0.0001 | 8.954 | <0.001 | *** | ~+21%/yr |
| classColpophyllia natans | −0.825 | 0.195 | — | <0.001 | *** | ⚠ low power |
| classOrbicella faveolata | +1.089 | 0.186 | — | <0.001 | *** | ⚠ low power |
| classPorites porites | −1.100 | 0.125 | — | <0.001 | *** | ⚠ low power |
| classPseudodiploria strigosa | −0.872 | 0.125 | — | <0.001 | *** | ⚠ low power |
| classDiploria labyrinthiformis | +0.703 | 0.183 | — | <0.001 | *** | ⚠ low power |
| classSiderastrea siderea | −0.346 | 0.136 | — | 0.011 | * | ⚠ low power |
| size_z | 1.091 | 0.027 | 40.803 | <0.001 | *** | ~3× area per +1 SD |
| time_days:size_z | −0.0005 | 0.0001 | −8.758 | <0.001 | *** | growth convergence |
| time_days:classSiderastrea siderea | +0.0009 | 0.0003 | 2.963 | 0.003 | ** | ⚠ low power |

> **Read:** PAST colonies grow ~21%/yr. Colony size is the dominant predictor: each +1 SD in initial size ≈ ×3 colony area (`size_z` β = 1.09). Larger initial colonies grow proportionally slower (`time_days:size_z` negative — ontogenetic growth slowdown, trajectories converge). SSID uniquely shows faster relative growth over time (positive `time_days:SSID` interaction). OFAV and DLAB start significantly larger than PAST; PPOR and PSTR start smaller.

---

### Colony growth — crossed model (species × size_z)

Extends the base model with `time_days × class × size_z` to test whether the size-growth relationship varies by species.

| Term | Estimate | SE | t | p | sig | note |
|------|----------|----|---|---|-----|------|
| time_days | 0.0006 | 0.0001 | — | <0.001 | *** | |
| size_z | ~1.07 | — | — | <0.001 | *** | |
| time_days:size_z | ~−0.0003–0.0005 | — | — | <0.001 | *** | |

> **Read:** Species × size interactions confirm that the size-growth gradient is not uniform. SSID shows pronounced size stratification; CNAT and SINT show weaker size scaling than PAST. Full coefficients mirror the base model main effects.

---

### Mortality occurrence — GLMM fixed effects

| Term | Estimate | SE | z | p | sig |
|------|----------|----|---|---|-----|
| (Intercept) | −0.717 | 0.113 | — | <0.001 | *** |

> **Read:** No species or size effect reached significance for mortality *occurrence*. All species experience similar rates of any partial mortality event at mean initial size.

---

### Mortality severity — LMM fixed effects

| Term | Estimate | SE | t | p | sig |
|------|----------|----|---|---|-----|
| (Intercept) | 3.061 | 0.202 | — | <0.001 | *** |
| size_z | −0.223 | 0.060 | — | <0.001 | *** |
| classPorites porites | +0.560 | 0.282 | — | 0.046 | * | ⚠ low power |

> **Read:** Larger initial colonies lose significantly less proportional tissue when mortality occurs (β = −0.22 per +1 SD; back-transformed ~20% less area lost). PPOR shows significantly higher severity than PAST.

---

### Colony fate — GLMM fixed effects

| Term | Estimate | SE | z | p | sig | note |
|------|----------|----|---|---|-----|------|
| (Intercept) | −1.108 | 0.235 | — | <0.001 | *** | PAST, S1P1 baseline |
| size_z | −1.007 | 0.201 | −5.017 | <0.001 | *** | OR = 0.37/SD |
| siteS2P2 | −1.517 | 0.360 | −4.209 | <0.001 | *** | OR = 0.22 |
| classPorites porites | +2.008 | 0.621 | +3.235 | 0.001 | ** | ⚠ low power |
| classColpophyllia natans | +2.128 | 1.038 | +2.049 | 0.040 | * | ⚠ low power |

> **Read:** Larger initial size strongly predicts survival — each +1 SD reduces odds of death by 63% (OR = 0.37). Colonies at S2P2 Farallon have 78% lower mortality odds than S1P1 Juanillo (OR = 0.22), indicating a strong site effect on colony survival. PPOR and CNAT have significantly higher complete mortality probability than PAST.

---

### Post-hoc pairwise — species fate (Tukey)

| Contrast | z | p | sig |
|----------|---|---|-----|
| Porites astreoides / Porites porites | −2.049 | 0.448 | ns |

> No significant pairwise fate contrasts after Tukey adjustment. Site effect dominates.

---

### Post-hoc pairwise — species mean size (Tukey, selected significant)

| Contrast | Estimate | SE | t | p | sig |
|----------|----------|----|---|---|-----|
| PAST − OFAV | −1.035 | 0.135 | — | <0.001 | *** |
| PAST − PPOR | +1.211 | 0.103 | — | <0.001 | *** |
| PAST − PSTR | +0.907 | 0.090 | — | <0.001 | *** |
| CNAT − DLAB | −1.266 | 0.201 | — | <0.001 | *** |
| CNAT − OFAV | −1.782 | 0.203 | — | <0.001 | *** |
| DLAB − PPOR | +1.730 | 0.161 | — | <0.001 | *** |
| DLAB − PSTR | +1.426 | 0.154 | — | <0.001 | *** |
| OFAV − PPOR | +2.246 | 0.168 | — | <0.001 | *** |
| OFAV − PSTR | +1.941 | 0.160 | — | <0.001 | *** |
| OFAV − SSID | +1.109 | 0.166 | — | <0.001 | *** |
| OFAV − SINT | +1.054 | 0.158 | — | <0.001 | *** |
| PPOR − SSID | −1.137 | 0.142 | — | <0.001 | *** |
| PPOR − SINT | −1.192 | 0.133 | — | <0.001 | *** |
| PSTR − SSID | −0.833 | 0.132 | — | <0.001 | *** |
| PSTR − SINT | −0.888 | 0.122 | — | <0.001 | *** |

> **Read:** OFAV is the largest species (significantly larger than all others). PPOR and PSTR are the smallest. CNAT, SSID, SINT, DLAB, and PAST occupy an intermediate cluster with multiple non-significant contrasts between them.

---

### Net turnover — one-sample t-tests (net change ≠ 0)

| Species | p | sig |
|---------|---|-----|
| PAST | 0.0001 | *** |
| SINT | 0.0001 | *** |
| PSTR | 0.027 | * |
| CNAT | 0.039 | * |
| PPOR | ns | |
| SSID | ns | |
| OFAV | ns | |
| DLAB | ns | |

> **Read:** PAST, SINT, PSTR, and CNAT show net tissue change significantly different from zero over the monitoring period. PPOR, SSID, OFAV, and DLAB do not — their gained and lost tissue are roughly in balance.

---

### Site cover comparison — LMM

| Term | Estimate | SE | t | p | sig |
|------|----------|----|---|---|-----|
| time_days | 0.0004 | — | 2.450 | 0.024 | * |

> **Read:** When the two sites are modelled separately (zone as fixed effect), total coral cover shows a significant positive trend over time (β = +0.0004 log-units/day, ~+15%/yr). The non-significant single-site pooled result above reflects greater residual variance when site structure is ignored.

---

### Species composition by site — Chi-square

| Statistic | Value |
|-----------|-------|
| χ² | 33.44 |
| df | 7 |
| p | < 0.001 *** |

> **Read:** Species composition differs significantly between sites (p < 0.001). Interpret with caution: chi-square approximation warning due to small cell counts for rare species (OFAV, DLAB, CNAT).

---

### Initial colony size by site — Wilcoxon rank-sum

| Statistic | Value |
|-----------|-------|
| W | 13051.5 |
| p | 0.035 * |

> **Read:** Colonies at the two sites have significantly different initial size distributions (W = 13051.5, p = 0.035). One site hosts larger colonies on average at the start of monitoring.

---

## Summary narrative

- **Site S2P2 Farallon outperforms S1P1 Juanillo in colony survival:** colonies at Farallon have 78% lower odds of complete mortality (OR = 0.22), the strongest effect in the fate model.
- **Total cover is increasing when site structure is accounted for:** the site-level LMM shows a significant positive trend (+15%/yr, p = 0.024), though the pooled single-site model is non-significant (p = 0.120).
- **Colony size is the dominant cross-model predictor:** larger initial size → larger area (×3 per SD), slower relative growth (ontogenetic convergence), lower partial mortality severity (~20% less tissue lost per event), and lower complete mortality odds (OR = 0.37 per SD).
- **SSID uniquely accelerates over time:** positive `time_days:SSID` interaction (β = +0.0009, p = 0.003) — SSID colonies gain a relative growth advantage as the monitoring period extends.
- **OFAV is the largest species** by a wide margin; PPOR and PSTR are the smallest.
- **PPOR and CNAT carry the highest mortality risk** (both p ≤ 0.040 above PAST in fate model).
- **Species composition and initial size differ between sites** (χ² p < 0.001; Wilcoxon p = 0.035), confirming that the two zones host structurally distinct communities.
- **Power caveat:** species with n < 26 flagged ⚠. Effect sizes are ecologically substantial and consistent across models; wider replication needed for confirmatory inference.

---

## Variable definitions

| Variable | Definition |
|----------|------------|
| `area` | Colony planar area (cm²) |
| `size_z` | Z-score of log(initial area) within species — 0 = species mean size, +1 = 1 SD above mean |
| `size_class` | Descriptive label only (small/medium/large, log-range thirds per species) — not used in models |
| `partial_mortality` | % area lost relative to previous survey (only when area decreased); complete death coded as 100% |
| `colony_fate` | `"died"` if colony absent before plot's last survey date |
| `time_days` | Days since first survey |
| `colony_id` | Unique ID: `site_plotnumber_genet_id` |
| `zone` | Site identifier (S1P1 = Juanillo, S2P2 = Farallon) |
| `PLOT_AREA_CM2` | 250,000 cm² (5 × 5 m plot) — used to convert area to % cover |
