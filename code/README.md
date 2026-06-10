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
| LMM: mortality severity by site | `lmer` | `log(% area lost)` | `time_sc × zone` | `(1\|plot_id) + (1\|colony_id)` |
| Chi-square: species composition | `chisq.test` | colony counts | by site × species | — |
| Wilcoxon: initial size by site | `wilcox.test` | `log(area)` | by site | — |

`time_sc` = `time_days / 100` — rescaled for numerical stability in Gamma GLMMs. Coefficients for `time_sc` describe change per 100 days; multiply by 3.65 for annual rate on log scale.
`size_z` = within-species z-score of log(initial area). Reference level: **PAST** (species), `size_z = 0` (average initial size).
Colonies with observation gaps (present → absent → present) were excluded prior to modelling.
The two **growth models** (colony growth GLMM and growth-by-size LMM) are fitted to follow-up observations only — each colony's baseline survey is dropped from the response so that the baseline-derived `size_z` is not also part of the outcome, avoiding a regression-to-the-mean inflation of the size terms (n = 947 obs, 264 colonies after exclusion).
Colony abundance (per species and per site) is modelled with a **Poisson GLMM** (log link) on a zero-filled plot × date × species grid; these are descriptive trend models and are not included in the inference table below.

**Trajectory visualisation models** (02a, 02d, 03c, 03d) use separate `glm(Gamma)` fitted on per-plot arithmetic means — these match the plotted points and are not used for inference. The partial-mortality severity models use partial events only (complete deaths excluded). The trend line in plot 07 (right panel) is one such visualisation GLM; its inferential counterpart is the **mortality severity by site** LMM below.

---

## Results

### Colony growth — GLMM fixed effects (Gamma, log link; follow-up obs only)

| Term | Estimate | SE | z | p | sig | note |
|------|----------|----|---|---|-----|------|
| (Intercept) | 3.640 | 0.054 | 67.61 | <0.001 | *** | PAST baseline |
| time_sc | 0.0597 | 0.0093 | 6.39 | <0.001 | *** | ~+24%/yr |
| size_z | 1.047 | 0.037 | 28.10 | <0.001 | *** | ~×2.8 area per +1 SD |
| time_sc:size_z | −0.0397 | 0.0081 | −4.89 | <0.001 | *** | growth convergence |
| siteS2P2 | +0.085 | 0.048 | 1.78 | 0.075 | . | n.s. |
| classOrbicella faveolata | +1.121 | 0.263 | 4.27 | <0.001 | *** | ⚠ very low power |
| classColpophyllia natans | −0.967 | 0.255 | −3.79 | <0.001 | *** | ⚠ very low power |
| classPorites porites | −1.415 | 0.251 | −5.65 | <0.001 | *** | ⚠ low power |
| classPseudodiploria strigosa | −1.029 | 0.202 | −5.08 | <0.001 | *** | ⚠ low power |
| classDiploria labyrinthiformis | −0.152 | 0.274 | −0.55 | 0.579 | | ⚠ n.s. |
| classSiderastrea siderea | +0.147 | 0.197 | 0.75 | 0.456 | | ⚠ n.s. |
| classStephanocoenia intersepta | +0.036 | 0.137 | 0.27 | 0.790 | | n.s. |
| time_sc:class (all species) | — | — | — | ≥0.32 | | none significant |

> **Read:** Colony size is the dominant predictor: each +1 SD in initial size ≈ ×2.8 colony area, and larger initial colonies grow proportionally slower (ontogenetic convergence; time_sc:size_z < 0). PAST grows ~+24%/yr (time_sc β = 0.0597 → ×3.65 → 0.218 log-units/yr → exp(0.218)−1 ≈ +24%). After excluding baseline observations, **no species differs significantly from PAST in growth *rate*** (all time × species interactions n.s.) — the earlier SSID-acceleration and DLAB-deceleration effects do not survive. OFAV colonies are larger and CNAT/PPOR/PSTR smaller than PAST at comparable timepoints; the S2P2 baseline-size difference is no longer significant (p = 0.075).

---

### Growth by size (LMM) — follow-up obs only; used for size-effect visualisation

| Term | Estimate | SE | t | p | sig |
|------|----------|----|---|---|-----|
| (Intercept) | 3.575 | 0.068 | 52.52 | <0.001 | *** |
| time_days | 0.00057 | 0.00009 | 6.38 | <0.001 | *** |
| size_z | 1.055 | 0.048 | 21.99 | <0.001 | *** |
| time_days:size_z | −0.00048 | 0.00009 | −5.52 | <0.001 | *** |

> **Read:** Consistent with the Gamma GLMM — the negative time × size term confirms larger colonies grow proportionally slower. Used to draw the size-effect prediction line in 02b only.

---

### Mortality severity — LMM fixed effects (partial events only)

| Term | Estimate | SE | t | p | sig |
|------|----------|----|---|---|-----|
| (Intercept) | 2.756 | 0.150 | 18.38 | <0.001 | *** |
| size_z | −0.059 | 0.064 | −0.91 | 0.362 | | |
| classOrbicella faveolata | −0.753 | 0.406 | −1.86 | 0.065 | . | ⚠ very low power |
| (all other species) | n.s. | | | ≥0.40 | | |

> **Read:** With complete deaths excluded, **no predictor is significant** — initial size does *not* predict how much tissue is lost in a partial-mortality event (size_z p = 0.36). The earlier "larger colonies lose less tissue" result was an artifact of including complete deaths (coded 100%), which cluster among small colonies. OFAV is marginally lower than PAST (p = 0.065).

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

### Site cover — GLMM (Gamma, log link; zone ref = Farallon)

| Term | Estimate | SE | z | p | sig |
|------|----------|----|---|---|-----|
| (Intercept) | 0.363 | 0.419 | 0.87 | 0.386 | |
| time_sc | 0.0371 | 0.0140 | 2.64 | 0.008 | ** |
| zoneJuanillo | −0.499 | 0.574 | −0.87 | 0.385 | |
| time_sc:zoneJuanillo | −0.041 | 0.021 | −1.95 | 0.052 | . |

> **Read:** Total cover increases significantly over time at Farallon, the reference zone (time_sc p = 0.008). Baseline cover does not differ significantly between sites (p = 0.385), and the difference in *trend* is only marginal (interaction p = 0.052), with Juanillo on a flatter trajectory. (This revises an earlier summary that reported Juanillo as increasing more steeply — that pattern does not hold in the current model.)

---

### Growth by site — GLMM (Gamma, log link; zone ref = Farallon)

| Term | Estimate | SE | z | p | sig |
|------|----------|----|---|---|-----|
| (Intercept) | 3.657 | 0.108 | 33.92 | <0.001 | *** |
| time_sc | 0.0644 | 0.0063 | 10.21 | <0.001 | *** |
| zoneJuanillo | −0.449 | 0.152 | −2.96 | 0.003 | ** |
| time_sc:zoneJuanillo | −0.019 | 0.011 | −1.75 | 0.080 | . |

> **Read:** Farallon colonies are significantly larger (Juanillo ~36% smaller, exp(−0.449) ≈ 0.64). Growth over time is positive; the zone × time interaction is not significant (p = 0.080), indicating similar growth trajectories at both sites after accounting for baseline size differences.

---

### Mortality severity by site — LMM (partial events only; zone ref = Farallon)

| Term | Estimate | SE | t | p | sig |
|------|----------|----|---|---|-----|
| (Intercept) | 2.321 | 0.250 | 9.30 | <0.001 | *** |
| zoneJuanillo | 1.041 | 0.426 | 2.45 | 0.022 | * |
| time_sc | 0.052 | 0.051 | 1.01 | 0.313 | |
| time_sc:zoneJuanillo | −0.149 | 0.092 | −1.62 | 0.106 | |

> **Read:** Inferential test of the site contrast shown in plot 07 (right panel). When a partial-mortality event occurs, Juanillo colonies lose significantly more tissue than Farallon colonies (exp(1.041) ≈ 2.8× on the % scale, p = 0.022). Severity does not change significantly over time (time_sc p = 0.31), and the trend does not differ between sites (interaction p = 0.11). Complete deaths are excluded, so this measures the magnitude of partial tissue loss, not the rate of colony death (see colony fate model).

---

### Species composition by site — Chi-square (Monte-Carlo)

| Statistic | Value |
|-----------|-------|
| χ² | 33.44 |
| p | < 0.001 *** (Monte-Carlo, B = 10⁵) |

> **Read:** Species composition differs significantly between sites. Because several rare species have expected cell counts below 5, a Monte-Carlo simulated _p_-value was used instead of the asymptotic chi-square approximation.

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
- **Total cover is increasing over time at Farallon** (time_sc p = 0.008); baseline cover does not differ significantly between sites and the site difference in trend is only marginal (interaction p = 0.052).
- **Colony size is the dominant cross-model predictor:** larger initial size → larger area (~×2.8 per SD), slower relative growth (ontogenetic convergence), and lower complete mortality odds (OR = 0.37 per SD).
- **Initial size does *not* predict partial-mortality severity** once complete deaths are excluded (size_z p = 0.36) — the earlier effect was driven by 100%-coded complete deaths.
- **Partial-mortality severity is significantly higher at Juanillo than Farallon** (≈2.8× more tissue lost per partial event, p = 0.022); the severity trend over time does not differ between sites (interaction p = 0.11) — the inferential test behind plot 07.
- **No species differs significantly from PAST in growth *rate*** (all time × species interactions n.s. on follow-up observations); the earlier SSID-acceleration and DLAB-deceleration findings do not survive removal of the baseline observation.
- **OFAV starts significantly larger than PAST; CNAT, PPOR and PSTR start smaller** (baseline-size contrasts in the growth model).
- **PPOR and CNAT carry the highest mortality risk** (both p ≤ 0.040 vs PAST in fate model).
- **Species composition and initial size differ between sites** (χ² Monte-Carlo p < 0.001; Wilcoxon p = 0.035).
- **Power caveat:** species with n < 26 flagged ⚠. Effects that survive (size across growth and fate; the site survival contrast) are ecologically substantial; wider replication is needed for confirmatory inference, and the site contrasts rest on only two sites (three plots each).

---

## Variable definitions

| Variable | Definition |
|----------|------------|
| `area` | Colony planar area (cm²) |
| `size_z` | Z-score of log(initial area) within species — 0 = species mean size, +1 = 1 SD above mean |
| `time_sc` | `time_days / 100` — used in Gamma GLMMs for numerical stability |
| `partial_mortality` | % area lost relative to previous survey (only when area decreased); complete death coded as 100% but excluded from the severity model and trajectory plots |
| `colony_fate` | `"died"` if colony absent before plot's last survey date |
| `time_days` | Days since first survey |
| `colony_id` | Unique ID: `site_plotnumber_genet_id` |
| `zone` | Site identifier (S1P1 = Juanillo, S2P2 = Farallon) |
| `PLOT_AREA_CM2` | 250,000 cm² (5 × 5 m plot) — used to convert area to % cover |
