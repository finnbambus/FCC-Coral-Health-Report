# FCC Coral Health Report — Statistical Results Summary

**Site:** Farallon (S2P2) · **Data:** `combined_data.csv`
**Species in analysis:** CNAT *(Colpophyllia natans)* · PAST *(Porites astreoides)* · PSTR *(Pseudodiploria strigosa)* · SSID *(Siderastrea siderea)* · SINT *(Stephanocoenia intersepta)*
**Note:** Species with fewer than 5 colonies excluded (hard floor). Species marked ⚠ are below the statistical power threshold (d=0.8, 1−β=0.80, n<26 colonies) — interpret with caution.
**Reference level:** *Porites astreoides* (PAST) — largest colony count, lowest SE.

---

## Model structure

| Model | Type | Response | Fixed effects | Random effects |
|-------|------|----------|---------------|----------------|
| LMM: total cover trend | `lmer` | `log(% cover)` | `time_days` | `(1\|plot_id)` |
| LMM: colony growth | `lmer` | `log(area)` | `time_days × class + time_days × size_z` | `(1\|plot_id) + (1\|colony_id)` |
| LMM: growth (species × size) | `lmer` | `log(area)` | `time_days × class × size_z` | same |
| GLMM: mortality occurrence | `glmer` (binomial) | P(any mortality) | `size_z + class` | same |
| LMM: mortality severity | `lmer` | `log(% area lost)` | `size_z + class` | same |
| GLMM: colony fate | `glmer` (binomial) | P(colony died) | `size_z + class` | `(1\|plot_id)` |

`size_z` = within-species z-score of log(initial area). Reference level: **PAST** (species), `size_z = 0` (average initial size).
Colonies with observation gaps (present → absent → present) were excluded prior to modelling.

---

## Significant results (p < 0.05)

### Total cover trend — LMM

| Term | Estimate | SE | t | p | sig |
|------|----------|----|---|---|-----|
| time_days | 0.0004 | 0.0001 | 3.154 | 0.009 | ** |

> **Read:** Total coral cover is increasing significantly over the monitoring period (β = +0.0004 log-units/day). Back-transformed: ~15% annual increase in % cover.

---

### Colony growth — LMM fixed effects

| Term | Estimate | SE | t | p | sig | note |
|------|----------|----|---|---|-----|------|
| (Intercept) | 3.768 | 0.090 | 41.816 | <0.001 | *** | PAST baseline |
| time_days | 0.0006 | 0.0001 | 10.319 | <0.001 | *** | |
| classColpophyllia natans | −0.857 | 0.171 | −5.028 | <0.001 | *** | ⚠ low power |
| classPseudodiploria strigosa | −0.728 | 0.162 | −4.488 | <0.001 | *** | ⚠ low power |
| size_z | 1.100 | 0.031 | 35.076 | <0.001 | *** | |
| time_days:size_z | −0.0003 | 0.0001 | −5.561 | <0.001 | *** | |

> **Read:** PAST colonies are significantly larger than CNAT and PSTR at baseline. `size_z` is strongly positive — each +1 SD in initial size corresponds to ~exp(1.10) ≈ 3× larger area. The `time_days:size_z` interaction is negative, meaning larger initial colonies grow more slowly in relative terms (ontogenetic growth slowdown — trajectories converge over time).

---

### Colony growth — crossed model (species × size_z)

| Term | Estimate | SE | t | p | sig | note |
|------|----------|----|---|---|-----|------|
| (Intercept) | 3.792 | 0.055 | 69.169 | <0.001 | *** | PAST baseline |
| time_days | 0.0006 | 0.0001 | 10.398 | <0.001 | *** | |
| classColpophyllia natans | −0.869 | 0.124 | −7.013 | <0.001 | *** | ⚠ low power |
| classPseudodiploria strigosa | −0.706 | 0.120 | −5.869 | <0.001 | *** | ⚠ low power |
| size_z | 1.068 | 0.025 | 42.410 | <0.001 | *** | |
| time_days:size_z | −0.0003 | 0.0001 | −4.675 | <0.001 | *** | |
| classColpophyllia natans:size_z | −0.510 | 0.138 | −3.689 | <0.001 | *** | ⚠ low power |
| classSiderastrea siderea:size_z | 1.213 | 0.100 | 12.115 | <0.001 | *** | ⚠ low power |
| classStephanocoenia intersepta:size_z | −0.206 | 0.072 | −2.876 | 0.004 | ** | ⚠ low power |

> **Read:** The species × size_z interactions reveal that the size effect is not uniform. SSID shows the most extreme size stratification (PAST size_z slope +1.068, SSID adds +1.213, combined ≈ 2.28 per SD). CNAT has a weaker size gradient (−0.51 below PAST's slope).

---

### Mortality occurrence — GLMM fixed effects

| Term | Estimate | SE | z | p | sig |
|------|----------|----|---|---|-----|
| (Intercept) | −0.888 | 0.088 | −10.039 | <0.001 | *** |

> **Read:** The intercept represents baseline mortality probability for PAST at mean size (logistic scale). No species or size_z effect reached significance in this model, meaning all species have similar occurrence rates of any mortality event.

---

### Mortality severity — LMM fixed effects

| Term | Estimate | SE | t | p | sig |
|------|----------|----|---|---|-----|
| (Intercept) | 2.681 | 0.227 | 11.788 | 0.014 | * |
| size_z | −0.227 | 0.088 | −2.583 | 0.011 | * |

> **Read:** When mortality occurs, larger colonies lose significantly less proportional area. Each +1 SD in initial size reduces log(% area lost) by −0.23 — back-transformed, ~20% less tissue lost per event. No significant species effect on severity.

---

### Colony fate — GLMM fixed effects

| Term | Estimate | SE | z | p | sig | note |
|------|----------|----|---|---|-----|------|
| (Intercept) | −3.153 | 0.501 | −6.297 | <0.001 | *** | PAST baseline |
| size_z | −1.467 | 0.416 | −3.528 | <0.001 | *** | |
| classColpophyllia natans | +2.553 | 1.140 | +2.240 | 0.025 | * | ⚠ low power |

> **Read:** Larger initial size strongly predicts survival — each +1 SD reduces log-odds of colony death by −1.47 (OR ≈ 0.23, ~77% reduction in odds). CNAT has significantly higher complete mortality probability than PAST (+2.55 log-odds). PAST is the most resilient species in the dataset.

---

### Post-hoc pairwise — species mean log(area) (Tukey)

| Contrast | Estimate | SE | t | p | sig | note |
|----------|----------|----|---|---|-----|------|
| PAST − CNAT | +0.814 | 0.155 | 5.244 | <0.001 | *** | ⚠ low power |
| PAST − PSTR | +0.689 | 0.143 | 4.816 | <0.001 | *** | ⚠ low power |
| CNAT − SINT | −0.730 | 0.172 | −4.236 | <0.001 | *** | ⚠ low power |
| CNAT − SSID | −0.761 | 0.191 | −3.987 | 0.001 | ** | ⚠ low power |
| PSTR − SINT | −0.605 | 0.163 | −3.717 | 0.003 | ** | ⚠ low power |
| PSTR − SSID | −0.636 | 0.180 | −3.533 | 0.005 | ** | ⚠ low power |

> **Read:** PAST and SINT and SSID are all significantly larger than CNAT and PSTR. SINT and SSID are the largest species; PAST and PSTR cluster together; CNAT is smallest. All contrasts involve at least one low-power species; effect sizes (|est| ≈ 0.6–0.8 log-units) are ecologically substantial.

---

## Summary narrative

- **Total coral cover is increasing** at Farallon S2P2 over the monitoring period (+15%/yr estimated from model), a positive signal for reef health.
- **Colony size is the strongest predictor** across all outcomes: larger initial size → larger area, slower relative growth, lower partial mortality severity, lower complete mortality risk (t = 35 in base model).
- **Growth trajectories converge:** larger initial colonies grow proportionally slower (`time_days:size_z` negative) — universal ontogenetic pattern in corals. SSID shows the most extreme size stratification.
- **PAST is the most resilient species:** lowest complete mortality probability, significantly lower than CNAT (OR ≈ 0.08 at mean size). This is ecologically consistent with PAST's known tolerance for marginal conditions.
- **CNAT shows elevated mortality risk:** +2.55 log-odds above PAST for complete colony death.
- **Power caveat:** most pairwise contrasts and several species fixed effects are flagged ⚠ low power. Effect sizes are ecologically substantial and consistent across models, suggesting the patterns are real, but wider replication is needed for confirmatory inference.

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
| `PLOT_AREA_CM2` | 250,000 cm² (5×5m plot) — used to convert area to % cover |
