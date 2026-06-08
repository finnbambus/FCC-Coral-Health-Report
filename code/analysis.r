# load packages
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(lme4, lmerTest, dplyr, tidyr, ggplot2, viridisLite, patchwork, emmeans, multcomp, multcompView, scales, pwr, here)

plots_dir <- here::here("plots")
dir.create(plots_dir, showWarnings = FALSE)

# read and combine both site data files
data <- dplyr::bind_rows(
  read.csv(here::here("data", "S1P1_combined_data.csv"), sep = ";"),
  read.csv(here::here("data", "S2P2_combined_data.csv"), sep = ";")
)

# unique IDs: genet_id and plot numbers are site-local, so prefix with site
data$plot_id   <- paste(data$site, gsub("Plot ", "", data$plot), sep = "_")
data$colony_id <- paste(data$plot_id, data$genet_id, sep = "_")

# convert date
data$date <- as.Date(data$date, format = "%Y-%m-%d")

# detect split and merge events:
# a colony that goes from 1 entry → 2+ entries per date has split
# a colony that goes from 2+ entries → 1 entry per date has merged
entry_counts <- data %>%
  group_by(colony_id, date) %>%
  summarise(n_entries = n(), .groups = "drop") %>%
  arrange(colony_id, date) %>%
  group_by(colony_id) %>%
  mutate(prev_n = lag(n_entries, default = first(n_entries))) %>%
  ungroup()

event_notes <- entry_counts %>%
  filter(n_entries != prev_n) %>%
  mutate(event_note = paste0(
    ifelse(n_entries > prev_n, "split", "merged"), " at ", date
  )) %>%
  group_by(colony_id) %>%
  summarise(event_note = paste(event_note, collapse = "; "), .groups = "drop")

# combine multiple entries per colony + date: sum areas, attach event note
data <- data %>%
  left_join(event_notes, by = "colony_id") %>%
  group_by(colony_id, date) %>%
  summarise(
    zone              = first(zone),
    site              = first(site),
    plot              = first(plot),
    plot_id           = first(plot_id),
    genet_id          = first(genet_id),
    class             = first(class),
    code              = first(code),
    area              = sum(area, na.rm = TRUE),
    partial_mortality = first(partial_mortality),
    notes             = ifelse(!is.na(first(event_note)), first(event_note), first(notes)),
    .groups           = "drop"
  )

# remove colonies with observation gaps (present → absent → present within their active span)
# a gap means the colony could not be identified at a survey; identity is unreliable
plot_survey_dates <- data %>%
  group_by(plot_id) %>%
  summarise(all_dates = list(sort(unique(date))), .groups = "drop")

colony_span <- data %>%
  group_by(colony_id, plot_id) %>%
  summarise(
    first_seen    = min(date),
    last_seen     = max(date),
    dates_present = list(sort(unique(date))),
    .groups       = "drop"
  ) %>%
  left_join(plot_survey_dates, by = "plot_id") %>%
  mutate(
    dates_expected = mapply(
      function(all, f, l) all[all >= f & all <= l],
      all_dates, first_seen, last_seen, SIMPLIFY = FALSE
    ),
    has_gap = mapply(
      function(present, expected) length(present) < length(expected),
      dates_present, dates_expected
    )
  )

n_gap <- sum(colony_span$has_gap)
if (n_gap > 0) {
  message(n_gap, " colonies removed: observation gap (present → absent → present)")
  data <- data %>% filter(!colony_id %in% colony_span$colony_id[colony_span$has_gap])
} else {
  message("No observation gaps detected")
}

# numeric time (recomputed after row-combining; dates unchanged so min is the same)
data$time_days <- as.numeric(data$date - min(data$date))

# partial mortality: percentage of area lost relative to previous time point
# only filled when area decreased; NA when area stayed the same or grew
data <- data %>%
  arrange(colony_id, date) %>%
  group_by(colony_id) %>%
  mutate(
    partial_mortality = ifelse(
      !is.na(lag(area)) & area < lag(area),
      round((lag(area) - area) / lag(area) * 100, 1),
      NA_real_
    )
  ) %>%
  ungroup()

# complete colony mortality: a colony that disappears before the last survey date
# in its plot is assumed dead (last seen alive = last observation date)
plot_max_date <- data %>%
  group_by(plot_id) %>%
  summarise(plot_max_date = max(date), .groups = "drop")

data <- data %>%
  group_by(colony_id) %>%
  mutate(colony_last_date = max(date)) %>%
  ungroup() %>%
  left_join(plot_max_date, by = "plot_id") %>%
  mutate(colony_fate = ifelse(colony_last_date < plot_max_date, "died", "alive"))

# initial size per colony: computed from each colony's first observation
# size_z   — continuous: z-score of log(area) within species (used in models)
#            captures relative size on the same scale across species
# size_class — categorical: log-range thirds within species (used for descriptive plots only)
#              single-colony / zero-variance species assigned z=0 and class="medium"
initial_size <- data %>%
  group_by(colony_id) %>%
  slice_min(date, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  group_by(class) %>%
  mutate(
    log_area   = log(area),
    size_z     = {
                   s <- sd(log_area)
                   if (n() < 2 || s == 0) rep(0, n())
                   else (log_area - mean(log_area)) / s
                 },
    size_class = if (n_distinct(area) == 1) {
      factor(rep("medium", n()), levels = c("small", "medium", "large"))
    } else {
      cut(log_area, breaks = seq(min(log_area), max(log_area), length.out = 4),
          labels = c("small", "medium", "large"), include.lowest = TRUE)
    }
  ) %>%
  ungroup() %>%
  dplyr::select(colony_id, size_z, size_class)

data <- left_join(data, initial_size, by = "colony_id")
data$size_class <- factor(data$size_class, levels = c("small", "medium", "large"))

# --- species sample size check (two-tier) ---
# Tier 1 — hard floor: species with fewer than n_hard_min colonies are removed;
#   below this threshold a species coefficient and random intercept cannot be estimated.
# Tier 2 — power flag: species above the floor but below the power-based n_power_min
#   remain in the analysis but are flagged as "low power" in results_sig.
#   This is the standard approach for field surveys with limited replication.
target_power <- 0.80
target_alpha <- 0.05
target_d     <- 0.8       # large Cohen's d (conservative LMM lower bound)
n_hard_min   <- 5         # absolute identifiability floor

n_power_min <- ceiling(pwr::pwr.t.test(
  d = target_d, sig.level = target_alpha, power = target_power, type = "two.sample"
)$n)

species_n <- data %>%
  distinct(colony_id, class) %>%
  count(class, name = "n_colonies") %>%
  mutate(
    n_power_min    = n_power_min,
    adequate_power = n_colonies >= n_power_min,
    low_power      = n_colonies >= 10 & n_colonies < n_power_min,
    very_low_power = n_colonies >= n_hard_min & n_colonies < 10,
    too_few        = n_colonies < n_hard_min
  )

message("Species n check  |  hard floor: ", n_hard_min,
        "  |  power threshold (d=", target_d, ", 1-β=", target_power, "): ", n_power_min, " colonies")
print(as.data.frame(species_n))

# apply only the hard floor; power-flagged species stay in but are marked in results
excluded_spp       <- species_n %>% filter(too_few)        %>% pull(class)
very_low_power_spp <- species_n %>% filter(very_low_power) %>% pull(class)
low_power_spp      <- species_n %>% filter(low_power)      %>% pull(class)

if (length(excluded_spp) > 0) {
  message("Removing (below hard floor n=", n_hard_min, "): ",
          paste(excluded_spp, collapse = ", "))
  data <- data %>% filter(!class %in% excluded_spp)
}
if (length(very_low_power_spp) > 0) {
  message("Very low power n=5-9 (flagged in results, not removed): ",
          paste(very_low_power_spp, collapse = ", "))
}
if (length(low_power_spp) > 0) {
  message("Low power n=10-", n_power_min - 1, " (flagged in results, not removed): ",
          paste(low_power_spp, collapse = ", "))
}

# reference species: PAST has most colonies → lowest SE on all pairwise comparisons
data$class <- relevel(droplevels(factor(data$class)), ref = "Porites astreoides")

# plot area constant for % cover conversion (5 × 5 m photoquadrats)
PLOT_AREA_CM2 <- 250000   # 25 m² = 250 000 cm²

# shared model flags
multi_site <- length(unique(data$site)) > 1
site_term  <- if (multi_site) "+ site" else ""

# significance star helper — used in plots and results summary
sig_stars <- function(p) {
  dplyr::case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    p < 0.1   ~ ".",
    TRUE      ~ ""
  )
}

# ── design system ──────────────────────────────────────────────────────────
accent_col    <- "#898C31"   # project accent  rgb(137, 140, 49)
accent_bg_col <- "#E5E6D2"   # accent lightened 78% (used for grid lines)
sp_font       <- "Times"     # Times New Roman for species names (italic)

theme_coral <- function(base_size = 9) {
  theme_minimal(base_size = base_size, base_family = "Helvetica") %+replace%
    theme(
      plot.title       = element_text(face = "bold", size = base_size + 1, hjust = 0,
                                      color = "grey10", margin = margin(b = 3)),
      plot.subtitle    = element_text(size = base_size - 1, hjust = 0, color = "grey45",
                                      margin = margin(b = 4)),
      axis.title       = element_text(size = base_size, color = "grey25"),
      axis.text        = element_text(size = base_size - 1, color = "grey35"),
      panel.grid.major = element_line(color = accent_bg_col, linewidth = 0.3),
      panel.grid.minor = element_blank(),
      legend.text      = element_text(size = base_size - 1),
      legend.title     = element_text(size = base_size),
      legend.key.size  = unit(0.75, "lines"),
      strip.text       = element_text(face = "italic", family = sp_font,
                                      size = base_size, color = "grey20"),
      plot.margin      = margin(6, 8, 6, 6)
    )
}

# canonical species palette — named by full species name so every scale_color/fill_manual
# call resolves to the same color regardless of which species appear in a given plot
sp_levels <- sort(unique(as.character(data$class)))
pal <- setNames(viridisLite::viridis(length(sp_levels), option = "D", begin = 0.05, end = 0.92),
                sp_levels)

# site palette: 2 viridis colors for site-level plots
n_sites  <- length(unique(data$site))
site_pal <- setNames(viridisLite::viridis(max(n_sites, 2), option = "D",
                                          begin = 0.15, end = 0.75)[seq_len(n_sites)],
                     sort(unique(data$zone)))

# species code lookup: full name → 4-letter code (SSID, DLAB, etc.)
sp_codes <- data %>% distinct(class, code) %>% { setNames(.$code, .$class) }

# species label with n colonies: plain two-line label for axis/strip
sp_n <- data %>% distinct(colony_id, class) %>% count(class, name = "n_col") %>%
  { setNames(.$n_col, .$class) }
sp_labels_n <- setNames(
  paste0(sp_codes[names(sp_n)], "
n = ", sp_n),
  names(sp_n)
)

# figure dimensions: two aspect ratios used throughout (300 DPI, px)
fig_sq_w  <- 945   # squarish  945 × 945 px  (half of document width)
fig_sq_h  <- 945
fig_ls_w  <- 1890  # landscape 1890 × 945 px (full document width, 2:1 ratio)
fig_ls_h  <- 945

# --- 0. Sample sizes ---
sample_sizes <- data %>%
  group_by(class) %>%
  summarise(
    n_colonies   = n_distinct(colony_id),
    n_obs        = n(),
    n_timepoints = n_distinct(date),
    .groups      = "drop"
  )
print(sample_sizes)


# ═══════════════════════════════════════════════════════════════════════════
# SECTION 0 — COMMUNITY TOTAL COVER
# primary reef health indicator: is total coral cover going up or down?
# ═══════════════════════════════════════════════════════════════════════════

cover <- data %>%
  group_by(site, zone, plot_id, date, time_days) %>%
  summarise(
    total_area = sum(area, na.rm = TRUE),
    n_colonies = n_distinct(colony_id),
    .groups    = "drop"
  ) %>%
  mutate(pct_cover = total_area / PLOT_AREA_CM2 * 100)

cover_sp <- data %>%
  group_by(plot_id, date, time_days, class) %>%
  summarise(sp_area = sum(area, na.rm = TRUE), .groups = "drop") %>%
  mutate(pct_cover_sp = sp_area / PLOT_AREA_CM2 * 100)

# LMM: is total coral cover changing significantly over time?
model_cover <- lmer(log(pct_cover) ~ time_days + (1 | plot_id), data = cover)
summary(model_cover)

cover_trend <- data.frame(
  time_days = seq(min(cover$time_days), max(cover$time_days), length.out = 80)
) %>% mutate(
  date      = min(cover$date) + time_days,
  pct_cover = exp(predict(model_cover, newdata = ., re.form = NA))
)

# 95% CI from fixed-effects vcov (population-level, no random effect uncertainty)
X_cover   <- model.matrix(~ time_days, data = cover_trend)
pred_se_c <- sqrt(rowSums((X_cover %*% as.matrix(vcov(model_cover))) * X_cover))
pred_log_c <- predict(model_cover, newdata = cover_trend, re.form = NA)
cover_trend$pct_cover_lo <- exp(pred_log_c - 1.96 * pred_se_c)
cover_trend$pct_cover_hi <- exp(pred_log_c + 1.96 * pred_se_c)

# Plot C1: stacked species cover — shows total cover AND composition shift
cover_sp_mean <- cover_sp %>%
  group_by(date, class) %>%
  summarise(mean_pct = mean(pct_cover_sp), .groups = "drop")

pC1 <- ggplot(cover_sp_mean, aes(x = date, y = mean_pct, fill = class)) +
  geom_area(alpha = 0.85, position = "stack") +
  scale_fill_manual(values = pal, labels = sp_codes) +
  scale_x_date(date_breaks = "6 months", date_labels = "%b '%y") +
  scale_y_continuous(labels = \(x) paste0(x, "%"),
                     expand = expansion(mult = c(0, 0.03))) +
  labs(x = "Survey date", y = "Mean cover (% of 25 m²)", fill = NULL) +
  theme_coral() +
  theme(legend.position = "top") +
  guides(fill = guide_legend(nrow = 1))

ggsave(file.path(plots_dir, "01_community_cover_composition.jpeg"), pC1,
       width = fig_ls_w, height = fig_ls_h, dpi = 300, units = "px")


# ═══════════════════════════════════════════════════════════════════════════
# SECTION 1 — SPECIES GROWTH TRAJECTORIES
# which species are growing, stable, or declining?
# size_z included as explanatory factor throughout
# ═══════════════════════════════════════════════════════════════════════════

model_base <- lmer(
  as.formula(paste("log(area) ~ time_days * class + time_days * size_z", site_term,
                   "+ (1 | plot_id) + (1 | colony_id)")),
  data = data
)
summary(model_base)

# --- 1a. LMM assumption checks ---

resid_df <- data.frame(
  fitted    = fitted(model_base),
  residuals = resid(model_base),
  sqrt_abs_resid = sqrt(abs(resid(model_base)))
)

# residuals vs fitted: checks linearity and homoscedasticity
# expect random scatter around zero with no trend or funnel shape
pa <- ggplot(resid_df, aes(x = fitted, y = residuals)) +
  geom_point(alpha = 0.20, size = 0.7, color = "grey40") +
  geom_hline(yintercept = 0, linetype = "dashed", color = accent_col, linewidth = 0.5) +
  geom_smooth(method = "loess", se = FALSE, color = accent_col, linewidth = 0.8) +
  labs(x = "Fitted values", y = "Residuals", title = "Residuals vs Fitted") +
  theme_coral()

pb <- ggplot(resid_df, aes(sample = residuals)) +
  stat_qq(alpha = 0.25, size = 0.7, color = "grey40") +
  stat_qq_line(color = accent_col, linewidth = 0.7) +
  labs(x = "Theoretical quantiles", y = "Sample quantiles",
       title = "Normal Q-Q — residuals") +
  theme_coral()

pc <- ggplot(resid_df, aes(x = fitted, y = sqrt_abs_resid)) +
  geom_point(alpha = 0.20, size = 0.7, color = "grey40") +
  geom_smooth(method = "loess", se = FALSE, color = accent_col, linewidth = 0.8) +
  labs(x = "Fitted values", y = "√|Residuals|", title = "Scale-Location") +
  theme_coral()

re_vals <- ranef(model_base)$colony_id[, 1]
pd <- ggplot(data.frame(re = re_vals), aes(sample = re)) +
  stat_qq(alpha = 0.3, size = 0.7, color = "grey40") +
  stat_qq_line(color = accent_col, linewidth = 0.7) +
  labs(x = "Theoretical quantiles", y = "Sample quantiles",
       title = "Q-Q: random effects (colony)") +
  theme_coral()

p_assumptions <- (pa | pb) / (pc | pd)
ggsave(file.path(plots_dir, "00_lmm_assumption_checks.jpeg"), p_assumptions,
       width = fig_ls_w, height = fig_ls_w, dpi = 300, units = "px")


# model 1b: species x size_z interaction — only run if size_z is significant in model_base
# tests whether the continuous size effect on growth differs by species
size_pvals <- coef(summary(model_base))
size_pvals <- size_pvals[grep("size_z", rownames(size_pvals)), "Pr(>|t|)"]

ran_crossed <- FALSE
if (any(size_pvals < 0.05, na.rm = TRUE)) {
  message("size_z significant — running crossed model (species x size_z)")
  model_crossed <- lmer(
    as.formula(paste("log(area) ~ time_days * class * size_z", site_term,
                     "+ (1 | plot_id) + (1 | colony_id)")),
    data = data
  )
  summary(model_crossed)
  ran_crossed <- TRUE
} else {
  message("size_z not significant — crossed model skipped")
}


# --- 1c. Post-hoc pairwise comparisons (Tukey-adjusted) ---
# emmeans computes estimated marginal means from the model for each group
# cld() assigns compact letter display: species sharing a letter are NOT significantly different

# species: all pairwise comparisons
em_class     <- emmeans::emmeans(model_base, ~ class)
pairs_class  <- pairs(em_class, adjust = "tukey")
cld_class    <- multcomp::cld(em_class, adjust = "tukey", Letters = letters, reversed = TRUE)
print(pairs_class)

# size_z is continuous — no categorical post-hoc; effect is read directly from the coefficient
# back-transform species means from log scale to cm²
cld_class_plot <- as.data.frame(cld_class) %>%
  mutate(
    mean_area = exp(emmean),
    lower     = exp(lower.CL),
    upper     = exp(upper.CL),
    .group    = trimws(.group)
  )


# Plot G1: per-species average colony area over time — observed mean±SE + LMM prediction
date_origin <- min(data$date)

traj_grid <- expand.grid(
  class     = levels(data$class),
  time_days = seq(min(data$time_days), max(data$time_days), length.out = 60),
  size_z    = 0   # average initial size
)
# predict() needs all fixed-effect variables; add reference site when multi-site
if (multi_site) traj_grid$site <- levels(factor(data$site))[1]
traj_grid <- traj_grid %>%
  mutate(
    date      = date_origin + time_days,
    pred_area = exp(predict(model_base, newdata = ., re.form = NA))
  )

obs_summary <- data %>%
  group_by(class, date) %>%
  summarise(mean_area = mean(area), se_area = sd(area) / sqrt(n()), .groups = "drop")

pG1 <- ggplot() +
  geom_line(data = traj_grid,
            aes(x = date, y = pred_area, color = class),
            linewidth = 1.0) +
  geom_pointrange(data = obs_summary,
                  aes(x = date, y = mean_area,
                      ymin = mean_area - se_area, ymax = mean_area + se_area,
                      color = class),
                  size = 0.3, alpha = 0.7) +
  scale_color_manual(values = pal, guide = "none") +
  scale_x_date(date_breaks = "1 year", date_labels = "'%y") +
  facet_wrap(~ class, scales = "free_y", nrow = 2,
             labeller = as_labeller(sp_labels_n)) +
  labs(x = "Survey date", y = "Colony area (cm²)") +
  theme_coral() +
  theme(strip.text = element_text(family = "Helvetica", face = "plain", size = 7))

# Plot G2: per-species log-area growth rate from model coefficients
# SE for non-reference species requires propagating variance of PAST slope + interaction:
# Var(a+b) = Var(a) + Var(b) + 2*Cov(a,b)
vcov_base  <- vcov(model_base)
past_var   <- vcov_base["time_days", "time_days"]

growth_coefs <- as.data.frame(coef(summary(model_base))) %>%
  mutate(term = rownames(.)) %>%
  filter(grepl("^time_days", term), !grepl("size_z", term)) %>%
  mutate(
    class = ifelse(term == "time_days", "Porites astreoides",
                   gsub("time_days:class", "", term)),
    p     = `Pr(>|t|)`,
    sig   = sig_stars(p)
  )

past_slope <- growth_coefs$Estimate[growth_coefs$class == "Porites astreoides"]

growth_coefs <- growth_coefs %>%
  mutate(
    abs_slope = ifelse(class == "Porites astreoides", Estimate, Estimate + past_slope),
    abs_se    = ifelse(
      class == "Porites astreoides",
      `Std. Error`,
      sqrt(past_var + `Std. Error`^2 + 2 * vcov_base[term, "time_days"])
    ),
    pct_yr    = (exp(abs_slope * 365) - 1) * 100,
    pct_yr_lo = (exp((abs_slope - abs_se) * 365) - 1) * 100,
    pct_yr_hi = (exp((abs_slope + abs_se) * 365) - 1) * 100
  )

pG2 <- ggplot(growth_coefs,
              aes(x = reorder(class, pct_yr), y = pct_yr, color = class)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = accent_col, linewidth = 0.5) +
  geom_errorbar(aes(ymin = pct_yr_lo, ymax = pct_yr_hi),
                width = 0.22, linewidth = 0.65) +
  geom_point(size = 3.2) +
  geom_text(aes(y = pct_yr_hi, label = sig),
            vjust = -0.5, size = 3.5, color = "grey30") +
  scale_color_manual(values = pal, guide = "none") +
  scale_x_discrete(labels = sp_labels_n) +
  coord_flip() +
  labs(x = NULL, y = "Annual growth rate (% per year)") +
  theme_coral() +
  theme(axis.text.y = element_text(size = 7))

ggsave(file.path(plots_dir, "02a_growth_trajectories.jpeg"), pG1,
       width = fig_ls_w, height = fig_ls_h, dpi = 300, units = "px")
ggsave(file.path(plots_dir, "02b_growth_rates.jpeg"), pG2,
       width = fig_sq_w, height = fig_sq_h, dpi = 300, units = "px")


# ═══════════════════════════════════════════════════════════════════════════
# SECTION 2 — TURNOVER: GROWTH vs TISSUE LOSS BALANCE
# ═══════════════════════════════════════════════════════════════════════════

turnover <- data %>%
  arrange(colony_id, date) %>%
  group_by(colony_id, class, plot_id) %>%
  mutate(
    prev_area    = lag(area),
    delta_area   = area - prev_area,
    days_elapsed = as.numeric(date - lag(date))
  ) %>%
  filter(!is.na(delta_area), days_elapsed > 0, prev_area > 0) %>%
  mutate(
    gained_pct = ifelse(delta_area > 0,  delta_area / prev_area * 100 / days_elapsed * 365, 0),
    lost_pct   = ifelse(delta_area < 0, -delta_area / prev_area * 100 / days_elapsed * 365, 0),
    net_pct    = delta_area / prev_area * 100 / days_elapsed * 365
  ) %>%
  ungroup()

# aggregate per-colony first so each colony gets equal weight regardless of survey count
turnover_colony <- turnover %>%
  group_by(colony_id, class) %>%
  summarise(
    mean_gained = mean(gained_pct),
    mean_lost   = mean(lost_pct),
    net_change  = mean(net_pct),
    .groups     = "drop"
  )

turnover_sp <- turnover_colony %>%
  group_by(class) %>%
  summarise(
    mean_gained = mean(mean_gained),
    mean_lost   = mean(mean_lost),
    net_change  = mean(net_change),
    se_net      = sd(net_change) / sqrt(n()),
    .groups     = "drop"
  )

# one-sample t-test per species: is net annual change significantly different from zero?
# H0: mean net change = 0 (no net growth or loss)
turnover_ttest <- turnover_colony %>%
  group_by(class) %>%
  summarise(
    t_p = tryCatch(t.test(net_change, mu = 0)$p.value, error = \(e) NA_real_),
    .groups = "drop"
  ) %>%
  mutate(t_sig = sig_stars(t_p))

turnover_sp <- turnover_sp %>% left_join(turnover_ttest, by = "class")

# Plot T1: growth vs loss diverging bar chart (% per year)
turnover_long <- turnover_sp %>%
  tidyr::pivot_longer(c(mean_gained, mean_lost),
               names_to = "direction", values_to = "pct_yr") %>%
  mutate(
    pct_yr    = ifelse(direction == "mean_lost", -pct_yr, pct_yr),
    direction = ifelse(direction == "mean_gained", "Gained", "Lost")
  )

pT1 <- ggplot(turnover_long,
              aes(x = reorder(class, net_change, mean), y = pct_yr, fill = direction)) +
  geom_col(width = 0.55, alpha = 0.85) +
  geom_point(data = turnover_sp,
             aes(x = reorder(class, net_change, mean), y = net_change),
             inherit.aes = FALSE, size = 2.8, color = "grey10") +
  geom_errorbar(data = turnover_sp,
                aes(x = reorder(class, net_change, mean),
                    ymin = net_change - se_net, ymax = net_change + se_net),
                inherit.aes = FALSE, width = 0.18, linewidth = 0.65, color = "grey10") +
  geom_hline(yintercept = 0, color = accent_col, linewidth = 0.55) +
  geom_text(data = turnover_sp,
            aes(x = reorder(class, net_change, mean),
                y = net_change + se_net, label = t_sig),
            inherit.aes = FALSE, hjust = -0.3, size = 3.2, color = "grey20") +
  scale_fill_manual(values = c(
    "Gained" = viridisLite::viridis(4, option = "D")[3],
    "Lost"   = viridisLite::viridis(4, option = "D")[1]
  )) +
  scale_x_discrete(labels = sp_codes) +
  coord_flip() +
  labs(x = NULL, y = "Mean annual change (% of colony area per year)", fill = NULL) +
  theme_coral() +
  theme(legend.position = "none")

ggsave(file.path(plots_dir, "03_turnover_balance.jpeg"), pT1,
       width = fig_sq_w, height = fig_sq_h, dpi = 300, units = "px")


# ═══════════════════════════════════════════════════════════════════════════
# SECTION 3 — INITIAL SIZE AS EXPLANATORY FACTOR
# ═══════════════════════════════════════════════════════════════════════════

pal_size <- setNames(viridisLite::viridis(3, option = "D", begin = 0.1, end = 0.85),
                     c("small", "medium", "large"))

# Plot S1: initial size vs empirical growth rate — visualises time_days:size_z interaction
colony_growth <- data %>%
  filter(!is.na(size_z)) %>%
  arrange(colony_id, date) %>%
  group_by(colony_id, class, size_z, size_class) %>%
  summarise(
    first_area = area[which.min(date)],
    last_area  = area[which.max(date)],
    total_days = as.numeric(max(date) - min(date)),
    n_obs      = n(),
    .groups    = "drop"
  ) %>%
  filter(n_obs >= 2, total_days > 0, first_area > 0) %>%
  mutate(pct_growth_yr = (exp((log(last_area) - log(first_area)) / total_days * 365) - 1) * 100)

time_coef      <- fixef(model_base)["time_days"]
time_size_coef <- fixef(model_base)["time_days:size_z"]
sz_range       <- range(colony_growth$size_z, na.rm = TRUE)

growth_pred <- data.frame(size_z = seq(sz_range[1], sz_range[2], length.out = 60)) %>%
  mutate(pred_pct = (exp((time_coef + time_size_coef * size_z) * 365) - 1) * 100)

pS1 <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = accent_col, linewidth = 0.5) +
  geom_point(data = colony_growth,
             aes(x = size_z, y = pct_growth_yr, color = class),
             alpha = 0.40, size = 1.6) +
  geom_line(data = growth_pred,
            aes(x = size_z, y = pred_pct),
            color = accent_col, linewidth = 1.1) +
  scale_color_manual(values = pal, name = NULL, labels = sp_codes) +
  labs(x = "Initial size z-score (within species)",
       y = "Growth rate (% per year)") +
  theme_coral() +
  theme(legend.position  = c(0.98, 0.98),
        legend.justification = c("right", "top"),
        legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
        legend.text = element_text(size = 7))

# Plot S2: size class composition per species — in-bar labels, no legend
pS2_dat <- data %>%
  filter(!is.na(size_class)) %>%
  mutate(size_class = factor(size_class, levels = c("small", "medium", "large"))) %>%
  count(class, size_class) %>%
  group_by(class) %>%
  arrange(size_class, .by_group = TRUE) %>%
  mutate(prop = n / sum(n), mid = cumsum(prop) - prop / 2) %>%
  ungroup()

pS2 <- ggplot(pS2_dat, aes(x = class, y = prop, fill = size_class)) +
  geom_col(alpha = 0.85, width = 0.6) +
  scale_fill_manual(values = pal_size) +
  scale_y_continuous(labels = scales::percent_format(),
                     breaks = c(0, 0.5, 1)) +
  scale_x_discrete(labels = sp_codes) +
  coord_flip() +
  labs(x = NULL, y = "Proportion", fill = NULL) +
  guides(fill = "none") +
  theme_coral()

ggsave(file.path(plots_dir, "04a_size_vs_growth.jpeg"), pS1,
       width = fig_sq_w, height = fig_sq_h, dpi = 300, units = "px")
ggsave(file.path(plots_dir, "04b_size_distribution.jpeg"), pS2,
       width = fig_sq_w, height = fig_sq_h, dpi = 300, units = "px")


# ═══════════════════════════════════════════════════════════════════════════
# SECTION 4 — MORTALITY RISK
# size_z and species as predictors; partial + complete mortality combined
# ═══════════════════════════════════════════════════════════════════════════

dead_terminal <- data %>%
  filter(colony_fate == "died", !is.na(size_z)) %>%
  group_by(colony_id) %>%
  slice_max(date, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(partial_mortality = 100, mortality_event = 1L, event_type = "complete death")

pm_data <- data %>%
  filter(!is.na(size_z)) %>%
  mutate(
    mortality_event = as.integer(!is.na(partial_mortality)),
    event_type      = ifelse(!is.na(partial_mortality), "partial mortality", "none")
  ) %>%
  bind_rows(dead_terminal)

# helper: warn if a (g)lmer model has a singular fit or optimizer convergence issues
check_convergence <- function(model, name) {
  if (isSingular(model))
    warning(name, ": singular fit — random effect estimates unreliable")
  msgs <- tryCatch(model@optinfo$conv$lme4$messages, error = function(e) NULL)
  if (!is.null(msgs))
    warning(name, ": convergence issue — ", paste(msgs, collapse = "; "))
}

# 4a. Occurrence
model_pm_occurrence <- glmer(
  mortality_event ~ size_z + class + (1 | colony_id) + (1 | plot_id),
  data = pm_data, family = binomial
)
summary(model_pm_occurrence)
check_convergence(model_pm_occurrence, "model_pm_occurrence")

# 4b. Severity
pm_events <- pm_data %>% filter(mortality_event == 1, partial_mortality > 0)

model_pm_severity <- lmer(
  log(partial_mortality) ~ size_z + class + (1 | colony_id) + (1 | plot_id),
  data = pm_events
)
summary(model_pm_severity)
check_convergence(model_pm_severity, "model_pm_severity")

# 4c. Colony fate
fate_data <- data %>%
  filter(!is.na(size_z)) %>%
  distinct(colony_id, class, site, plot_id, size_z, size_class, colony_fate) %>%
  mutate(fate_binary = as.integer(colony_fate == "died"))

model_fate <- glm(
  as.formula(paste("fate_binary ~ size_z + class", site_term)),
  data = fate_data, family = binomial
)
summary(model_fate)

em_fate_class  <- emmeans::emmeans(model_fate, ~ class, type = "response")
pairs_fate_cl  <- pairs(em_fate_class, adjust = "tukey")
cld_fate_class <- multcomp::cld(em_fate_class, adjust = "tukey",
                                Letters = letters, reversed = TRUE)
print(pairs_fate_cl)

cld_fate_class_plot <- as.data.frame(cld_fate_class) %>% mutate(.group = trimws(.group))

# Plot M1: survival probability by species (Tukey CLD)
pM1 <- ggplot(cld_fate_class_plot,
              aes(x = reorder(class, prob), y = prob, color = class)) +
  geom_point(size = 3.0) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.25, linewidth = 0.65) +
  geom_text(aes(y = asymp.UCL, label = .group),
            hjust = -0.4, size = 4.0, color = "grey30", fontface = "bold") +
  scale_color_manual(values = pal, guide = "none") +
  scale_x_discrete(labels = sp_labels_n) +
  scale_y_continuous(labels = scales::percent_format(),
                     expand = expansion(mult = c(0.05, 0.2))) +
  coord_flip() +
  labs(x = NULL, y = "P(complete colony death)") +
  theme_coral() +
  theme(axis.text.y = element_text(size = 7))

# Plot M2: size_z → % area lost (partial + complete mortality combined)
sev_grid <- data.frame(
  size_z = seq(min(pm_events$size_z, na.rm = TRUE),
               max(pm_events$size_z, na.rm = TRUE), length.out = 60),
  class  = "Porites astreoides"
) %>% mutate(pred_pct = exp(predict(model_pm_severity, newdata = ., re.form = NA)))

pM2 <- ggplot() +
  geom_hline(yintercept = 100, linetype = "dashed", color = accent_col, linewidth = 0.5) +
  geom_jitter(data = pm_events,
              aes(x = size_z, y = partial_mortality, color = event_type),
              alpha = 0.40, width = 0.05, height = 0, size = 1.4) +
  geom_line(data = sev_grid,
            aes(x = size_z, y = pred_pct),
            color = accent_col, linewidth = 1.1) +
  scale_color_manual(
    name   = NULL,
    values = c("partial mortality" = accent_bg_col, "complete death" = "#3d2b1f"),
    labels = c("partial mortality" = "Partial mortality", "complete death" = "Complete death (100%)")
  ) +
  scale_y_continuous(labels = \(x) paste0(x, "%"), limits = c(0, 108)) +
  labs(x = "Initial size z-score (within species)", y = "Partial mortality (%)") +
  theme_coral() +
  theme(legend.position      = c(0.98, 0.5),
        legend.justification = c("right", "center"),
        legend.background    = element_rect(fill = alpha("white", 0.7), color = NA),
        legend.text          = element_text(size = 7))

ggsave(file.path(plots_dir, "05a_mortality_risk_by_species.jpeg"), pM1,
       width = fig_sq_w, height = fig_sq_h, dpi = 300, units = "px")
ggsave(file.path(plots_dir, "05b_mortality_severity_by_size.jpeg"), pM2,
       width = fig_sq_w, height = fig_sq_h, dpi = 300, units = "px")


# ═══════════════════════════════════════════════════════════════════════════
# SECTION 6 — SITE COMPARISON
# only produced when data contains more than one site
# ═══════════════════════════════════════════════════════════════════════════

if (multi_site) {

  # per-plot cover for site comparison (individual plot observations as points)
  cover_site_plots <- data %>%
    group_by(site, zone, plot_id, date, time_days) %>%
    summarise(total_area = sum(area, na.rm = TRUE), .groups = "drop") %>%
    mutate(pct_cover = total_area / PLOT_AREA_CM2 * 100)

  pSC1 <- ggplot(cover_site_plots, aes(x = date, y = pct_cover, color = zone, fill = zone)) +
    geom_point(alpha = 0.55, size = 2.0) +
    geom_smooth(method = "loess", formula = y ~ x, se = TRUE, alpha = 0.18,
                linewidth = 1.0) +
    scale_color_manual(values = site_pal) +
    scale_fill_manual(values  = site_pal, guide = "none") +
    scale_x_date(date_breaks = "6 months", date_labels = "%b '%y") +
    scale_y_continuous(labels = \(x) paste0(round(x, 1), "%")) +
    labs(x = "Survey date", y = "Total cover (% of 25 m²)", color = NULL) +
    theme_coral() +
    theme(legend.position        = c(0.02, 0.98),
          legend.justification   = c("left", "top"),
          legend.background      = element_rect(fill = alpha("white", 0.7), color = NA))

  # species composition per site: proportion of colonies
  sp_site <- data %>%
    distinct(colony_id, zone, class) %>%
    count(zone, class, name = "n_colonies") %>%
    group_by(zone) %>%
    mutate(prop = n_colonies / sum(n_colonies)) %>%
    ungroup()

  pSC2 <- ggplot(sp_site, aes(x = zone, y = prop, fill = class)) +
    geom_col(width = 0.5, alpha = 0.88) +
    scale_fill_manual(values = pal, labels = sp_codes) +
    scale_y_continuous(labels = scales::percent_format()) +
    coord_flip() +
    labs(x = NULL, y = "Proportion of colonies", fill = NULL) +
    theme_coral() +
    theme(legend.text     = element_text(size = 7),
          legend.position = "bottom") +
    guides(fill = guide_legend(nrow = 2))

  # initial size distribution per site
  init_size_site <- data %>%
    group_by(colony_id, zone) %>%
    slice_min(date, n = 1, with_ties = FALSE) %>%
    ungroup()

  pSC3 <- ggplot(init_size_site, aes(x = zone, y = log(area), color = zone, fill = zone)) +
    geom_violin(alpha = 0.18, linewidth = 0.55, trim = FALSE) +
    geom_boxplot(width = 0.12, alpha = 0.85, outlier.shape = NA, linewidth = 0.5) +
    scale_color_manual(values = site_pal, guide = "none") +
    scale_fill_manual(values  = site_pal, guide = "none") +
    labs(x = NULL, y = "log(colony area, cm²)") +
    theme_coral()

  # --- statistical tests for site comparison ---

  # 7a: LMM testing whether sites differ in cover level and rate of change
  # zone main effect = baseline cover difference; time_days:zone = diverging trajectories
  model_cover_site <- lmer(
    log(pct_cover) ~ time_days * zone + (1 | plot_id),
    data = cover_site_plots
  )
  summary(model_cover_site)

  # 7b: chi-square test of species composition independence across sites
  sp_comp_mat <- sp_site %>%
    tidyr::pivot_wider(id_cols = zone, names_from = class,
                       values_from = n_colonies, values_fill = 0L) %>%
    tibble::column_to_rownames("zone") %>%
    as.matrix()
  chisq_sp_comp <- chisq.test(sp_comp_mat)
  message("Species composition chi-square: X2=", round(chisq_sp_comp$statistic, 2),
          "  df=", chisq_sp_comp$parameter,
          "  p=", round(chisq_sp_comp$p.value, 4))

  # 7c: Wilcoxon rank-sum test comparing initial colony size between sites
  # non-parametric: size distributions are right-skewed even after log-transform
  wilcox_size_site <- wilcox.test(log(area) ~ zone, data = init_size_site)
  message("Initial size Wilcoxon: W=", wilcox_size_site$statistic,
          "  p=", round(wilcox_size_site$p.value, 4))

  ggsave(file.path(plots_dir, "07a_site_cover_trajectory.jpeg"), pSC1,
         width = fig_ls_w, height = fig_ls_h, dpi = 300, units = "px")
  ggsave(file.path(plots_dir, "07b_site_species_composition.jpeg"), pSC2,
         width = fig_sq_w, height = fig_sq_h, dpi = 300, units = "px")
  ggsave(file.path(plots_dir, "07c_site_initial_size.jpeg"), pSC3,
         width = fig_sq_w, height = fig_sq_h, dpi = 300, units = "px")

} else {
  message("Single site — site comparison plot skipped")
}


# ═══════════════════════════════════════════════════════════════════════════
# APPENDIX — SPECIES BIOLOGY (mean size, Tukey CLD)
# ═══════════════════════════════════════════════════════════════════════════

# sig_mean_size: derive significance vs PAST from pairs already computed in Section 1c
sig_mean_size <- as.data.frame(pairs_class) %>%
  filter(grepl("Porites astreoides", contrast)) %>%
  mutate(
    other_sp = trimws(gsub("Porites astreoides - | - Porites astreoides", "", contrast)),
    stars    = sig_stars(p.value)
  ) %>%
  filter(stars != "") %>%
  { setNames(.$stars, .$other_sp) }

# per-cell: mean area and n colonies
size_means <- data %>%
  filter(!is.na(size_class)) %>%
  group_by(class, size_class) %>%
  summarise(
    mean_area  = mean(area, na.rm = TRUE),
    n_colonies = n_distinct(colony_id),
    .groups    = "drop"
  )

# species total n for y-axis labels
species_n_total <- size_means %>%
  group_by(class) %>%
  summarise(n_sp = sum(n_colonies), .groups = "drop")

# size-class total row: grand mean from raw observations (not mean of species means)
sizeclass_totals <- data %>%
  filter(!is.na(size_class)) %>%
  group_by(size_class) %>%
  summarise(
    mean_area  = mean(area, na.rm = TRUE),
    n_colonies = n_distinct(colony_id),
    class      = "All species",
    .groups    = "drop"
  )

# assemble species rows with enriched labels
size_plot_data <- size_means %>%
  left_join(species_n_total, by = "class") %>%
  mutate(
    class_label = paste0(
      sp_codes[class],
      ifelse(class %in% names(sig_mean_size), paste0(" ", sig_mean_size[class]), ""),
      "
n = ", n_sp
    )
  ) %>%
  bind_rows(
    sizeclass_totals %>% mutate(n_sp = NA_integer_, class_label = "All species")
  )

# factor: species sorted by mean area ascending, total row at top
sp_label_order <- size_plot_data %>%
  filter(class != "All species") %>%
  group_by(class_label) %>%
  summarise(m = mean(mean_area), .groups = "drop") %>%
  arrange(m) %>%
  pull(class_label)

size_plot_data$class_label <- factor(
  size_plot_data$class_label,
  levels = c(sp_label_order, "All species")
)

# species use the canonical pal; "All species" aggregate row gets grey
dot_colors <- c(pal, "All species" = "grey45")

# separator y-position: just above the last species row
sep_y <- length(sp_label_order) + 0.5

p3 <- ggplot(size_plot_data,
             aes(x = size_class, y = class_label, size = mean_area, color = class)) +
  geom_hline(yintercept = sep_y, linetype = "dashed",
             color = accent_bg_col, linewidth = 0.5) +
  geom_point(alpha = 0.85) +
  geom_text(aes(label = n_colonies), size = 2.5, color = "black", fontface = "bold") +
  geom_text(aes(label = ifelse(!is.na(mean_area), paste0(round(mean_area), " cm² →"), "")),
            hjust = 1.15, vjust = 0.5, size = 2.2, color = "grey30") +
  scale_size_area(max_size = 22, guide = "none") +
  scale_color_manual(values = dot_colors, guide = "none") +
  labs(x = "Initial size class", y = NULL) +
  theme_coral() +
  theme(axis.text.y = element_text(size = 7, colour = "grey20"))

# Appendix plot A2: Tukey CLD for species mean size
pA2 <- ggplot(cld_class_plot,
              aes(x = reorder(class, mean_area), y = mean_area, color = class)) +
  geom_point(size = 3.0) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.25, linewidth = 0.65) +
  geom_text(aes(y = upper, label = .group),
            hjust = -0.4, size = 4.0, color = "grey30", fontface = "bold") +
  scale_color_manual(values = pal, guide = "none") +
  scale_x_discrete(labels = sp_labels_n) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  coord_flip() +
  labs(x = NULL, y = "Estimated mean area (cm²)") +
  theme_coral() +
  theme(axis.text.y = element_text(size = 7))

ggsave(file.path(plots_dir, "06a_species_size_bubble.jpeg"), p3,
       width = fig_ls_w, height = fig_ls_h, dpi = 300, units = "px")
ggsave(file.path(plots_dir, "06b_species_size_tukey.jpeg"), pA2,
       width = fig_sq_w, height = fig_sq_h, dpi = 300, units = "px")


# ═══════════════════════════════════════════════════════════════════════════
# SECTION 5 — RESULTS SUMMARY DATAFRAME
# ═══════════════════════════════════════════════════════════════════════════

# helper: extract fixed effects from any lmer / glmer model
extract_fixed <- function(model, source, response) {
  s     <- coef(summary(model))
  pcol  <- grep("Pr\\(", colnames(s), value = TRUE)[1]
  # match "t value" (lmer) or "z value" (glmer); exclude "Std. Error"
  tcol  <- grep("value$", colnames(s), value = TRUE)
  tcol  <- tcol[grepl("^[tz]", tcol)][1]
  data.frame(
    source    = source,
    response  = response,
    type      = "fixed effect",
    term      = rownames(s),
    estimate  = round(as.numeric(s[, "Estimate"]),   4),
    se        = round(as.numeric(s[, "Std. Error"]), 4),
    statistic = round(as.numeric(s[, tcol]),         3),
    p_value   = round(as.numeric(s[, pcol]),         4),
    row.names = NULL
  )
}

# helper: extract emmeans pairwise contrasts
extract_pairs <- function(pairs_obj, source, response) {
  df   <- as.data.frame(pairs_obj)
  # estimate column: "estimate" (log scale) or "odds.ratio"/"ratio" (response scale)
  ecol <- grep("^estimate$|^odds\\.ratio$|^ratio$", colnames(df), value = TRUE)[1]
  # statistic column: "t.ratio" or "z.ratio" — exclude "odds.ratio"
  rcol <- grep("^[tz]\\.ratio$", colnames(df), value = TRUE)[1]
  data.frame(
    source    = source,
    response  = response,
    type      = "pairwise (Tukey)",
    term      = as.character(df$contrast),
    estimate  = round(as.numeric(df[, ecol]),   4),
    se        = round(as.numeric(df$SE),        4),
    statistic = round(as.numeric(df[, rcol]),   3),
    p_value   = round(as.numeric(df$p.value),   4),
    row.names = NULL
  )
}

# turnover t-test results (one row per species)
turnover_ttest_results <- turnover_ttest %>%
  transmute(
    source    = "t-test: net turnover vs zero",
    response  = "net annual change (% per year)",
    type      = "one-sample t-test",
    term      = as.character(class),
    estimate  = round(turnover_sp$net_change[match(class, turnover_sp$class)], 4),
    se        = round(turnover_sp$se_net[match(class,    turnover_sp$class)], 4),
    statistic = NA_real_,
    p_value   = round(t_p, 4),
    sig       = t_sig
  )

# collect all results
results_all <- bind_rows(
  extract_fixed(model_cover,         "LMM: total cover trend",    "log(% cover)"),
  extract_fixed(model_base,          "LMM: colony growth",        "log(area)"),
  extract_fixed(model_pm_occurrence, "GLMM: mortality occurrence","P(any mortality)"),
  extract_fixed(model_pm_severity,   "LMM: mortality severity",   "log(% area lost)"),
  extract_fixed(model_fate,          "GLMM: colony fate",         "P(colony died)"),
  extract_pairs(pairs_class,         "post-hoc: species size",    "log(area)"),
  extract_pairs(pairs_fate_cl,       "post-hoc: colony fate",     "P(colony died)"),
  turnover_ttest_results
) %>%
  mutate(sig = sig_stars(p_value)) %>%
  dplyr::arrange(source, p_value)

# also add crossed model if it was fitted this run
if (ran_crossed) {
  results_all <- bind_rows(
    results_all,
    extract_fixed(model_crossed, "LMM: growth (species x size)", "log(area)") %>%
      mutate(sig = sig_stars(p_value))
  )
}

# add site-comparison results if multi-site run was performed
if (multi_site) {
  results_all <- bind_rows(
    results_all,
    extract_fixed(model_cover_site, "LMM: site cover comparison", "log(% cover)") %>%
      mutate(sig = sig_stars(p_value)),
    data.frame(
      source    = "chi-square: species composition by site",
      response  = "colony counts",
      type      = "chi-square test of independence",
      term      = "all species × site",
      estimate  = NA_real_,
      se        = NA_real_,
      statistic = round(as.numeric(chisq_sp_comp$statistic), 3),
      p_value   = round(chisq_sp_comp$p.value, 4),
      sig       = sig_stars(chisq_sp_comp$p.value)
    ),
    data.frame(
      source    = "Wilcoxon: initial colony size by site",
      response  = "log(colony area)",
      type      = "Wilcoxon rank-sum test",
      term      = "zone",
      estimate  = NA_real_,
      se        = NA_real_,
      statistic = round(as.numeric(wilcox_size_site$statistic), 3),
      p_value   = round(wilcox_size_site$p.value, 4),
      sig       = sig_stars(wilcox_size_site$p.value)
    )
  )
}

# significant results only (p < 0.05), annotated with low-power species flag
results_sig <- results_all %>%
  filter(p_value < 0.05) %>%
  dplyr::select(source, response, type, term, estimate, se, statistic, p_value, sig) %>%
  mutate(power_note = dplyr::case_when(
    sapply(term, function(t) any(sapply(very_low_power_spp, function(sp) grepl(sp, t, fixed = TRUE)))) ~ "very low power (n 5-9)",
    sapply(term, function(t) any(sapply(low_power_spp,      function(sp) grepl(sp, t, fixed = TRUE)))) ~ "low power (n 10-25)",
    TRUE ~ ""
  ))

print(results_sig)
