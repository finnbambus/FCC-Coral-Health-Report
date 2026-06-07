# load packages
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(lme4, lmerTest, dplyr, tidyr, ggplot2, wesanderson, patchwork, emmeans, multcomp, multcompView, scales, pwr)

# set wd to project folder
this_file <- rstudioapi::getSourceEditorContext()$path
setwd(dirname(this_file))
setwd("..")

plots_dir <- "plots"
dir.create(plots_dir, showWarnings = FALSE)

# read data
data <- read.csv("./data/combined_data.csv", sep = ";")

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
  group_by(colony_id, class) %>%
  slice_min(date, n = 1, with_ties = FALSE) %>%
  group_by(class) %>%
  mutate(
    log_area   = log(area),
    size_z     = if (n() < 2 || sd(log_area) == 0) 0
                 else (log_area - mean(log_area)) / sd(log_area),
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
    n_power_min  = n_power_min,
    adequate_power = n_colonies >= n_power_min,
    low_power      = n_colonies >= n_hard_min & !adequate_power,
    too_few        = n_colonies < n_hard_min
  )

message("Species n check  |  hard floor: ", n_hard_min,
        "  |  power threshold (d=", target_d, ", 1-β=", target_power, "): ", n_power_min, " colonies")
print(as.data.frame(species_n))

# apply only the hard floor; power-flagged species stay in but are marked in results
excluded_spp  <- species_n %>% filter(too_few) %>% pull(class)
low_power_spp <- species_n %>% filter(low_power) %>% pull(class)

if (length(excluded_spp) > 0) {
  message("Removing (below hard floor n=", n_hard_min, "): ",
          paste(excluded_spp, collapse = ", "))
  data <- data %>% filter(!class %in% excluded_spp)
}
if (length(low_power_spp) > 0) {
  message("Low-power species (flagged in results, not removed): ",
          paste(low_power_spp, collapse = ", "))
}

# reference species: PAST has most colonies → lowest SE on all pairwise comparisons
data$class <- relevel(factor(data$class), ref = "Porites astreoides")

# plot area constant for % cover conversion (5 × 5 m photoquadrats)
PLOT_AREA_CM2 <- 250000   # 25 m² = 250 000 cm²

# shared model flags
multi_site <- length(unique(data$site)) > 1
site_term  <- if (multi_site) "+ site" else ""

# color palette: one color per species, interpolated from Moonrise3
n_species <- length(unique(data$class))
pal <- wesanderson::wes_palette("Moonrise3", n_species, type = "continuous")
names(pal) <- sort(unique(data$class))


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
  group_by(plot_id, date, time_days) %>%
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

# Plot C1: stacked species cover — shows total cover AND composition shift
cover_sp_mean <- cover_sp %>%
  group_by(date, class) %>%
  summarise(mean_pct = mean(pct_cover_sp), .groups = "drop")

pC1 <- ggplot(cover_sp_mean, aes(x = date, y = mean_pct, fill = class)) +
  geom_area(alpha = 0.88, color = "white", linewidth = 0.25) +
  scale_fill_manual(values = pal) +
  scale_y_continuous(labels = \(x) paste0(x, "%"),
                     expand = expansion(mult = c(0, 0.03))) +
  labs(x = NULL, y = "Mean coral cover (% of 25 m²)", fill = NULL,
       title = "Species composition of total coral cover over time") +
  theme_minimal(base_size = 11) +
  theme(legend.position  = "bottom",
        legend.text      = element_text(face = "italic", size = 8),
        panel.grid.minor = element_blank())

# Plot C2: per-plot cover trajectories with model trend line
pC2 <- ggplot(cover, aes(x = date, y = pct_cover)) +
  geom_line(aes(group = plot_id, color = plot_id), alpha = 0.55, linewidth = 0.9) +
  geom_point(aes(color = plot_id), size = 2.2, alpha = 0.8) +
  geom_ribbon(data = cover_trend,
              aes(x = date, ymin = pct_cover * 0.85, ymax = pct_cover * 1.15),
              inherit.aes = FALSE, alpha = 0.12, fill = "grey30") +
  geom_line(data = cover_trend, aes(x = date, y = pct_cover),
            inherit.aes = FALSE, color = "grey15", linewidth = 1.5) +
  scale_y_continuous(labels = \(x) paste0(round(x, 1), "%")) +
  labs(x = NULL, y = "Total coral cover (% of 25 m²)", color = "Plot",
       title = "Total cover trajectory per plot",
       subtitle = "Dark line = LMM-predicted trend  |  Band = ±15% envelope") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())

p_community <- pC1 / pC2 +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(
    title = "Community coral cover — Farallon site (S2P2)",
    theme = theme(plot.title = element_text(face = "bold", size = 13))
  )
print(p_community)
ggsave(file.path(plots_dir, "01_community_cover.jpeg"), p_community,
       width = 11.69, height = 8, dpi = 300, units = "in")


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
  geom_point(alpha = 0.25, size = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_smooth(method = "loess", se = FALSE, color = pal[[3]], linewidth = 0.8) +
  labs(x = "Fitted values", y = "Residuals", title = "Residuals vs Fitted") +
  theme_minimal()

# QQ plot of residuals: checks normality of residuals
# expect points to fall along the line
pb <- ggplot(resid_df, aes(sample = residuals)) +
  stat_qq(alpha = 0.3, size = 0.8) +
  stat_qq_line(color = pal[[3]]) +
  labs(x = "Theoretical quantiles", y = "Sample quantiles",
       title = "Normal Q-Q — residuals") +
  theme_minimal()

# scale-location plot: second check for homoscedasticity
# expect flat smoother (constant spread across fitted range)
pc <- ggplot(resid_df, aes(x = fitted, y = sqrt_abs_resid)) +
  geom_point(alpha = 0.25, size = 0.8) +
  geom_smooth(method = "loess", se = FALSE, color = pal[[3]], linewidth = 0.8) +
  labs(x = "Fitted values", y = "√|Residuals|", title = "Scale-Location") +
  theme_minimal()

# QQ plot of colony random effects: checks normality assumption for random intercepts
re_vals <- ranef(model_base)$colony_id[, 1]
pd <- ggplot(data.frame(re = re_vals), aes(sample = re)) +
  stat_qq(alpha = 0.4, size = 0.8) +
  stat_qq_line(color = pal[[3]]) +
  labs(x = "Theoretical quantiles", y = "Sample quantiles",
       title = "Normal Q-Q — random effects (colony)") +
  theme_minimal()

print((pa | pb) / (pc | pd) +
  plot_annotation(title = "LMM assumption checks",
                  theme = theme(plot.title = element_text(face = "bold"))))


# model 1b: species x size_z interaction — only run if size_z is significant in model_base
# tests whether the continuous size effect on growth differs by species
size_pvals <- coef(summary(model_base))
size_pvals <- size_pvals[grep("size_z", rownames(size_pvals)), "Pr(>|t|)"]

if (any(size_pvals < 0.05, na.rm = TRUE)) {
  message("size_z significant — running crossed model (species x size_z)")
  model_crossed <- lmer(
    as.formula(paste("log(area) ~ time_days * class * size_z", site_term,
                     "+ (1 | plot_id) + (1 | colony_id)")),
    data = data
  )
  summary(model_crossed)
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
) %>%
  mutate(
    date      = date_origin + time_days,
    pred_area = exp(predict(model_base, newdata = ., re.form = NA))
  )

obs_summary <- data %>%
  group_by(class, date) %>%
  summarise(mean_area = mean(area), se_area = sd(area) / sqrt(n()), .groups = "drop")

pG1 <- ggplot() +
  geom_ribbon(data = traj_grid,
              aes(x = date, ymin = pred_area * 0.7, ymax = pred_area * 1.3, fill = class),
              alpha = 0.15) +
  geom_line(data = traj_grid,
            aes(x = date, y = pred_area, color = class),
            linewidth = 1.1) +
  geom_pointrange(data = obs_summary,
                  aes(x = date, y = mean_area,
                      ymin = mean_area - se_area, ymax = mean_area + se_area,
                      color = class),
                  size = 0.35, alpha = 0.75) +
  scale_color_manual(values = pal, guide = "none") +
  scale_fill_manual(values = pal, guide = "none") +
  facet_wrap(~ class, scales = "free_y", ncol = 3) +
  labs(x = NULL, y = "Mean colony area (cm²)",
       title = "Species growth trajectories",
       subtitle = "Line = LMM prediction at mean initial size  |  Points = observed mean ± SE  |  Band = ±30%") +
  theme_minimal(base_size = 10) +
  theme(strip.text = element_text(face = "italic", size = 8))

# Plot G2: per-species log-area growth rate from model coefficients
growth_coefs <- as.data.frame(coef(summary(model_base))) %>%
  mutate(term = rownames(.)) %>%
  filter(grepl("^time_days", term), !grepl("size_z", term)) %>%
  mutate(
    class     = ifelse(term == "time_days", "Porites astreoides",
                       gsub("time_days:class", "", term)),
    se        = `Std. Error`,
    p         = `Pr(>|t|)`,
    sig       = dplyr::case_when(p < 0.001 ~ "***", p < 0.01 ~ "**",
                                 p < 0.05 ~ "*", TRUE ~ "ns")
  )

past_slope <- growth_coefs$Estimate[growth_coefs$class == "Porites astreoides"]
growth_coefs <- growth_coefs %>%
  mutate(abs_slope = ifelse(class == "Porites astreoides", Estimate,
                            Estimate + past_slope))

growth_coefs <- growth_coefs %>%
  mutate(
    pct_yr     = (exp(abs_slope * 365) - 1) * 100,
    pct_yr_lo  = (exp((abs_slope - se) * 365) - 1) * 100,
    pct_yr_hi  = (exp((abs_slope + se) * 365) - 1) * 100
  )

pG2 <- ggplot(growth_coefs,
              aes(x = reorder(class, pct_yr), y = pct_yr, color = class)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbar(aes(ymin = pct_yr_lo, ymax = pct_yr_hi),
                width = 0.25, linewidth = 0.7) +
  geom_point(size = 3.5) +
  geom_text(aes(y = pct_yr_hi, label = sig),
            vjust = -0.5, size = 3.8, color = "grey20") +
  scale_color_manual(values = pal, guide = "none") +
  coord_flip() +
  labs(x = NULL, y = "Annual growth rate (% per year)",
       title = "Per-species annual growth rate",
       subtitle = "Back-transformed from LMM: (exp(slope × 365) − 1) × 100  |  Positive = growing") +
  theme_minimal() +
  theme(axis.text.y = element_text(face = "italic"))

p_growth <- pG1 / pG2 +
  plot_layout(heights = c(2, 1)) +
  plot_annotation(title = "Species growth trajectories",
                  theme = theme(plot.title = element_text(face = "bold", size = 13)))
print(p_growth)
ggsave(file.path(plots_dir, "02_species_growth_trajectories.jpeg"), p_growth,
       width = 11.69, height = 9, dpi = 300, units = "in")


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

turnover_sp <- turnover %>%
  group_by(class) %>%
  summarise(
    mean_gained = mean(gained_pct),
    mean_lost   = mean(lost_pct),
    net_change  = mean(net_pct),
    se_net      = sd(net_pct) / sqrt(n()),
    .groups     = "drop"
  )

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
  geom_col(width = 0.6, alpha = 0.88) +
  geom_point(data = turnover_sp,
             aes(x = reorder(class, net_change, mean), y = net_change),
             inherit.aes = FALSE, size = 3.2, color = "grey10") +
  geom_errorbar(data = turnover_sp,
                aes(x = reorder(class, net_change, mean),
                    ymin = net_change - se_net, ymax = net_change + se_net),
                inherit.aes = FALSE, width = 0.2, linewidth = 0.7, color = "grey10") +
  geom_hline(yintercept = 0, linewidth = 0.6) +
  scale_fill_manual(values = c("Gained" = "#7fba84", "Lost" = "#c8706b")) +
  coord_flip() +
  labs(x = NULL, y = "Mean annual change (% of colony area per year)", fill = NULL,
       title = "Growth vs. tissue loss balance per species",
       subtitle = "Bars = mean % area gained/lost per year  |  Dot = net change ± SE") +
  theme_minimal() +
  theme(axis.text.y = element_text(face = "italic"), legend.position = "top")

print(pT1)
ggsave(file.path(plots_dir, "03_turnover_balance.jpeg"), pT1,
       width = 11.69, height = 5, dpi = 300, units = "in")


# ═══════════════════════════════════════════════════════════════════════════
# SECTION 3 — INITIAL SIZE AS EXPLANATORY FACTOR
# ═══════════════════════════════════════════════════════════════════════════

pal_size <- setNames(wesanderson::wes_palette("Moonrise3")[c(1, 3, 5)],
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
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey55") +
  geom_point(data = colony_growth,
             aes(x = size_z, y = pct_growth_yr, color = class),
             alpha = 0.45, size = 1.8) +
  geom_line(data = growth_pred,
            aes(x = size_z, y = pred_pct),
            color = "grey10", linewidth = 1.2) +
  scale_color_manual(values = pal, name = NULL) +
  labs(x = "Initial size z-score (within species)",
       y = "Growth rate (% per year)",
       title = "Initial colony size vs. subsequent growth rate",
       subtitle = "Larger colonies grow proportionally slower (ontogenetic slowdown)  |  Line = model (PAST reference)") +
  theme_minimal() +
  theme(legend.position = "top")

# Plot S2: size class composition per species
pS2 <- ggplot(filter(data, !is.na(size_class)), aes(x = class, fill = size_class)) +
  geom_bar(position = "fill", alpha = 0.88, width = 0.65) +
  scale_fill_manual(values = pal_size, labels = c("Small", "Medium", "Large")) +
  scale_y_continuous(labels = scales::percent_format()) +
  coord_flip() +
  labs(x = NULL, y = "Proportion of observations", fill = "Initial size class",
       title = "Initial size class composition per species") +
  theme_minimal() +
  theme(axis.text.y = element_text(face = "italic"), legend.position = "top")

p_size <- pS1 + pS2 +
  plot_annotation(title = "Initial colony size as explanatory factor",
                  theme = theme(plot.title = element_text(face = "bold")))
print(p_size)
ggsave(file.path(plots_dir, "04_initial_size_effects.jpeg"), p_size,
       width = 11.69, height = 5.5, dpi = 300, units = "in")


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

# 4a. Occurrence
model_pm_occurrence <- glmer(
  mortality_event ~ size_z + class + (1 | colony_id) + (1 | plot_id),
  data = pm_data, family = binomial
)
summary(model_pm_occurrence)

# 4b. Severity
pm_events <- pm_data %>% filter(mortality_event == 1, partial_mortality > 0)

model_pm_severity <- lmer(
  log(partial_mortality) ~ size_z + class + (1 | colony_id) + (1 | plot_id),
  data = pm_events
)
summary(model_pm_severity)

# 4c. Colony fate
fate_summary <- data %>%
  filter(!is.na(size_class)) %>%
  distinct(colony_id, size_class, class, colony_fate) %>%
  group_by(class, size_class, colony_fate) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(class, size_class) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

fate_data <- data %>%
  filter(!is.na(size_z)) %>%
  distinct(colony_id, class, site, plot_id, size_z, size_class, colony_fate) %>%
  mutate(fate_binary = as.integer(colony_fate == "died"))

model_fate <- glmer(
  as.formula(paste("fate_binary ~ size_z + class", site_term, "+ (1 | plot_id)")),
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
  geom_point(size = 3.5) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.3, linewidth = 0.7) +
  geom_text(aes(y = asymp.UCL, label = .group),
            hjust = -0.4, size = 4.5, color = "grey25", fontface = "bold") +
  scale_color_manual(values = pal, guide = "none") +
  scale_y_continuous(labels = scales::percent_format(),
                     expand = expansion(mult = c(0.05, 0.2))) +
  coord_flip() +
  labs(x = NULL, y = "P(complete colony death)",
       title = "Complete mortality risk by species (Tukey post-hoc)") +
  theme_minimal() +
  theme(axis.text.y = element_text(face = "italic"))

# Plot M2: size_z → % area lost (partial + complete mortality combined)
sev_grid <- data.frame(
  size_z = seq(min(pm_events$size_z, na.rm = TRUE),
               max(pm_events$size_z, na.rm = TRUE), length.out = 60),
  class  = "Porites astreoides"
) %>% mutate(pred_pct = exp(predict(model_pm_severity, newdata = ., re.form = NA)))

pM2 <- ggplot() +
  geom_hline(yintercept = 100, linetype = "dashed", color = "grey55", linewidth = 0.5) +
  geom_jitter(data = pm_events,
              aes(x = size_z, y = partial_mortality, color = event_type),
              alpha = 0.45, width = 0.05, height = 0, size = 1.6) +
  geom_line(data = sev_grid,
            aes(x = size_z, y = pred_pct),
            color = "grey15", linewidth = 1.2) +
  scale_color_manual(
    values = c("partial mortality" = "#c8706b", "complete death" = "#3d2b1f"),
    labels = c("partial mortality" = "Partial mortality", "complete death" = "Complete death (100%)")
  ) +
  scale_y_continuous(labels = \(x) paste0(x, "%"), limits = c(0, 108)) +
  labs(x = "Initial size z-score (within species)", y = "Area lost (%)",
       color = NULL,
       title = "Size effect on mortality severity",
       subtitle = "Line = LMM prediction (reference: PAST)  |  Dashed = complete colony death threshold") +
  theme_minimal() +
  theme(legend.position = "top")

p_mortality <- (pM1 | pM2) +
  plot_annotation(title = "Colony mortality risk — species and initial size",
                  theme = theme(plot.title = element_text(face = "bold", size = 13)))
print(p_mortality)
ggsave(file.path(plots_dir, "05_mortality_risk.jpeg"), p_mortality,
       width = 11.69, height = 5.5, dpi = 300, units = "in")


# ═══════════════════════════════════════════════════════════════════════════
# APPENDIX — SPECIES BIOLOGY (mean size, Tukey CLD)
# ═══════════════════════════════════════════════════════════════════════════

# post-hoc species mean size
em_class    <- emmeans::emmeans(model_base, ~ class)
pairs_class <- pairs(em_class, adjust = "tukey")
cld_class   <- multcomp::cld(em_class, adjust = "tukey", Letters = letters, reversed = TRUE)
print(pairs_class)

cld_class_plot <- as.data.frame(cld_class) %>%
  mutate(mean_area = exp(emmean), lower = exp(lower.CL), upper = exp(upper.CL),
         .group = trimws(.group))

sig_mean_size <- c(
  "Pseudodiploria strigosa"    = "**",
  "Siderastrea siderea"       = "***",
  "Stephanocoenia intersepta" = "***"
)

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

# size-class total row: mean of species means + total n per size class
sizeclass_totals <- size_means %>%
  group_by(size_class) %>%
  summarise(
    mean_area  = mean(mean_area, na.rm = TRUE),
    n_colonies = sum(n_colonies),
    class      = "All species",
    .groups    = "drop"
  )

# assemble species rows with enriched labels
size_plot_data <- size_means %>%
  left_join(species_n_total, by = "class") %>%
  mutate(
    class_label = paste0(
      class,
      ifelse(class %in% names(sig_mean_size), paste0("  ", sig_mean_size[class]), ""),
      "  (n=", n_sp, ")"
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

# colors: species use pal; total row grey
dot_colors <- c(pal[names(pal) %in% unique(size_means$class)],
                "All species" = "grey45")

# separator y-position: just above the last species row
sep_y <- length(sp_label_order) + 0.5

p3 <- ggplot(size_plot_data,
             aes(x = size_class, y = class_label, size = mean_area, color = class)) +
  geom_hline(yintercept = sep_y, linetype = "dashed",
             color = "grey65", linewidth = 0.4) +
  geom_point(alpha = 0.85) +
  geom_text(aes(label = n_colonies), size = 2.6, color = "white", fontface = "bold") +
  geom_text(aes(label = ifelse(!is.na(mean_area), paste0(round(mean_area), " cm²"), "")),
            vjust = 3.2, size = 2.3, color = "grey35") +
  scale_size_area(max_size = 22, guide = "none") +
  scale_color_manual(values = dot_colors, guide = "none") +
  labs(x = "Initial size class", y = NULL,
       title = "Mean colony area by species and initial size class",
       subtitle = "Circle labels: n colonies per cell  |  Below circles: mean area  |  Top row: size-class totals") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(
      face   = c(rep("italic", length(sp_label_order)), "bold"),
      colour = c(rep("grey20", length(sp_label_order)), "grey40")
    )
  )

# Appendix plot A2: Tukey CLD for species mean size
pA2 <- ggplot(cld_class_plot,
              aes(x = reorder(class, mean_area), y = mean_area, color = class)) +
  geom_point(size = 3.5) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.3, linewidth = 0.7) +
  geom_text(aes(y = upper, label = .group),
            hjust = -0.4, size = 4.5, color = "grey25", fontface = "bold") +
  scale_color_manual(values = pal, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  coord_flip() +
  labs(x = NULL, y = "Estimated mean area (cm²)",
       title = "Appendix: pairwise species size comparison (Tukey)",
       subtitle = "Shared letter = not significantly different  |  Geometric mean ± 95% CI") +
  theme_minimal() +
  theme(axis.text.y = element_text(face = "italic"))

p_appendix <- p3 + pA2 +
  plot_annotation(title = "Species biology reference",
                  theme = theme(plot.title = element_text(face = "bold")))
print(p_appendix)
ggsave(file.path(plots_dir, "06_appendix_species_size.jpeg"), p_appendix,
       width = 11.69, height = 6.5, dpi = 300, units = "in")


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

# significance labels
sig_stars <- function(p) {
  dplyr::case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    p < 0.1   ~ ".",
    TRUE      ~ ""
  )
}

# collect all results
results_all <- bind_rows(
  extract_fixed(model_cover,         "LMM: total cover trend",    "log(% cover)"),
  extract_fixed(model_base,          "LMM: colony growth",        "log(area)"),
  extract_fixed(model_pm_occurrence, "GLMM: mortality occurrence","P(any mortality)"),
  extract_fixed(model_pm_severity,   "LMM: mortality severity",   "log(% area lost)"),
  extract_fixed(model_fate,          "GLMM: colony fate",         "P(colony died)"),
  extract_pairs(pairs_class,         "post-hoc: species size",    "log(area)"),
  extract_pairs(pairs_fate_cl,       "post-hoc: colony fate",     "P(colony died)")
) %>%
  mutate(sig = sig_stars(p_value)) %>%
  dplyr::arrange(source, p_value)

# also add crossed model if it was fitted
if (exists("model_crossed")) {
  results_all <- bind_rows(
    results_all,
    extract_fixed(model_crossed, "LMM: growth (species x size)", "log(area)") %>%
      mutate(sig = sig_stars(p_value))
  )
}

# significant results only (p < 0.05), annotated with low-power species flag
results_sig <- results_all %>%
  filter(p_value < 0.05) %>%
  dplyr::select(source, response, type, term, estimate, se, statistic, p_value, sig) %>%
  mutate(power_note = ifelse(
    sapply(term, function(t) any(sapply(low_power_spp, function(sp) grepl(sp, t, fixed = TRUE)))),
    "low power", ""
  ))

print(results_sig)
