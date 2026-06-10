# load packages
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(lme4, lmerTest, dplyr, tidyr, ggplot2, viridisLite, patchwork, scales, pwr, here, DHARMa)

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

# detect split and merge events
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

# numeric time
data$time_days <- as.numeric(data$date - min(data$date))

# partial mortality: % area lost since last survey (only when area decreased)
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

# colony fate: died = absent before plot's last survey date
plot_max_date <- data %>%
  group_by(plot_id) %>%
  summarise(plot_max_date = max(date), .groups = "drop")

data <- data %>%
  group_by(colony_id) %>%
  mutate(colony_last_date = max(date)) %>%
  ungroup() %>%
  left_join(plot_max_date, by = "plot_id") %>%
  mutate(colony_fate = ifelse(colony_last_date < plot_max_date, "died", "alive"))

# initial size per colony: z-score of log(area) within species
initial_size <- data %>%
  group_by(colony_id) %>%
  slice_min(date, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  group_by(class) %>%
  mutate(
    log_area = log(area),
    size_z   = {
      s <- sd(log_area)
      if (n() < 2 || s == 0) rep(0, n()) else (log_area - mean(log_area)) / s
    }
  ) %>%
  ungroup() %>%
  dplyr::select(colony_id, size_z)

data <- left_join(data, initial_size, by = "colony_id")
time_scale <- 100  # used to rescale time_days for numerical stability in glmer
data       <- data %>% mutate(time_sc = time_days / time_scale)

# --- species sample size check ---
target_power <- 0.80; target_alpha <- 0.05; target_d <- 0.8; n_hard_min <- 5
n_power_min  <- ceiling(pwr::pwr.t.test(
  d = target_d, sig.level = target_alpha, power = target_power, type = "two.sample"
)$n)

species_n <- data %>%
  distinct(colony_id, class) %>%
  count(class, name = "n_colonies") %>%
  mutate(
    adequate_power = n_colonies >= n_power_min,
    low_power      = n_colonies >= 10 & n_colonies < n_power_min,
    very_low_power = n_colonies >= n_hard_min & n_colonies < 10,
    too_few        = n_colonies < n_hard_min
  )

message("Species n check  |  hard floor: ", n_hard_min,
        "  |  power threshold (d=", target_d, ", 1-β=", target_power, "): ", n_power_min)
print(as.data.frame(species_n))

excluded_spp       <- species_n %>% filter(too_few)        %>% pull(class)
very_low_power_spp <- species_n %>% filter(very_low_power) %>% pull(class)
low_power_spp      <- species_n %>% filter(low_power)      %>% pull(class)

if (length(excluded_spp) > 0) {
  message("Removing (below hard floor n=", n_hard_min, "): ", paste(excluded_spp, collapse = ", "))
  data <- data %>% filter(!class %in% excluded_spp)
}

data$class <- relevel(droplevels(factor(data$class)), ref = "Porites astreoides")

PLOT_AREA_CM2 <- 250000
multi_site    <- length(unique(data$site)) > 1

sig_stars <- function(p) {
  dplyr::case_when(
    p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", p < 0.1 ~ ".", TRUE ~ ""
  )
}

# ── design system ──────────────────────────────────────────────────────────
accent_col    <- "#898C31"
accent_bg_col <- "#E5E6D2"
sp_font       <- "Times"

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

sp_levels <- sort(unique(as.character(data$class)))
pal <- setNames(viridisLite::viridis(length(sp_levels), option = "D", begin = 0.05, end = 0.92),
                sp_levels)

n_sites  <- length(unique(data$site))
site_pal <- setNames(
  viridisLite::viridis(max(n_sites, 2), option = "D", begin = 0.15, end = 0.75)[seq_len(n_sites)],
  sort(unique(data$zone))
)

sp_codes <- data %>% distinct(class, code) %>% { setNames(.$code, .$class) }

sp_n <- data %>% distinct(colony_id, class) %>% count(class, name = "n_col") %>%
  { setNames(.$n_col, .$class) }
sp_labels_n <- setNames(
  paste0(sp_codes[names(sp_n)], "\nn = ", sp_n),
  names(sp_n)
)

fig_sq_w <- 945;  fig_sq_h <- 945
fig_ls_w <- 1890; fig_ls_h <- 945


# ═══════════════════════════════════════════════════════════════════════════
# SECTION 1 — SPECIES OVERVIEW: COVER AND ABUNDANCE
# ═══════════════════════════════════════════════════════════════════════════

cover_sp <- data %>%
  group_by(plot_id, date, time_days, class) %>%
  summarise(sp_area = sum(area, na.rm = TRUE), .groups = "drop") %>%
  mutate(pct_cover_sp = sp_area / PLOT_AREA_CM2 * 100)

cover_sp_mean <- cover_sp %>%
  group_by(date, class) %>%
  summarise(mean_pct = mean(pct_cover_sp), .groups = "drop")

# smooth each species with smooth.spline on absolute cover %, predict on weekly grid
dates_grid <- seq(min(cover_sp_mean$date), max(cover_sp_mean$date), by = "week")

cover_sp_smooth <- cover_sp_mean %>%
  group_by(class) %>%
  group_modify(~ {
    x        <- as.numeric(.x$date)
    y        <- .x$mean_pct
    xout     <- as.numeric(dates_grid)
    n_unique <- length(unique(x))
    pred <- if (n_unique >= 4) {
      fit <- smooth.spline(x, y, df = 5)
      pmax(predict(fit, xout)$y, 0)
    } else if (n_unique >= 2) {
      pmax(predict(lm(y ~ x), newdata = data.frame(x = xout)), 0)
    } else {
      rep(y[1], length(xout))
    }
    data.frame(date = dates_grid, mean_pct = pred)
  }) %>%
  ungroup()

# order species by overall mean cover ascending — smallest at bottom of stack
sp_cover_order <- cover_sp_mean %>%
  group_by(class) %>%
  summarise(overall_mean = mean(mean_pct), .groups = "drop") %>%
  arrange(desc(overall_mean)) %>%
  pull(class)

cover_sp_smooth$class <- factor(cover_sp_smooth$class, levels = sp_cover_order)

# 01a: smoothed stacked area — stack height = total cover, relative area = composition
p01a <- ggplot(cover_sp_smooth, aes(x = date, y = mean_pct, fill = class)) +
  geom_area(position = "stack", alpha = 0.88) +
  scale_fill_manual(values = pal, labels = sp_codes) +
  scale_x_date(date_breaks = "6 months", date_labels = "%b '%y") +
  scale_y_continuous(labels = \(x) paste0(x, "%"),
                     expand = expansion(mult = c(0, 0.03))) +
  labs(x = "Survey date", y = "Mean cover (% of 25 m²)", fill = NULL) +
  theme_coral() +
  theme(legend.position = "top") +
  guides(fill = guide_legend(nrow = 1))

date_origin <- min(data$date)

# 01b: colony abundance per species over time — Poisson GLMM trend per species
# per-plot counts as the modelling unit; zero-fill species absent from a surveyed
# plot-date so the trend reflects true abundance (including structural zeros), not
# abundance conditional on presence. Counts modelled with a log link → non-negative
# predictions (a Gaussian LMM on counts could predict negative abundance).
abundance_sp_plot <- data %>%
  group_by(date, time_days, time_sc, plot_id, class) %>%
  summarise(n_colonies = n_distinct(colony_id), .groups = "drop") %>%
  tidyr::complete(
    tidyr::nesting(date, time_days, time_sc, plot_id), class,
    fill = list(n_colonies = 0)
  )

model_abundance <- glmer(
  n_colonies ~ time_sc * class + (1 | plot_id),
  data = abundance_sp_plot, family = poisson(link = "log"),
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

# GLMM predictions on a dense date grid per species (response scale)
abund_grid <- expand.grid(
  class     = levels(data$class),
  time_days = seq(min(data$time_days), max(data$time_days), length.out = 60)
) %>%
  mutate(time_sc = time_days / time_scale, date = date_origin + time_days) %>%
  mutate(pred_abund = predict(model_abundance, newdata = ., re.form = NA, type = "response"))

p01b <- ggplot() +
  geom_point(data = abundance_sp_plot,
             aes(x = date, y = n_colonies, color = class),
             alpha = 0.25, size = 1.4) +
  geom_line(data = abund_grid,
            aes(x = date, y = pred_abund, color = class),
            linewidth = 0.8) +
  scale_color_manual(values = pal, labels = sp_codes) +
  scale_x_date(date_breaks = "6 months", date_labels = "%b '%y") +
  labs(x = "Survey date", y = "Mean colonies per plot (25 m²)", color = NULL) +
  theme_coral() +
  theme(legend.position = "top") +
  guides(color = guide_legend(nrow = 1))

p01_combined <- guide_area() /
  ((p01a + theme(legend.position = "none")) |
   (p01b + guides(color = guide_legend(nrow = 1)))) +
  plot_layout(guides = "collect", heights = c(0.12, 1))
ggsave(file.path(plots_dir, "01_species_overview.jpeg"), p01_combined,
       width = fig_ls_w, height = fig_sq_h, dpi = 300, units = "px")


# ═══════════════════════════════════════════════════════════════════════════
# SECTION 2 — SPECIES GROWTH TRAJECTORIES AND SIZE EFFECTS
# ═══════════════════════════════════════════════════════════════════════════

# follow-up observations only: drop each colony's baseline survey before fitting the
# size models. size_z is the z-scored log of the *first* area, so leaving the baseline
# in the response makes the colony predict its own initial value — a mechanical
# coupling that inflates the size_z effect and the size×time interaction
# (regression to the mean). Fitting on follow-ups separates predictor from response.
# NOTE: residual RTM remains because size_z is measured with error; a size-transition
# (IPM-style) model on consecutive (area_t, area_t+1) pairs would remove it fully.
data_fu <- data %>%
  group_by(colony_id) %>%
  filter(date > min(date)) %>%
  ungroup()

# full model with size_z — used for coefficient tables and size-effect plots
model_growth <- glmer(
  area ~ time_sc * class + time_sc * size_z + site +
    (1 | plot_id) + (1 | colony_id),
  data = data_fu, family = Gamma(link = "log"),
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)
summary(model_growth)

# per-plot means for trajectory visualisation (averaging unit = plot × species × date)
obs_summary_traj <- data %>%
  group_by(class, zone, plot_id, date, time_days, time_sc) %>%
  summarise(mean_area = mean(area), .groups = "drop")

# GLM on per-plot means — predictions match the arithmetic mean shown in the plot
model_growth_traj <- glm(
  mean_area ~ time_sc * class,
  data = obs_summary_traj, family = Gamma(link = "log")
)

# simplified size-only model — used for growth-by-size prediction line
# (follow-up observations only, same RTM rationale as model_growth)
model_growth_size <- lmer(
  log(area) ~ time_days * size_z + (1 | plot_id) + (1 | colony_id),
  data = data_fu
)
summary(model_growth_size)

# --- GLMM assumption checks via DHARMa simulation-based residuals ---
set.seed(42)
sim_res <- simulateResiduals(model_growth, n = 500, plot = FALSE)

# panel a: simulated residuals vs fitted (should be uniform, flat)
dharma_df <- data.frame(
  fitted   = fitted(model_growth),
  sim_resid = sim_res$scaledResiduals
)

pa <- ggplot(dharma_df, aes(x = fitted, y = sim_resid)) +
  geom_point(alpha = 0.15, size = 0.7, color = "grey40") +
  geom_hline(yintercept = c(0.25, 0.5, 0.75), linetype = "dashed",
             color = "grey70", linewidth = 0.4) +
  geom_smooth(method = "loess", se = FALSE, color = accent_col, linewidth = 0.8) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Fitted values (cm²)", y = "Simulated residual quantile",
       title = "Residuals vs Fitted (DHARMa)") +
  theme_coral()

# panel b: QQ of simulated residuals (should follow diagonal)
sim_q   <- sort(sim_res$scaledResiduals)
unif_q  <- qunif(ppoints(length(sim_q)))
pb <- ggplot(data.frame(x = unif_q, y = sim_q), aes(x = x, y = y)) +
  geom_abline(slope = 1, intercept = 0, color = accent_col, linewidth = 0.7) +
  geom_point(alpha = 0.20, size = 0.7, color = "grey40") +
  labs(x = "Expected (uniform)", y = "Observed quantile",
       title = "Uniform Q-Q — simulated residuals") +
  theme_coral()

# panel c: Pearson residuals vs fitted (correct scale for Gamma)
pearson_df <- data.frame(
  fitted        = fitted(model_growth),
  pearson_resid = residuals(model_growth, type = "pearson")
)

pc <- ggplot(pearson_df, aes(x = fitted, y = pearson_resid)) +
  geom_point(alpha = 0.15, size = 0.7, color = "grey40") +
  geom_hline(yintercept = 0, linetype = "dashed", color = accent_col, linewidth = 0.5) +
  geom_smooth(method = "loess", se = FALSE, color = accent_col, linewidth = 0.8) +
  labs(x = "Fitted values (cm²)", y = "Pearson residual",
       title = "Pearson Residuals vs Fitted") +
  theme_coral()

# panel d: random effects QQ (still normally distributed by assumption)
re_vals <- ranef(model_growth)$colony_id[[1]]
pd <- ggplot(data.frame(re = re_vals), aes(sample = re)) +
  stat_qq(alpha = 0.3, size = 0.7, color = "grey40") +
  stat_qq_line(color = accent_col, linewidth = 0.7) +
  labs(x = "Theoretical quantiles", y = "Sample quantiles",
       title = "Q-Q: random effects (colony)") +
  theme_coral()

p_assumptions <- (pa | pb) / (pc | pd) +
  plot_annotation(
    title   = "Assumption checks — Gamma GLMM: colony growth (model_growth)",
    subtitle = "Visualisation GLMs (trajectory plots) are not used for inference and are not checked here.",
    theme   = theme(
      plot.title    = element_text(family = "Helvetica", size = 9, face = "bold", color = "grey20"),
      plot.subtitle = element_text(family = "Helvetica", size = 7, color = "grey50")
    )
  )

ggsave(file.path(plots_dir, "00_assumption_checks.jpeg"), p_assumptions,
       width = fig_ls_w, height = fig_ls_h, dpi = 300, units = "px")

# 02a: per-species colony area trajectories — observed mean±SE + model prediction
traj_grid <- expand.grid(
  class     = levels(data$class),
  time_days = seq(min(data$time_days), max(data$time_days), length.out = 60)
) %>%
  mutate(time_sc = time_days / time_scale, date = date_origin + time_days) %>%
  mutate(pred_area = predict(model_growth_traj, newdata = ., type = "response"))

# per-plot mean ± SE per species per date, colored by site
obs_summary <- data %>%
  group_by(class, zone, plot_id, date) %>%
  summarise(
    mean_area = mean(area),
    se_area   = sd(area) / sqrt(n()),
    .groups   = "drop"
  )

p02a <- ggplot() +
  geom_line(data = traj_grid,
            aes(x = date, y = pred_area, color = class), linewidth = 1.0) +
  geom_pointrange(data = obs_summary,
                  aes(x = date, y = mean_area,
                      ymin = mean_area - se_area, ymax = mean_area + se_area,
                      color = zone),
                  size = 0.2, linewidth = 0.4, alpha = 0.35) +
  scale_color_manual(values = c(pal, site_pal), guide = "none") +
  scale_x_date(date_breaks = "1 year", date_labels = "'%y") +
  facet_wrap(~ class, scales = "free_y", nrow = 2,
             labeller = as_labeller(sp_labels_n)) +
  labs(x = "Survey date", y = "Colony area (cm²)") +
  theme_coral() +
  theme(strip.text = element_text(family = "Helvetica", face = "plain", size = 7))

ggsave(file.path(plots_dir, "02_growth_trajectories.jpeg"), p02a,
       width = fig_ls_w, height = fig_ls_h, dpi = 300, units = "px")

# partial mortality data — used by 02d, 02c, and 03d
dead_terminal <- data %>%
  filter(colony_fate == "died", !is.na(size_z)) %>%
  group_by(colony_id) %>%
  slice_max(date, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(partial_mortality = 100, event_type = "complete death")

pm_data <- data %>%
  filter(!is.na(size_z)) %>%
  mutate(event_type = ifelse(!is.na(partial_mortality), "partial mortality", "none")) %>%
  bind_rows(dead_terminal)

pm_events  <- pm_data %>% filter(!is.na(partial_mortality), partial_mortality > 0)
pm_partial <- pm_events %>% filter(event_type != "complete death")

# 02d: partial mortality trajectories over time per species
pm_obs_summary <- pm_partial %>%
  group_by(class, zone, plot_id, date, time_days, time_sc) %>%
  summarise(
    mean_pm = mean(partial_mortality),
    se_pm   = sd(partial_mortality) / sqrt(n()),
    .groups = "drop"
  )

# GLM on per-plot means so predictions match the arithmetic mean shown in points
model_pm_traj <- glm(
  mean_pm ~ time_sc + class,
  data = pm_obs_summary, family = Gamma(link = "log")
)

pm_traj_grid <- expand.grid(
  class     = levels(data$class),
  time_days = seq(min(pm_partial$time_days, na.rm = TRUE),
                  max(pm_partial$time_days, na.rm = TRUE), length.out = 60)
) %>%
  mutate(time_sc = time_days / time_scale, date = date_origin + time_days) %>%
  mutate(pred_pm = predict(model_pm_traj, newdata = ., type = "response"))

p02d <- ggplot() +
  geom_line(data = pm_traj_grid,
            aes(x = date, y = pred_pm, color = class), linewidth = 1.0) +
  geom_pointrange(data = pm_obs_summary,
                  aes(x = date, y = mean_pm,
                      ymin = mean_pm - se_pm, ymax = mean_pm + se_pm,
                      color = zone),
                  size = 0.2, linewidth = 0.4, alpha = 0.35) +
  scale_color_manual(values = c(pal, site_pal), guide = "none") +
  scale_x_date(date_breaks = "1 year", date_labels = "'%y") +
  scale_y_continuous(labels = \(x) paste0(x, "%"),
                     expand = expansion(mult = c(0.02, 0.05))) +
  coord_cartesian(ylim = c(0, NA)) +
  facet_wrap(~ class, nrow = 2, labeller = as_labeller(sp_labels_n)) +
  labs(x = "Survey date", y = "Partial mortality (%)") +
  theme_coral() +
  theme(strip.text = element_text(family = "Helvetica", face = "plain", size = 7))

ggsave(file.path(plots_dir, "03_mortality_trajectories.jpeg"), p02d,
       width = fig_ls_w, height = fig_ls_h, dpi = 300, units = "px")

# 02b: growth rate by initial size — empirical points + model prediction line
colony_growth <- data %>%
  filter(!is.na(size_z)) %>%
  arrange(colony_id, date) %>%
  group_by(colony_id, class, zone, size_z) %>%
  summarise(
    first_area = area[which.min(date)],
    last_area  = area[which.max(date)],
    total_days = as.numeric(max(date) - min(date)),
    n_obs      = n(),
    .groups    = "drop"
  ) %>%
  filter(n_obs >= 2, total_days > 0, first_area > 0) %>%
  mutate(pct_growth_yr = (exp((log(last_area) - log(first_area)) / total_days * 365) - 1) * 100)

time_coef      <- fixef(model_growth_size)["time_days"]
time_size_coef <- fixef(model_growth_size)["time_days:size_z"]
sz_range       <- range(colony_growth$size_z, na.rm = TRUE)

growth_pred <- data.frame(size_z = seq(sz_range[1], sz_range[2], length.out = 60)) %>%
  mutate(pred_pct = (exp((time_coef + time_size_coef * size_z) * 365) - 1) * 100)

p02b <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = accent_col, linewidth = 0.5) +
  geom_point(data = colony_growth,
             aes(x = size_z, y = pct_growth_yr, color = class),
             alpha = 0.40, size = 1.6) +
  geom_line(data = growth_pred,
            aes(x = size_z, y = pred_pct),
            color = accent_col, linewidth = 1.1) +
  scale_color_manual(values = pal, name = NULL, labels = sp_codes) +
  labs(x = "Initial size z-score (within species)", y = "Growth rate (% per year)") +
  theme_coral() +
  theme(legend.position      = c(0.98, 0.98),
        legend.justification = c("right", "top"),
        legend.background    = element_rect(fill = alpha("white", 0.7), color = NA),
        legend.text          = element_text(size = 7))

ggsave(file.path(plots_dir, "04_growth_by_size.jpeg"), p02b,
       width = fig_sq_w, height = fig_sq_h, dpi = 300, units = "px")

# 02c: mortality severity by initial size
# fit on partial-mortality events only — complete deaths (forced to 100%) are a
# different process and a ceiling spike, so including them would bias the severity
# slope. The scatter below still shows both event types for context.
model_pm_severity <- lmer(
  log(partial_mortality) ~ size_z + class + (1 | colony_id) + (1 | plot_id),
  data = pm_partial
)
summary(model_pm_severity)

sev_grid <- data.frame(
  size_z = seq(min(pm_partial$size_z, na.rm = TRUE),
               max(pm_partial$size_z, na.rm = TRUE), length.out = 60),
  class  = "Porites astreoides"
) %>% mutate(pred_pct = exp(predict(model_pm_severity, newdata = ., re.form = NA)))

p02c <- ggplot() +
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

ggsave(file.path(plots_dir, "05_mortality_by_size.jpeg"), p02c,
       width = fig_sq_w, height = fig_sq_h, dpi = 300, units = "px")



# ═══════════════════════════════════════════════════════════════════════════
# SECTION 3 — SITE COMPARISON
# only produced when data contains more than one site
# ═══════════════════════════════════════════════════════════════════════════

if (multi_site) {

  zones_sorted <- sort(unique(data$zone))  # time_sc already set globally (see time_scale)

  # per-plot cover for site-level LMM
  cover_site_plots <- data %>%
    group_by(site, zone, plot_id, date, time_days) %>%
    summarise(total_area = sum(area, na.rm = TRUE), .groups = "drop") %>%
    mutate(pct_cover = total_area / PLOT_AREA_CM2 * 100,
           time_sc   = time_days / time_scale)

  # Gamma GLMM: site cover trend over time
  model_cover_site <- glmer(
    pct_cover ~ time_sc * zone + (1 | plot_id),
    data = cover_site_plots, family = Gamma(link = "log"),
    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
  )
  summary(model_cover_site)

  # Gamma GLMM predictions for 03a on a dense date grid per zone
  cover_site_grid <- expand.grid(
    zone      = zones_sorted,
    time_days = seq(min(cover_site_plots$time_days), max(cover_site_plots$time_days), length.out = 60)
  ) %>%
    mutate(time_sc = time_days / time_scale, date = min(data$date) + time_days) %>%
    mutate(pred_cover = predict(model_cover_site, newdata = ., re.form = NA, type = "response"))

  # 03a: site cover — per-plot points (background) + LMM trend lines
  p03a <- ggplot() +
    geom_point(data = cover_site_plots,
               aes(x = date, y = pct_cover, color = zone),
               alpha = 0.25, size = 1.4) +
    geom_line(data = cover_site_grid,
              aes(x = date, y = pred_cover, color = zone),
              linewidth = 1.0) +
    scale_color_manual(values = site_pal) +
    scale_x_date(date_breaks = "6 months", date_labels = "%b '%y") +
    scale_y_continuous(labels = \(x) paste0(round(x, 1), "%")) +
    labs(x = "Survey date", y = "Total cover (% of 25 m²)", color = NULL) +
    theme_coral() +
    theme(legend.position      = c(0.02, 0.98),
          legend.justification = c("left", "top"),
          legend.background    = element_rect(fill = alpha("white", 0.7), color = NA))


  # 03b: colony abundance per site — per-plot counts (background) + LMM trend lines
  abundance_site_plot <- data %>%
    group_by(date, time_days, time_sc, plot_id, zone) %>%
    summarise(n_colonies = n_distinct(colony_id), .groups = "drop")

  # Poisson GLMM (log link) — counts, non-negative predictions
  model_abundance_site <- glmer(
    n_colonies ~ time_sc * zone + (1 | plot_id),
    data = abundance_site_plot, family = poisson(link = "log"),
    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
  )

  abund_site_grid <- expand.grid(
    zone      = zones_sorted,
    time_days = seq(min(abundance_site_plot$time_days), max(abundance_site_plot$time_days), length.out = 60)
  ) %>%
    mutate(time_sc = time_days / time_scale, date = min(data$date) + time_days) %>%
    mutate(pred_abund = predict(model_abundance_site, newdata = ., re.form = NA, type = "response"))

  p03b <- ggplot() +
    geom_point(data = abundance_site_plot,
               aes(x = date, y = n_colonies, color = zone),
               alpha = 0.25, size = 1.4) +
    geom_line(data = abund_site_grid,
              aes(x = date, y = pred_abund, color = zone),
              linewidth = 1.0) +
    scale_color_manual(values = site_pal) +
    scale_x_date(date_breaks = "6 months", date_labels = "%b '%y") +
    labs(x = "Survey date", y = "Colonies per plot (25 m²)", color = NULL) +
    theme_coral() +
    theme(legend.position      = c(0.02, 0.98),
          legend.justification = c("left", "top"),
          legend.background    = element_rect(fill = alpha("white", 0.7), color = NA))

  # LMM: growth by site — site shifts intercept and overall growth rate (time_days:zone)
  # size_z interaction still assumed equal across sites
  model_growth_site <- glmer(
    area ~ time_sc * zone + (1 | colony_id),
    data = data, family = Gamma(link = "log"),
    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
  )
  summary(model_growth_site)


  # 03c: colony area trajectories over time by site — per-plot mean + model prediction
  obs_site <- data %>%
    group_by(zone, plot_id, date, time_days, time_sc) %>%
    summarise(
      mean_area = mean(area),
      se_area   = sd(area) / sqrt(n()),
      .groups   = "drop"
    )

  # GLM on per-plot means so predictions match the arithmetic mean shown in points
  model_growth_site_traj <- glm(
    mean_area ~ time_sc * zone,
    data = obs_site, family = Gamma(link = "log")
  )

  traj_site_grid <- expand.grid(
    zone      = zones_sorted,
    time_days = seq(min(data$time_days), max(data$time_days), length.out = 60)
  ) %>%
    mutate(time_sc = time_days / time_scale, date = date_origin + time_days) %>%
    mutate(pred_area = predict(model_growth_site_traj, newdata = ., type = "response"))

  p04_site <- (p03a | p03b) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  ggsave(file.path(plots_dir, "06_site_cover_abundance.jpeg"), p04_site,
         width = fig_ls_w, height = fig_ls_h, dpi = 300, units = "px")

  p03c_single <- ggplot() +
    geom_line(data = traj_site_grid,
              aes(x = date, y = pred_area, color = zone), linewidth = 1.0) +
    geom_pointrange(data = obs_site,
                    aes(x = date, y = mean_area,
                        ymin = mean_area - se_area, ymax = mean_area + se_area,
                        color = zone),
                    size = 0.2, linewidth = 0.4, alpha = 0.35) +
    scale_color_manual(values = site_pal, labels = zones_sorted) +
    scale_x_date(date_breaks = "1 year", date_labels = "'%y") +
    labs(x = "Survey date", y = "Mean colony area (cm²)", color = NULL) +
    theme_coral()

  # 03d: partial mortality % over time by site
  # per-plot mean ± SE of partial mortality per date (complete deaths excluded)
  mort_obs_site <- pm_partial %>%
    group_by(zone, plot_id, date, time_days, time_sc) %>%
    summarise(
      mean_pm = mean(partial_mortality),
      se_pm   = sd(partial_mortality) / sqrt(n()),
      .groups = "drop"
    )

  # GLM on per-plot means so predictions match the arithmetic mean shown in points
  model_mortality_site_traj <- glm(
    mean_pm ~ time_sc * zone,
    data = mort_obs_site, family = Gamma(link = "log")
  )

  mort_site_grid <- expand.grid(
    zone      = zones_sorted,
    time_days = seq(min(pm_partial$time_days, na.rm = TRUE),
                    max(pm_partial$time_days, na.rm = TRUE), length.out = 60)
  ) %>%
    mutate(time_sc = time_days / time_scale, date = min(data$date) + time_days) %>%
    mutate(pred_pct = predict(model_mortality_site_traj, newdata = ., type = "response"))

  p03d_single <- ggplot() +
    geom_line(data = mort_site_grid,
              aes(x = date, y = pred_pct, color = zone), linewidth = 1.0) +
    geom_pointrange(data = mort_obs_site,
                    aes(x = date, y = mean_pm,
                        ymin = mean_pm - se_pm, ymax = mean_pm + se_pm,
                        color = zone),
                    size = 0.2, linewidth = 0.4, alpha = 0.35) +
    scale_color_manual(values = site_pal, labels = zones_sorted) +
    scale_x_date(date_breaks = "1 year", date_labels = "'%y") +
    scale_y_continuous(labels = \(x) paste0(x, "%"),
                       expand = expansion(mult = c(0.02, 0.05))) +
    coord_cartesian(ylim = c(0, NA)) +
    labs(x = "Survey date", y = "Partial mortality (%)", color = NULL) +
    theme_coral()

  p07_site_traj <- (p03c_single | p03d_single) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  ggsave(file.path(plots_dir, "07_site_trajectories.jpeg"), p07_site_traj,
         width = fig_ls_w, height = fig_ls_h, dpi = 300, units = "px")

  # 03e: species composition by site — proportion of colonies per species
  sp_site <- data %>%
    distinct(colony_id, zone, class) %>%
    count(zone, class, name = "n_colonies") %>%
    group_by(zone) %>%
    mutate(prop = n_colonies / sum(n_colonies)) %>%
    ungroup()

  p03e <- ggplot(sp_site, aes(x = zone, y = prop, fill = class)) +
    geom_col(width = 0.5, alpha = 0.88) +
    scale_fill_manual(values = pal, labels = sp_codes) +
    scale_y_continuous(labels = scales::percent_format()) +
    coord_flip() +
    labs(x = NULL, y = "Proportion of colonies", fill = NULL) +
    theme_coral() +
    theme(legend.text     = element_text(size = 7),
          legend.position = "bottom") +
    guides(fill = guide_legend(nrow = 2))

  ggsave(file.path(plots_dir, "08_species_by_site.jpeg"), p03e,
         width = fig_sq_w, height = fig_sq_h, dpi = 300, units = "px")

  # chi-square: species composition independence by site
  sp_comp_mat <- sp_site %>%
    tidyr::pivot_wider(id_cols = zone, names_from = class,
                       values_from = n_colonies, values_fill = 0L) %>%
    tibble::column_to_rownames("zone") %>% as.matrix()
  # several species have <5 colonies, so the asymptotic chi-square approximation is
  # unreliable; use a Monte-Carlo simulated p-value (no df reported in that mode)
  set.seed(42)
  chisq_sp_comp <- chisq.test(sp_comp_mat, simulate.p.value = TRUE, B = 1e5)
  message("Species composition chi-square (Monte-Carlo): X2=", round(chisq_sp_comp$statistic, 2),
          "  p=", round(chisq_sp_comp$p.value, 4))

  # Wilcoxon: initial colony size by site
  init_size_site <- data %>%
    group_by(colony_id, zone) %>%
    slice_min(date, n = 1, with_ties = FALSE) %>%
    ungroup()
  wilcox_size_site <- wilcox.test(log(area) ~ zone, data = init_size_site)
  message("Initial size Wilcoxon: W=", wilcox_size_site$statistic,
          "  p=", round(wilcox_size_site$p.value, 4))

  # colony fate model (site effect on survival probability)
  fate_data <- data %>%
    filter(!is.na(size_z)) %>%
    distinct(colony_id, class, site, zone, plot_id, size_z, colony_fate) %>%
    mutate(fate_binary = as.integer(colony_fate == "died"))

  model_fate <- glm(
    fate_binary ~ size_z + class + site,
    data = fate_data, family = binomial
  )
  summary(model_fate)

} else {
  message("Single site — site comparison skipped")
}


# ═══════════════════════════════════════════════════════════════════════════
# RESULTS SUMMARY
# ═══════════════════════════════════════════════════════════════════════════

extract_fixed <- function(model, source, response) {
  s    <- coef(summary(model))
  pcol <- grep("Pr\\(", colnames(s), value = TRUE)[1]
  tcol <- grep("value$", colnames(s), value = TRUE)
  tcol <- tcol[grepl("^[tz]", tcol)][1]
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

results_all <- bind_rows(
  extract_fixed(model_growth,      "GLMM: colony growth (Gamma)",      "area [log link]"),
  extract_fixed(model_growth_size, "LMM: growth by size",              "log(area)"),
  extract_fixed(model_pm_severity, "LMM: mortality severity",          "log(% area lost)")
)

if (multi_site) {
  results_all <- bind_rows(
    results_all,
    extract_fixed(model_cover_site,  "GLMM: site cover (Gamma)",         "pct_cover [log link]"),
    extract_fixed(model_growth_site, "GLMM: growth by site (Gamma)",     "area [log link]"),
    extract_fixed(model_fate,        "GLMM: colony fate",   "P(colony died)"),
    data.frame(
      source    = "chi-square: species composition by site",
      response  = "colony counts",
      type      = "chi-square",
      term      = "all species x site",
      estimate  = NA_real_, se = NA_real_,
      statistic = round(as.numeric(chisq_sp_comp$statistic), 3),
      p_value   = round(chisq_sp_comp$p.value, 4)
    ),
    data.frame(
      source    = "Wilcoxon: initial size by site",
      response  = "log(colony area)",
      type      = "Wilcoxon rank-sum",
      term      = "zone",
      estimate  = NA_real_, se = NA_real_,
      statistic = round(as.numeric(wilcox_size_site$statistic), 3),
      p_value   = round(wilcox_size_site$p.value, 4)
    )
  )
}

results_all <- results_all %>%
  mutate(sig = sig_stars(p_value)) %>%
  arrange(source, p_value)

results_sig <- results_all %>%
  filter(p_value < 0.05) %>%
  mutate(power_note = dplyr::case_when(
    sapply(term, function(t) any(sapply(very_low_power_spp, function(sp) grepl(sp, t, fixed = TRUE)))) ~ "very low power (n 5-9)",
    sapply(term, function(t) any(sapply(low_power_spp,      function(sp) grepl(sp, t, fixed = TRUE)))) ~ "low power (n 10-25)",
    TRUE ~ ""
  ))

print(results_sig)
