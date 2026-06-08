# load packages
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(lme4, lmerTest, dplyr, tidyr, ggplot2, viridisLite, patchwork, scales, pwr, here)

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

# 01a: stacked species cover over time
p01a <- ggplot(cover_sp_mean, aes(x = date, y = mean_pct, fill = class)) +
  geom_area(alpha = 0.85, position = "stack") +
  scale_fill_manual(values = pal, labels = sp_codes) +
  scale_x_date(date_breaks = "6 months", date_labels = "%b '%y") +
  scale_y_continuous(labels = \(x) paste0(x, "%"),
                     expand = expansion(mult = c(0, 0.03))) +
  labs(x = "Survey date", y = "Mean cover (% of 25 m²)", fill = NULL) +
  theme_coral() +
  theme(legend.position = "top") +
  guides(fill = guide_legend(nrow = 1))

ggsave(file.path(plots_dir, "01a_species_cover.jpeg"), p01a,
       width = fig_ls_w, height = fig_ls_h, dpi = 300, units = "px")

# 01b: colony abundance per species over time — mean colonies per 25 m² plot
abundance_sp <- data %>%
  group_by(date, class) %>%
  summarise(
    n_colonies    = n_distinct(colony_id),
    n_plots       = n_distinct(plot_id),
    mean_per_plot = n_colonies / n_plots,
    .groups       = "drop"
  )

p01b <- ggplot(abundance_sp, aes(x = date, y = mean_per_plot, color = class)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.8, alpha = 0.8) +
  scale_color_manual(values = pal, labels = sp_codes) +
  scale_x_date(date_breaks = "6 months", date_labels = "%b '%y") +
  labs(x = "Survey date", y = "Mean colonies per plot (25 m²)", color = NULL) +
  theme_coral() +
  theme(legend.position = "top") +
  guides(color = guide_legend(nrow = 1))

ggsave(file.path(plots_dir, "01b_species_abundance.jpeg"), p01b,
       width = fig_ls_w, height = fig_ls_h, dpi = 300, units = "px")


# ═══════════════════════════════════════════════════════════════════════════
# SECTION 2 — SPECIES GROWTH TRAJECTORIES AND SIZE EFFECTS
# ═══════════════════════════════════════════════════════════════════════════

# model for species trajectory predictions (includes species to allow per-species lines)
model_growth <- lmer(
  log(area) ~ time_days * class + time_days * size_z + site +
    (1 | plot_id) + (1 | colony_id),
  data = data
)
summary(model_growth)

# simplified size-only model — used for growth-by-size prediction line
model_growth_size <- lmer(
  log(area) ~ time_days * size_z + (1 | plot_id) + (1 | colony_id),
  data = data
)
summary(model_growth_size)

# --- LMM assumption checks ---
resid_df <- data.frame(
  fitted         = fitted(model_growth),
  residuals      = resid(model_growth),
  sqrt_abs_resid = sqrt(abs(resid(model_growth)))
)

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

re_vals <- ranef(model_growth)$colony_id[, 1]
pd <- ggplot(data.frame(re = re_vals), aes(sample = re)) +
  stat_qq(alpha = 0.3, size = 0.7, color = "grey40") +
  stat_qq_line(color = accent_col, linewidth = 0.7) +
  labs(x = "Theoretical quantiles", y = "Sample quantiles",
       title = "Q-Q: random effects (colony)") +
  theme_coral()

p_assumptions <- (pa | pb) / (pc | pd)
ggsave(file.path(plots_dir, "00_lmm_assumption_checks.jpeg"), p_assumptions,
       width = fig_ls_w, height = fig_ls_w, dpi = 300, units = "px")

# 02a: per-species colony area trajectories — observed mean±SE + model prediction
date_origin <- min(data$date)
ref_site    <- levels(factor(data$site))[1]

traj_grid <- expand.grid(
  class     = levels(data$class),
  time_days = seq(min(data$time_days), max(data$time_days), length.out = 60),
  size_z    = 0,
  site      = ref_site
) %>%
  mutate(
    date      = date_origin + time_days,
    pred_area = exp(predict(model_growth, newdata = ., re.form = NA))
  )

obs_summary <- data %>%
  group_by(class, date) %>%
  summarise(mean_area = mean(area), se_area = sd(area) / sqrt(n()), .groups = "drop")

p02a <- ggplot() +
  geom_line(data = traj_grid,
            aes(x = date, y = pred_area, color = class), linewidth = 1.0) +
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

ggsave(file.path(plots_dir, "02a_growth_trajectories.jpeg"), p02a,
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

ggsave(file.path(plots_dir, "02b_growth_by_size.jpeg"), p02b,
       width = fig_sq_w, height = fig_sq_h, dpi = 300, units = "px")

# 02c: mortality severity by initial size
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

pm_events <- pm_data %>% filter(!is.na(partial_mortality), partial_mortality > 0)

model_pm_severity <- lmer(
  log(partial_mortality) ~ size_z + class + (1 | colony_id) + (1 | plot_id),
  data = pm_events
)
summary(model_pm_severity)

sev_grid <- data.frame(
  size_z = seq(min(pm_events$size_z, na.rm = TRUE),
               max(pm_events$size_z, na.rm = TRUE), length.out = 60),
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

ggsave(file.path(plots_dir, "02c_mortality_by_size.jpeg"), p02c,
       width = fig_sq_w, height = fig_sq_h, dpi = 300, units = "px")


# ═══════════════════════════════════════════════════════════════════════════
# SECTION 3 — SITE COMPARISON
# only produced when data contains more than one site
# ═══════════════════════════════════════════════════════════════════════════

if (multi_site) {

  # per-plot cover for site-level LMM
  cover_site_plots <- data %>%
    group_by(site, zone, plot_id, date, time_days) %>%
    summarise(total_area = sum(area, na.rm = TRUE), .groups = "drop") %>%
    mutate(pct_cover = total_area / PLOT_AREA_CM2 * 100)

  # LMM: site cover trend over time
  model_cover_site <- lmer(
    log(pct_cover) ~ time_days * zone + (1 | plot_id),
    data = cover_site_plots
  )
  summary(model_cover_site)

  # 03a: site cover trajectories
  p03a <- ggplot(cover_site_plots, aes(x = date, y = pct_cover, color = zone, fill = zone)) +
    geom_point(alpha = 0.55, size = 2.0) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, alpha = 0.18, linewidth = 1.0) +
    scale_color_manual(values = site_pal) +
    scale_fill_manual(values  = site_pal, guide = "none") +
    scale_x_date(date_breaks = "6 months", date_labels = "%b '%y") +
    scale_y_continuous(labels = \(x) paste0(round(x, 1), "%")) +
    labs(x = "Survey date", y = "Total cover (% of 25 m²)", color = NULL) +
    theme_coral() +
    theme(legend.position      = c(0.02, 0.98),
          legend.justification = c("left", "top"),
          legend.background    = element_rect(fill = alpha("white", 0.7), color = NA))

  ggsave(file.path(plots_dir, "03a_site_cover.jpeg"), p03a,
         width = fig_ls_w, height = fig_ls_h, dpi = 300, units = "px")

  # 03b: colony abundance per site over time — mean colonies per 25 m² plot
  abundance_site <- data %>%
    group_by(date, zone) %>%
    summarise(
      n_colonies    = n_distinct(colony_id),
      n_plots       = n_distinct(plot_id),
      mean_per_plot = n_colonies / n_plots,
      .groups       = "drop"
    )

  p03b <- ggplot(abundance_site, aes(x = date, y = mean_per_plot, color = zone, fill = zone)) +
    geom_point(alpha = 0.65, size = 2.0) +
    geom_smooth(method = "loess", formula = y ~ x, se = TRUE, alpha = 0.18, linewidth = 1.0) +
    scale_color_manual(values = site_pal) +
    scale_fill_manual(values  = site_pal, guide = "none") +
    scale_x_date(date_breaks = "6 months", date_labels = "%b '%y") +
    labs(x = "Survey date", y = "Mean colonies per plot (25 m²)", color = NULL) +
    theme_coral() +
    theme(legend.position      = c(0.02, 0.98),
          legend.justification = c("left", "top"),
          legend.background    = element_rect(fill = alpha("white", 0.7), color = NA))

  ggsave(file.path(plots_dir, "03b_site_abundance.jpeg"), p03b,
         width = fig_ls_w, height = fig_ls_h, dpi = 300, units = "px")

  # LMM: growth by site — site shifts intercept and overall growth rate (time_days:zone)
  # size_z interaction still assumed equal across sites
  model_growth_site <- lmer(
    log(area) ~ time_days * size_z + time_days * zone + (1 | plot_id) + (1 | colony_id),
    data = data
  )
  summary(model_growth_site)

  # derive model prediction lines per site
  coefs_gs     <- fixef(model_growth_site)
  zones_sorted <- sort(unique(data$zone))

  growth_pred_site <- bind_rows(lapply(zones_sorted, function(z) {
    key        <- paste0("time_days:zone", z)
    zone_extra <- if (key %in% names(coefs_gs)) coefs_gs[key] else 0
    data.frame(
      zone   = z,
      size_z = seq(sz_range[1], sz_range[2], length.out = 60)
    ) %>%
      mutate(pred_pct = (exp((coefs_gs["time_days"] + zone_extra +
                                coefs_gs["time_days:size_z"] * size_z) * 365) - 1) * 100)
  }))

  # 03c: colony area trajectories over time by site — observed mean±SE + model prediction
  obs_site <- data %>%
    group_by(zone, date) %>%
    summarise(mean_area = mean(area), se_area = sd(area) / sqrt(n()), .groups = "drop")

  traj_site_grid <- expand.grid(
    zone      = zones_sorted,
    time_days = seq(min(data$time_days), max(data$time_days), length.out = 60),
    size_z    = 0
  ) %>%
    mutate(
      date      = date_origin + time_days,
      pred_area = exp(predict(model_growth_site, newdata = ., re.form = NA))
    )

  p03c <- ggplot() +
    geom_line(data = traj_site_grid,
              aes(x = date, y = pred_area, color = zone), linewidth = 1.0) +
    geom_pointrange(data = obs_site,
                    aes(x = date, y = mean_area,
                        ymin = mean_area - se_area, ymax = mean_area + se_area,
                        color = zone),
                    size = 0.3, alpha = 0.7) +
    scale_color_manual(values = site_pal, name = NULL) +
    scale_x_date(date_breaks = "1 year", date_labels = "'%y") +
    labs(x = "Survey date", y = "Mean colony area (cm²)") +
    theme_coral() +
    theme(legend.position      = c(0.02, 0.98),
          legend.justification = c("left", "top"),
          legend.background    = element_rect(fill = alpha("white", 0.7), color = NA))

  ggsave(file.path(plots_dir, "03c_growth_by_site.jpeg"), p03c,
         width = fig_ls_w, height = fig_ls_h, dpi = 300, units = "px")

  # 03d: mortality rate over time by site
  # proportion of colonies experiencing any mortality event per survey interval
  mortality_site_time <- data %>%
    arrange(colony_id, date) %>%
    group_by(colony_id, zone) %>%
    mutate(had_mortality = as.integer(!is.na(partial_mortality))) %>%
    ungroup() %>%
    group_by(zone, date) %>%
    summarise(
      prop_mortality = mean(had_mortality, na.rm = TRUE),
      n              = n(),
      se             = sqrt(prop_mortality * (1 - prop_mortality) / n),
      .groups        = "drop"
    )

  p03d <- ggplot(mortality_site_time,
                 aes(x = date, y = prop_mortality, color = zone, fill = zone)) +
    geom_ribbon(aes(ymin = pmax(prop_mortality - se, 0),
                    ymax = pmin(prop_mortality + se, 1)),
                alpha = 0.18, color = NA) +
    geom_line(linewidth = 1.0) +
    geom_point(size = 2.0, alpha = 0.8) +
    scale_color_manual(values = site_pal) +
    scale_fill_manual(values  = site_pal, guide = "none") +
    scale_x_date(date_breaks = "6 months", date_labels = "%b '%y") +
    scale_y_continuous(labels = scales::percent_format(), limits = c(0, NA),
                       expand = expansion(mult = c(0, 0.05))) +
    labs(x = "Survey date", y = "Proportion of colonies with mortality", color = NULL) +
    theme_coral() +
    theme(legend.position      = c(0.02, 0.98),
          legend.justification = c("left", "top"),
          legend.background    = element_rect(fill = alpha("white", 0.7), color = NA))

  ggsave(file.path(plots_dir, "03d_mortality_by_site.jpeg"), p03d,
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

  ggsave(file.path(plots_dir, "03e_species_by_site.jpeg"), p03e,
         width = fig_sq_w, height = fig_sq_h, dpi = 300, units = "px")

  # chi-square: species composition independence by site
  sp_comp_mat <- sp_site %>%
    tidyr::pivot_wider(id_cols = zone, names_from = class,
                       values_from = n_colonies, values_fill = 0L) %>%
    tibble::column_to_rownames("zone") %>% as.matrix()
  chisq_sp_comp <- chisq.test(sp_comp_mat)
  message("Species composition chi-square: X2=", round(chisq_sp_comp$statistic, 2),
          "  df=", chisq_sp_comp$parameter,
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
  extract_fixed(model_growth,      "LMM: colony growth",        "log(area)"),
  extract_fixed(model_growth_size, "LMM: growth (size only)",   "log(area)"),
  extract_fixed(model_pm_severity, "LMM: mortality severity",   "log(% area lost)")
)

if (multi_site) {
  results_all <- bind_rows(
    results_all,
    extract_fixed(model_cover_site,  "LMM: site cover",  "log(% cover)"),
    extract_fixed(model_growth_site, "LMM: growth by site", "log(area)"),
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
