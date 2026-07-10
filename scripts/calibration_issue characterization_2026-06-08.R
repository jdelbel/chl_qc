# ─────────────────────────────────────────────────────────────────────────────
# Back-calibration % difference panel
# Reproduces: figures/backcali_perc_diff_panel_test.png
#
# Estimates/visualizes error introduced by variability in calibration slopes
# and acid ratios (fm), assuming instrument stability and a linear backward
# drift from a benchmark (DFO 2018-05-04) calibration.
#
# Instruments: 0982 (left column), 1154 (right column).
# Rows: calibration slope | acid ratio (fm) | % difference in chl.
# ─────────────────────────────────────────────────────────────────────────────

library(tidyverse)
library(readxl)
library(here)
library(patchwork)
library(hakaiApi)

# ── Data ─────────────────────────────────────────────────────────────────────
# Only the Hakai chl sheet is needed for this figure.
chl <- read_xlsx(here("files", "hakai_chl_all.xlsx"), sheet = "Hakai Data")

# ── Figure 2 data: calibration summary ────────────────────────────────────────
cali <- read_xlsx(here("files", "cali_summary_2022-04-06.xlsx")) %>%
  mutate_at(vars(r2), funs(round(., 4)))   # round r2 for plotting labels

#Download data using the API

# Initialize the client
client <- Client$new()

#Setting up Query
endpoint <- "/eims/views/output/chlorophyll"
endpoint_hplc <- "/eims/views/output/hplc"
filter <- "site_id=QU39"
chl_url <- paste0("https://portal.hakai.org/api", endpoint,"?limit=-1&", filter)
hplc_url <- paste0("https://portal.hakai.org/api", endpoint_hplc,"?limit=-1&", filter)
data_chl <- client$get(chl_url)
data_hplc <- client$get(hplc_url)

# ── Calvert acid-contamination analysis (Issue 3) ─────────────────────────────
filter_cal <- "work_area=CALVERT"
chl_cal_url <- paste0("https://portal.hakai.org/api", endpoint, "?limit=-1&", filter_cal)
data_chl_cal <- client$get(chl_cal_url)

# ── Calvert HPLC pull (Issue 3 corroboration) ─────────────────────────────────
hplc_cal_url <- paste0("https://portal.hakai.org/api", endpoint_hplc,
                       "?limit=-1&", filter_cal)
data_hplc_cal <- client$get(hplc_cal_url)

# ── Instrument 1154: back-calibrate to DFO 2018-05-04, with linear drift ──────
# DFO benchmark slope/fm held constant; a graduated backward drift
# (1.5% per step) is applied to model slope error over successive calibrations.
chl_corrected_1154 <- chl %>%
  filter(filter_type == "Bulk GF/F" & flurometer_serial_no == "720001154" &
           !is.na(calibration) & !is.na(before_acid) & !is.na(after_acid)) %>%
  select(date, volume, acetone_volume_ml, flurometer_serial_no:after_acid, chla) %>%
  mutate(chl_dfo = 0.0005026532 *
           (1.717642 / (1.717642 - 1)) *
           (before_acid - after_acid) *
           (acetone_volume_ml / volume),
         rel_diff = ((chl_dfo - chla) / chla) * 100) %>%
  mutate_at(vars(rel_diff), funs(round(., 2)))

chl_corrected_1154 <- chl_corrected_1154 %>%
  distinct(calibration, .keep_all = TRUE) %>%
  filter(!as.character(calibration) == "2014-11-21" &
           !as.character(calibration) == "2019-05-09") %>%
  mutate(slope_dfo = 0.0005026532,
         fm_dfo = 1.717642) %>%
  arrange(desc(calibration)) %>%
  mutate(perc_slope = seq(0.015, 0.075, by = 0.015),
         slope_dfo_corr = slope_dfo - (slope_dfo * perc_slope),
         chl_dfo_new = slope_dfo_corr *
           (1.717642 / (1.717642 - 1)) *
           (before_acid - after_acid) *
           (acetone_volume_ml / volume),
         rel_diff_corr = ((chl_dfo_new - chla) / chla) * 100) %>%
  mutate_at(vars(rel_diff_corr), funs(round(., 2))) %>%
  arrange(calibration)

chl_corr1154_plot <- chl_corrected_1154 %>%
  select(flurometer_serial_no:calibration_slope, slope_dfo:slope_dfo_corr,
         rel_diff_corr) %>%
  add_row(flurometer_serial_no = "720001154",
          calibration = as.Date("2018-05-04", tz = "UTC"),
          acid_coefficient = NA, calibration_slope = 0.0005026532,
          slope_dfo = 0.0005026532, fm_dfo = 1.717642, perc_slope = 0,
          slope_dfo_corr = 0.0005026532, rel_diff_corr = 0)

# ── Instrument 0982: back-correct using DFO 2018 acid ratio ───────────────────
chl_corrected_0982 <- chl %>%
  filter(filter_type == "Bulk GF/F" & flurometer_serial_no == "720000982" &
           !is.na(calibration) & !is.na(before_acid) & !is.na(after_acid)) %>%
  select(date, volume, acetone_volume_ml, flurometer_serial_no:after_acid, chla) %>%
  mutate(chl_dfo = calibration_slope *
           (1.831313 / (1.831313 - 1)) *
           (before_acid - after_acid) *
           (acetone_volume_ml / volume),
         rel_diff = ((chl_dfo - chla) / chla) * 100) %>%
  mutate_at(vars(rel_diff), funs(round(., 2)))

chl_corrected_0982 <- chl_corrected_0982 %>%
  distinct(calibration, .keep_all = TRUE) %>%
  filter(!as.character(calibration) == "2019-05-09") %>%
  select(flurometer_serial_no:calibration_slope, rel_diff) %>%
  mutate(fm_dfo = 1.831313) %>%
  mutate(acid_coefficient = case_when(acid_coefficient > 1.8 ~ 0,
                                      TRUE ~ as.numeric(acid_coefficient)))

# ── Panels ────────────────────────────────────────────────────────────────────
# Left column: instrument 0982
p1 <- chl_corrected_0982 %>%
  select(calibration, calibration_slope) %>%
  mutate(lab = if_else(as.character(calibration) == "2018-05-04",
                       "DFO", "Hakai")) %>%
  ggplot(aes(x = factor(calibration), y = calibration_slope * 10000, fill = lab)) +
  geom_point(size = 7, pch = 21, color = "black", stroke = 2) +
  theme_bw() +
  ylim(4, 7) +
  labs(title = "Instrument #0982",
       subtitle = "Stable slope; error reflects fm difference only",
       y = expression("Calibration Slope (×" * 10^4 * ")")) +
  scale_fill_manual(values = c(Hakai = "#E41A1C", DFO = "#377EB8")) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        plot.subtitle = element_text(size = 22),
        text = element_text(size = 30),
        axis.text = element_text(colour = "black"))

p2 <- chl_corrected_0982 %>%
  select(calibration, acid_coefficient, fm_dfo) %>%
  pivot_longer(c(acid_coefficient, fm_dfo), names_to = "type", values_to = "fm") %>%
  ggplot(aes(x = factor(calibration), y = fm, fill = type)) +
  geom_bar(stat = "identity", position = "dodge", color = "Black", size = 1,
           width = 0.6) +
  theme_bw() +
  coord_cartesian(ylim = c(1.6, 1.9)) +
  labs(y = "Acid Ratio (fm)") +
  scale_fill_brewer(palette = "Set1", name = element_blank(),
                    labels = c("Hakai", "DFO 2018")) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        text = element_text(size = 30),
        axis.text = element_text(colour = "black"))

p3 <- chl_corrected_0982 %>%
  mutate(pos = rel_diff > 0) %>%
  ggplot(aes(x = factor(calibration), y = rel_diff, fill = pos)) +
  geom_col(position = "identity", colour = "black", width = 0.6, size = 1) +
  geom_hline(yintercept = 0) +
  geom_text(aes(label = rel_diff), vjust = -0.5, size = 7) +
  scale_fill_manual(values = c("#007FFF", "#FF6347"), guide = FALSE) +
  lims(y = c(-12, 12)) +
  theme_bw() +
  labs(x = "Calibration Date", y = "% Diff. Chl") +
  theme(text = element_text(size = 30),
        axis.text = element_text(colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1))

# Right column: instrument 1154
p4 <- chl_corr1154_plot %>%
  select(calibration, calibration_slope, slope_dfo_corr) %>%
  pivot_longer(c(calibration_slope, slope_dfo_corr),
               names_to = "type", values_to = "slope") %>%
  ggplot(aes(x = factor(calibration), y = slope * 10000, fill = type)) +
  geom_point(size = 7, pch = 21, color = "black", stroke = 2) +
  theme_bw() +
  ylim(4, 7) +
  labs(title = "Instrument #1154",
       subtitle = "Slope drift modeled from 0982 stability",
       y = "Calibration Slope") +
  scale_fill_brewer(palette = "Set1", name = element_blank(),
                    labels = c("Hakai", "DFO 2018 Drift Corrected (-1.5%)")) +
  theme(legend.position = c(0.6, 0.9),
        legend.background = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.subtitle = element_text(size = 22),
        text = element_text(size = 30),
        axis.text = element_text(colour = "black"))

p5 <- chl_corr1154_plot %>%
  select(calibration, acid_coefficient, fm_dfo) %>%
  pivot_longer(c(acid_coefficient, fm_dfo), names_to = "type", values_to = "fm") %>%
  ggplot(aes(x = factor(calibration), y = fm, fill = type)) +
  geom_bar(stat = "identity", position = "dodge", color = "Black", size = 1) +
  theme_bw() +
  coord_cartesian(ylim = c(1.6, 1.9)) +
  labs(y = "Acid Ratio (fm)") +
  scale_fill_brewer(palette = "Set1", name = element_blank(),
                    labels = c("Hakai", "DFO 2018")) +
  theme(legend.position = c(0.85, 0.93),
        legend.background = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        text = element_text(size = 30),
        axis.text = element_text(colour = "black"))

p6 <- chl_corr1154_plot %>%
  mutate(pos = rel_diff_corr > 0) %>%
  ggplot(aes(x = factor(calibration), y = rel_diff_corr, fill = pos)) +
  geom_col(position = "identity", colour = "black", width = 0.75, size = 1) +
  geom_hline(yintercept = 0) +
  geom_text(aes(label = rel_diff_corr), vjust = -0.5, size = 7) +
  scale_fill_manual(values = c("#007FFF", "#FF6347"), guide = FALSE) +
  lims(y = c(-12, 12)) +
  theme_bw() +
  labs(x = "Calibration Date", y = "% Diff. Chl") +
  theme(text = element_text(size = 30),
        axis.text = element_text(colour = "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

# ── Assemble & save ───────────────────────────────────────────────────────────
fig <- p1 + p4 + p2 + p5 + p3 + p6 + plot_layout(ncol = 2) +
  plot_annotation(tag_levels = "A",
                  theme = theme(plot.tag = element_text(size = 30)))

ggsave(here("figures", "backcali_perc_diff_panel_test.png"), fig,
       width = 16, height = 17, dpi = 300)


# ─────────────────────────────────────────────────────────────────────────────
# HPLC vs fluorometric chlorophyll comparison — ADL vs in-range samples
#
# Left:  scatter of fluorometric Chl vs HPLC TChla, 1:1 reference line,
#        points coloured by group (in-range [AV/NA] vs over-range [ADL]).
# Right: boxplot of Chl/TChla ratio by group, reference line at 1.0.
#
# Assumptions (adjust as needed):
#   - grouping comes from the fluorometric chla_flag: ADL = over-range,
#     AV/NA collapsed into a single "In range" group
#   - inner join: only samples with both fluorometric and HPLC values
#   - ratio = chla (fluorometric) / all_chl_a (HPLC)
#   - HPLC analyzing-lab column is `analyzing_lab`  <-- CONFIRM
# ─────────────────────────────────────────────────────────────────────────────

# ── Filter fluorometric data ──────────────────────────────────────────────────
chl_f <- data_chl %>%
  filter(filter_type == "Bulk GF/F",
         is.na(chla_flag) | chla_flag %in% c("AV", "ADL")) %>%
  select(collected, line_out_depth, site_id, chla, chla_flag)

# ── Filter HPLC data ──────────────────────────────────────────────────────────
hplc_f <- data_hplc %>%
  filter(analyzing_lab == "USC",
         is.na(all_chl_a_flag) | all_chl_a_flag == "AV") %>%
  select(collected, line_out_depth, site_id, all_chl_a)

# ── Join and build grouping / ratio ───────────────────────────────────────────
comp <- chl_f %>%
  inner_join(hplc_f, by = c("collected", "line_out_depth", "site_id")) %>%
  mutate(group = case_when(
    chla_flag %in% "ADL"            ~ "Over range (ADL)",
    all_chl_a > 5                   ~ "In range, >5 \u00b5g/L",
    TRUE                            ~ "In range"),
    group = factor(group, levels = c("In range",
                                     "In range, >5 \u00b5g/L",
                                     "Over range (ADL)")),
    ratio = chla / all_chl_a) %>%
  filter(is.finite(ratio))   # drops any divide-by-zero / NA ratios

grp_cols <- c("In range"      = "#1B9E77",
              "In range, >5 \u00b5g/L" = "#7570B3",
              "Over range (ADL)"      = "#D95F02")

# ── Panel A: scatter vs 1:1 line ──────────────────────────────────────────────
ax_max <- max(c(comp$all_chl_a, comp$chla), na.rm = TRUE)

p_scatter <- comp %>%
  ggplot(aes(all_chl_a, chla, fill = group)) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.8) +   # 1:1 reference
  geom_point(size = 3, pch = 21, colour = "black", stroke = 0.6, alpha = 0.85) +
  coord_cartesian(xlim = c(0, ax_max), ylim = c(0, ax_max)) +
  scale_fill_manual(values = grp_cols, name = NULL) +
  labs(x = expression("HPLC TChl " * italic(a) * " (" * mu * "g/L)"),
       y = expression("Fluorometric Chl (" * mu * "g/L)")) +
  theme_bw() +
  theme(legend.position = c(0.7, 0.12),
        aspect.ratio = 1,
        plot.margin = margin(1, 1, 1, 1),
        legend.background = element_blank(),
        text = element_text(size = 22),
        axis.text = element_text(colour = "black"))
# ── Panel B: ratio boxplot ────────────────────────────────────────────────────
p_box <- comp %>%
  ggplot(aes(group, ratio, fill = group)) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.6) + # agreement
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.85) +
  geom_jitter(width = 0.12, size = 1.6, alpha = 0.5) +
  scale_fill_manual(values = grp_cols, guide = "none") +
  labs(x = NULL,
       y = expression("Fluorometric Chl / HPLC TChl " * italic(a))) +
  theme_bw() +
  theme(text = element_text(size = 22),
        # aspect.ratio = 1,
        plot.margin = margin(1, 1, 1, 1),
        axis.text = element_text(colour = "black"),
        axis.text.x = element_text(angle = 15, hjust = 1))
# ── Assemble & save ───────────────────────────────────────────────────────────
fig <- p_scatter + p_box +
  plot_annotation(tag_levels = "A",
                  theme = theme(plot.tag = element_text(size = 18)))

ggsave(here("figures", "hplc_adl_comparison.png"),
       fig, width = 13, height = 6.5, dpi = 300)

# Quick numeric summary for the writeup
comp %>%
  group_by(group) %>%
  summarise(n = n(),
            median_ratio = median(ratio),
            mean_ratio = mean(ratio),
            iqr = IQR(ratio))


# ═════════════════════════════════════════════════════════════════════════════
# Issue 3 — Acid contamination (Calvert station)
#
# Acid contamination acidifies the before-acid fluorometric reading, driving the
# acid ratio (fm = before/after) down and suppressing calculated chlorophyll.
# It is an analysis-SESSION effect: a contaminated acetone batch or cuvette
# affects most samples run that day, so contamination is identified at the
# ANALYSIS-DAY level, not the individual sample level.
#
# Approach:
#   1. Compute fm per sample (all filter types, all depths, full record).
#   2. Classify each analysis day by the proportion of its samples with fm < 1.2.
#      A day is "Contaminated" if >= 80% of its samples fall below 1.2.
#   3. Figure 1 (fm): per-day proportion distribution + fm by day type.
#   4. Figure 2 (HPLC): fluorometric Chl vs HPLC TChla, Bulk GF/F matched,
#      showing that low-fm / contaminated-day samples are suppressed vs HPLC.
#
# Cutoff justification: the 0.80 threshold was chosen using the HPLC comparison
# (Figure 2 / evidence table) — analysis days above ~0.80 showed clear
# fluorometric suppression relative to HPLC, while lower-proportion days
# reflected isolated, individually-flagged low-fm samples rather than
# session-wide contamination. 0.80 also sits in an empty gap in the per-day
# proportion distribution (no days between ~0.69 and ~0.87).
#
# Inputs: data_chl_cal (work_area=CALVERT chlorophyll), data_hplc_cal
#         (work_area=CALVERT hplc).
# ═════════════════════════════════════════════════════════════════════════════

library(tidyverse)
library(here)
library(patchwork)

# ── Parameters ────────────────────────────────────────────────────────────────
FM_THRESH     <- 1.2     # low-fm (contamination / high-phaeopigment) threshold
PROP_THRESH   <- 0.80    # fraction of a day's samples below FM_THRESH to flag it
MIN_N         <- 4       # minimum samples per analysis day to classify it
DISPLAY_FLOOR <- 0.01    # floor for <=0 fluorometric values on the log axis

grp_cols  <- c("Normal" = "#1B9E77", "Contaminated" = "#D95F02")
flag_cols <- c("fm \u2265 1.2" = "#1B9E77", "fm < 1.2" = "#D95F02")

# ═════════════════════════════════════════════════════════════════════════════
# 1. Sample-level fm and analysis-day classification
# ═════════════════════════════════════════════════════════════════════════════

# All filter types / depths / years: the day-level proportion rule is robust to
# individual low-fm samples (deep, size-fractionated), so no pre-filtering needed.
fm_dat <- data_chl_cal %>%
  filter(!is.na(before_acid), !is.na(after_acid), after_acid > 0) %>%
  mutate(analysis_date = as.Date(substr(analyzed, 1, 10)),  # tz-proof date
         fm = before_acid / after_acid)

day_class <- fm_dat %>%
  group_by(analysis_date) %>%
  summarise(n = n(),
            prop_low = mean(fm < FM_THRESH),
            median_fm = median(fm),
            .groups = "drop") %>%
  filter(n >= MIN_N) %>%
  mutate(day_type = if_else(prop_low >= PROP_THRESH, "Contaminated", "Normal"))

# Pool for the fm comparison: Normal = AV-flagged (accepted) data only;
# Contaminated = all samples, to show the full extent of the low-fm effect.
fm_grouped <- fm_dat %>%
  inner_join(day_class %>% select(analysis_date, day_type), by = "analysis_date") %>%
  filter(day_type == "Contaminated" | (day_type == "Normal" & chla_flag == "AV")) %>%
  mutate(day_type = factor(day_type, levels = c("Normal", "Contaminated")))

# ═════════════════════════════════════════════════════════════════════════════
# 2. Figure 1 — fm evidence (per-day proportion + fm by day type)
# ═════════════════════════════════════════════════════════════════════════════

# Panel A: distribution of per-day low-fm proportion (clustering by day)
p_prop <- day_class %>%
  ggplot(aes(prop_low, fill = day_type)) +
  geom_histogram(breaks = seq(0, 1, by = 0.1), colour = "black",
                 linewidth = 0.3, alpha = 0.85) +
  geom_vline(xintercept = PROP_THRESH, linetype = "dashed", linewidth = 0.6) +
  scale_fill_manual(values = grp_cols, name = NULL) +
  labs(x = "Proportion of day's samples with fm < 1.2",
       y = "Number of analysis days") +
  theme_bw() +
  theme(legend.position = c(0.5, 0.85),
        legend.background = element_blank(),
        text = element_text(size = 22),
        axis.text = element_text(colour = "black"))

# Panel B: fm by day type
p_fm <- fm_grouped %>%
  ggplot(aes(day_type, fm, fill = day_type)) +
  geom_hline(yintercept = FM_THRESH, linetype = "dashed", linewidth = 0.6) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.85) +
  geom_jitter(width = 0.12, size = 1.4, alpha = 0.35) +
  scale_fill_manual(values = grp_cols, guide = "none") +
  labs(x = NULL, y = "Acid Ratio (fm)") +
  theme_bw() +
  theme(text = element_text(size = 22),
        axis.text = element_text(colour = "black"))

fig_fm <- p_prop + p_fm +
  plot_annotation(tag_levels = "A",
                  theme = theme(plot.tag = element_text(size = 22)))

ggsave(here("figures", "acid_contamination_fm.png"),
       fig_fm, width = 14, height = 7, dpi = 300)

# ═════════════════════════════════════════════════════════════════════════════
# 3. Figure 2 — HPLC corroboration (Bulk GF/F matched, log-log)
# ═════════════════════════════════════════════════════════════════════════════

# HPLC side: USC lab, valid TChla
hplc_cal_f <- data_hplc_cal %>%
  filter(analyzing_lab == "USC",
         is.na(all_chl_a_flag) | all_chl_a_flag == "AV") %>%
  select(collected, line_out_depth, site_id, all_chl_a)

# Fluorometric side: Bulk GF/F only (HPLC TChla is whole-water; avoids
# size-fraction fan-out), carrying day_type and prop_low.
fluor_bulk <- fm_dat %>%
  filter(filter_type == "Bulk GF/F") %>%
  inner_join(day_class %>% select(analysis_date, day_type, prop_low),
             by = "analysis_date") %>%
  select(collected, line_out_depth, site_id, chla, fm,
         day_type, prop_low, analysis_date)

matched <- fluor_bulk %>%
  inner_join(hplc_cal_f, by = c("collected", "line_out_depth", "site_id")) %>%
  mutate(fm_flag = if_else(fm < FM_THRESH, "fm < 1.2", "fm \u2265 1.2"),
         fm_flag = factor(fm_flag, levels = c("fm \u2265 1.2", "fm < 1.2")),
         chla_disp = pmax(chla, DISPLAY_FLOOR),
         year = factor(format(analysis_date, "%Y")))

# Join sanity check: should be zero (one bulk fluorometric per HPLC)
dup_check <- matched %>%
  count(collected, line_out_depth, site_id) %>%
  filter(n > 1)
cat("Duplicated join keys (should be 0):", nrow(dup_check), "\n")
cat("Fluorometric values <= 0 (floored for display):",
    sum(matched$chla <= 0, na.rm = TRUE), "\n")

ax_min <- DISPLAY_FLOOR
ax_max <- max(c(matched$all_chl_a, matched$chla_disp), na.rm = TRUE)

# Shared scatter scaffold to avoid repetition
base_scatter <- function(df, fill_aes) {
  ggplot(df, aes(all_chl_a, chla_disp, fill = {{ fill_aes }})) +
    geom_abline(slope = 1, intercept = 0, linewidth = 0.8) +
    geom_point(size = 4, pch = 21, colour = "black", stroke = 0.6, alpha = 0.85) +
    scale_x_log10(limits = c(ax_min, ax_max)) +
    scale_y_log10(limits = c(ax_min, ax_max)) +
    annotation_logticks(sides = "bl") +
    coord_equal() +
    labs(x = expression("HPLC TChl " * italic(a) * " (" * mu * "g/L)"),
         y = expression("Fluorometric Chl (" * mu * "g/L)")) +
    theme_bw() +
    theme(legend.position = c(0.02, 0.98),
          legend.justification = c(0, 1),
          legend.background = element_blank(),
          text = element_text(size = 24),
          axis.text = element_text(colour = "black"))
}

# Panel A: sample-level fm flag
pA <- base_scatter(arrange(matched, fm_flag), fm_flag) +
  scale_fill_manual(values = flag_cols, name = NULL)

# Panel B: day-level prop_low (continuous)
pB <- base_scatter(arrange(matched, prop_low), prop_low) +
  scale_fill_viridis_c(name = "Day prop.\nfm < 1.2", option = "inferno",
                       direction = -1, limits = c(0, 1)) +
  theme(legend.key.size = unit(0.5, "cm"))

# Panel C: analysis year
pC <- base_scatter(arrange(matched, year), year) +
  scale_fill_viridis_d(name = "Analysis\nyear", option = "turbo") +
  theme(legend.key.size = unit(0.4, "cm"))

fig_hplc <- pA + pB + pC +
  plot_annotation(tag_levels = "A",
                  theme = theme(plot.tag = element_text(size = 22)))

ggsave(here("figures", "acid_contamination_hplc.png"),
       fig_hplc, width = 21, height = 8, dpi = 300)

# ═════════════════════════════════════════════════════════════════════════════
# 4. Summaries for the writeup
# ═════════════════════════════════════════════════════════════════════════════

# Days and samples per group
day_class %>% count(day_type, name = "n_days")
fm_grouped %>%
  group_by(day_type) %>%
  summarise(n_samples = n(),
            median_fm = median(fm),
            mean_fm = mean(fm),
            .groups = "drop")

# Contaminated days and their dates
day_class %>%
  filter(day_type == "Contaminated") %>%
  arrange(analysis_date) %>%
  select(analysis_date, n, prop_low, median_fm) %>%
  print(n = Inf)

# HPLC matched-pair summary by fm flag
matched %>%
  group_by(fm_flag) %>%
  summarise(n = n(),
            median_hplc  = median(all_chl_a, na.rm = TRUE),
            median_fluor = median(chla, na.rm = TRUE),
            median_ratio = median(all_chl_a / chla, na.rm = TRUE),
            n_fluor_le0  = sum(chla <= 0, na.rm = TRUE),
            .groups = "drop")

# Evidence table — low-fm matched pairs with prop_low and day_type
matched %>%
  filter(fm_flag == "fm < 1.2") %>%
  select(analysis_date, site_id, line_out_depth,
         all_chl_a, chla, fm, prop_low, day_type) %>%
  arrange(desc(prop_low), fm) %>%
  print(n = Inf)

# Day-classification table — confirm 0.80 sits in an empty gap
day_class %>% arrange(desc(prop_low)) %>% print(n = Inf)




# All samples from that analysis day, full detail
day_raw <- fm_dat %>%
  filter(analysis_date == as.Date("2024-01-29")) %>%
  select(analysis_date, collected, site_id, line_out_depth, filter_type,
         chla_flag, before_acid, after_acid, fm, chla) %>%
  arrange(fm)

print(day_raw, n = Inf)

# Quick diagnostics on what kind of samples these are
day_raw %>% count(filter_type)              # all size-fractionated?
day_raw %>% count(site_id)                  # one site / cast?
day_raw %>% summarise(min_depth = min(line_out_depth),
                      max_depth = max(line_out_depth))   # all deep?
day_raw %>% count(chla_flag)                # already flagged?

# ─────────────────────────────────────────────────────────────────────────────
# Figure 2 — Calibration slope and acid ratio (fm) by calibration date
# Instruments #0982 (left, A/C) and #1154 (right, B/D)
# Output: figures/slope_fm_compare_summary3.png
# ─────────────────────────────────────────────────────────────────────────────

# ── Data ──────────────────────────────────────────────────────────────────────
cali <- read_xlsx(here("files", "cali_summary_2022-04-06.xlsx")) %>%
  mutate(across(r2, ~ round(., 4)))   # round r2 for plotting labels

# ── Panel A: #0982 slope ──────────────────────────────────────────────────────
p1 <- cali %>%
  filter(Fluorometer == 720000982) %>%
  ggplot(aes(factor(Date), slope * 10000, fill = Lab)) +
  geom_point(size = 6, pch = 21, color = "black") +
  geom_text(aes(label = r2), vjust = -0.8, size = 7) +
  theme_bw() +
  coord_cartesian(ylim = c(4, 7)) +
  annotate("text", x = 3.6, y = 4.7,
           label = "Pre-DFO Test calibration", size = 6, angle = 90) +
  annotate("rect", xmin = 3.5, xmax = 4.5, ymin = 3, ymax = 8,
           alpha = 0.2, fill = "red") +
  labs(title = "Instrument #0982",
       y = expression("Calibration Slope (×" * 10^4 * ")"),
       x = "Cali. Date") +
  scale_fill_brewer(palette = "Set1") +
  theme(legend.position = c(0.21, 0.18),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        text = element_text(size = 30),
        axis.text = element_text(colour = "black"))

# ── Panel B: #1154 slope ──────────────────────────────────────────────────────
p2 <- cali %>%
  filter(Fluorometer == 720001154) %>%
  ggplot(aes(factor(Date), slope * 10000, fill = Lab)) +
  geom_point(size = 6, pch = 21, color = "black") +
  geom_text(aes(label = r2), vjust = -0.8, size = 7) +
  scale_fill_manual(breaks = c("Hakai-Quadra", "DFO"),
                    values = c("#4DAF4A", "#E41A1C")) +
  theme_bw() +
  coord_cartesian(ylim = c(4, 7)) +
  annotate("rect", xmin = 5.5, xmax = 6.5, ymin = 3, ymax = 8,
           alpha = 0.2, fill = "red") +
  annotate("rect", xmin = 8.5, xmax = 9.5, ymin = 3, ymax = 8,
           alpha = 0.2, fill = "blue") +
  annotate("text", x = 8.6, y = 6.4,
           label = "Hakai DFO methods", size = 6, angle = 90) +
  labs(title = "Instrument #1154",
       y = expression("Calibration Slope (×" * 10^4 * ")"),
       x = "Cali. Date") +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        text = element_text(size = 30),
        axis.text = element_text(colour = "black"))

# ── Panel C: #0982 fm ─────────────────────────────────────────────────────────
p3 <- cali %>%
  filter(Fluorometer == 720000982) %>%
  ggplot(aes(factor(Date), fm, fill = Lab)) +
  geom_bar(stat = "identity", color = "Black", size = 1, width = 0.6) +
  theme_bw() +
  coord_cartesian(ylim = c(1.65, 1.90)) +
  annotate("rect", xmin = 3.5, xmax = 4.5, ymin = 1, ymax = 2,
           alpha = 0.2, fill = "red") +
  labs(y = "Acid Ratio (fm)",
       x = "Cali. Date") +
  scale_fill_brewer(palette = "Set1") +
  theme(legend.position = "none",
        text = element_text(size = 30),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(colour = "black"))

# ── Panel D: #1154 fm ─────────────────────────────────────────────────────────
p4 <- cali %>%
  filter(Fluorometer == 720001154) %>%
  ggplot(aes(factor(Date), fm, fill = Lab)) +
  geom_bar(stat = "identity", color = "Black", size = 1, width = 0.6) +
  theme_bw() +
  coord_cartesian(ylim = c(1.65, 1.90)) +
  annotate("rect", xmin = 8.5, xmax = 9.5, ymin = 1, ymax = 2,
           alpha = 0.2, fill = "blue") +
  annotate("rect", xmin = 5.5, xmax = 6.5, ymin = 1, ymax = 2,
           alpha = 0.2, fill = "red") +
  labs(y = "Acid Ratio (fm)",
       x = "Cali. Date") +
  scale_fill_manual(breaks = c("Hakai-Quadra", "DFO"),
                    values = c("#4DAF4A", "#E41A1C")) +
  theme(legend.position = "none",
        text = element_text(size = 30),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text = element_text(colour = "black"))

# ── Assemble & save ───────────────────────────────────────────────────────────
fig <- p1 + p2 + p3 + p4 +
  plot_annotation(tag_levels = "A",
                  theme = theme(plot.tag = element_text(size = 30)))

ggsave(here("figures", "slope_fm_compare_2026-07-10.png"),
       fig, width = 16, height = 15, dpi = 300)


