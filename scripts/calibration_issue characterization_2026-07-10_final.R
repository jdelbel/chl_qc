# ═════════════════════════════════════════════════════════════════════════════
# Benchtop chlorophyll fluorometer — data-quality review figures
#
# Reproduces all figures for the QC review document, in document order:
#   Figure 1  Back-calibration % difference panel        (Issue 1)
#   Figure 2  Calibration slope & acid ratio by date     (Issue 1)
#   Figure 3  HPLC vs fluorometric, ADL comparison        (Issue 2)
#   Figure 4  Acid contamination — fm evidence            (Issue 3)
#   Figure 5  Acid contamination — HPLC corroboration     (Issue 3)
#
# Plus supporting numeric summaries and a diagnostics section at the end.
#
# Data sources:
#   files/hakai_chl_all.xlsx            historical chl (back-cali, Fig 1)
#   files/cali_summary_2022-04-06.xlsx  calibration summary (Fig 2)
#   Hakai API                           QU39 chl + hplc (Fig 3),
#                                       CALVERT chl + hplc (Fig 4, 5)
# ═════════════════════════════════════════════════════════════════════════════

# ── Libraries ─────────────────────────────────────────────────────────────────
library(tidyverse)
library(readxl)
library(here)
library(patchwork)
library(hakaiApi)

# ── Shared parameters / palettes ──────────────────────────────────────────────
FM_THRESH     <- 1.2     # low-fm (contamination / high-phaeopigment) threshold
PROP_THRESH   <- 0.80    # fraction of a day's samples below FM_THRESH to flag it
MIN_N         <- 4       # minimum samples per analysis day to classify it
DISPLAY_FLOOR <- 0.01    # floor for <=0 fluorometric values on the log axis

# Issue 2 (ADL groups)
adl_cols  <- c("In range"                = "#1B9E77",
               "In range, >5 \u00b5g/L"  = "#7570B3",
               "Over range (ADL)"        = "#D95F02")
# Issue 3 (day type / fm flag)
grp_cols  <- c("Normal" = "#1B9E77", "Contaminated" = "#D95F02")
flag_cols <- c("fm \u2265 1.2" = "#1B9E77", "fm < 1.2" = "#D95F02")

# ═════════════════════════════════════════════════════════════════════════════
# DATA
# ═════════════════════════════════════════════════════════════════════════════

# ── Local files ───────────────────────────────────────────────────────────────
chl  <- read_xlsx(here("files", "hakai_chl_all.xlsx"), sheet = "Hakai Data")
cali <- read_xlsx(here("files", "cali_summary_2022-04-06.xlsx")) %>%
  mutate(across(r2, ~ round(., 4)))     # round r2 for plotting labels

# ── Hakai API ─────────────────────────────────────────────────────────────────
client <- Client$new()

endpoint      <- "/eims/views/output/chlorophyll"
endpoint_hplc <- "/eims/views/output/hplc"

# QU39 (Issue 2)
data_chl  <- client$get(paste0("https://portal.hakai.org/api", endpoint,
                               "?limit=-1&site_id=QU39"))
data_hplc <- client$get(paste0("https://portal.hakai.org/api", endpoint_hplc,
                               "?limit=-1&site_id=QU39"))

# CALVERT (Issue 3)
data_chl_cal  <- client$get(paste0("https://portal.hakai.org/api", endpoint,
                                   "?limit=-1&work_area=CALVERT"))
data_hplc_cal <- client$get(paste0("https://portal.hakai.org/api", endpoint_hplc,
                                   "?limit=-1&work_area=CALVERT"))

# ═════════════════════════════════════════════════════════════════════════════
# FIGURE 1 — Back-calibration % difference panel (Issue 1)
#
# Estimates error introduced by variability in calibration slopes and acid ratios
# (fm), assuming instrument stability and a linear backward drift from the DFO
# 2018-05-04 benchmark calibration.
#   #0982 (left column):  stable slope; error reflects fm difference only
#   #1154 (right column): slope drift modelled from #0982 stability
#   Rows: calibration slope | acid ratio (fm) | % difference in chl
# ═════════════════════════════════════════════════════════════════════════════

# ── #1154: back-calibrate to DFO 2018-05-04 with graduated linear drift ───────
chl_corrected_1154 <- chl %>%
  filter(filter_type == "Bulk GF/F" & flurometer_serial_no == "720001154" &
           !is.na(calibration) & !is.na(before_acid) & !is.na(after_acid)) %>%
  select(date, volume, acetone_volume_ml, flurometer_serial_no:after_acid, chla) %>%
  mutate(chl_dfo = 0.0005026532 *
           (1.717642 / (1.717642 - 1)) *
           (before_acid - after_acid) *
           (acetone_volume_ml / volume),
         rel_diff = ((chl_dfo - chla) / chla) * 100) %>%
  mutate(across(rel_diff, ~ round(., 2)))

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
  mutate(across(rel_diff_corr, ~ round(., 2))) %>%
  arrange(calibration)

chl_corr1154_plot <- chl_corrected_1154 %>%
  select(flurometer_serial_no:calibration_slope, slope_dfo:slope_dfo_corr,
         rel_diff_corr) %>%
  add_row(flurometer_serial_no = "720001154",
          calibration = as.Date("2018-05-04", tz = "UTC"),
          acid_coefficient = NA, calibration_slope = 0.0005026532,
          slope_dfo = 0.0005026532, fm_dfo = 1.717642, perc_slope = 0,
          slope_dfo_corr = 0.0005026532, rel_diff_corr = 0)

# ── #0982: back-correct using DFO 2018 acid ratio (slope stable, no drift) ────
chl_corrected_0982 <- chl %>%
  filter(filter_type == "Bulk GF/F" & flurometer_serial_no == "720000982" &
           !is.na(calibration) & !is.na(before_acid) & !is.na(after_acid)) %>%
  select(date, volume, acetone_volume_ml, flurometer_serial_no:after_acid, chla) %>%
  mutate(chl_dfo = calibration_slope *
           (1.831313 / (1.831313 - 1)) *
           (before_acid - after_acid) *
           (acetone_volume_ml / volume),
         rel_diff = ((chl_dfo - chla) / chla) * 100) %>%
  mutate(across(rel_diff, ~ round(., 2)))

chl_corrected_0982 <- chl_corrected_0982 %>%
  distinct(calibration, .keep_all = TRUE) %>%
  filter(!as.character(calibration) == "2019-05-09") %>%
  select(flurometer_serial_no:calibration_slope, rel_diff) %>%
  mutate(fm_dfo = 1.831313) %>%
  mutate(acid_coefficient = case_when(acid_coefficient > 1.8 ~ 0,
                                      TRUE ~ as.numeric(acid_coefficient)))

# ── Left column: #0982 ────────────────────────────────────────────────────────
f1_p1 <- chl_corrected_0982 %>%
  select(calibration, calibration_slope) %>%
  mutate(lab = if_else(as.character(calibration) == "2018-05-04", "DFO", "Hakai")) %>%
  ggplot(aes(factor(calibration), calibration_slope * 10000, fill = lab)) +
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

f1_p2 <- chl_corrected_0982 %>%
  select(calibration, acid_coefficient, fm_dfo) %>%
  pivot_longer(c(acid_coefficient, fm_dfo), names_to = "type", values_to = "fm") %>%
  ggplot(aes(factor(calibration), fm, fill = type)) +
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

f1_p3 <- chl_corrected_0982 %>%
  mutate(pos = rel_diff > 0) %>%
  ggplot(aes(factor(calibration), rel_diff, fill = pos)) +
  geom_col(position = "identity", colour = "black", width = 0.6, size = 1) +
  geom_hline(yintercept = 0) +
  geom_text(aes(label = rel_diff), vjust = -0.5, size = 7) +
  scale_fill_manual(values = c("#007FFF", "#FF6347"), guide = "none") +
  lims(y = c(-12, 12)) +
  theme_bw() +
  labs(x = "Calibration Date", y = "% Diff. Chl") +
  theme(text = element_text(size = 30),
        axis.text = element_text(colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1))

# ── Right column: #1154 ───────────────────────────────────────────────────────
f1_p4 <- chl_corr1154_plot %>%
  select(calibration, calibration_slope, slope_dfo_corr) %>%
  pivot_longer(c(calibration_slope, slope_dfo_corr),
               names_to = "type", values_to = "slope") %>%
  ggplot(aes(factor(calibration), slope * 10000, fill = type)) +
  geom_point(size = 7, pch = 21, color = "black", stroke = 2) +
  theme_bw() +
  ylim(4, 7) +
  labs(title = "Instrument #1154",
       subtitle = "Slope drift modeled from 0982 stability",
       y = expression("Calibration Slope (×" * 10^4 * ")")) +
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

f1_p5 <- chl_corr1154_plot %>%
  select(calibration, acid_coefficient, fm_dfo) %>%
  pivot_longer(c(acid_coefficient, fm_dfo), names_to = "type", values_to = "fm") %>%
  ggplot(aes(factor(calibration), fm, fill = type)) +
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

f1_p6 <- chl_corr1154_plot %>%
  mutate(pos = rel_diff_corr > 0) %>%
  ggplot(aes(factor(calibration), rel_diff_corr, fill = pos)) +
  geom_col(position = "identity", colour = "black", width = 0.75, size = 1) +
  geom_hline(yintercept = 0) +
  geom_text(aes(label = rel_diff_corr), vjust = -0.5, size = 7) +
  scale_fill_manual(values = c("#007FFF", "#FF6347"), guide = "none") +
  lims(y = c(-12, 12)) +
  theme_bw() +
  labs(x = "Calibration Date", y = "% Diff. Chl") +
  theme(text = element_text(size = 30),
        axis.text = element_text(colour = "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

fig1 <- f1_p1 + f1_p4 + f1_p2 + f1_p5 + f1_p3 + f1_p6 + plot_layout(ncol = 2) +
  plot_annotation(tag_levels = "A",
                  theme = theme(plot.tag = element_text(size = 30)))

ggsave(here("figures", "fig1_backcali_perc_diff_2.png"), fig1,
       width = 16, height = 17, dpi = 300)

# ═════════════════════════════════════════════════════════════════════════════
# FIGURE 2 — Calibration slope & acid ratio (fm) by calibration date (Issue 1)
#   #0982 (A/C) and #1154 (B/D); coloured by calibrating lab.
# ═════════════════════════════════════════════════════════════════════════════

f2_p1 <- cali %>%
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

f2_p2 <- cali %>%
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

f2_p3 <- cali %>%
  filter(Fluorometer == 720000982) %>%
  ggplot(aes(factor(Date), fm, fill = Lab)) +
  geom_bar(stat = "identity", color = "Black", size = 1, width = 0.6) +
  theme_bw() +
  coord_cartesian(ylim = c(1.65, 1.90)) +
  annotate("rect", xmin = 3.5, xmax = 4.5, ymin = 1, ymax = 2,
           alpha = 0.2, fill = "red") +
  labs(y = "Acid Ratio (fm)", x = "Cali. Date") +
  scale_fill_brewer(palette = "Set1") +
  theme(legend.position = "none",
        text = element_text(size = 30),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(colour = "black"))

f2_p4 <- cali %>%
  filter(Fluorometer == 720001154) %>%
  ggplot(aes(factor(Date), fm, fill = Lab)) +
  geom_bar(stat = "identity", color = "Black", size = 1, width = 0.6) +
  theme_bw() +
  coord_cartesian(ylim = c(1.65, 1.90)) +
  annotate("rect", xmin = 8.5, xmax = 9.5, ymin = 1, ymax = 2,
           alpha = 0.2, fill = "blue") +
  annotate("rect", xmin = 5.5, xmax = 6.5, ymin = 1, ymax = 2,
           alpha = 0.2, fill = "red") +
  labs(y = "Acid Ratio (fm)", x = "Cali. Date") +
  scale_fill_manual(breaks = c("Hakai-Quadra", "DFO"),
                    values = c("#4DAF4A", "#E41A1C")) +
  theme(legend.position = "none",
        text = element_text(size = 30),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text = element_text(colour = "black"))

fig2 <- f2_p1 + f2_p2 + f2_p3 + f2_p4 +
  plot_annotation(tag_levels = "A",
                  theme = theme(plot.tag = element_text(size = 30)))

ggsave(here("figures", "fig2_slope_fm_compare_2.png"), fig2,
       width = 16, height = 15, dpi = 300)

# ═════════════════════════════════════════════════════════════════════════════
# FIGURE 3 — HPLC vs fluorometric chlorophyll, ADL comparison (Issue 2)
#   A: scatter vs 1:1 line, coloured by group
#   B: Chl/TChla ratio boxplot by group (reference line at 1.0)
#   Groups: in-range, in-range >5 ug/L (post-2018 diluted high-biomass control),
#           over-range (ADL).
# ═════════════════════════════════════════════════════════════════════════════

chl_f <- data_chl %>%
  filter(filter_type == "Bulk GF/F",
         is.na(chla_flag) | chla_flag %in% c("AV", "ADL")) %>%
  select(collected, line_out_depth, site_id, chla, chla_flag)

hplc_f <- data_hplc %>%
  filter(analyzing_lab == "USC",
         is.na(all_chl_a_flag) | all_chl_a_flag == "AV") %>%
  select(collected, line_out_depth, site_id, all_chl_a)

comp <- chl_f %>%
  inner_join(hplc_f, by = c("collected", "line_out_depth", "site_id")) %>%
  mutate(group = case_when(
    chla_flag %in% "ADL" ~ "Over range (ADL)",
    all_chl_a > 5        ~ "In range, >5 \u00b5g/L",
    TRUE                 ~ "In range"),
    group = factor(group, levels = c("In range",
                                     "In range, >5 \u00b5g/L",
                                     "Over range (ADL)")),
    ratio = chla / all_chl_a) %>%
  filter(is.finite(ratio))

ax_max <- max(c(comp$all_chl_a, comp$chla), na.rm = TRUE)

f3_scatter <- comp %>%
  ggplot(aes(all_chl_a, chla, fill = group)) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.8) +
  geom_point(size = 3, pch = 21, colour = "black", stroke = 0.6, alpha = 0.85) +
  coord_cartesian(xlim = c(0, ax_max), ylim = c(0, ax_max)) +
  scale_fill_manual(values = adl_cols, name = NULL) +
  labs(x = expression("HPLC TChl " * italic(a) * " (" * mu * "g/L)"),
       y = expression("Fluorometric Chl (" * mu * "g/L)")) +
  theme_bw() +
  theme(legend.position = c(0.7, 0.12),
        aspect.ratio = 1,
        plot.margin = margin(1, 1, 1, 1),
        legend.background = element_blank(),
        text = element_text(size = 22),
        axis.text = element_text(colour = "black"))

f3_box <- comp %>%
  ggplot(aes(group, ratio, fill = group)) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.6) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.85) +
  geom_jitter(width = 0.12, size = 1.6, alpha = 0.5) +
  scale_fill_manual(values = adl_cols, guide = "none") +
  labs(x = NULL,
       y = expression("Fluorometric Chl / HPLC TChl " * italic(a))) +
  theme_bw() +
  theme(text = element_text(size = 22),
        plot.margin = margin(1, 1, 1, 1),
        axis.text = element_text(colour = "black"),
        axis.text.x = element_text(angle = 15, hjust = 1))

fig3 <- f3_scatter + f3_box +
  plot_annotation(tag_levels = "A",
                  theme = theme(plot.tag = element_text(size = 18)))

ggsave(here("figures", "fig3_hplc_adl_comparison_2.png"), fig3,
       width = 13, height = 6.5, dpi = 300)

# ═════════════════════════════════════════════════════════════════════════════
# FIGURE 4 & 5 — Acid contamination (Issue 3), Calvert station
#
# Contamination acidifies the before-acid reading, driving fm (= before/after)
# toward 1 and suppressing calculated chlorophyll. It is an analysis-SESSION
# effect, so days are classified by the proportion of samples with fm < 1.2;
# a day is "Contaminated" if >= 80% of its samples fall below the threshold.
#
# Cutoff (0.80) justification: the HPLC comparison (Fig 5 / evidence table)
# shows days above ~0.80 with clear fluorometric suppression, while
# lower-proportion days reflect isolated, individually-flagged samples. 0.80
# also sits in an empty gap in the per-day proportion distribution.
# ═════════════════════════════════════════════════════════════════════════════

# ── Sample-level fm and analysis-day classification ───────────────────────────
# All filter types / depths / years: the day-level proportion rule is robust to
# individual low-fm samples (deep, size-fractionated), so no pre-filtering here.
fm_dat <- data_chl_cal %>%
  filter(!is.na(before_acid), !is.na(after_acid), after_acid > 0) %>%
  mutate(analysis_date = as.Date(substr(analyzed, 1, 10)),   # tz-proof date
         fm = before_acid / after_acid)

day_class <- fm_dat %>%
  group_by(analysis_date) %>%
  summarise(n = n(),
            prop_low = mean(fm < FM_THRESH),
            median_fm = median(fm),
            .groups = "drop") %>%
  filter(n >= MIN_N) %>%
  mutate(day_type = if_else(prop_low >= PROP_THRESH, "Contaminated", "Normal"))

# Pool for the fm comparison: Normal = AV-flagged data only; Contaminated = all
# samples, to show the full extent of the low-fm effect.
fm_grouped <- fm_dat %>%
  inner_join(day_class %>% select(analysis_date, day_type), by = "analysis_date") %>%
  filter(day_type == "Contaminated" | (day_type == "Normal" & chla_flag == "AV")) %>%
  mutate(day_type = factor(day_type, levels = c("Normal", "Contaminated")))

# ── Figure 4: fm evidence (per-day proportion + fm by day type) ───────────────
f4_prop <- day_class %>%
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

f4_fm <- fm_grouped %>%
  ggplot(aes(day_type, fm, fill = day_type)) +
  geom_hline(yintercept = FM_THRESH, linetype = "dashed", linewidth = 0.6) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.85) +
  geom_jitter(width = 0.12, size = 1.4, alpha = 0.35) +
  scale_fill_manual(values = grp_cols, guide = "none") +
  labs(x = NULL, y = "Acid Ratio (fm)") +
  theme_bw() +
  theme(text = element_text(size = 22),
        axis.text = element_text(colour = "black"))

fig4 <- f4_prop + f4_fm +
  plot_annotation(tag_levels = "A",
                  theme = theme(plot.tag = element_text(size = 22)))

ggsave(here("figures", "fig4_acid_contamination_fm.png"), fig4,
       width = 14, height = 7, dpi = 300)

# ── Figure 5: HPLC corroboration (Bulk GF/F matched, log-log) ─────────────────
hplc_cal_f <- data_hplc_cal %>%
  filter(analyzing_lab == "USC",
         is.na(all_chl_a_flag) | all_chl_a_flag == "AV") %>%
  select(collected, line_out_depth, site_id, all_chl_a)

# Bulk GF/F only (HPLC TChla is whole-water; avoids size-fraction fan-out).
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

# Join sanity check
cat("Duplicated join keys (should be 0):",
    nrow(matched %>% count(collected, line_out_depth, site_id) %>% filter(n > 1)), "\n")
cat("Fluorometric values <= 0 (floored for display):",
    sum(matched$chla <= 0, na.rm = TRUE), "\n")

f5_min <- DISPLAY_FLOOR
f5_max <- max(c(matched$all_chl_a, matched$chla_disp), na.rm = TRUE)

f5_base <- function(df, fill_aes) {
  ggplot(df, aes(all_chl_a, chla_disp, fill = {{ fill_aes }})) +
    geom_abline(slope = 1, intercept = 0, linewidth = 0.8) +
    geom_point(size = 4, pch = 21, colour = "black", stroke = 0.6, alpha = 0.85) +
    scale_x_log10(limits = c(f5_min, f5_max)) +
    scale_y_log10(limits = c(f5_min, f5_max)) +
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

f5_pA <- f5_base(arrange(matched, fm_flag), fm_flag) +
  scale_fill_manual(values = flag_cols, name = NULL)

f5_pB <- f5_base(arrange(matched, prop_low), prop_low) +
  scale_fill_viridis_c(name = "Day prop.\nfm < 1.2", option = "inferno",
                       direction = -1, limits = c(0, 1)) +
  theme(legend.key.size = unit(0.5, "cm"))

f5_pC <- f5_base(arrange(matched, year), year) +
  scale_fill_viridis_d(name = "Analysis\nyear", option = "turbo") +
  theme(legend.key.size = unit(0.4, "cm"))

fig5 <- f5_pA + f5_pB + f5_pC +
  plot_annotation(tag_levels = "A",
                  theme = theme(plot.tag = element_text(size = 22)))

ggsave(here("figures", "fig5_acid_contamination_hplc_2.png"), fig5,
       width = 21, height = 8, dpi = 300)

# ═════════════════════════════════════════════════════════════════════════════
# NUMERIC SUMMARIES FOR THE WRITEUP
# ═════════════════════════════════════════════════════════════════════════════

# Issue 2 — ADL group summary
comp %>%
  group_by(group) %>%
  summarise(n = n(),
            median_ratio = median(ratio),
            mean_ratio = mean(ratio),
            iqr = IQR(ratio),
            .groups = "drop")

# Issue 3 — days and samples per group
day_class %>% count(day_type, name = "n_days")
fm_grouped %>%
  group_by(day_type) %>%
  summarise(n_samples = n(),
            median_fm = median(fm),
            mean_fm = mean(fm),
            .groups = "drop")

# Issue 3 — contaminated days and their dates
day_class %>%
  filter(day_type == "Contaminated") %>%
  arrange(analysis_date) %>%
  select(analysis_date, n, prop_low, median_fm) %>%
  print(n = Inf)

# Issue 3 — HPLC matched-pair summary by fm flag
matched %>%
  group_by(fm_flag) %>%
  summarise(n = n(),
            median_hplc  = median(all_chl_a, na.rm = TRUE),
            median_fluor = median(chla, na.rm = TRUE),
            median_ratio = median(all_chl_a / chla, na.rm = TRUE),
            n_fluor_le0  = sum(chla <= 0, na.rm = TRUE),
            .groups = "drop")

# Issue 3 — evidence table: low-fm matched pairs with prop_low and day_type
matched %>%
  filter(fm_flag == "fm < 1.2") %>%
  select(analysis_date, site_id, line_out_depth,
         all_chl_a, chla, fm, prop_low, day_type) %>%
  arrange(desc(prop_low), fm) %>%
  print(n = Inf)

# Issue 3 — full day-classification table (confirm 0.80 sits in an empty gap)
day_class %>% arrange(desc(prop_low)) %>% print(n = Inf)

# ═════════════════════════════════════════════════════════════════════════════
# DIAGNOSTICS (exploratory — not part of the final figures)
# ═════════════════════════════════════════════════════════════════════════════

# 2024-01-29: a "contaminated"-flagged day that is actually a benign batch of
# deep samples (false positive). Kept here to document the check.
day_raw <- fm_dat %>%
  filter(analysis_date == as.Date("2024-01-29")) %>%
  select(analysis_date, collected, site_id, line_out_depth, filter_type,
         chla_flag, before_acid, after_acid, fm, chla) %>%
  arrange(fm)

print(day_raw, n = Inf)
day_raw %>% count(filter_type)                                   # size fractions?
day_raw %>% count(site_id)                                       # one cast?
day_raw %>% summarise(min_depth = min(line_out_depth),
                      max_depth = max(line_out_depth))           # all deep?
day_raw %>% count(chla_flag)                                     # already flagged?