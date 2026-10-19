# ============================================================
# CE plane: 3 specificity groups (50/500/1000) + 1 split point
# Fixes + small clarity improvements (no framework changes):
# - Make WTP line dashed (no linetype scale -> no warning)
# - Add spec label to midpoint labels: "All Cancers (0.xx)\n$ICER/QALY"
# ============================================================

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(scales)
library(ggrepel)

fmt_icer <- function(x) paste0("$", format(round(x), big.mark = ","), "/QALY")
wtp <- 100000

# -----------------------------
# 1) Hard filter to your run set
# -----------------------------
keep_runs <- c(
  "spec_0.95_cost_50",   "spec_0.95_cost_500",   "spec_0.95_cost_1000",
  "spec_0.98_cost_50",   "spec_0.98_cost_500",   "spec_0.98_cost_1000",
  "spec_0.9999_cost_50", "spec_0.9999_cost_500", "spec_0.9999_cost_1000",
  "screen_0.98_nonscreen_0.99_cost_195"
)

df <- df_all %>%
  filter(run_name %in% keep_runs) %>%
  mutate(
    test_cost = as.numeric(str_extract(run_name, "(?<=cost_)[0-9]+")),
    key       = str_remove(run_name, "_cost_[0-9]+$")
  )

# -----------------------------
# 2) Split into: spec-series vs split-single
# -----------------------------
df_spec <- df %>%
  filter(str_detect(run_name, "^spec_")) %>%
  mutate(
    cancer_group = case_when(
      str_detect(run_name, "^spec_0\\.95_cost_")   ~ "Specificity 0.95",
      str_detect(run_name, "^spec_0\\.98_cost_")   ~ "Specificity 0.98",
      str_detect(run_name, "^spec_0\\.9999_cost_") ~ "Specificity 0.9999",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(cancer_group)) %>%
  filter(test_cost %in% c(50, 500, 1000))

df_split <- df %>%
  filter(run_name == "screen_0.98_nonscreen_0.99_cost_195") %>%
  mutate(
    point_label = paste0("With screening (0.98) / without (0.99)", "\n", "MCD Test Cost = $195"),
    icer_label  = fmt_icer(icer),
    ce_flag = case_when(
      is.na(icer)    ~ "ICER unavailable",
      inc_qalys <= 0 ~ "Not cost-effective",
      icer <= wtp    ~ "Cost-effective",
      TRUE           ~ "Not cost-effective"
    ),
    ce_label = paste0(icer_label)
  )

# -----------------------------
# 3) Pivot spec-series wide to get 50/500/1000 endpoints
# -----------------------------
plot_df <- df_spec %>%
  select(key, cancer_group, test_cost, inc_cost_M, inc_qalys, icer) %>%
  pivot_wider(
    id_cols     = c(key, cancer_group),
    names_from  = test_cost,
    values_from = c(inc_cost_M, inc_qalys, icer),
    names_glue  = "{.value}_{test_cost}"
  ) %>%
  filter(
    !is.na(inc_cost_M_50),   !is.na(inc_cost_M_500),   !is.na(inc_cost_M_1000),
    !is.na(inc_qalys_50),    !is.na(inc_qalys_500),    !is.na(inc_qalys_1000)
  ) %>%
  mutate(
    icer_label = fmt_icer(icer_500),
    spec_label = case_when(
      cancer_group == "Specificity 0.95"   ~ "All Cancers (Specificity = 0.95)",
      cancer_group == "Specificity 0.98"   ~ "*Base Case* All Cancers (Specificity = 0.98)",
      cancer_group == "Specificity 0.9999" ~ "All Cancers (Specificity = 0.9999)",
      TRUE ~ "All Cancers"
    ),
    pt_label   = paste0(spec_label, "\n", "MCD Test Cost = $500", "\n", icer_label),
    label_low  = "MCD Test Cost = $50",
    label_high = "MCD Test Cost = $1000"
  )

# -----------------------------
# 4) Aesthetics
# -----------------------------
col_map <- c(
  "Specificity 0.95"   = "#6BAED6",
  "Specificity 0.98"   = "#3182BD",
  "Specificity 0.9999" = "#08519C"
)

shape_map <- c(
  "Specificity midpoint ($500)" = 16,
  "With screening (0.98) / without (0.99)" = 15
)

# WTP label placement (stable)
x_wtp <- 8000
y_wtp <- 500

# -----------------------------
# 5) Plot
# -----------------------------
p <- ggplot() +
  geom_hline(yintercept = 0, linewidth = 0.7, color = "black") +
  geom_vline(xintercept = 0, linewidth = 0.7, color = "black") +
  
  # WTP line (force dashed; no linetype mapping -> no warning)
  geom_abline(
    slope = wtp / 1e6,
    linewidth = 1.1,
    color = "black",
    linetype = "dashed"
  ) +
  annotate(
    "text",
    x = x_wtp, y = y_wtp,
    label = "WTP = $100,000/QALY",
    hjust = 0, vjust = -0.75, size = 4
  ) +
  
  # --- SPEC RANGES ---
  geom_errorbar(
    data = plot_df,
    aes(
      x = inc_qalys_500,
      ymin = inc_cost_M_50,
      ymax = inc_cost_M_1000,
      color = cancer_group
    ),
    linewidth = 0.8,
    width = 8000,
    show.legend = FALSE
  ) +
  geom_segment(
    data = plot_df,
    aes(
      x = inc_qalys_50, xend = inc_qalys_1000,
      y = inc_cost_M_500, yend = inc_cost_M_500,
      color = cancer_group
    ),
    linewidth = 0.8,
    show.legend = FALSE
  ) +
  geom_point(
    data = plot_df,
    aes(
      x = inc_qalys_500, y = inc_cost_M_500,
      color = cancer_group,
      shape = "Specificity midpoint ($500)"
    ),
    size = 3.2
  ) +
  
  # Midpoint labels: spec label + ICER
  geom_label_repel(
    data = plot_df,
    aes(
      x = inc_qalys_500, y = inc_cost_M_500,
      label = pt_label,
      color = cancer_group
    ),
    fill = "white",
    label.size = 0,
    size = 4.5,
    box.padding = 0.45,
    point.padding = 0.35,
    min.segment.length = 0,
    seed = 1,
    show.legend = FALSE
  ) +
  
  # $50 / $1000 endpoint labels (vertical only)
  geom_text(
    data = plot_df,
    aes(x = inc_qalys_500, y = inc_cost_M_50, label = label_low, color = cancer_group),
    vjust = 1.2, hjust = 0.8, nudge_x = -10000, size = 4,
    show.legend = FALSE
  ) +
  geom_text(
    data = plot_df,
    aes(x = inc_qalys_500, y = inc_cost_M_1000, label = label_high, color = cancer_group),
    vjust = -0.5, hjust = 0.8, nudge_x = -10000, size = 4,
    show.legend = FALSE
  ) +
  
  # --- SPLIT POINT ---
  geom_point(
    data = df_split,
    aes(
      x = inc_qalys, y = inc_cost_M,
      shape = "With screening (0.98) / without (0.99)"
    ),
    size = 3.2,
    color = "black"
  ) +
  geom_label_repel(
    data = df_split,
    aes(
      x = inc_qalys, y = inc_cost_M,
      label = paste0(point_label, "\n", ce_label)
    ),
    nudge_x = 10000,
    nudge_y = -1000,
    fill = "white",
    label.size = 0,
    size = 4.5,
    box.padding = 0.45,
    point.padding = 0.35,
    min.segment.length = 0,
    seed = 1,
    color = "black",
    show.legend = FALSE
  ) +
  
  # Legends / scales
  scale_color_manual(values = col_map, name = NULL) +
  scale_shape_manual(values = shape_map, name = NULL) +
  
  # Formatting
  scale_x_continuous(labels = label_comma(), expand = expansion(mult = 0.15)) +
  scale_y_continuous(labels = label_comma(), expand = expansion(mult = 0.06)) +
  labs(
    x = "Incremental QALYs",
    y = "Incremental Cost (Millions)"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "none"
  )

print(p)
ggsave("ce_plane_plot_millions.png", p, width = 10, height = 6, dpi = 300)
