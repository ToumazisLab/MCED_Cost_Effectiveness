# ---- libraries -------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(scales)      # nice comma labels (optional)

# ---- 1. key inputs ---------------------------------------------------------
baseline_icer <- -4.660230e+04

# ---- 2. order parameters by overall ICER swing ----------------------------
# (largest absolute range appears at the top of the chart)
ordered <- dsa_results %>% 
  group_by(plot_label) %>% 
  summarise(span = max(Incremental_Utility) - min(Incremental_Utility), .groups = "drop") %>% 
  arrange(desc(span))

dsa_results <- dsa_results %>% 
  mutate(plot_label = factor(plot_label, levels = ordered$plot_label))

# ---- 3. main plot ----------------------------------------------------------
p <- ggplot(dsa_results,
            aes(x = Incremental_Utility, y = plot_label, group = plot_label)) +
  
  # bar for each ± pair (segment between the two ICERs)
  geom_segment(
    data = dsa_results %>% 
      group_by(plot_label) %>% 
      summarise(xmin = min(Incremental_Utility),
                xmax = max(Incremental_Utility),
                .groups = "drop"),
    aes(x = xmin, xend = xmax, yend = plot_label),
    linewidth = 8,    # bar thickness
    color = "#69b3a2"
  ) +
  
  # endpoints (optional – helps when bars overlap the baseline line)
  geom_point(size = 1, colour = "grey20") +
  
  # annotate the endpoints with the associated new_value
  geom_text(aes(label = new_value,
                hjust = ifelse(change_type == "plus", -0.3, 1.3)),   # <- was -0.1 / 1.1
            size = 3.2) +
  
  # baseline ICER
  geom_vline(xintercept = baseline_icer,
             linetype = "dashed", colour = "black") +
  
  # axis & theme tweaks
  scale_x_continuous(labels = comma) +
  labs(x = "Incremental Utility",
       y = NULL) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank())

print(p)
