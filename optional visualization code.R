# optional visualization
# Complex visual
# Custom missingness plot: color-coded by study group
# Each row = one participant, each column = one neuro test
# Colors show which group the missing value belongs to
# Grey = data present, Red = Acute IPV-BI missing,Blue = Acute Non-IPV-BI missing, Yellow = Chronic IPV-BI missing

# Step 1: reshape the neuro test columns to long format
# pivot_longer() converts wide data (one column per test) to long data
# (one row per participant-test combination), which ggplot2 needs to
# map missingness and group onto the same plot

miss_plot_data <- cohort %>%
  select(participant_id, study_group, all_of(neuro_test_cols)) %>%
  pivot_longer(
    cols      = all_of(neuro_test_cols),
    names_to  = "test",
    values_to = "value"
  ) %>%
  mutate(
    # Create a fill variable that captures both presence/absence and group
    fill_group = case_when(
      !is.na(value)                        ~ "Present",
      study_group == "Acute IPV-BI"        ~ "Acute IPV-BI missing",
      study_group == "Acute Non-IPV-BI"    ~ "Acute Non-IPV-BI missing",
      study_group == "Chronic IPV-BI"      ~ "Chronic IPV-BI missing"
    ),
    fill_group = factor(fill_group,
                        levels = c("Present",
                                   "Acute IPV-BI missing",
                                   "Acute Non-IPV-BI missing",
                                   "Chronic IPV-BI missing")),
    # Clean up test names for display by removing prefixes
    test = str_remove(test, "WAIS_|SDMT_|RAVLT_|TMT_|TOPF_|COWA_"),
    test = factor(test)
  )

# Step 2: define participant order by group so dashed dividers fall
# cleanly between groups
participant_order <- cohort %>%
  arrange(study_group) %>%
  pull(participant_id)

miss_plot_data <- miss_plot_data %>%
  mutate(participant_id = factor(participant_id, levels = participant_order))

# Step 3: calculate y-axis positions for the dashed group divider lines
# These sit between the last participant of one group and the first of the next
n_acute_ipv    <- sum(cohort$study_group == "Acute IPV-BI")
n_acute_nonipv <- sum(cohort$study_group == "Acute Non-IPV-BI")

divider_1 <- n_acute_ipv + 0.5           # between Acute IPV-BI and Acute Non-IPV-BI
divider_2 <- n_acute_ipv + n_acute_nonipv + 0.5  # between Acute Non-IPV-BI and Chronic IPV-BI

# Step 4: build the plot
ggplot(miss_plot_data,
       aes(x = test, y = participant_id, fill = fill_group)) +
  geom_tile(color = "white", linewidth = 0.3) +
  # Dashed horizontal lines separating the three groups
  geom_hline(yintercept = divider_1, linetype = "dashed",
             color = "black", linewidth = 0.7) +
  geom_hline(yintercept = divider_2, linetype = "dashed",
             color = "black", linewidth = 0.7) +
  scale_fill_manual(
    values = c(
      "Present"                 = "grey85",
      "Acute IPV-BI missing"    = "#D62728",   # red
      "Acute Non-IPV-BI missing"= "#1F77B4",   # blue
      "Chronic IPV-BI missing"  = "#FFBF00"    # yellow
    ),
    name = "Data status"
  ) +
  # Group labels on the right side of the plot
  annotate("text", x = Inf, y = n_acute_ipv / 2,
           label = "Acute\nIPV-BI", hjust = -0.1, size = 3, fontface = "bold") +
  annotate("text", x = Inf, y = n_acute_ipv + n_acute_nonipv / 2,
           label = "Acute\nNon-IPV-BI", hjust = -0.1, size = 3, fontface = "bold") +
  annotate("text", x = Inf, y = n_acute_ipv + n_acute_nonipv + 10 / 2,
           label = "Chronic\nIPV-BI", hjust = -0.1, size = 3, fontface = "bold") +
  labs(
    title    = "Missing data in neuropsychological test scores by study group",
    subtitle = "Each row = one participant | Each column = one test | ~3-4% MCAR expected",
    x        = "Neuropsychological test",
    y        = "Participant"
  ) +
  theme_minimal() +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y      = element_blank(),   # participant IDs too cluttered to show
    axis.ticks.y     = element_blank(),
    panel.grid       = element_blank(),
    legend.position  = "bottom",
    plot.margin      = margin(5, 60, 5, 5)  # extra right margin for group labels
  )