# Growth Rates Sensitivity Analysis - Vary the bandwidths from 7, 14, 21, 28 and produce growth rate plots of each
library(epipredict)
library(tidyverse)

ww_adj_df <- readRDS("quantiles_ww_adj_df.RDS")

# Pre-process the data in one go
x_base <- ww_adj_df %>%
  select(geo_value, time_value, infections = adj_infect_q50, cases = case_count_7d_av) %>%
  mutate(cases = if_else(cases == 0, 1, cases)) %>% 
  arrange(geo_value, time_value) %>%
  as_epi_df()

# Thresholds for growth rate highlights
upper <- 0.01
lower <- -0.01

# Bandwidths to loop through
bandwidths <- c(7, 14, 21, 28)

# Date breaks for plotting
manual_breaks <- as.Date(c("2022-07-01", "2022-10-01", "2023-01-01"))
for (h in bandwidths) {
  x <- x_base %>%
    group_by(geo_value) %>%
    mutate(
      infections_gr = growth_rate(time_value, infections, method = "rel_change", log_scale = TRUE, h = h),
      cases_gr = growth_rate(time_value, cases, method = "rel_change", log_scale = TRUE, h = h)
    )
  
  p <- ggplot(x, aes(x = time_value)) +
    geom_line(aes(y = infections_gr, color = "Infection-based"), size = 1) +
    geom_line(aes(y = cases_gr, color = "Case-based"), size = 1) +
    geom_hline(yintercept = upper, linetype = 2, col = "darkgrey", size = 0.8) +
    geom_hline(yintercept = lower, linetype = 2, col = "darkgrey", size = 0.8) +
    scale_color_manual(values = c("Infection-based" = "midnightblue", "Case-based" = "darkorange2"),
                       guide = guide_legend(title = NULL)) +
    facet_wrap(vars(geo_value), scales = "free_y", ncol = 2) +
    scale_x_date(name = "", breaks = manual_breaks, date_labels = "%b %Y", expand = c(0, 0)) +
    labs(
      x = "Date", y = "Growth rate", col = "State"
    ) +
    coord_cartesian(ylim = c(-0.1, 0.1)) +
    theme_bw() +
    theme(legend.position = c(0.75, 0.08))
  
  # Save the plot
  ggsave(
    filename = paste0("growth_rate_plot_h", h, ".pdf"),
    plot = p
  )
}
