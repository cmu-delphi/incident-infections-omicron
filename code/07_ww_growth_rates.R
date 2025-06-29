# Epiprocess Growth Rates Comparison
library(epiprocess)
library(tidyverse)
library(RColorBrewer)
library(grid)
library(gtable)

ww_adj_df = readRDS("quantiles_ww_adj_df.RDS")

# Format ww_adj_df_list into dataframe that is compatible with epiprocess Estimate Growth Signals vignette
# https://cmu-delphi.github.io/epiprocess/articles/growth_rate.html
x <- ww_adj_df %>%
  select(geo_value, time_value, infections = adj_infect_q50, cases = case_count_7d_av) %>%
  mutate(cases = if_else(cases == 0, 1, cases)) %>% # If 0 for cases for a day, pad with 1 so don't get -Inf or Inf growth rate estimates around it
  arrange(geo_value, time_value) %>%
  as_epi_df()

# Relative change
# The default method is “rel_change”, which is the simplest way to estimate growth rates.
# The default bandwidth is h = 7, which for daily data, considers the relative change in a signal over adjacent weeks.
x <- x %>%
  group_by(geo_value) %>%
  mutate(infections_gr1 = growth_rate(time_value, infections, method = "rel_change", log_scale = TRUE, h = 14),
         cases_gr1 = growth_rate(time_value, cases, method = "rel_change", log_scale = TRUE, h = 14))

head(x, 10)

# We can visualize these growth rate estimates by plotting the signal values
# and highlighting the periods in time for which the relative change is above 1% (in red) and below -1% (in blue),
# faceting by geo value.
theme_set(theme_bw())

upper <- 0.01
lower <- -0.01

# As a more direct visualization, we plot the estimated growth rates themselves,
# faceting the curves for the states.
manual_breaks <- as.Date(c("2022-07-01", "2022-10-01", "2023-01-01"))

ggplot(x, aes(x = time_value)) +
  geom_hline(yintercept = 0, linetype = 1, col = "darkgrey", size = 0.8) +
  geom_line(aes(y = infections_gr1, color = "Infection-based"), size = 0.8) +
  # Plot cases growth rate line (solid)
  geom_line(aes(y = cases_gr1, color = "Case-based"), size = 0.8) +
  # Customize color scale for Infections and Cases
  scale_color_manual(values = c("Infection-based" = "midnightblue", "Case-based" = "darkorange2"), 
                     guide = guide_legend(title = NULL)) + # Remove the legend title here
  facet_wrap(vars(geo_value), scales = "fixed", ncol = 2) +
  scale_x_date(name = "", breaks = manual_breaks, date_labels = "%b %Y", expand = c(0, 0)) +
  labs(x = "Date", y = "Growth rate", col = "State") +
  # Use coord_cartesian to set y-axis limits
  coord_cartesian(ylim = c(-0.1, 0.1)) + 
  theme(legend.position = c(0.75, 0.08)) # Legend in the bottom right corner

####################################################################################################################################
#  Make plots with subvariant proportions in the background
var_props_omic_sub = readr::read_csv("seq_df_post_decon_omicron.csv")

`%nin%` <- Negate(`%in%`)

df_long <- var_props_omic_sub %>%
  gather(key = "Variant", value = "Prop", -State, -Date) %>%  # gather variant columns, keeping State and Date
  # Only need Omicron sub-variants
  filter(Variant %nin% c("Alpha", "Beta", "Epsilon", "Iota", "Gamma", "Delta", "Other")) %>%
  # Only need dates >= first date for the adjusted infections
  filter(Date >= ww_adj_df$time_value[1]) %>% #
  filter(State %in% c("CA", "NY", "NC", "MN", "HI", "WA", "NV")) %>%
  arrange(State, Date) %>%  # Arrange by state and date
  filter(Date <= max(ww_adj_df$time_value)) %>%
  mutate(geo_value = State) # For keeping this the same in below plotting using geom_area

# Rename some columns of the df_long for easy plotting
df_long$time_value = df_long$Date

# Remove (Omicron) from after subvariant name
df_long$Variant <- gsub(" \\(Omicron\\)", "", df_long$Variant)

####################################################################################################################################
# Now plot the variant proportions in the background of the growth rate plot for cases and infections

variant_colors <- viridis::viridis(14, begin = 0.3, end = 1)

# Only have the variants in the legend that have a Prop of > 10%
geom_area_breaks <- unique(df_long$Variant)[which(unique(df_long$Variant) %in% unique(df_long$Variant[df_long$Prop > 0.1]))] 

####################################################################################################################################
# As a more direct visualization, we plot the estimated growth rates themselves,
# faceting the curves for the states.
manual_breaks <- as.Date(c("2022-07-01", "2022-10-01", "2023-01-01"))

# 1. Background: variant proportions
p_background <- ggplot(df_long, aes(Date, y = Prop, fill = Variant, group = interaction(geo_value, Variant))) +
  geom_area(position = "fill", alpha = 0.7) +
  scale_fill_manual(values = setNames(variant_colors, unique(df_long$Variant)), breaks = geom_area_breaks) +
  facet_wrap(~ geo_value, ncol = 2) +
  scale_x_date(name = "", date_breaks = "3 month", date_labels = "%b %Y", expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  theme_bw() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.justification = c(0, 1), 
    legend.position = "top",
    axis.text.x = element_text(size = 9.5),
    axis.title.x = element_text(size = 11.5),
    strip.text = element_text(size = 8, face = "bold")
  )

# 2. Overlay plot: growth rate lines (y-axis on left and no background or x-axis)
p_lines <- ggplot(x, aes(x = time_value)) +
  geom_hline(yintercept = 0, linetype = 1, col = "darkgrey", size = 0.8) +
  geom_line(aes(y = infections_gr1, color = "Infection-based"), size = 0.8) +
  geom_line(aes(y = cases_gr1, color = "Case-based"), size = 0.8) +
  scale_color_manual(values = c("Infection-based" = "midnightblue", "Case-based" = "darkorange2"), 
                     guide = guide_legend(title = NULL)) + 
  facet_wrap(vars(geo_value), scales = "fixed", ncol = 2) +
  scale_x_date(name = "", date_breaks = "3 month", date_labels = "%b %Y", expand = c(0, 0)) +
  labs(x = "Date", y = "Growth rate", col = "State") +
  # Use coord_cartesian to set y-axis limits
  coord_cartesian(ylim = c(-0.1, 0.1)) + 
  theme(panel.background = element_rect(fill = "transparent", color = NA),  # transparent background
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid = element_blank(), # Remove all grid lines 
        axis.title.x = element_text(margin = margin(t = 10)), # extra space at bottom of plot
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.justification = c(1, 1), 
        legend.position = "top",
        legend.box.margin = margin(10, 10, 10, 10),
        strip.text = element_text(size = 8, face = "bold"))

# Convert both ggplots to grobs
g1 <- ggplotGrob(p_background)
g2 <- ggplotGrob(p_lines)

# Ensure same dimensions: copy widths and heights
g1$widths <- g2$widths
g1$heights <- g2$heights

# Set the output file and size
pdf("rel-growth-rates.pdf", width = 10, height = 7)  

# Overlay the grobs using grid
grid.newpage()
grid.draw(g1)        # draw base 
grid.draw(g2)        # overlay lines

dev.off()