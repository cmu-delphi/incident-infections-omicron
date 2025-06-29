# Rt estimation

#remotes::install_github("dajmcdon/rtestim")
library(rtestim)
library(tidyverse)
library(grid)
library(gtable)
library(RColorBrewer)
library(cowplot)

set.seed(12345)

ww_adj_df = readRDS("quantiles_ww_adj_df.RDS") %>% mutate(adj_infect = adj_infect_q50)
# Split the data frame by the geo_value column
ww_adj_df_list <- split(ww_adj_df, ww_adj_df$geo_value) 
delay_distns_by_var <- readRDS("delay-distns-byvar.rds") 

GT_omicron_meta = delay_distns_by_var %>% dplyr::filter(para == "GT" & type == "Omicron")
SI_omicron_meta = delay_distns_by_var %>% dplyr::filter(para == "SI" & type == "Omicron")
IP_omicron_meta = delay_distns_by_var %>% dplyr::filter(para == "IP" & type == "Omicron")

###############################################################################################################################################
# Define the start and end date for Rt estimation
start_date <- as.Date("2022-05-01")
end_date <- as.Date("2023-02-01")
###############################################################################################################################################

# Use lapply to apply the filtering to each data frame in the list
ww_filtered_list <- lapply(ww_adj_df_list, function(df) {
  df %>%
    filter(time_value >= start_date & time_value <= end_date)
})

# Rt estimation on infections
rt_list_res = lapply(ww_filtered_list, function(x) cv_estimate_rt(x$adj_infect, nsol = 20, dist_gamma = c(GT_omicron_meta$shape, GT_omicron_meta$scale)))
# Plot Rt for a state like NY
plot(rt_list_res$NY, which_lambda = "lambda.1se")

# Rt estimation on cases
rt_list_res_cases = lapply(ww_filtered_list, function(x) cv_estimate_rt(x$case_count_7d_av, nsol = 20, dist_gamma = c(SI_omicron_meta$shape, SI_omicron_meta$scale)))
# Plot Rt for a state like NY
plot(rt_list_res_cases$NY, which_lambda = "lambda.1se")


###############################################################################################################################################
# Function to extract Rt values for lambda.min for a given state
extract_rt_for_lambda_1se <- function(state_data) {
  # Extract Rt matrix and corresponding lambda values
  rt_matrix <- state_data$full_fit$Rt  
  lambda_values <- state_data$lambda  
  lambda_1se <- state_data$lambda.1se  # Best lambda value according to lambda.1se 

  # Find the index of the lambda.1se
  lambda_1se_index <- which.min(abs(lambda_values - lambda_1se))

  # Extract the Rt values for lambda.1se (the column corresponding to lambda_1se_index)
  rt_values <- rt_matrix[, lambda_1se_index]

  return(rt_values)
}

###############################################################################################################################################
# For infections
# Create a data frame for all states
state_names <- names(rt_list_res)
rt_data_list <- lapply(state_names, function(state) {
  # Extract Rt for lambda.1se for each state
  rt_values <- extract_rt_for_lambda_1se(rt_list_res[[state]])

  # Create a data frame with dates and the extracted Rt values
  data.frame(
    time_value = seq.Date(from = ww_filtered_list$CA$time_value[1], to = ww_filtered_list$CA$time_value[length(ww_filtered_list$CA$time_value)], by = "day"),
    Rt = rt_values,
    state = state
  )
})

# For cases (same as above) except with date shift by mean incubation period
# Create a data frame for all states
state_names <- names(rt_list_res_cases)
rt_data_list_cases <- lapply(state_names, function(state) {
  # Extract Rt for lambda.min for each state
  rt_values <- extract_rt_for_lambda_1se(rt_list_res_cases[[state]])
  
  # Dates before incubation period shift:
  dates_before_inc_shift <- seq.Date(from = ww_filtered_list$CA$time_value[1], to = ww_filtered_list$CA$time_value[length(ww_filtered_list$CA$time_value)], by = "day")
  dates_before_inc_shift <- dates_before_inc_shift - round(IP_omicron_meta$mean)
  
  # Create a data frame with dates and the extracted Rt values
  data.frame(
    time_value = dates_before_inc_shift,
    Rt = rt_values,
    state = state
  )
})

# Combine all the state data into a single data frame
rt_df <- bind_rows(rt_data_list) # infections Rt df
rt_df_cases <- bind_rows(rt_data_list_cases) #cases Rt df
# Save RDS files of each
saveRDS(rt_df, file = "rt_df.RDS")
saveRDS(rt_df_cases, file = "rt_df_cases.RDS")

# Plot Rt for lambda.min for all 7 states
# Just infections
ggplot(rt_df, aes(x = time_value, y = Rt, color = state)) +
  geom_hline(yintercept = 1, color = "grey") +  # Grey dashed line at y = 1
  geom_line(size = 1) +
  facet_wrap(~ state, scales = "free_y", ncol = 2) +  # Facet by state with free y-axis scaling
  scale_x_date(name = "Date", date_breaks = "3 month", date_labels = "%b %Y") +
  theme_bw() +
  ylab("Estimated Rt") +
  theme(
    axis.text.x = element_text(size = 9.5),
    axis.text.y = element_text(size = 9.5),
    axis.title.x = element_text(size = 11.5),
    axis.title.y = element_text(size = 11.5),
    legend.position = "none",
    strip.text = element_text(size = 8, face = "bold")
  )

# Rt from infections and cases on the same plot
ggplot() +
  geom_hline(yintercept = 1, color = "grey") +  # Grey dashed line at y = 1
  geom_line(data = rt_df, aes(x = time_value, y = Rt, color = "Infection Rt"), size = 0.8) +  # Infections
  geom_line(data = rt_df_cases, aes(x = time_value, y = Rt, color = "Case Rt"), size = 0.8) +  # Cases
  facet_wrap(~ state, ncol = 2) + # scales = "free_y", 
  scale_x_date(name = "Date", date_breaks = "3 month", date_labels = "%b %Y") +
  scale_color_manual(values = c("Infection Rt" = "midnightblue", "Case Rt" = "darkorange2")) +  
  coord_cartesian(ylim = c(0.8, 1.2)) +  
  theme_bw() +
  ylab("Estimated Rt") +
  theme(
    axis.text.x = element_text(size = 9.5),
    axis.text.y = element_text(size = 9.5),
    axis.title.x = element_text(size = 11.5),
    axis.title.y = element_text(size = 11.5),
    legend.position = "top",  # Legend to differentiate Infections vs Cases
    strip.text = element_text(size = 8, face = "bold")
  )

####################################################################################################################################
#  Make plots with variant proportions in the background
var_props_omic_sub = readr::read_csv("seq_df_post_decon_omicron.csv")

`%nin%` <- Negate(`%in%`)

df_long <- var_props_omic_sub %>%
  gather(key = "Variant", value = "Prop", -State, -Date) %>%  # gather variant columns, keeping State and Date
  # Only need Omicron sub-variants
  filter(Variant %nin% c("Alpha", "Beta", "Epsilon", "Iota", "Gamma", "Delta", "Other")) %>%
  # Only need dates >= first date for the adjusted infections
  filter(Date >= rt_df_cases$time_value[1]) %>% # use cases first time value (not infections) because those dates have been shifted back by incubation period
  filter(State %in% c("CA", "NY", "NC", "MN", "HI", "WA", "NV")) %>%
  arrange(State, Date) %>%  # Arrange by state and date
  filter(Date <= max(ww_adj_df$time_value)) %>%
  mutate(state = State) # For keeping this the same in below plotting using geom_area
####################################################################################################################################
# Transformations
# Define transformation functions
trans <- function(x, from_range, to_range) {
  (x - from_range[1]) / diff(from_range) * diff(to_range) + to_range[1]
}
inv_trans <- function(x, from_range, to_range) {
  (x - to_range[1]) / diff(to_range) * diff(from_range) + from_range[1]
}

# Choose a sensible range for Rt lines
range_rt <- range(c(min(rt_df$Rt), max(rt_df$Rt)))

# Add a normalized Rt variable to plot on the fill-scale (0–1)
rt_df$Rt_scaled <- trans(rt_df$Rt, from_range = range_rt, to_range = c(0, 1))
rt_df_cases$Rt_scaled <- trans(rt_df_cases$Rt, from_range = range_rt, to_range = c(0, 1))
####################################################################################################################################
# Rename some columns of the df_long for easy plotting
df_long$time_value = df_long$Date

# Remove (Omicron) from after subvariant name
df_long$Variant <- gsub(" \\(Omicron\\)", "", df_long$Variant)
####################################################################################################################################
# Now plot the variant proportions in the background of the Rt plot for cases and infections

variant_colors <- viridis::viridis(14, begin = 0.3, end = 1)

# Only have the variants in the legend that have a Prop of > 10%
geom_area_breaks <- unique(df_long$Variant)[which(unique(df_long$Variant) %in% unique(df_long$Variant[df_long$Prop > 0.1]))] 

# 1. Background: variant proportions
p_background <- ggplot(df_long, aes(Date, y = Prop, fill = Variant, group = interaction(state, Variant))) +
  geom_area(position = "fill", alpha = 0.7) +
  scale_fill_manual(values = setNames(variant_colors, unique(df_long$Variant)), breaks = geom_area_breaks) +
  facet_wrap(~ state, ncol = 2) +
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
    strip.text = element_text(size = 8, face = "bold")
  )

# 2. Overlay plot: growth rate lines (y-axis on left and no background or x-axis)
p_lines <- ggplot() +
  geom_hline(yintercept = 1, color = "grey") +  # Grey dashed line at y = 1
  geom_line(data = rt_df, aes(x = time_value, y = Rt, color = "Infection-based Rt"), size = 0.8) +  # Infections
  geom_line(data = rt_df_cases, aes(x = time_value, y = Rt, color = "Case-based Rt"), size = 0.8) +  # Cases
  facet_wrap(~ state, ncol = 2) + # scales = "free_y", 
  scale_x_date(name = "", date_breaks = "3 month", date_labels = "%b %Y", expand = c(0, 0)) +
  scale_color_manual(name = "", values = c("Infection-based Rt" = "midnightblue", "Case-based Rt" = "darkorange2")) + 
  coord_cartesian(ylim = c(0.8, 1.2)) +  
  theme_bw() +
  ylab("Estimated Rt") +
  theme(
    panel.background = element_rect(fill = "transparent", color = NA),  # transparent background
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.grid = element_blank(), # Remove all grid lines 
    axis.title.x = element_text(margin = margin(t = 10), size = 11.5),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.justification = c(1, 1), 
    legend.position = "top",
    legend.box.margin = margin(10, 10, 10, 10), 
    strip.text = element_text(size = 8, face = "bold")
  )

# Convert both ggplots to grobs
g1 <- ggplotGrob(p_background)
g2 <- ggplotGrob(p_lines)

# Ensure same dimensions: copy widths and heights
g1$widths <- g2$widths
g1$heights <- g2$heights

# Set the output file and size
pdf("rt-lambdamin.pdf", width = 10, height = 7)  

# Overlay the grobs using grid
grid.newpage()
grid.draw(g1)        # draw base (variant area)
grid.draw(g2)        # overlay lines exactly

dev.off()

#########################################################################################################################
# Function for plotting all lambda values for each state
plot_all_lambda_rt <- function(rt_list_res, states = c("CA", "WA", "NY", "NV", "NC", "MN", "HI"), cases = FALSE) {

  # Create a list to store data frames for each state
  all_data <- list()

  # Loop through each state and create the data frame for each state
  for (state in states) {

    # Extract Rt matrix and lambda values
    Rt_matrix <- rt_list_res[[state]]$full_fit$Rt  
    lambda_values <- rt_list_res[[state]]$lambda    
    time_points <- rt_list_res[[state]]$full_fit$x  

    # Check if Rt_matrix and lambda_values have the expected dimensions
    if (ncol(Rt_matrix) != length(lambda_values)) {
      stop(paste("Error: The number of lambda values does not match the number of columns in Rt_matrix for state", state))
    }

    # Get the number of lambdas (should be same as the number of columns in Rt_matrix)
    k <- length(lambda_values)
    
    # Dates before incubation period shift:
    if(cases == TRUE){
      days_dates <- seq.Date(from = start_date, to = end_date, by = "day")
      days_dates <- days_dates - round(IP_omicron_meta$mean)
    } else{
      days_dates <- seq.Date(from = start_date, to = end_date, by = "day")
    }
    
    # Reshape Rt_matrix into a long format
    df <- data.frame(
      Rt = as.vector(Rt_matrix),   
      lambda = rep(lambda_values, each = nrow(Rt_matrix)), 
      time_value = rep(days_dates, k)   
    )

    # Add a column for the state name
    df$State <- state

    # Add the data frame to the list
    all_data[[state]] <- df
  }

  # Combine all the state data into one data frame
  all_data_df <- do.call(rbind, all_data)

  # Plot using ggplot2
  ggplot(all_data_df, aes(x = .data$time_value, y = .data$Rt, colour = .data$lambda, group = .data$lambda)) +
    geom_line() +
    geom_hline(yintercept = 1, color = "grey") +  # Add a horizontal line at y = 1
    ylab("Estimated Rt") +
    xlab("Date") +
    theme_bw() +
    scale_colour_viridis_c(trans = "log10", name = "Lambda") +  # Ensure lambda is treated numerically
    facet_wrap(~ State, scales = "free_y", ncol = 2) + 
    theme(
      axis.text.x = element_text(size = 9.5),
      axis.title.y = element_text(size = 11.5),
      axis.text.y = element_text(size = 9.5),
      strip.text = element_text(size = 10, face = "bold")
    )
}

# Use the function to plot Rt for each lambda value for infections
plot_all_lambda_rt(rt_list_res)
# And then for cases
plot_all_lambda_rt(rt_list_res_cases, cases = TRUE)


