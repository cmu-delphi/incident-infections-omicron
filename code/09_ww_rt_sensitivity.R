# Rt sensitivity analysis using highest and lowest generation time
# and serial intervals (the same incubation period as before)

#remotes::install_github("dajmcdon/rtestim")
library(rtestim)
library(tidyverse)

set.seed(12345)

ww_adj_df = readRDS("quantiles_ww_adj_df.RDS") %>% mutate(adj_infect = adj_infect_q50)
# Split the data frame by the geo_value column
ww_adj_df_list <- split(ww_adj_df, ww_adj_df$geo_value) 
delay_distns_by_var <- readRDS("delay-distns-byvar.rds") 

GT_omicron_meta = delay_distns_by_var %>% dplyr::filter(para == "GT" & type == "Omicron")
SI_omicron_meta = delay_distns_by_var %>% dplyr::filter(para == "SI" & type == "Omicron")
IP_omicron_meta = delay_distns_by_var %>% dplyr::filter(para == "IP" & type == "Omicron")

# Following extraction of xu-etal-DATA_RAW.xlsx is from https://dajmcdon.github.io/rtestim/articles/delay-distributions.html
data_raw <- readxl::read_excel("xu-etal-DATA_RAW.xlsx") |>
  select(type, para, n = Sample_size, mean, sd, se, median) |>
  filter(!is.na(type)) |>
  mutate(across(-c(type, para), as.numeric)) |> 
  filter(para %in% c("GT", "IP", "SI") & type == "Omicron")

# Test lowest and highest GT and SIs impact on final Rt estimates

# GT - min and max
GT_omicron_meta_all = data_raw %>% dplyr::filter(para == "GT" & type == "Omicron" & mean != 6.84) # Remove intrinsic generation time (only consider realized)

GT_omicron_meta_min = GT_omicron_meta_all |> 
  dplyr::filter(mean == min(mean, na.rm = TRUE)) |> 
  dplyr::mutate(sd = ifelse(is.na(sd), GT_omicron_meta$sd, sd)) |> 
  mutate(shape = mean^2 / sd^2, scale = mean / shape) 

GT_omicron_meta_max = GT_omicron_meta_all |> 
  dplyr::filter(mean == max(mean, na.rm = TRUE)) |> 
  dplyr::mutate(sd = ifelse(is.na(sd), GT_omicron_meta$sd, sd)) |> 
  mutate(shape = mean^2 / sd^2, scale = mean / shape) 

# SI - min and max
SI_omicron_meta_all = data_raw %>% dplyr::filter(para == "SI" & type == "Omicron")

SI_omicron_meta_min = SI_omicron_meta_all |> 
  dplyr::filter(mean == min(mean, na.rm = TRUE)) |> 
  dplyr::mutate(sd = ifelse(is.na(sd), SI_omicron_meta$sd, sd)) |> 
  mutate(shape = mean^2 / sd^2, scale = mean / shape) 

SI_omicron_meta_max = SI_omicron_meta_all |> 
  dplyr::filter(mean == max(mean, na.rm = TRUE)) |> 
  dplyr::mutate(sd = ifelse(is.na(sd), SI_omicron_meta$sd, sd)) |> 
  mutate(shape = mean^2 / sd^2, scale = mean / shape) 

# IP, as before
IP_omicron_meta = delay_distns_by_var %>% dplyr::filter(para == "IP" & type == "Omicron")

#########################################################################################################################################
# Settings for min and max
# To external user: Use the below line to turn on which want to produce plot for
min_res <- "TRUE"

if(min_res == "TRUE"){
  GT_omicron_meta_to_use = GT_omicron_meta_min
  SI_omicron_meta_to_use = SI_omicron_meta_min
}else{
  GT_omicron_meta_to_use = GT_omicron_meta_max
  SI_omicron_meta_to_use = SI_omicron_meta_max
}
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
rt_list_res = lapply(ww_filtered_list, function(x) cv_estimate_rt(x$adj_infect, nsol = 20, dist_gamma = c(GT_omicron_meta_to_use$shape, GT_omicron_meta_to_use$scale)))
# Plot Rt for a state like NY
plot(rt_list_res$NY, which_lambda = "lambda.1se")

# Rt estimation on cases
rt_list_res_cases = lapply(ww_filtered_list, function(x) cv_estimate_rt(x$case_count_7d_av, nsol = 20, dist_gamma = c(SI_omicron_meta_to_use$shape, SI_omicron_meta_to_use$scale)))
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

# Plot Rt for lambda.min for all 7 states
# Just infections
ggplot(rt_df, aes(x = time_value, y = Rt, color = state)) +
  geom_hline(yintercept = 1, color = "grey") +  # Grey dashed line at y = 1
  geom_line(size = 1) +
  facet_wrap(~ state, scales = "free_y", ncol = 2) +
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
  geom_line(data = rt_df, aes(x = time_value, y = Rt, color = "Infection-based"), size = 0.8) +  # Infections
  geom_line(data = rt_df_cases, aes(x = time_value, y = Rt, color = "Case-based"), size = 0.8) +  # Cases
  facet_wrap(~ state, ncol = 2) + # scales = "free_y", # Facet by state with free y-axis scaling 
  scale_x_date(name = "", date_breaks = "3 month", date_labels = "%b %Y", expand = c(0, 0)) +
  scale_color_manual(values = c("Infection-based" = "midnightblue", "Case-based" = "darkorange2"), name = NULL) +  
  coord_cartesian(ylim = c(0.8, 1.2)) +  # Set y-axis limit without clipping data
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