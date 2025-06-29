library(tidyverse)
library(mgcv)
# pak::pkg_install("cmu-delphi/epipredict@main")
# library(epipredict)
# remotes::install_github("cmu-delphi/epiprocess@v0.7.0", force = T) # Older version where no .window_size
library(epiprocess)
library(rtestim) 
library(tidyr)
library(covidcast)
library(Matrix)
library(lubridate)
library(RColorBrewer)
library(grid)
library(gtable)

# Load data

start_date = as.Date("2020-06-01") 
end_date = as.Date("2023-02-01") 
decon_start_date = as.Date("2020-03-01")

ww_data <- readr::read_csv("wastewater_by_county.csv")
ww_data = ww_data %>% filter(sampling_week <= end_date)
# Filter out all non-us states
ww_data = ww_data %>% filter(state %in% state.abb) 

abvs = read.table("us_states_abbr_list.txt", header = FALSE)
length(unique(abvs$V1)) # 56, only get 50 states
abvs = abvs %>% filter(V1 %in% state.abb) 

fips_tab = readr::read_csv("fips_table.txt")

# %Change US code to 0
impute_us <- function(x){
  if(x == 'US'){
    return(0)
  } else{
    return(x)
  }
}

fips_tab$location = sapply(fips_tab$location, function(x) impute_us(x))

ww_data$fipscode = as.character(ww_data$fipscode)

ww_data_fips = as.numeric(ww_data$fipscode)

second_grp <- as.numeric(fips_tab$location)

fips_idx = c()
for(w in 1:length(ww_data_fips)){
  if(ww_data_fips[w] %in% second_grp){
    fips_idx = c(fips_idx, which(second_grp == ww_data_fips[w])) # extract element of second_grp (fips_codes) where ww_data_fips[w] shows up
  } else{
    fips_idx = c(fips_idx, 0)
  }
}

# Get indices of states according to abvs
abvs_idx = c()
for(k in ww_data$state){
  if(k %in% as.vector(abvs$V1)){
    abvs_idx = c(abvs_idx, which(abvs == k)) # where is k in abvs (index of k in abvs)
  }
}

# ww_ts matrix
ww_ts = matrix(0, nrow = length(sort(unique(ww_data$sampling_week))), ncol = length(unique(as.vector(abvs$V1))))
ww_ts = data.frame(ww_ts)
rownames(ww_ts) = as.character(sort(unique(ww_data$sampling_week)))
colnames(ww_ts) = abvs$V1

# total county pop matrix for the populations used to make each element of ww_ts
tot_county_pop = matrix(0, nrow = length(sort(unique(ww_data$sampling_week))), ncol = length(unique(as.vector(abvs$V1))))
tot_county_pop = data.frame(tot_county_pop)
rownames(tot_county_pop) = as.character(sort(unique(ww_data$sampling_week)))
colnames(tot_county_pop) = abvs$V1


for(i in 1:length(ww_data$sampling_week)){
  # fill in ww_ts
  ww_ts[as.character(ww_data$sampling_week[i]), abvs_idx[i]] = ww_ts[as.character(ww_data$sampling_week[i]), abvs_idx[i]]  +
    ww_data$effective_concentration_rolling_average[i] * fips_tab$population[fips_idx[i]]

  # fill in tot_county_pop
  tot_county_pop[as.character(ww_data$sampling_week[i]), abvs_idx[i]] = tot_county_pop[as.character(ww_data$sampling_week[i]), abvs_idx[i]] +
    fips_tab$population[fips_idx[i]]
}

# Divide each element in ww_ts by corresponding tot_county_pop
ww_ts = ww_ts / tot_county_pop

# Set any ww_data data to NA if <=0
ww_ts[ww_ts <= 0] = NA

# Convert dataframe wide to long (ie. each state's ww values stacked on top of each other from first week to last week
ww_ts_long <- gather(ww_ts, state, ww, WA:WV)
ww_ts_long <- cbind(date = rep(rownames(ww_ts), times = 50), ww_ts_long)
ww_ts_long$date = as.Date(ww_ts_long$date)

ww_ts_long = ww_ts_long %>% rename(geo_value = state, time_value = date) %>% filter(!is.nan(ww) & !is.na(ww))


# Convert tot_county_pop to long
tot_county_pop_long <- gather(tot_county_pop, state, tot_c_pop, WA:WV)
tot_county_pop_long <- cbind(date = rep(rownames(tot_county_pop), times = 50), tot_county_pop_long)
tot_county_pop_long$date = as.Date(tot_county_pop_long$date)

############################################################################################################################################
# Investigation into ww_ts_long

# Plot of raw ww estimates for state
# Check just one state
state_abb = "CA"
ggplot(ww_ts_long %>% filter(geo_value == state_abb), aes(time_value, y = ww)) +
  geom_point() +
  ggtitle(paste0("Raw ww estimates for ", state_abb, " according to Biobot WW data")) +
  ylab("Average effective genome copies per mL of ww\n across the observed state counties")

# Convert to epi_df and rolling centered median to smooth
ww_ts_long_final = ww_ts_long %>%
  as_epi_df() %>%
  group_by(geo_value) %>%
  epi_slide(~ median(.x$ww, na.rm = TRUE), .window_size = as.difftime(9, units = "weeks"), .align = "center", .new_col_name = "ww_sm") %>%
  ungroup()

state_abb = "GA"
ggplot(ww_ts_long_final %>% filter(geo_value == state_abb), aes(time_value, y = ww_sm)) +
  geom_point() +
  ggtitle(paste0("Adj ww estimates for ", state_abb, " according to Biobot WW data")) +
  ylab("Average effective genome copies per mL of wastewater in the state")

###################################################################################################################################
# Shedding rate cutoff date
sr_cutoff_date = as.Date("2022-02-16") 

# Obtain df with our raw unadjusted infections for each sample state from start_date to end_date
state_unadj_list = vector(mode = "list", length = length(samp_states))
names(state_unadj_list) <- samp_states
vmix_list = vector(mode = "list", length = length(samp_states))
names(vmix_list) <- samp_states
for(state_abb in samp_states){
  
  setwd(paste0("/data/state-deconvolved-case-data/", state_abb))
  
  unadj_infect_df = read_rds("final-thetas-df.rds")

  unadj_infect_df_day = unadj_infect_df %>% group_by(time_value) %>% summarise(infect = sum(infect))
  unadj_infect = unadj_infect_df_day$infect

  unadj_infect = unadj_infect[(start_date - decon_start_date + 1):(end_date - decon_start_date + 1)]

  state_unadj_list[[state_abb]] = data.frame(geo_value = state_abb, unadj_infect = unadj_infect, time_value = seq(from = start_date, to = end_date, by = "1 day"))

  # Separate application
  # Convert unadj_infect_df to wider to esimate prop of variants in circulation
  unadj_infect_df_wider = pivot_wider(unadj_infect_df, names_from = "variant", values_from = "infect")
 
  # Convert frequencies of deconvolved cases to proportions of total number of sequences (not cases),
  # over time, that fall into defined variant groups (so now new_row should sum to 1)
  unadj_infect_df_wider = unadj_infect_df_wider %>%
    ungroup() %>% 
    rowwise() %>%
    mutate(total = sum(c_across(-c(time_value, geo_value)))) %>%
    mutate(across(Other:Iota, function(x) x/total)) %>%
    select(-total) %>%
    ungroup()

  # List the variant columns want to zero out (excluding Omicron)
  variant_cols <- c("Other", "Alpha", "Delta", "Beta", "Epsilon", "Gamma", "Iota")
  
  # Set Omicron = 1 and others variant cols = 0 for each row from cutoff_date and beyond
  # To ensure that Omicron dominates (else there could be possible boundary issues with allotting proportions when reported case counts are close to 0).
  unadj_infect_df_wider <- unadj_infect_df_wider %>%
    mutate(across(-c(time_value, geo_value), ~ if_else(time_value >= sr_cutoff_date, 0, .))) %>%
    mutate(Omicron = if_else(time_value >= sr_cutoff_date, 1, Omicron))  
  
  vmix_list[[state_abb]] = unadj_infect_df_wider
}
vmix <- do.call("rbind", vmix_list) # combine all dfs into a data.frame

# Unadjusted infections string of 0s very long (more than 7 days)
# then 7dav doesn't work very well on infections

# Centered 7 day rolling average for each unadj_infect
unadj_window = 7
unadj_df = bind_rows(state_unadj_list) %>%
  as_epi_df() %>%
  group_by(geo_value) %>%
  epi_slide(~ mean(.x$unadj_infect, na.rm = TRUE), .window_size = unadj_window, .align = "center", .new_col_name = "unadj_infect_7dav")

# Use linear interpolation where 0 values in unadj_infect
unadj_df$unadj_infect_no_0 = ifelse(unadj_df$unadj_infect == 0, NA, unadj_df$unadj_infect)
unadj_df$unadj_infect_no_0 = zoo::na.approx(unadj_df$unadj_infect_no_0, rule = 2)

# Now centered 7 day rolling average for each unadj_infect_no_0
unadj_df = unadj_df %>%
  group_by(geo_value) %>%
  epi_slide(~ mean(.x$unadj_infect_no_0, na.rm = TRUE), .window_size = unadj_window, .align = "center", .new_col_name = "unadj_infect_no_0_7dav")

#############################################################################################################################################
# Total infectiousness at each time point - Set-up gamma delay distribution for viral shedding into ww 

# Set gamma shedding (days) parameters
# Source for peak viral shedding mean and sd: https://pmc.ncbi.nlm.nih.gov/articles/PMC11250817/
# Aligns well with CFA and is for Omicron

# Okada and Nishiura (2024)
Okada = "TRUE"
if(Okada == "TRUE"){
  # Use mean and sd given in paper
  mean_days <- 3.46
  sd_days <- 3.0
  var_days <- sd_days^2
  
}else{ 
# Huisman et al (2022)
  mean_days <- 5
  sd_days <- 0.5
  var_days <- sd_days^2
}
# Create file suffix using above mean sd pair, e.g., "5.2-2.9"
suffix_d <- sprintf("sens-%.1f-%.1f", mean_days, sd_days)
suffix_u <- gsub("-", "_", suffix_d)

# Parameters for gamma distribution (MOM)
shape_gam <- mean_days^2 / var_days
scale_gam <- var_days / mean_days

# Discretize from day 0 to 17
days_gam <- 0:17

# Calculate discrete probabilities: P(t) = P(t ≤ X < t+1)
disc_probs <- pgamma(days_gam + 1, shape = shape_gam, scale = scale_gam) - pgamma(days_gam, shape = shape_gam, scale = scale_gam)  
# plot(disc_probs) # Check plot looks reasonable

# Normalize to ensure the total sums to 1
ga_probs <- disc_probs / sum(disc_probs)
# sum(ga_probs) # Check sums to 1
######################################################################################################################################################
# lm() regression using region (not state) ww, prop circ, and infection data
# w(t) = case(t)*a(t)*\sum_{i\geq 0} prev_var_i(t)*Cd_i
# w(t) = infection(t)*\sum_{i\geq 0} prev_var_i(t)*Cd_i
# Goal: Estimate Cd_i for delta and omicron (during period of delta and omicron mixing)

###################################################################################################################################
# Estimate shedding rate coefficients
# During omicron and delta mixing
`%ni%` = Negate(`%in%`)

# Load data

# Load pop_df
pop_df = readRDS("pop_df.RDS")

# Reshape pop_df from wide to long format
pop_long <- pop_df %>%
  tidyr::pivot_longer(
    cols = starts_with("population_"),
    names_to = "year",
    values_to = "state_pop"
  ) %>%
  mutate(year = as.numeric(gsub("population_", "", year))) 

# pop_used = "population_2021" 

# Load adjusted and unadjusted infections df
adj_df_list = readRDS("gp_adj_df_list_Mar20.RDS")
adj_df = bind_rows(adj_df_list) 

# Load region ww data

ww_data_region <- readr::read_csv("wastewater_by_region.csv")
ww_data_region = ww_data_region %>% filter(sampling_week <= end_date)
# Filter out all non-US states
ww_data_region = ww_data_region %>% filter(region %in% c("Northeast", "Midwest", "South", "West"))

# Make a df that says what region each state belongs to
us_west_states <- data.frame(region = "West", geo_value = c('AK', 'AZ', 'CA', 'CO', 'HI', 'ID', 'MT', 'NM', 'NV', 'OR', 'UT', 'WA', 'WY'))
us_south_states <- data.frame(region = "South", geo_value = c('AL', 'AR', 'DE', 'FL', 'GA', 'KY', 'LA', 'MD', 'MS', 'NC', 'OK', 'SC', 'TN', 'TX', 'VA', 'WV'))
us_midwest_states <- data.frame(region = "Midwest", geo_value = c('IA', 'IL', 'IN', 'KS', 'MI', 'MN', 'MO', 'ND', 'NE', 'OH', 'SD', 'WI'))
us_northeast_states <- data.frame(region = "Northeast", geo_value = c('CT', 'MA', 'ME', 'NH', 'NJ', 'NY', 'PA', 'RI', 'VT'))

us_region_state <- rbind(us_west_states, us_south_states, us_midwest_states, us_northeast_states)

unique_regions = unique(ww_data_region$region)

############################################################################################################################################

# Fix sample states and cutoff_date
samp_states = unique(adj_df$geo_value)

# Regional and state ratio and final infection estimate cutoff date
cutoff_date = as.Date("2022-05-01") 

infect_est_to_use = names(adj_df[,-c(1:2)]) 
lm_res_list = vector(mode = "list", length = length(infect_est_to_use))

# To save off shrinkage state ww data 
ww_state_region_list <- vector("list", length = length(samp_states))
names(ww_state_region_list) <- samp_states

for(z in 1:length(infect_est_to_use)){
  states_df = data.frame()
  for(j in 1:length(samp_states)){ # Change to state samples

    # Shrink state ww data towards region
    ww_data_one_state = ww_ts_long_final %>%
      rename(ww_state = ww) %>%  
      filter(geo_value == samp_states[j] & time_value <= end_date) %>%
      select(geo_value, time_value, ww_state)

    state_abb = samp_states[j]

    region_abb = us_region_state %>% filter(geo_value == state_abb) %>% pull(region)

    ww_data_one_region = ww_data_region %>%
      rename(time_value = sampling_week, ww_region = effective_concentration_rolling_average) %>%
      select(time_value, ww_region, region) %>%
      filter(region == region_abb & time_value <= end_date)

    ww_state_region = left_join(ww_data_one_region, ww_data_one_state, by = "time_value")

    # Compute weights (fraction of the state population accounted for in the state ww estimate - time-varying so differs for each time)
    tot_county_pop_state = tot_county_pop_long %>%
      filter(state == state_abb) %>%
      filter(date %in% ww_state_region$time_value) %>%
      rename(time_value = date) %>%
      select(time_value, tot_c_pop)

    ww_state_region = ww_state_region %>% left_join(tot_county_pop_state, by = "time_value")

    # Get state population
    
    # Reshape pop_df from wide to long format
    pop_long_state <- pop_long %>% # remove population_ in front of ex. population_2020
      filter(geo_value == state_abb) %>%
      select(-geo_value)
    
    ww_state_region$year <- as.numeric(format(ww_state_region$time_value, "%Y"))
    ww_state_region$year <- ifelse(ww_state_region$year == 2023, 2022, ww_state_region$year) # Deal with having only 2022 pop (no 2023)
    
    # Join the data based on year and geo_value/state abbreviation
    ww_state_region <- ww_state_region %>%
      left_join(pop_long_state, by = c("year" = "year"))
    
    # Shrink ww towards regional and store as ww
    ww_state_region = ww_state_region %>%
      mutate(weight = ifelse(!is.na(ww_state), tot_c_pop/state_pop, 0),
             ww = ifelse(!is.na(ww_state), weight*ww_state + (1-weight)*(ww_region), ww_region)) 
    
    # Mean and sd summary of weight by state
    # Roughly assesses coverage of wastewater concentration data of respective state populations on average
    ww_state_region %>%
      summarise(
        Count = n(), # number of observations in each group
        Mean = mean(weight, na.rm = TRUE) * 100,
        Sd = sd(weight, na.rm = TRUE) * 100,
        Median = median(weight, na.rm = TRUE) * 100,
        quantile_0.05 = quantile(weight, 0.05, na.rm = TRUE) * 100,
        quantile_0.95 = quantile(weight, 0.95, na.rm = TRUE) * 100
      ) # %>% print()
    
      # Fill missing geo_value entries for the state by carrying the last non-NA value upward
      ww_state_region <- ww_state_region %>%
        fill(geo_value, .direction = "up")
      
      # View one estimate of ww per week (complete set of weekly values)
      # ggplot(ww_state_region %>% filter(time_value <= end_date), aes(time_value, y = ww)) +
      #  geom_point(color = "#1f77b4")
      
      # Linearly interpolate ww values from weekly to daily
      ww_daily <- ww_state_region %>%
        group_by(region, geo_value) %>%
        complete(time_value = seq(min(time_value), max(time_value), by = "day")) %>%
        arrange(time_value) %>%
        mutate(ww = approx(time_value, ww, time_value, method = "linear", rule = 2)$y) %>%
        ungroup()
      
      # View daily interpolated estimates of ww
      # ggplot(ww_daily %>% filter(time_value <= end_date), aes(time_value, y = ww)) +
      #  geom_point(color = "#1f77b4")
      
    if (z == 1){
      ww_state_region_list[[state_abb]] <- ww_daily %>%
        filter(time_value >= cutoff_date) %>% # starting date for infection estimation
        select(time_value, geo_value, region, weight, ww_state, ww_region, ww)
    }
    
    prop_circ_df_state_capped = vmix %>%
      filter(geo_value == samp_states[j]) %>%
      filter((Delta > 0.01) & (Omicron > 0.01)|(Delta > 0.01)|(Omicron > 0.01)) %>% 
      filter(time_value %in% ww_daily$time_value & time_value <= sr_cutoff_date & time_value >= as.Date("2021-06-01")) 

    prop_df_first_row = prop_circ_df_state_capped %>% arrange(time_value) %>% slice(1)

    adj_df_subset_states = adj_df %>% filter(geo_value == samp_states[j])
    adj_infects_state_df = adj_df_subset_states %>%
      filter(date %in% prop_circ_df_state_capped$time_value)

    adj_infects_state_df = adj_infects_state_df %>%
      mutate(year = as.numeric(format(adj_infects_state_df$date, "%Y"))) %>%
      left_join(pop_long_state, by = c("year" = "year"))
    
    adj_infects_state = adj_infects_state_df %>%
      mutate(adj_infects_state_per100k = (!! rlang::sym(infect_est_to_use[z]) / state_pop) * 100000) %>% 
      select(geo_value, date, adj_infects_state_per100k)

    ww_state_within_adj_dates = ww_daily %>% filter(time_value %in% prop_circ_df_state_capped$time_value)

    # w(t) = \sum_{i\geq 0} infection(t)*prev_var_i(t) * Cd_i 
    adj_prop_circ_df_state = prop_circ_df_state_capped
    for(i in 1:nrow(adj_prop_circ_df_state)) adj_prop_circ_df_state[3:10][i, ] = adj_infects_state$adj_infects_state_per100k[i] * prop_circ_df_state_capped[3:10][i, ]
    
    total_z = rtestim::delay_calculator(adj_prop_circ_df_state$Delta, prop_circ_df_state_capped$time_value, delay_distn = ga_probs)
    total_v = rtestim::delay_calculator(adj_prop_circ_df_state$Omicron, prop_circ_df_state_capped$time_value, delay_distn = ga_probs)
    
    tmp_df = cbind(ww_state_within_adj_dates, data.frame(Delta = total_z, Omicron = total_v), first_prop_dt = prop_df_first_row$time_value, adj_inf_100k = adj_infects_state$adj_infects_state_per100k)
    states_df = rbind(states_df, tmp_df)
  }

  lm_res_list[[z]] = lm(ww ~ Delta + Omicron - 1, data = states_df) 
  # summary(lm_res_list[[z]])

}

# Delta and Omicron coefficients
coef_imp = lapply(lm_res_list, FUN = coef)
coef_imp

# matrix of coefficients
coef_mat <- str2str::lv2m(coef_imp, along = 1)
saveRDS(coef_mat, paste0("delta_omicron_coef_mat_", suffix_u, ".RDS"))

# vcov matrices list
vcov_imp <- lapply(X = lm_res_list, FUN = vcov)
vcov_imp

##############################################################################################################################
# Side-by-side boxplots of the 100 Delta and Omicron coefficients
# Reshape the data to long format
coef_mat_long <- as.data.frame(coef_mat) %>%
  pivot_longer(cols = c("Delta", "Omicron"), names_to = "Variant", values_to = "Coefficient")

# Custom colors for box and points 
box_colors <- c("Delta" = "#1f77b4", "Omicron" = "#DAA520")
point_colors <- c("Delta" = "#1f77b4", "Omicron" = "#8B6914")

# Create the plot
ggplot(coef_mat_long, aes(x = Variant, y = Coefficient, fill = Variant)) +
  geom_boxplot(outlier.shape = NA, color = "black", width = 0.6, alpha = 0.8) +
  geom_jitter(aes(color = Variant), size = 2, width = 0.1, alpha = 0.7) +
  scale_fill_manual(values = box_colors) + 
  scale_color_manual(values = box_colors) + 
  theme_minimal() +
  labs(y = "Shedding rate coefficient", x = "Variant") +
  theme(legend.position = "none") # Remove legend


# Plot of the 100 Omicron coefficients only (for flowchart)
# Extract Omicron values and create a data frame for ggplot2
omicron_data <- data.frame(
  Set = 1:nrow(coef_mat),         
  Omicron = coef_mat[,"Omicron"]    
)

ggplot(omicron_data, aes(x = 1, y = Omicron)) +  
  geom_boxplot(fill = "#DAA520", color = "black", alpha = 0.8) +  
  geom_jitter(color = "#8B6914", alpha = 0.7, width = 0.1, height = 0) +  
  coord_flip() +  # Make the plot horizontal
  theme_minimal() +
  labs(
    title = "",
    y = "Omicron shedding rate coefficients",
    x = ""
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

##############################################################################################################################################
# Source necessary files to run tf cv
Rcpp::sourceCpp("/code/supporting-files-for-deconvolutions/src/estim_path.cpp")
source("utils-arg.R")
source("utils-enlist.R")
source("cv_estimate_tf.R") 
##############################################################################################################################################
pop_used = "population_2022"

# Define a list to store results for all 100 sets
unadj_at_df_list <- vector("list", length = length(infect_est_to_use))
for(z in 1:length(infect_est_to_use)){
  # Store df of results for the sample of seven states for each z in a list
  all_samp_states_df <- vector(mode = "list", length = length(samp_states))
  names(all_samp_states_df) <- samp_states
  for(state in samp_states){
    # Get variant-specific part
    unadj_prop_circ_df_state = vmix %>%
      filter(geo_value == state) %>%
      filter(time_value >= cutoff_date & time_value <= end_date)
    
    # For E(aX + bY) est 
    for(i in 1:nrow(unadj_prop_circ_df_state)) unadj_prop_circ_df_state[c("Delta", "Omicron")][i, ] = coef_mat[z,] * unadj_prop_circ_df_state[c("Delta", "Omicron")][i, ]
    
    # aX + bY mean and var est
    variant_specific_mean_est = rowSums(unadj_prop_circ_df_state[c("Delta", "Omicron")]) 
    
    # Next, get deconvolved cases per 100k 
    state_pop = pop_df %>%
      filter(geo_value == state) %>%
      pull(pop_used)
    
    ww_data_one_state <- ww_state_region_list[[state]]
    
    state_unadj_df = unadj_df %>% filter(geo_value == state & time_value %in% ww_data_one_state$time_value) 
    unadj_infects_state = (state_unadj_df$unadj_infect / state_pop) * 100000
    
    # Construct diagonal matrix where the diagonal entries are unadjusted_infects * variant specific component for each day
    diag_coef <- unadj_infects_state * variant_specific_mean_est
    D <- sparseMatrix(i = 1:length(diag_coef), j = 1:length(diag_coef), x = diag_coef)
    
    # Construct lower-triangular Gamma/L matrix w/row sums = 1
    
    # Row and column indices
    row_idx <- row(matrix(0, length(diag_coef), length(diag_coef)))
    col_idx <- col(matrix(0, length(diag_coef), length(diag_coef)))
    
    # Calculate index into ga_probs vec to input in lower triangle
    ga_idx <- row_idx - col_idx + 1
    
    # Logical mask for lower triangle
    mask <- ga_idx >= 1
    
    # Build final matrix using ga_probs
    L <- matrix(0, length(diag_coef), length(diag_coef))
    L[mask] <- ga_probs[ga_idx[mask]]
    
    # Fill NA values in lower triangle (beyond the length of ga_probs) to be 0
    L[is.na(L)] <- 0
    
    # View L
    # L
    # Check all upper triangular elements are 0
    #all(L[upper.tri(L)] == 0)
    
    # Re-normalize L so each row-sum = 1:
    L <- L / rowSums(L)
    
    # Now multiply L %*% D to get C matrix
    C = L %*% D
    
    # Get a dense matrix as a result, so convert to sparse (compressed column format) before use in cv_estimate_tf 
    C <- as(C, "dgCMatrix")
    
    # CV estimate using tf 
    res <- cv_estimate_tf(ww_data_one_state$ww, 
                          x = 1:ncol(C),
                          cmat = C, 
                          korder = 3, 
                          error_measure = "mse") 
    
    # Plot of the thetas for all lambdas
    # matplot(res$full_fit$thetas, type = "l")
    # Plot line of thetas for lambda_min
    # matplot(res$full_fit$thetas[,which(res$lambda == res$lambda.min)], type = "l")
    
    final_thetas = res$full_fit$thetas[,which(res$lambda == res$lambda.min)]
    
    # Join pop temporarily to compute adjusted infection rates
    pop_df_sub = pop_df %>%
      filter(geo_value %in% state) %>%
      select(geo_value, pop_used) %>%
      rename(pop = pop_used)
    
    # Assign column for at_val in state_unadj_df
    state_unadj_df$at_val = final_thetas
    
    # Calculate infections: at_val * unadj_infect and pop-adjust rate
    state_unadj_df <- state_unadj_df %>% left_join(pop_df_sub, by = "geo_value") %>%
      mutate(adj_infect = as.double(at_val*unadj_infect),
             adj_inf_rate = (adj_infect / pop) * 100000)
    
    all_samp_states_df[[state]] <- state_unadj_df
  }
  # Combine all rows for states into a single data frame
  unadj_at_df_list[[z]] <- dplyr::bind_rows(all_samp_states_df)
}

# Save inv_ratios_all_states as RDS
saveRDS(unadj_at_df_list, paste0("unadj_at_df_list_Mar20_", suffix_u, ".RDS"))

##########################################################################################################################################
# Compute time-varying quantiles for the adj_infect column by first combining the list of data frames into a single data frame

unadj_at_df_list <- readRDS(paste0("unadj_at_df_list_Mar20_", suffix_u, ".RDS")) 
combined_unadj_at_df <- bind_rows(unadj_at_df_list)
                                                         
# Compute the 5th and 95th percentiles for the adj_infect col
quantiles_unadj_at_df <- combined_unadj_at_df %>%
  group_by(geo_value, time_value) %>%
  summarize(
    adj_infect_q05 = quantile(adj_infect, 0.05),
    adj_infect_q50 = quantile(adj_infect, 0.50),
    adj_infect_q95 = quantile(adj_infect, 0.95),
    .groups = "drop"
  )

# View the resulting data frame
head(quantiles_unadj_at_df)
#############################################################################################################################################
# Now, load pop and case data to produce plots 

# Load pop data
pop_df = readRDS("pop_df.RDS")

quantiles_unadj_at_df = left_join(quantiles_unadj_at_df, pop_df %>% select(geo_value, pop_used))
names(quantiles_unadj_at_df)[ncol(quantiles_unadj_at_df)] <- "pop"
##############################################################################################################################################
# Load covidcast finalized case data (for plotting later on)

cases <- covidcast_signal("jhu-csse", "confirmed_7dav_incidence_num",
                          start_day = start_date, end_day = end_date,
                          as_of = as.Date("2023-07-06"),
                          geo_type = "state") %>%
  select(geo_value, time_value, case_count_7d_av = value) %>% 
  mutate(geo_value = toupper(geo_value))

quantiles_unadj_at_df = left_join(quantiles_unadj_at_df, cases, by = c("geo_value", "time_value"))

########################################################################################################################
# Compute rates per 100k for infection quantiles and cases

quantiles_unadj_at_df = quantiles_unadj_at_df %>% mutate(adj_inf_q50_rate = (adj_infect_q50 / pop) * 100000,
                                                         cases_rate_7d_av = (case_count_7d_av / pop) * 100000,
                                                         adj_infect_q05_rate = (adj_infect_q05 / pop) * 100000,
                                                         adj_infect_q95_rate = (adj_infect_q95 / pop) * 100000)


##############################################################################################################################################
# Save resulting df of case and infection rates
saveRDS(quantiles_unadj_at_df, file = paste0("quantiles_ww_adj_df_", suffix_u, ".RDS"))

# Read in file and store as ww_adj_df
ww_adj_df = readRDS(paste0("quantiles_ww_adj_df_", suffix_u, ".RDS"))
#############################################################################################################################
# Plots
# Side-by-side plots of infections with CI

# Deal with multiple consecutive values of cases being the same (ex. one measurement per week or every couple of days) for a state
# To make the plotting lines for cases smooth-ish
# Omit multiple consecutive values of cases that are the same for each state
# Keep the first of each consecutive values & the first for each state
filtered_df_cases <- ww_adj_df %>%
  select(time_value, geo_value, cases_rate_7d_av) %>%
  group_by(geo_value) %>%  
  mutate(change_flag = cases_rate_7d_av != lag(cases_rate_7d_av, default = NA)) %>%  # Flag if the value changes
  filter(change_flag | row_number() == 1) %>%  # Keep the first row, even if no change
  select(time_value, geo_value, cases_rate_7d_av)

# Plot specifications
cases_to_plot_df = filtered_df_cases 
manual_breaks <- as.Date(c("2022-07-01", "2022-10-01", "2023-01-01")) # Define specific date breaks manually

# Adj. infect. / 100k plot and cases / 100k with all infection trajectories
# Combine all the data frames in unadj_at_df_list with ids for each of the 100 dfs/trajectories under source
combined_unadj_at_df_wsource <- bind_rows(unadj_at_df_list, .id = "source")

ggplot(ww_adj_df) +
  geom_line(data = cases_to_plot_df, aes(time_value, cases_rate_7d_av, color = "Cases / 100k"), size = 0.4) +
  # Add lines for every trajectory in the combined_unadj_at_df_wsource
  geom_line(data = combined_unadj_at_df_wsource, aes(time_value, adj_inf_rate, group = source, color = "Inf. trajectories"), size = 0.2, alpha = 0.7) +
  # Overlay median
  geom_line(aes(time_value, adj_inf_q50_rate, color = "Adj. inf. / 100k"), size = 0.65) + 
  scale_x_date(name = "", breaks = manual_breaks, date_labels = "%b %Y", expand = c(0, 0)) +
  scale_y_continuous(labels = scales::comma) +
  scale_color_manual(name = 'Estimates', 
                     breaks = c('Adj. inf. / 100k', 'Cases / 100k', 'Inf. trajectories'),
                     values = c('Adj. inf. / 100k' = "midnightblue", 
                                'Cases / 100k' = "darkorange2", 
                                'Inf. trajectories' = "#78C1F2")) +
  facet_wrap(. ~ geo_value, ncol = 2, scales = "free_y") +
  theme_bw(16) +
  ylab("Incidence per 100k population") +
  theme(axis.text.x = element_text(size = 9.5),
        axis.title.y = element_text(size = 11.5),
        axis.text.y = element_text(size = 9.5),
        title = element_blank(),
        legend.title = element_text(size = 10.5),
        legend.text = element_text(size = 9.5),
        legend.position = "bottomright",
        legend.box.spacing = unit(-15, "pt"),
        strip.text = element_text(size = 7, face = "bold", margin = margin(1, 0, 1, 0, "mm")))
ggsave(filename = paste0("fig-infect-cases-omicron-no-dual-", suffix_d, ".pdf"), width = 10, height = 7)


# Now, Adj. infect. / 100k plot and cases / 100k plot
# But with two axes for cases and infections (may help with cases being so much lower than infections)
# Source code for transforms from https://github.com/cmu-delphi/insightnet-workshop-2024/blob/main/_code/ca_cases_deaths.R
# and https://github.com/cmu-delphi/insightnet-workshop-2024/blob/main/slides/day1-morning.qmd

# Handy function to produce a transformation from one range to another
trans = function(x, from_range, to_range) {
  (x - from_range[1]) / (from_range[2] - from_range[1]) *
    (to_range[2] - to_range[1]) + to_range[1]
}

# Compute ranges of the two signals, and transformations in b/w them
range1 = range(ww_adj_df$adj_inf_q50_rate) 
range1[1] <- 0  # Force y-axis to start at 0
range2 = range(ww_adj_df$cases_rate_7d_av)
trans12 = function(x) trans(x, range1, range2)
trans21 = function(x) trans(x, range2, range1)

ggplot(ww_adj_df) +
  # Add lines for every trajectory in the combined_unadj_at_df_wsource
  geom_line(data = combined_unadj_at_df_wsource, aes(time_value, adj_inf_rate, group = source, color = "Inf. trajectories"), size = 0.2, alpha = 0.7) +
  
  geom_line(aes(time_value,  adj_inf_q50_rate, color = "Adj. inf. / 100k"), size = 0.6) +
  geom_line(data = cases_to_plot_df %>% mutate(cases_rate_7d_av = trans21(cases_rate_7d_av)), aes(time_value, cases_rate_7d_av, color = "Cases / 100k"), size = 0.4) +
  scale_x_date(name = "", breaks = manual_breaks, date_labels = "%b %Y", expand = c(0, 0)) +
  scale_y_continuous(labels = scales::comma, 
                     expand = expansion(c(0, 0.05)),
                     name = "Incident infections per 100k", 
                     #limits = range1,
                     sec.axis = sec_axis(
                       transform = trans12,  # Secondary axis is purely based on the predefined transformation (trans12), which was set up using the range of both variables beforehand. 
                       name = "Incident cases per 100k"))  + 
  coord_cartesian(ylim = range1) + 
  scale_color_manual(name = 'Estimates', # Legend
                     breaks=c('Adj. inf. / 100k',
                              'Cases / 100k',
                              'Inf. trajectories'),
                     values=c('Adj. inf. / 100k' = "midnightblue",
                              'Cases / 100k' = "darkorange2",
                              'Inf. trajectories' = "#78C1F2")) +
  scale_fill_manual(name = 'CI', values = c("skyblue3")) + 
  facet_wrap(. ~ geo_value, ncol = 2, scales = "free_y") +
  theme_bw(16) +
  theme(axis.text.x = element_text(size = 11.5),
        axis.title.y = element_text(size = 13.5),
        axis.text.y = element_text(size = 11.5),
        title = element_blank(),
        legend.title = element_text(size = 13.5),
        legend.text = element_text(size = 11.5),
        legend.position = "bottomright",
        legend.box.spacing = unit(-15, "pt"),
        strip.text = element_text(size = 14, face = "bold", margin = margin(1, 0, 1, 0, "mm")))
ggsave(filename = paste0("fig-infect-cases-omicron-", suffix_d, ".pdf"), width = 10, height = 7) 

# Now, Adj. infect. / 100k plot and cases / 100k plot (same plot as previous)
# but let's colour background by Omicron subvariant 

`%nin%` <- Negate(`%in%`)

var_props_omic_sub = readr::read_csv("seq_df_post_decon_omicron.csv")

df_long <- var_props_omic_sub %>%
  gather(key = "Variant", value = "Prop", -State, -Date) %>%  # gather variant columns, keeping State and Date
  # Only need Omicron sub-variants
  filter(Variant %nin% c("Alpha", "Beta", "Epsilon", "Iota", "Gamma", "Delta", "Other")) %>%
  # Only need dates >= first date for the adjusted infections
  filter(Date >= ww_adj_df$time_value[1]) %>% 
  filter(State %in% c("CA", "NY", "NC", "MN", "HI", "WA", "NV")) %>%
  arrange(State, Date) %>%  # Arrange by state and date
  filter(Date <= max(ww_adj_df$time_value)) %>%
  mutate(geo_value = State)

# Remove (Omicron) from after subvariant name
df_long$Variant <- gsub(" \\(Omicron\\)", "", df_long$Variant)

variant_colors <- viridis::viridis(14, begin = 0.3, end = 1)

geom_area_breaks <- unique(df_long$Variant)[which(unique(df_long$Variant) %in% unique(df_long$Variant[df_long$Prop > 0.1]))] 

# Base plot of variant proportions 
p_background <- ggplot(df_long %>% mutate(geo_value = State), aes(Date, y = Prop, fill = Variant, group = interaction(geo_value, Variant))) +
  geom_area(position = "fill", alpha = 0.7) +
  scale_fill_manual(name = "Variant", values = setNames(variant_colors, unique(df_long$Variant)), breaks = geom_area_breaks) +  # To ensure legend colors are fully opaque
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

p_lines <- ggplot(ww_adj_df) +
  # Add lines for every trajectory in the combined_unadj_at_df_wsource
  geom_line(data = combined_unadj_at_df_wsource, aes(time_value, adj_inf_rate, group = source, color = "Inf. trajectories"), size = 0.2, alpha = 0.7) +
  geom_line(aes(time_value, adj_inf_q50_rate, color = "Median inf. / 100k"), size = 0.6) +
  geom_line(data = cases_to_plot_df %>% mutate(cases_rate_7d_av = trans21(cases_rate_7d_av)), aes(time_value, cases_rate_7d_av, color = "Cases / 100k"), size = 0.4) +
  scale_x_date(name = "", date_breaks = "3 month", date_labels = "%b %Y", expand = c(0, 0)) +
  scale_y_continuous(labels = scales::comma, 
                     expand = expansion(c(0, 0.05)),
                     name = "Incident infections per 100k", 
                     sec.axis = sec_axis(
                       transform = trans12,  # Secondary axis is purely based on the predefined transformation (trans12), which was set up using the range of both variables beforehand. 
                       name = "Incident cases per 100k"))  + 
  coord_cartesian(ylim = range1) + 
  scale_color_manual(name = '', # Legend
                     breaks=c('Median inf. / 100k',
                              'Cases / 100k',
                              'Inf. trajectories'),
                     values=c('Median inf. / 100k' = "midnightblue",
                              'Cases / 100k' = "darkorange2",
                              'Inf. trajectories' = "cyan2")) +
  labs(fill = "", color = "") +  
  facet_wrap(. ~ geo_value, ncol = 2, scales = "free_y") +
  theme_bw(16) +
  theme(panel.background = element_rect(fill = "transparent", color = NA),  # transparent background
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid = element_blank(), # Remove all grid lines 
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 11.5),
        axis.text.y = element_text(size = 9.5),
        title = element_blank(),
        legend.justification = c(1, 1), 
        legend.position = "top",
        legend.box.margin = margin(5, 5, 5, 5),  
        legend.title = element_text(size = 10.5),
        legend.text = element_text(size = 9.5),
        strip.text = element_text(size = 7, face = "bold", margin = margin(1, 0, 1, 0, "mm")))

# Convert both ggplots to grobs
g1 <- ggplotGrob(p_background)
g2 <- ggplotGrob(p_lines)

# Ensure same dimensions: copy widths and heights
g1$widths <- g2$widths
g1$heights <- g2$heights

# Set the output file and size
pdf(paste0("fig-infect-cases-omicron-with-shading-", suffix_d, ".pdf"), width = 10, height = 7)  

# Overlay the grobs using grid
grid.newpage()
grid.draw(g1)        # draw base
grid.draw(g2)        # overlay

dev.off()


# Plot inverse reporting ratios for all infection trajectories
ggplot(combined_unadj_at_df_wsource) +
  geom_line(data = combined_unadj_at_df_wsource, aes(time_value, at_val, group = source), color = "maroon", size = 0.2, alpha = 0.5) +
  scale_x_date(name = "", breaks = manual_breaks, date_labels = "%b %Y", expand = c(0, 0)) +
  ylab("Inverse reporting ratio") +
  facet_wrap(. ~ geo_value, ncol = 2) +
  theme_bw(16) +
  theme(axis.text.x = element_text(size = 9.5),
        axis.title.y = element_text(size = 11.5),
        axis.text.y = element_text(size = 9.5),
        title = element_blank(),
        legend.title = element_text(size = 10.5),
        legend.text = element_text(size = 9.5),
        legend.position = "top",
        strip.text = element_text(size = 7, face = "bold", margin = margin(1, 0, 1, 0, "mm")))
ggsave(filename = paste0("inverse-reporting-ratios-", suffix_d, ".pdf"), width = 10, height = 7)

