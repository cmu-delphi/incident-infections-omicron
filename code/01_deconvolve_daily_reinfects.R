library(tidyverse)

##############################################################################################################################################
# Some initial settings
start_date = as.Date("2020-03-01") 
end_date =  as.Date("2023-03-01") 
length_start_end = end_date - start_date + 1

##############################################################################################################################################
# Source necessary files to run tf cv
Rcpp::sourceCpp("/code/supporting-files-for-deconvolutions/src/estim_path.cpp") # Change this as needed.
source("utils-arg.R")
source("utils-enlist.R")
source("cv_estimate_tf.R") 

##############################################################################################################################################
# Load reinfection data

states = c("NC", "MN")

for(state in states){
  # Load reinfection data
  setwd('/data/state-reinfection-csvs')
  dat <- read.csv(paste0(state, "_reinfect.csv"))
  dat$time_value <- as.Date(dat$time_value)
  
  setwd(paste0("/data/state-deconvolved-case-data/", state)) 
  params <- read_rds("convolution-mat-list.rds")
  
  # Deconvolve each of reinfections and cases to infection onset date by each variant, then sum them up
  
  # Start with y = reinfections * prop of variant in circulation
  # Get cmat
  
  final_thetas_op_list = vector(mode = "list", length = length(params$Variant))
  names(final_thetas_op_list) <- params$Variant
  for(var in params$Variant){
    idx_subtract_from_start = length_start_end + 1 # + 1 to length_start_end to ensure that index is beyond cmat_var and y indices
    idx_subtract_from_end = length_start_end + 1
    if(dat$time_value[1] != start_date) idx_subtract_from_start = 1:(dat$time_value[1] - start_date)
    if(dat$time_value[length(dat$time_value)] != end_date) idx_subtract_from_end = (length_start_end - (end_date - dat$time_value[length(dat$time_value)]) + 1):length_start_end
    
    cmat_var = params[params$Variant == var,]$Cmat[[1]][-c(idx_subtract_from_start, idx_subtract_from_end), -c(idx_subtract_from_start, idx_subtract_from_end)]  
    y_reinfect_var = dat$reinfect * params[params$Variant == var, ]$probs[[1]][-c(idx_subtract_from_start, idx_subtract_from_end)]
    y_cases_var = dat$cases * params[params$Variant == var, ]$probs[[1]][-c(idx_subtract_from_start, idx_subtract_from_end)] 
    
    # Checks that dimensions of y_var and cmat are reasonable
    if(nrow(cmat_var) != ncol(cmat_var)) cli::cli_abort("`cmat` is not square. Check the dimensions.")
    if(nrow(cmat_var) != length(y_reinfect_var)) cli::cli_abort("Nrow of `cmat` is not the same as the length of `y_reinfect_var`. Check the dimensions and dates used of each.")
    if(ncol(cmat_var) != length(y_reinfect_var)) cli::cli_abort("Ncol of `cmat` is not the same as the length of `y_reinfect_var`. Check the dimensions and dates used of each.")
    if(nrow(cmat_var) != length(y_cases_var)) cli::cli_abort("Nrow of `cmat` is not the same as the length of `y_cases_var`. Check the dimensions and dates used of each.")
    if(ncol(cmat_var) != length(y_cases_var)) cli::cli_abort("Ncol of `cmat` is not the same as the length of `y_cases_var`. Check the dimensions and dates used of each.")
    
    # CV estimate using tf 
    res_op_reinfect <- cv_estimate_tf(y_reinfect_var, 
                             cmat = cmat_var,
                             error_measure = "mse")
    
    res_op_cases <- cv_estimate_tf(y_cases_var, 
                                   cmat = cmat_var,
                                   error_measure = "mse")

    final_thetas_op_reinfect = res_op_reinfect$full_fit$thetas[,which(res_op_reinfect$lambda == res_op_reinfect$lambda.min)]
    final_thetas_op_cases = res_op_cases$full_fit$thetas[,which(res_op_cases$lambda == res_op_cases$lambda.min)]
    final_thetas_op_list[[var]] = data.frame(time_value = seq(dat$time_value[1], dat$time_value[length(dat$time_value)], by = "day"), geo_value = state, reinfect = final_thetas_op_reinfect, cases = final_thetas_op_cases) 
  }
  final_thetas_op_df_state = bind_rows(final_thetas_op_list, .id="variant")
  
  setwd(paste0("data/state-deconvolved-case-data/", state)) # Where the results and plots are to be saved for the state
  write_rds(final_thetas_op_df_state, "reinfect-final-thetas-df.rds")
}