library(tidyverse)
library(Matrix)

# Load inverse ratio data (alpha to be used as mean in gp)
inv_ratios_all_states = readRDS("ww_inv_ratios_all_states_M20.RDS")
########################################################################################################################################
# Helper functions

one_thresh <- function(x) pmax(x, 1)

companion_mat <- function(k, direction = c("forward", "backward")) {
  mat <- matrix(0, k, k)
  mat[1, 1:k] <- dspline::d_mat(k, 1:(k + 1), row_idx = 1)[1, -1] *
    (-1)^(k + 1)
  if (k > 1) for (it in 1:(k - 1)) mat[1 + it, it] <- 1
  direction <- match.arg(direction)
  if (direction == "backward") {
    matold <- mat
    for (i in 1:k) for (j in 1:k) mat[i, j] <- matold[k - i + 1, k - j + 1]
  }
  mat
}

fwd_bkd_extrap <- function(am_matrix_row, start_date, end_date, mon_obs_weekly, k, svar, multiplier){
  # Construct trend1, trend2, trend3 and fill in the missing trend components
  sample_m = data.frame(trend1 = am_matrix_row,
                        trend2 = dplyr::lag(am_matrix_row, n = 1),
                        trend3 = dplyr::lag(am_matrix_row, n = 2))
  sample_m[1,2] = sample_m[2,3] = sample_m[3,1] - (3*sample_m[2,1] - 3*sample_m[2,2])
  sample_m[1,3] = sample_m[2,1] - (3*sample_m[1,1] - 3*sample_m[1,2])


  all_dates <- seq(from = start_date, to = end_date, by = "days")
  N <- length(which(lubridate::wday(all_dates) == 2)) # count number of M between start_date and end_date
  ndays_diff_M = as.Date(mon_obs_weekly[1]) - start_date + 1 # how many days between 1st Monday relative to start date?
  first_rdd_mon_num = sum(lubridate::wday(all_dates[1:ndays_diff_M]) == 2) # where is yd$rdd[1] # It is the 7th Monday.
  ktrend <- k # allowed to be <= k if we want lower-smoothness extrapolation
  RQR <- matrix(0, k, k)

  Ttrendf <- matrix(0, k, k)
  Ttrendb <- matrix(0, k, k)
  Ttrendf[1:ktrend, 1:ktrend] <- companion_mat(ktrend)
  Ttrendb[1:ktrend, 1:ktrend] <- companion_mat(ktrend, "backward")

  a <- matrix(0, k, N)
  a[,first_rdd_mon_num:(length(mon_obs_weekly) + first_rdd_mon_num-1)] <- t(sample_m) 
  # backward
  RQR[1, 1] <- svar * multiplier 
  for (i in first_rdd_mon_num:2) {
    a[,i-1] <- Ttrendb %*% a[,i]
  }
  # forward
  for (i in (length(mon_obs_weekly) + first_rdd_mon_num-1):(N - 1)) {
    a[,i+1] <- Ttrendf %*% a[,i]
  }
  alpha <- a[1, ]

  # Linearly interpolate from M to all days from start_date to end_date

  alpha_df = data.frame(date = all_dates[lubridate::wday(all_dates) == 2], alpha = alpha)
  alpha_li_df <- alpha_df %>%
    complete(date = full_seq(date, 1)) %>%
    mutate(alpha = approx(date[!is.na(.$alpha)], y = alpha[!is.na(.$alpha)], xout = date, method = "linear")$y)

  # Threshold so alpha are at least 1
  alpha_li_df = alpha_li_df %>% mutate(alpha = one_thresh(alpha))

  alpha_li_df$alpha
}
########################################################################################################################################
# For loop to get covariance matrix to be used in gp
states_gp_res <- vector("list", length(names(inv_ratios_all_states)))
names(states_gp_res) <- names(inv_ratios_all_states)
for(state in names(inv_ratios_all_states)){

  # State ss and fs
  ss = inv_ratios_all_states[[state]]$ss
  fs = inv_ratios_all_states[[state]]$fs

  # Create U and B diagonal matrices
  U <- bdiag(purrr::array_branch(ss$H, 3))
  dim(U)


  y_long = as.vector(t(ss$y))
  z_df = do.call(rbind, rep(list(ss$Z[,,1]), length(fs$alphahat[,2])))
  z_mat = as.matrix(z_df)
  for(i in 1:nrow(z_mat)){
    if(is.na(y_long[i])){
      z_mat[i,1] = 0
    }
  }

  z_df2 = data.frame(id = rep(1:length(fs$alphahat[,2]), each = 2), z_mat)
  rownames(z_df2) = NULL
  z_list = split(z_df2 %>% select(-id), f = z_df2$id)
  z_list = lapply(z_list, function(x) as.matrix(x))

  names(z_list) = NULL
  rownames(z_list) = NULL

  B <- bdiag(z_list)
  dim(B)


  # Create G matrix
  G <- bdiag(c(list(fs$P[,,1]), rep(list(ss$R[,,1] %*% ss$Q[,,1] %*% t(ss$R[,,1])), length(fs$alphahat[,2])-1)))
  dim(G)


  # Construct block lower triangular A matrix
  n = dim(ss$T)[3]
  res = diag(4*n)
  for (i in 1:(n - 1)) {
    cidx <- (i - 1) * 4 + 1
    Tt <- ss$T[,,i]

    for (j in (i + 1):(n)) {
      ridx <- (j - 1) * 4 + 1
      res[ridx:(ridx + 3), cidx:(cidx + 3)] <- Tt
      Tt <- ss$T[,,j] %*% Tt
    }
  }
  A = res


  sigma_aa = A %*% G %*% t(A)

  sigma_aa_given_y = sigma_aa - sigma_aa %*% t(B) %*% solve(B %*% sigma_aa %*% t(B) + U) %*% B %*% sigma_aa

  # Extract correct part of the res matrix (corresponding to the cov between the state elements of interest)
  cond_cov_mat = sigma_aa_given_y[seq(2, 4*n, by = 4), seq(2, 4*n, by = 4)]

  # Add some variance in all directions
  cond_cov_mat <- cond_cov_mat + diag(ncol(cond_cov_mat))*0.05 
  det(cond_cov_mat)

  # Simulate a_m vector from mvnorm
  am_matrix <- MASS::mvrnorm(n = 100, mu = as.vector(fs$alphahat[,2]), Sigma = cond_cov_mat) 
  dim(am_matrix) # each row is a sample of a_m


  # Each column is a set of results pertaining to a row of am_matrix.
  state_res = apply(am_matrix, MARGIN = 1, function(x) fwd_bkd_extrap(x,
                                                                      inv_ratios_all_states[[state]]$date[1], # start_date
                                                                      inv_ratios_all_states[[state]]$date[length(inv_ratios_all_states[[state]]$date)], # end_date
                                                                      inv_ratios_all_states[[state]]$mon_obs_weekly,
                                                                      3L,
                                                                      inv_ratios_all_states[[state]]$var_sero_est,
                                                                      inv_ratios_all_states[[state]]$mult_est))

  states_gp_res[[state]] <- list(gp_alphas = state_res,
                                 unadj_infect = inv_ratios_all_states[[state]]$unadj_infect)
}

# Get adjusted infection estimates using each col of inverse ratio estimtates from states_gp_res

gp_adjusted_infections_fun <- function(state, states_gp_res){

  res_state = states_gp_res[[state]]

  state_alpha <- res_state$gp_alphas

  dates <- seq(from = as.Date(inv_ratios_all_states[[state]]$date[1], format = "%d-%m-%Y"),
               to = as.Date(inv_ratios_all_states[[state]]$date[length(inv_ratios_all_states[[state]]$date)], format = "%d-%m-%Y"),
               by = "days")

  ######################################################################################################################
  # k(t) = \hat{a}(t) * I(t)
  # Note only get adjusted infections up to last available sero
  # So must get subset of unadj_infect up to that date
  res_state$unadj_infect = ifelse(res_state$unadj_infect == 0, NA, res_state$unadj_infect)
  res_state$unadj_infect = zoo::na.approx(res_state$unadj_infect, rule = 2)

  unadj_infect_sub = res_state$unadj_infect
  adj_infect = unadj_infect_sub * state_alpha

  ######################################################################################################################
  # Adjusted infection results

  state_df_res = data.frame(date = dates,
                             geo_value = state,
                             adj_infect = adj_infect)

  state_df_res
}


gp_adj_df_list <- lapply(names(inv_ratios_all_states), function(x) gp_adjusted_infections_fun(x, states_gp_res))
names(gp_adj_df_list) <- names(inv_ratios_all_states)

# Save file of adj_df_list
saveRDS(gp_adj_df_list, file = "gp_adj_df_list_Mar20.RDS")








