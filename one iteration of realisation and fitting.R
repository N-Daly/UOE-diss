# this script runs expecting that all the appropriate variables have been
# defined in the environment already - double check this as each fit takes a while
source("function Definitions.R")
iteration_start <- proc.time()

# simulate a ground truth
sim_info <- simulate_lcgp_dual_obs_HR_thinning(
  true_sigmaA = true_sigmaA, true_sigmaB = true_sigmaB,
  true_gammaA = true_gammaA, true_gammaB = true_gammaB,
  true_beta0 = -4,
  true_rho = true_range, true_sigma_GRF = true_sigma_grf
)

lets_have_a_look_at_you(sim_info, detect_df = true_detect_prob_df)


#the sampled points are in observed_points
observed_points <- sim_info$samples_df

true_loglambda_at_ip <- fm_evaluate(
  proj = proj_from_DGP_mesh_to_intensity_pred_locs,
  field = sim_info$log_lambda
)

list_of_models = list()

# need to save space and these model fits are large
# so we delete the fitted object before fitting another
rm(fit); invisible(gc())
# # Model what A and B saw as a combined observer with a hn detection function
catt("fitting single HN")
start <- proc.time()
fit <- model_one_observer_hn(
  observed_points, matern_prior, ips=ips, bru_verbose=how_verbose
)
end <- proc.time()
print((end-start)[3])

start <- proc.time()
list_of_models$one_obs_hn <- get_preds_from_one_observer_hn_fit(fit)
list_of_models$one_obs_hn$name <- "one_obs_hn"
end <- proc.time()
print((end-start)[3])

plot_detect_pred(
  pred_df = list_of_models$one_obs_hn$pred_detect,
  main = list_of_models$one_obs_hn$name
)

rm(fit); invisible(gc())
# # Model what A and B saw as a two observer likelihood with hn detection functions
catt("fitting two HN ")
start <- proc.time()
fit <- model_two_observers_hn(
  observed_points,
  prior_on_sigma = bm_marginal(qexp, pexp, dexp, rate = 2/8),
  mtrn_prior = matern_prior,
  ips = ips_with_detection_states,
  bru_inits = list(
    sigmaA = qexp(pnorm(2), rate = 2/8),
    sigmaB = qexp(pnorm(2), rate = 2/8)
  ),
  bru_verbose=how_verbose
)
end <- proc.time()
print((end-start)[3])

start <- proc.time()
list_of_models$two_obs_hn <- get_preds_from_two_observers_hn_fit(fit)
list_of_models$two_obs_hn$name <- "two_obs_hn"
end <- proc.time()
print((end-start)[3])

plot_detect_pred(
  pred_df = list_of_models$two_obs_hn$pred_detect,
  main = list_of_models$two_obs_hn$name
)

rm(fit); invisible(gc())
# Model what A and B saw as a combined observer with a  hazard-rate detection function
catt("fitting merged HR")
start <- proc.time()
fit <- model_one_observer_hr(
  observed_points, mtrn_prior =  matern_prior, ips=ips,
  prior_on_gamma = bm_marginal(qgamma, pgamma, dgamma, shape=2, rate=1),
  bru_initial_params = list(
    sigma = qnorm(pexp(2, rate= 1/8)),
    gamma = qnorm(pgamma(2, shape=2, rate=1))
  ),
  bru_verbose = how_verbose
)
end <- proc.time()
print((end-start)[3])


start <- proc.time()
list_of_models$one_obs_HR <- get_preds_from_one_observer_hr_fit(fit)
list_of_models$one_obs_HR$name <- "one_obs_HR"
end <- proc.time()
print((end-start)[3])

plot_detect_pred(
  pred_df = list_of_models$one_obs_HR$pred_detect,
  main = list_of_models$one_obs_HR$name
)


rm(fit); invisible(gc())
catt("fitting dual HR")
start <- proc.time()
fit <- model_two_observers_hr(
  observed_points,
  prior_on_gamma = bm_marginal(qunif, punif, dunif, min = 0.01, max = 10),
  prior_on_sigma = bm_marginal(qexp, pexp, dexp, rate = 1/2),
  bru_initial_params = list(
    gammaA = qnorm(punif(1, min = 0.01, max = 10)),
    gammaB = qnorm(punif(4, min = 0.01, max = 10)),
    sigmaA = qnorm(pexp(2, rate = 1/2)),
    sigmaB = qnorm(pexp(2, rate = 1/2))
  ),
  mtrn_prior = matern_prior,
  ips=ips_with_detection_states, 
  bru_verbose=how_verbose
)
end <- proc.time()
print((end-start)[3])

start <- proc.time()
list_of_models$two_obs_HR <- get_preds_from_two_observers_hr_fit(fit)
list_of_models$two_obs_HR$name <- "two_obs_HR"
end <- proc.time()
print((end-start)[3])

plot_detect_pred(
  pred_df = list_of_models$two_obs_HR$pred_detect,
  main = list_of_models$two_obs_HR$name
)


rm(fit); invisible(gc())
# Model what A and B saw as a combined observer with a spline(like) detection function
catt("fitting spline")
start <- proc.time()
fit <- model_one_observer_spline(
  observed_points, matern_prior, detect_matern, ips=ips, bru_verbose=how_verbose
)
end <- proc.time()
print((end-start)[3])

start <- proc.time()
list_of_models$one_obs_spline <- get_preds_from_one_observer_spline_fit(fit)
list_of_models$one_obs_spline$name <- "one_obs_spline"
end <- proc.time()
print((end-start)[3])

plot_detect_pred(
  pred_df = list_of_models$one_obs_spline$pred_detect,
  main = list_of_models$one_obs_spline$name
)


rm(fit); invisible(gc())
# scoring
some_results <- get_scoring_differences(
  list_of_models, 
  true_detect = true_detect_prob,
  true_loglambda = true_loglambda_at_ip, 
  distance_ips = distance_ips
)

iteration_end <- proc.time()
iteration_duration <- (iteration_end-iteration_start)[[3]]
catt("Iteration took ", pretty_print_seconds(iteration_duration))

