different_configs <- list(
  list(
    name = "Gamma Gamma(2,1) Sigma Exp(1/8)",
    gamma_prior = bm_marginal(qgamma, pgamma, dgamma, shape=2, rate=1),
    sigma_prior = bm_marginal(qexp, pexp, dexp, rate = 1/8),
    bru_inits = list(
      gammaA = qnorm(pgamma(2, shape=2, rate=1)),
      gammaB = qnorm(pgamma(2, shape=2, rate=1)),
      sigmaA = qnorm(pexp(2, rate = 1/8)),
      sigmaB = qnorm(pexp(2, rate = 1/8))
    )
  ),
  list(
    name = "Gamma Gamma(2,1) Sigma Exp(1/2)",
    gamma_prior = bm_marginal(qgamma, pgamma, dgamma, shape=2, rate=1),
    sigma_prior = bm_marginal(qexp, pexp, dexp, rate = 1/2),
    bru_inits = list(
      gammaA = qnorm(pgamma(2, shape=2, rate=1)),
      gammaB = qnorm(pgamma(3, shape=2, rate=1)),
      sigmaA = qnorm(pexp(2, rate = 1/2)),
      sigmaB = qnorm(pexp(2, rate = 1/2))
    )
  ),
  list(
    name = "Gamma Unif(0.01, 10) Sigma Exp(1/2)",
    gamma_prior = bm_marginal(qunif, punif, dunif, min = 0.01, max = 10),
    sigma_prior = bm_marginal(qexp, pexp, dexp, rate = 1/2),
    bru_inits = list(
      gammaA = qnorm(punif(2, min = 0.01, max = 10)),
      gammaB = qnorm(punif(3, min = 0.01, max = 10)),
      sigmaA = qnorm(pexp(2, rate = 1/2)),
      sigmaB = qnorm(pexp(2, rate = 1/2))
    )
  ),
  list(
    name = "Gamma Unif(0.01, 10) Sigma Exp(1/8)",
    gamma_prior = bm_marginal(qunif, punif, dunif, min = 0.01, max = 10),
    sigma_prior = bm_marginal(qexp, pexp, dexp, rate = 1/8),
    bru_inits = list(
      gammaA = qnorm(punif(2, min = 0.01, max = 10)),
      gammaB = qnorm(punif(3, min = 0.01, max = 10)),
      sigmaA = qnorm(pexp(2, rate = 1/8)),
      sigmaB = qnorm(pexp(2, rate = 1/8))
    )
  )
)

for (config in different_configs){
  
  
  rm(fit); invisible(gc())
  catt("fitting ", config$name)
  fit_start <- proc.time()
  fit <-  model_two_observers_hr(
    observed_points,
    ips = ips_with_detection_states,
    bru_verbose = 1,
    mtrn_prior = matern_prior,
    prior_on_gamma = config$gamma_prior,
    prior_on_sigma = config$sigma_prior,
    bru_initial_params = config$bru_inits
  )
  fit_end <- proc.time()
  fit_time <- (fit_end-fit_start)[[3]]
  print(fit_time)
  
  # # a picture
  # print(bru_convergence_plot(fit))
  
  start <- proc.time()
  list_of_models[[config$name]] <- get_preds_from_two_observers_hr_fit(fit)
  list_of_models[[config$name]]$name <- config$name
  list_of_models[[config$name]]$fit_time <- round(fit_time)
  end <- proc.time()
  print((end-start)[3])
  
  # plot some predictions
  plot_detect_pred(
    pred_df = list_of_models[[config$name]]$pred_detect,
    main = config$name
  )
  
}
