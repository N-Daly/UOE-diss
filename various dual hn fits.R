different_configs <- list(
  list(
    name = "Gamma(2,1)",
    sigma_prior = bm_marginal(qgamma, pgamma, dgamma, shape=2, rate=1),
    bru_inits = list(
      sigmaA = qnorm(pgamma(2, shape=2, rate=1)),
      sigmaB = qnorm(pgamma(2, shape=2, rate=1))
    )
  ),
  list(
    name = "Exp(1/8)",
    sigma_prior = bm_marginal(qexp, pexp, dexp, rate = 1/8),
    bru_inits = list(
      sigmaA = qnorm(pexp(2, rate = 1/8)),
      sigmaB = qnorm(pexp(2, rate = 1/8))
    )
  ),
  list(
    name = "Exp(1/2)",
    sigma_prior = bm_marginal(qexp, pexp, dexp, rate = 1/2),
    bru_inits = list(
      sigmaA = qnorm(pexp(2, rate = 1/2)),
      sigmaB = qnorm(pexp(2, rate = 1/2))
    )
  ),
  list(
    name = "Unif(0.001, 10)",
    sigma_prior = bm_marginal(qunif, punif, dunif, min=0.001, max=10),
    bru_inits = list(
      sigmaA = qnorm(punif(2, min=0.001, max=10)),
      sigmaB = qnorm(punif(2, min=0.001, max=10))
    )
  )
)

for (config in different_configs){
  
  
  rm(fit); invisible(gc())
  catt("fitting dual HN", config$name)
  fit_start <- proc.time()
  fit <-  model_two_observers_hn(
    observed_points,
    ips = ips_with_detection_states,
    bru_verbose = how_verbose,
    mtrn_prior = matern_prior,
    prior_on_sigma = config$sigma_prior,
    bru_initial_params = config$bru_inits
  )
  fit_end <- proc.time()
  fit_time <- (fit_end-fit_start)[[3]]
  print(fit_time)
  
  # # a picture
  # print(bru_convergence_plot(fit))
  
  start <- proc.time()
  list_of_models[[config$name]] <- get_preds_from_two_observers_hn_fit(fit)
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
