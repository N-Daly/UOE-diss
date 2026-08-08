model_one_observer_hn <- function(
    dd, 
    ips,
    mtrn_prior =matern_prior,
    half_width=8,
    bru_verbose = 0
){

  #assuming hn, log_hn, dist_to_nearest_line_transect are defined globally
  
  #save compute and space if a previous model already made it
  matern_prior = mtrn_prior
  
  lambda_sigma <- 1/half_width
  
  cmp <- ~ Intercept(1) +
    sigma(1,
           prec.linear = 1,
           marginal = bm_marginal(qexp, pexp, dexp, rate = lambda_sigma) #need to think more about this - effect on detection at end of halfwidth
    ) +
    spde(main=geometry, model = matern_prior) 
  
  
  form <- geometry  ~ Intercept  +
    log_hn(distance, sigma) + log(2) + 
    spde
  
  
  fit <- lgcp(
    components = cmp,
    formula = form,
    data =  dd,
    ips = ips,
    options = list(
      bru_verbose= bru_verbose, 
      bru_initial = list(sigma = qnorm(pexp(2, lambda_sigma)) )
    )
  )
  
  fit
}

get_preds_from_one_observer_hn_fit <- function(
    fit,
    ips_for_pred = fm_int(list(geometry = mexdolphin_sf$mesh), samplers=mexdolphin_sf$ppoly),
    desired = c("detect", "lambda", "avg_prob"),
    dists_ips_for_pred = distance_ips,
    m = 750, # MC samples
    halfwidth = 8
){
  # assumes the detection functions are defined already
  
  lst <- list()
  # log lambda and lambda
  if ( "lambda" %in% desired){
    
    lst <- predict(
      fit,
      data.frame(geometry = ips_for_pred$geometry),
      formula = ~ list(
          pred_lambda = exp(Intercept + spde),
          pred_loglambda = Intercept + spde
      ),
      n.samples = m
    )
  }
  
  # detection prob  
  if( "detect" %in% desired ){
    lst$pred_detect <- predict(
      fit,
     data.frame(distance = dists_ips_for_pred$distance, weight = dists_ips_for_pred$weight),
      formula = ~ hn(distance, sigma),
      n.samples = m
    )
  }
  # average prob of detection
  if( "avg_prob" %in% desired ){
    lst$pred_avg_prob <- predict(
      fit, 
      data.frame(distance = dists_ips_for_pred$distance, weight = dists_ips_for_pred$weight),
      formula = ~ sum(weight*hn(distance, sigma))/halfwidth,
      n.samples = m
    )
  }
  
  lst$name <- deparse(substitute(fit))
  lst
}

model_one_observer_spline <- function(
    dd,
    ips,
    mtrn_prior = matern_prior,
    mtrn_detect_prior = detect_matern,
    half_width = 8,
    bru_verbose = 0
){
  
  
  matern_prior = mtrn_prior
  detect_matern = mtrn_detect_prior
  
  cmp <- ~ Intercept(1) +
    spline_spde(main=distance, model=detect_matern) +
    spde(main=geometry, model=matern_prior)
  
  form <- geometry  ~ Intercept  + spde + log(2) + 
    ln_spline_detect_func(spline_spde)
  
  
  fit <- lgcp(
    components = cmp,
    formula = form,
    data =  dd,
    ips=ips,
    options = list(bru_verbose= bru_verbose)
  )
  fit
  
}

get_preds_from_one_observer_spline_fit <- function(
    fit,
    ips_for_pred = fm_int(list(geometry = mexdolphin_sf$mesh), samplers=mexdolphin_sf$ppoly),
    desired = c("detect", "lambda", "avg_prob"),
    dists_ips_for_pred = distance_ips,
    m = 750, # MC samples
    halfwidth = 8
){
  # assumes the detection functions are defined already
  
  lst <- list()
  
  # log lambda and lambda
  if ( "lambda" %in% desired){
    
    lst <- predict(
      fit,
      data.frame(geometry = ips_for_pred$geometry),
      formula = ~ list(
        pred_lambda = exp(Intercept  + spde),
        pred_loglambda = Intercept  + spde
      ),
      n.samples = m
    )
  }
  
  # detection prob
  if ( "detect" %in% desired){
    
    lst$pred_detect <- predict(
      fit,
      data.frame(distance = dists_ips_for_pred$distance, weight = dists_ips_for_pred$weight),
      formula = ~ spline_detect_func(spline_spde),
      n.samples = m
    )
  }
  
  if ( "avg_prob" %in% desired ){
    
    lst$pred_avg_prob <- predict(
      fit,
      data.frame(distance = dists_ips_for_pred$distance, weight = dists_ips_for_pred$weight),
      formula = ~ sum( weight * spline_detect_func(spline_spde) )/halfwidth,
      n.samples = m
    )
  }
  
  
  lst$name <- deparse(substitute(fit))
  lst
}


model_two_observers_hn <- function(
  dd, 
  ips,
  prior_on_sigma,
  mtrn_prior = matern_prior,
  half_width = 8,
  bru_verbose = 0,
  bru_inits = list()
){

  matern_prior <- mtrn_prior
  
  if (identical(bru_inits, list())){
    lambda_sigma <- 2/half_width
    bru_inits <- list(
      sigmaA = qnorm(pexp(2, lambda_sigma)),
      sigmaB = qnorm(pexp(2, lambda_sigma))
    ) 
  }
  
  cmp <- ~ Intercept(1) +
    sigmaA(1,
           prec.linear = 1,
           marginal = prior_on_sigma 
    ) +
    sigmaB(1,
           prec.linear = 1,
           marginal = prior_on_sigma
    )+
    spde(main=geometry, model = matern_prior)

  
  form <- geometry + detected  ~ Intercept  +
    log_g_2observer_hn(
      distance,
      detected,
      sigmaA, sigmaB
    ) +log(2)+ spde
  
  
  fit_two_observers <- lgcp(
    components = cmp,
    formula = form,
    data =  dd,
    ips = ips,
    options = list(
      bru_verbose= bru_verbose,
      bru_initial = bru_inits
    )
  )
  
  fit_two_observers
}


get_preds_from_two_observers_hn_fit <- function(
    fit,
    ips_for_pred = fm_int(list(geometry = mexdolphin_sf$mesh), samplers=mexdolphin_sf$ppoly),
    desired = c("detect", "lambda", "avg_prob"),
    dists_ips_for_pred = distance_ips,
    m = 750, # MC samples
    halfwidth = 8
){
  # assumes the detection functions are defined already
  
  lst <- list()
  # log lambda and lambda
  if ( "lambda" %in% desired){
    
    lst <- predict(
      fit,
      data.frame(geometry = ips_for_pred$geometry),
      formula = ~ list(
        pred_lambda = exp(Intercept + spde),
        pred_loglambda = Intercept + spde
      ),
      n.samples = m
    )
  }
  
  # detection prob
  if ( "detect" %in% desired){
    
    lst$pred_detect <- predict(
      fit,
      data.frame(distance = dists_ips_for_pred$distance, weight = dists_ips_for_pred$weight),
      formula = ~ detect_func_2observer_hn(distance, sigmaA, sigmaB),
      n.samples = m
    )
  }
  
  # average prob of detection
  if( "avg_prob" %in% desired ){
    lst$pred_avg_prob <- predict(
      fit, 
      data.frame(distance = dists_ips_for_pred$distance, weight = dists_ips_for_pred$weight),
      formula = ~ sum(weight*detect_func_2observer_hn(distance, sigmaA, sigmaB))/halfwidth,
      n.samples = m
    )
  }
  
  lst$name <- deparse(substitute(fit))
  lst
}


model_one_observer_hr_fixed_gamma <- function(
    dd, 
    ips,
    mtrn_prior = matern_prior, # save memory
    half_width = 8,
    bru_verbose = 0,
    fixed_gamma_val = 1
){
  matern_prior <- mtrn_prior
  
  cmp <- ~ Intercept(1) +
    sigma(1,
           prec.linear = 1,
           marginal = bm_marginal(qexp, pexp, dexp, rate = 1/half_width) 
    ) +
    spde(main=geometry, model = matern_prior)
  
  form <- geometry ~ Intercept  +
    log_hr(distance, sigma, gamma=fixed_gamma_val) + 
    log(2) + spde
  
  
  fit <- lgcp(
    components = cmp, 
    formula = form,
    data = dd,
    ips = ips,
    options = list(
      bru_verbose = bru_verbose,
      bru_initial = list(sigma = half_width/4), # need to review this
      bru_max_iter = 20
    )
  )
}

model_one_observer_hr <- function(
    dd, 
    ips,
    prior_on_gamma,
    bru_initial_params = list(),
    mtrn_prior = matern_prior, # save memory
    half_width = 8,
    bru_verbose = 0
){
  matern_prior <- mtrn_prior
  
  cmp <- ~ Intercept(1) +
    sigma(1,
          prec.linear = 1,
          marginal = bm_marginal(qexp, pexp, dexp, rate = 1/half_width) 
    ) +
    gamma(
      1,
      prec.linear = 1,
      marginal = prior_on_gamma
    ) +
    spde(main=geometry, model = mtrn_prior)
  
  form <- geometry ~ Intercept  +
    log_hr(distance, sigma, gamma) + 
    log(2) + spde
  
  
  fit <- lgcp(
    components = cmp, 
    formula = form,
    data = dd,
    ips = ips,
    options = list(
      bru_verbose = bru_verbose,
      bru_initial = bru_initial_params,
      bru_max_iter = 20
    )
  )
}

get_preds_from_one_observer_hr_fit <- function(
    fit,
    ips_for_pred = fm_int(list(geometry = mexdolphin_sf$mesh), samplers=mexdolphin_sf$ppoly),
    fixed_gamma_val = NULL, # if NULL, gamma is assumed to be a r.v
    desired = c("detect", "lambda", "avg_prob"),
    dists_ips_for_pred = distance_ips,
    m = 750, # MC samples
    halfwidth = 8
){
  # assumes the detection functions are defined already
  
  lst <- list()
  
  # use the gamma value if given otherwise it should be defined as a model parameter
  if ( !is.null(fixed_gamma_val) ){
    gamma <- fixed_gamma_val
  }
  
  
  # log lambda and lambda
  if ( "lambda" %in% desired){
    
    lst <- predict(
      fit,
      data.frame(geometry = ips_for_pred$geometry),
      formula = ~ list(
        pred_lambda = exp(Intercept + spde),
        pred_loglambda = Intercept + spde
      ),
      n.samples = m
    )
  }
  
  # detection prob
  if ( "detect" %in% desired){
    
    lst$pred_detect <- predict(
      fit,
      data.frame(distance = dists_ips_for_pred$distance, weight = dists_ips_for_pred$weight),
      formula = ~ hr(distance, sigma, gamma),
      n.samples = m
    )
  }
  
  # average prob of detection
  if( "avg_prob" %in% desired ){
    lst$pred_avg_prob <- predict(
      fit, 
      data.frame(distance = dists_ips_for_pred$distance, weight = dists_ips_for_pred$weight),
      formula = ~ sum( weight * hr(distance, sigma, gamma) )/halfwidth,
      n.samples = m
    )
  }
  
  lst$name <- deparse(substitute(fit))
  lst
}


model_two_observers_hr_fixed_gamma <- function(
    dd, 
    ips,
    mtrn_prior = matern_prior, # save memory
    half_width = 8,
    bru_verbose = 0,
    fixed_gamma_val = 1
){
  matern_prior <- mtrn_prior
  lambda_sigma <- 2/half_width
  
  cmp <- ~ Intercept(1) +
    sigmaA(1,
           prec.linear = 1,
           marginal = bm_marginal(qexp, pexp, dexp, rate = lambda_sigma) 
    ) +
    sigmaB(1,
           prec.linear = 1,
           marginal = bm_marginal(qexp, pexp, dexp, rate = lambda_sigma)
    )+
    spde(main=geometry, model = matern_prior)
  
  form <- geometry + detected  ~ Intercept  +
    log_g_2observer_hr(
      distance,
      detected,
      sigmaA, sigmaB,
      gammaA = fixed_gamma_val, gammaB = fixed_gamma_val
    ) + log(2) + spde
  
  
  fit <- lgcp(
    components = cmp, 
    formula = form,
    data = dd,
    ips = ips,
    options = list(
      bru_verbose = bru_verbose,
      bru_initial = list( # need to review this
        sigmaA = qnorm(pexp(2, lambda_sigma)), 
        sigmaB = qnorm(pexp(2, lambda_sigma))
        ), 
      bru_max_iter = 20
      )
  )
  fit
}

model_two_observers_hr <- function(
    dd, 
    ips,
    prior_on_gamma,
    prior_on_sigma,
    bru_initial_params = list(),
    mtrn_prior = matern_prior, # save memory
    half_width = 8,
    bru_verbose = 0
){
  matern_prior <- mtrn_prior
  lambda_sigma <- 1/half_width
  
  # bru_initial_params$sigmaA <- qnorm(pexp(2, lambda_sigma))
  # bru_initial_params$sigmaB <- qnorm(pexp(2, lambda_sigma))
  
  cmp <- ~ Intercept(1) +
    sigmaA(1,
           prec.linear = 1,
           marginal = prior_on_sigma# bm_marginal(qexp, pexp, dexp, rate = lambda_sigma) 
    ) +
    sigmaB(1,
           prec.linear = 1,
           marginal = prior_on_sigma # bm_marginal(qexp, pexp, dexp, rate = lambda_sigma)
    )+ 
    gammaA(1,
           prec.linear = 1,
           marginal = prior_on_gamma
    )+ 
    gammaB(1,
           prec.linear = 1,
           marginal = prior_on_gamma
    )+
    spde(main=geometry, model = matern_prior)
  
  form <- geometry + detected  ~ Intercept  +
    log_g_2observer_hr(
      distance,
      detected,
      sigmaA, sigmaB,
      gammaA, gammaB
    ) + 
    log(2) + spde
  
  
  fit <- lgcp(
    components = cmp, 
    formula = form,
    data = dd,
    ips = ips,
    options = list(
      verbose = F,
      bru_verbose = bru_verbose,
      bru_initial = bru_initial_params,
      bru_max_iter = 20
    )
  )
  fit
}

get_preds_from_two_observers_hr_fit <- function(
    fit,
    ips_for_pred = fm_int(list(geometry = mexdolphin_sf$mesh), samplers=mexdolphin_sf$ppoly),
    fixed_gamma_val = NULL, # if NULL, gamma is assumed to be a r.v
    desired = c("detect", "lambda", "avg_prob"),
    dists_ips_for_pred = distance_ips,
    m = 750, # MC samples
    halfwidth = 8
){
  # if NULL, gamma is assumed to be a r.v
  if ( !is.null(fixed_gamma_val)){
    gammaA <- fixed_gamma_val
    gammaB <- fixed_gamma_val
  }
  
  lst <- list()
  # log lambda and lambda
  if ( "lambda" %in% desired){
    
    lst <- predict(
      fit,
      data.frame(geometry = ips_for_pred$geometry),
      formula = ~ list(
        pred_lambda = exp(Intercept + spde),
        pred_loglambda = Intercept + spde
      ),
      n.samples = m
    )
  }
  
  # detection prob
  if ( "detect" %in% desired){

      lst$pred_detect <- predict(
      fit,
      data.frame(distance = dists_ips_for_pred$distance, weight = dists_ips_for_pred$weight),
      formula = ~ detect_func_2_observer_hr(distance, sigmaA, sigmaB, gammaA, gammaB),
      n.samples = m
    )
  }
  
  # average prob of detection
  if( "avg_prob" %in% desired ){
    lst$pred_avg_prob <- predict(
      fit, 
      data.frame(distance = dists_ips_for_pred$distance, weight = dists_ips_for_pred$weight),
      formula = ~ sum( weight * detect_func_2_observer_hr(distance, sigmaA, sigmaB, gammaA, gammaB) )/halfwidth,
      n.samples = m
    )
  }
 
  lst$name <- deparse(substitute(fit))
  lst
}


