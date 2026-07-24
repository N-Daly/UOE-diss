make_spatially_varying_mesh <- function(num_vertices, sim_info){
  
  # this determines the target max edge length from each initially seeded vertex
  qual_loc <- function(locs){
    dist <- dist_to_nearest_line_transect(locs, sim_info$line_transects)
    #round the distance so that there are less unique edge lengths initially
    pmin(20, pmax(0.5, floor(dist) ))
  }
  
  # a grid of points across the study region, 40km apart
  seed_points <- fm_hexagon_lattice(sim_info$boundary, edge_len = 40)
  
  # also include points on the line transects themselves
  # so that in the following delauney triangulation the transects will be better captured
  pts_near_line <- st_cast(sim_info$line_transects$geometry, "POINT")
  seed_points <- c(seed_points, pts_near_line)
  
  
  # make a mesh with:
  # an extension outside the study region, around 10% farther
  # no triangle angles smaller than 23 degrees
  # there is some user set number of triangles, once reached the algorithm starts resizing edge lengths
  mesh <- fm_rcdt_2d(
    loc = seed_points
    , extend = list(offset=-.1),
    refine = list(
      min.angle = 23,
      max.n = num_vertices*1000
    ),
    quality.spec = list(
      loc = qual_loc(seed_points)
    ),
    crs=sim_info$the_mesh$crs
  )
  mesh
}

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
    ips,
    desired = c("detect", "lambda"),
    dists = seq(0,8, length.out=1000),
    m = 1000 # MC samples
){
  # assumes the detection functions are defined already
  
  lst <- list()
  # log lambda and lambda
  if ( "lambda" %in% desired){
    
    lst <- predict(
      fit,
      data.frame(geometry = ips$geometry),
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
      data.frame(distance=dists),
      formula = ~ hn(distance, sigma),
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
    typical_spde(main=geometry, model=matern_prior)
  
  form <- geometry  ~ Intercept  + typical_spde +
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
    ips,
    desired = c("detect", "lambda"),
    dists = seq(0,8, length.out=1000),
    m = 1000 # MC samples
){
  # assumes the detection functions are defined already
  
  lst <- list()
  
  # log lambda and lambda
  if ( "lambda" %in% desired){
    
    lst <- predict(
      fit,
      data.frame(geometry = ips$geometry),
      formula = ~ list(
        pred_lambda = exp(Intercept  + typical_spde),
        pred_loglambda = Intercept  + typical_spde
      ),
      n.samples = m
    )
  }
  
  # detection prob
  if ( "detect" %in% desired){
    
    lst$pred_detect <- predict(
      fit,
      data.frame(distance = dists),
      formula = ~ spline_detect_func(spline_spde),
      n.samples = m
    )
  }
  
  lst$name <- deparse(substitute(fit))
  lst
}


model_two_observers_hn <- function(
  dd, 
  ips,
  mtrn_prior = matern_prior, # save memory
  half_width = 8,
  bru_verbose = 0
){

  matern_prior <- mtrn_prior
  
  lambda_sigma <- 1/half_width
  
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
      bru_initial = list(
        sigmaA = qnorm(pexp(2, lambda_sigma)),
        sigmaB = qnorm(pexp(2, lambda_sigma))
      ) 
    )
  )
  
  fit_two_observers
}


get_preds_from_two_observers_hn_fit <- function(
    fit,
    ips,
    desired = c("detect", "lambda"),
    dists = seq(0,8, length.out=1000),
    m = 1000 # MC samples
){
  # assumes the detection functions are defined already
  
  lst <- list()
  # log lambda and lambda
  if ( "lambda" %in% desired){
    
    lst <- predict(
      fit,
      data.frame(geometry = ips$geometry),
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
      data.frame(distance = dists),
      formula = ~ detect_func_2observer_sigma(distance, sigmaA, sigmaB),
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
  
  cmp <- ~ Intercept(1) +
    sigma(1,
           prec.linear = 1,
           marginal = bm_marginal(qexp, pexp, dexp, rate = 1/half_width) 
    ) +
    spde(main=geometry, model = matern_prior)
  
  form <- geometry + detected  ~ Intercept  +
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

model_two_observers_hr_fixed_gamma <- function(
    dd, 
    ips,
    mtrn_prior = matern_prior, # save memory
    half_width = 8,
    bru_verbose = 0,
    fixed_gamma_val = 1
){
  
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
    bru_initial_params = list(),
    mtrn_prior = matern_prior, # save memory
    half_width = 8,
    bru_verbose = 0
){
  
  lambda_sigma <- 1/half_width
  
  bru_initial_params$sigmaA <- qnorm(pexp(2, lambda_sigma))
  bru_initial_params$sigmaB <- qnorm(pexp(2, lambda_sigma))
  
  cmp <- ~ Intercept(1) +
    sigmaA(1,
           prec.linear = 1,
           marginal = bm_marginal(qexp, pexp, dexp, rate = lambda_sigma) 
    ) +
    sigmaB(1,
           prec.linear = 1,
           marginal = bm_marginal(qexp, pexp, dexp, rate = lambda_sigma)
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
      gammaA, gammaB,
      eps=1e-9
    ) + 
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
  fit
}

