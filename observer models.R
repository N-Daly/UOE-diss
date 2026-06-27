model_one_observer_hn <- function(
    dd, 
    ground_truth, # for the mesh, buffered transects and so on
    mtrn_prior =matern_prior,
    half_width=8
){

  #assuming hn, log_hn, dist_to_nearest_line_transect are defined globally
  
  mesh1 <- ground_truth$the_mesh
  
  #save compute and space if a previous model already made it
  matern_prior = mtrn_prior
  # matern_prior <- inla.spde2.pcmatern(
  #   mesh1,
  #   prior.range = c(600, 0.1), # true rho is 500
  #   prior.sigma = c(.5, 0.5)   # true sigma is 1
  # )
  
  cmp <- ~ Intercept(1) +
    sigma(1,
           prec.linear = 1,
           marginal = bm_marginal(qexp, pexp, dexp, rate = 1/half_width) #need to think more about this - effect on detection at end of halfwidth
    ) +spde(main=geometry, model = matern_prior)
  
  
  form <- geometry  ~ Intercept  +
    log_hn(
      dist_to_nearest_line_transect(geometry, ground_truth$line_transects),
      sigma
    ) +log(2)+ spde
  
  
  fit <- lgcp(
    components = cmp,
    formula = form,
    data =  dd,
    domain = list(geometry=fm_subdivide(mesh1, 3) ),
    # the observed regions are line transects expanded outwards by the transect widths, right?
    samplers = ground_truth$buffered_transects,
    options = list(bru_verbose= 1,
                   # verbose = 4,
                   bru_initial = list(sigma = half_width/4) # need to review this 
    )
  )
  
  fit
}


model_two_observers_hn <- function(
  dd, 
  ground_truth, # for the mesh, buffered transects and other info
  mtrn_prior = matern_prior, # save memory
  half_width = 8
){
  mesh1 <- ground_truth$the_mesh
  # the space is now the 2D map and detection state
  the_whole_domain = list(geometry=fm_subdivide(mesh1, 3), detected=1:3)
  
  # this can be reused by models
  matern_prior <- mtrn_prior
  # matern_prior <- inla.spde2.pcmatern(
  #   mesh1,
  #   prior.range = c(600, 0.1), # true rho is 500
  #   prior.sigma = c(.5, 0.5)   # true sigma is 1
  # )
  
  
  cmp <- ~ Intercept(1) +
    sigmaA(1,
           prec.linear = 1,
           marginal = bm_marginal(qexp, pexp, dexp, rate = 1/8) 
    ) +
    sigmaB(1,
           prec.linear = 1,
           marginal = bm_marginal(qexp, pexp, dexp, rate = 1/8)
    )+
    spde(main=geometry, model = matern_prior)
  
  
  form <- geometry + detected  ~ Intercept  +
    log_g_2observer(
      dist_to_nearest_line_transect(geometry, ground_truth$line_transects),
      detected,
      sigmaA, sigmaB
    ) +log(2)+ spde
  
  
  fit_two_observers <- lgcp(
    components = cmp,
    formula = form,
    data =  dd,
    domain = the_whole_domain,
    # the observed regions are line transects expanded outwards by the transect widths, right?
    # not sure what'll happen with the detected states
    samplers = ground_truth$buffered_transects,
    options = list(bru_verbose= 1,
                   # verbose = 4,
                   bru_initial = list(sigmaA = half_width/4, sigmaB = half_width/4) # need to review this
    )
  )
  
  fit_two_observers
}
