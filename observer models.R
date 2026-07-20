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
    construction_info, # for the mesh, buffered transects and so on
    ips,
    mtrn_prior =matern_prior,
    half_width=8,
    bru_verbose = 0
){

  #assuming hn, log_hn, dist_to_nearest_line_transect are defined globally
  
  #save compute and space if a previous model already made it
  matern_prior = mtrn_prior
  # matern_prior <- inla.spde2.pcmatern(
  #    construction_info$the_mesh,
  #   prior.range = c(600, 0.1), # true rho is 500
  #   prior.sigma = c(.5, 0.5)   # true sigma is 1
  # )
  
  cmp <- ~ Intercept(1) +
    sigma(1,
           prec.linear = 1,
           marginal = bm_marginal(qexp, pexp, dexp, rate = 1/half_width) #need to think more about this - effect on detection at end of halfwidth
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
    options = list(bru_verbose= bru_verbose,
                   # verbose = 4,
                   bru_initial = list(sigma = half_width/4) # need to review this 
    )
  )
  
  fit
}

get_preds_from_one_observer_hn_fit <- function(
    fit,
    ips,
    dists = seq(0,8, length.out=1000),
    m = 1000 # MC samples
){
  # assumes the detection functions are defined already
  
  # log lambda and lambda
  lst <- predict(
    fit,
    data.frame(geometry = ips$geometry),
    formula = ~ list(
        pred_lambda = exp(Intercept + spde),
        pred_loglambda = Intercept + spde
    ),
    n.samples = m
  )
  
  # detection prob
  lst$pred_detect <- predict(
    fit,
    data.frame(distance=dists),
    formula = ~ hn(distance, sigma),
    n.samples = m
  )
  
  lst$name <- deparse(substitute(fit))
  lst
}

model_one_observer_spline <- function(
    dd,
    construction_info,
    ips,
    mtrn_prior = matern_prior,
    mtrn_detect_prior = detect_matern,
    half_width = 8,
    bru_verbose = 0
){
  
  
  matern_prior = mtrn_prior
  # matern_prior <- inla.spde2.pcmatern(
  #   sim_info$the_mesh,
  #   prior.range = c(600, 0.1), # true rho is 500
  #   prior.sigma = c(.5, 0.5)   # true sigma is 1
  # )
  
  detect_matern = mtrn_detect_prior
  # # alpha=2, rho=2, sigma=.25
  # detect_matern <- inla.spde2.pcmatern(
  #   detect_mesh,
  #   alpha = 2,
  #   prior.range = c(2, 0.99), # P(rho < val) = alpha
  #   prior.sigma = c(0.25, 0.01) # P(sigma > val) = alpha
  # )
  
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
    dists = seq(0,8, length.out=1000),
    m = 1000 # MC samples
){
  # assumes the detection functions are defined already
  
  # log lambda and lambda
  lst <- predict(
    fit,
    data.frame(geometry = ips$geometry),
    formula = ~ list(
      pred_lambda = exp(Intercept  + typical_spde),
      pred_loglambda = Intercept  + typical_spde
    ),
    n.samples = m
  )
  
  # detection prob
  lst$pred_detect <- predict(
    fit,
    data.frame(distance = dists),
    formula = ~ spline_detect_func(spline_spde),
    n.samples = m
  )
  
  lst$name <- deparse(substitute(fit))
  lst
}


model_two_observers_hn <- function(
  dd, 
  construction_info, # for the mesh, buffered transects and other info
  ips,
  mtrn_prior = matern_prior, # save memory
  half_width = 8,
  bru_verbose = 0
){

  matern_prior <- mtrn_prior
  # matern_prior <- inla.spde2.pcmatern(
  #   construction_info$the_mesh,
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
      distance,
      detected,
      sigmaA, sigmaB
    ) +log(2)+ spde
  
  
  fit_two_observers <- lgcp(
    components = cmp,
    formula = form,
    data =  dd,
    ips = ips,
    options = list(bru_verbose= bru_verbose,
                   # verbose = 4,
                   bru_initial = list(sigmaA = half_width/4, sigmaB = half_width/4) # need to review this
    )
  )
  
  fit_two_observers
}


get_preds_from_two_observers_hn_fit <- function(
    fit,
    ips,
    dists = seq(0,8, length.out=1000),
    m = 1000 # MC samples
){
  # assumes the detection functions are defined already
  
  # log lambda and lambda
  lst <- predict(
    fit,
    data.frame(geometry = ips$geometry),
    formula = ~ list(
      pred_lambda = exp(Intercept + spde),
      pred_loglambda = Intercept + spde
    ),
    n.samples = m
  )
  
  # detection prob
  lst$pred_detect <- predict(
    fit,
    data.frame(distance = dists),
    formula = ~ detect_func_2observer_sigma(distance, sigmaA, sigmaB),
    n.samples = m
  )
  
  lst$name <- deparse(substitute(fit))
  lst
}



