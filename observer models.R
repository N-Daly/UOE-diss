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
    mtrn_prior =matern_prior,
    half_width=8
){

  #assuming hn, log_hn, dist_to_nearest_line_transect are defined globally
  
  mesh1 <- construction_info$the_mesh
  
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
    ) +
    spde(main=geometry, model = matern_prior) +
    dist_to_transect(
      dist_to_nearest_line_transect(geometry, construction_info$line_transects),
      model = "const"
    )
  
  
  
  form <- geometry  ~ Intercept  +
    log_hn(
      dist_to_transect,
      sigma
    ) +log(2)+ spde
  
  
  fit <- lgcp(
    components = cmp,
    formula = form,
    data =  dd,
    domain = list(geometry=fm_subdivide(mesh1, 3) ),
    # the observed regions are line transects expanded outwards by the transect widths, right?
    samplers = construction_info$buffered_transects,
    options = list(bru_verbose= 0,
                   # verbose = 4,
                   bru_initial = list(sigma = half_width/4) # need to review this 
    )
  )
  
  fit
}


model_two_observers_hn <- function(
  dd, 
  construction_info, # for the mesh, buffered transects and other info
  mtrn_prior = matern_prior, # save memory
  half_width = 8
){
  mesh1 <- construction_info$the_mesh
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
    spde(main=geometry, model = matern_prior) +
    dist_to_transect(
      dist_to_nearest_line_transect(geometry, construction_info$line_transects),
      model = "const"
    )
  
  
  form <- geometry + detected  ~ Intercept  +
    log_g_2observer(
      dist_to_transect,
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
    samplers = construction_info$buffered_transects,
    options = list(bru_verbose= 0,
                   # verbose = 4,
                   bru_initial = list(sigmaA = half_width/4, sigmaB = half_width/4) # need to review this
    )
  )
  
  fit_two_observers
}
