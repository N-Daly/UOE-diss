# rm(list=ls());while(dev.cur()>1){dev.off()};old_par<- par(no.readonly = T, pch=19);options(digits=3)

library(INLA)
library(inlabru)
library(sf)
library(fmesher)
library(ggplot2)
library(patchwork)
library(ggspatial)

my_dir <- r"(C:\Users\ND\OneDrive - University of Edinburgh\Dissertation\UOE-diss)"
setwd(my_dir)
source("simulate 2D ground truth.R")
source("observer models.R")
source("function Definitions.R")
source("observer models.R")


get_nonoverlapping_samplers <- function(ltransects=mexdolphin_sf$samplers){

  nonoverlapping <- c(1)
  b_trans <- st_buffer(ltransects, dist = 8, endCapStyle = "ROUND")
  
  for (i in 2:nrow(ltransects)){
    l <- b_trans[i,]
    
    overlaps <- st_intersects(
      b_trans[nonoverlapping,],
      b_trans[i, ] ,
      sparse = F
    )
    
    if ( !any(overlaps) ){
      nonoverlapping <- rbind(nonoverlapping, i)
    }
  }
  ltransects[nonoverlapping,]
}

# reconstruct the points' perpendicular distance using the mesh
eval_mesh <- function(a_mesh){

  field <- dist_to_nearest_line_transect(fm_vertices(a_mesh))

   approximated_distance <- fm_evaluate(
    a_mesh, field = field, loc= pts$geometry
  )

  approximated_distance
}

# reconstruct the point's

# plot the reconstructed distances against actual ones
# and get the mean absolute deviation, as a percentage of the true value
plot_results <- function(reconstructed, actual = actual_distance){

  mean_percent_error <-  mean( abs(actual-reconstructed)/actual )
  
  plot(
    actual, reconstructed,
    xlim = range(0, actual),
    ylim = range(0, actual, reconstructed),
    main = paste(
      "mean abs % error", round(mean_percent_error, 2)
    )
  )
  abline(0:1, col ="red", lwd=2)
  
  mean_percent_error
}

# fit some model to see how long it takes with a given ips from a mesh
fit_model_with_it <- function(mesh, ips){
  
  
  coarse_hex <- fm_hexagon_lattice(mexdolphin_sf$ppoly, edge_len = 30)
  coarse_mesh <- fm_mesh_2d(coarse_hex, crs= fm_crs(mexdolphin_sf$mesh))

  matern_prior <- inla.spde2.pcmatern(
    coarse_mesh, # the integration scheme is on a finer resolution than the model
    prior.range = c(600, 0.1), # true rho is 500
    prior.sigma = c(.5, 0.5)   # true sigma is 1
  )

  fit <- model_one_observer_hr(
    sim_info$samples_df,
    mtrn_prior =  matern_prior,
    ips=ips,
    prior_on_gamma = bm_marginal(qgamma, pgamma, dgamma, shape=2, rate=1),
    bru_initial_params = list(
      sigma = qnorm(pexp(2, rate= 1/8)),
      gamma = qnorm(pgamma(2, shape=2, rate=1))
    ),
    bru_verbose = 1
  )
  
  fit
}

# it was checked if there was any variation in model fitting times for the same 
# mesh/ips. there is not.
repeated_model_fitting <- function(..., nreps=5){
  
  fit_times <- numeric(nreps)
  
  for (i in 1:nreps){ 
    start <- proc.time()
    fit <- fit_model_with_it(...)
    end <- proc.time()
    
    fit_times[i] <- (end-start)[[3]]
  }
  
  list(
    last_model_fit = fit,
    fit_times = fit_times
  )
}
################ the actual function ############

make_spatially_varying_mesh1 <- function(num_vertices){

  # this determines the target max edge length from each initially seeded vertex
  qual_loc <- function(locs){
    dist <- dist_to_nearest_line_transect(locs, mexdolphin_sf$samplers)
    #round the distance so that there are less unique edge lengths initially
    pmin(20, pmax(0.5, floor(dist) ))
  }

  # a grid of points across the study region, 40km apart
  seed_points <- fm_hexagon_lattice(st_buffer(mexdolphin_sf$ppoly, 100), edge_len = 40)

  # also include points on the line transects themselves
  # so that in the following delauney triangulation the transects will be better captured
  pts_near_line <- st_cast(mexdolphin_sf$samplers$geometry, "POINT")
  seed_points <- c(seed_points, pts_near_line)


  # make a mesh with:
  # an extension outside the study region, around 10% farther
  # no triangle angles smaller than 23 degrees
  # there is some user set number of triangles, once reached the algorithm starts resizing edge lengths
  mesh <- fm_rcdt_2d(
    loc = seed_points,
    ext= fm_nonconvex_hull(seed_points, convex=100),
    refine = list(
      min.angle = 23,
      max.n = num_vertices*1000
    ),
    quality.spec = list(
      loc = qual_loc(seed_points)
    ),
    crs=fm_crs(mexdolphin_sf$mesh)
  )
  mesh
}

make_spatially_varying_mesh2 <- function(param, nonoverlapping_line_transects=nonoverlapping){

  
  hex_points <- fm_hexagon_lattice(mexdolphin_sf$ppoly, edge_len = 20)
  # transect_points <- st_cast(nonoverlapping$geometry, "POINT")
  seed_points <- hex_points # c(transect_points, hex_points)
  
  extension_amount = 130 # seems to work well in practice
  
  ext <- fm_extensions(
    mexdolphin_sf$ppoly,
    convex = extension_amount,
    crs = fm_crs(mexdolphin_sf$mesh)
  )
  
  mesh1 <- fm_rcdt_2d(
    loc = seed_points,
    boundary = ext[[1]],
    refine=F,
    # extend = list(offset = -0.2),
    crs = fm_crs(mexdolphin_sf$mesh)
  )

  ## Save the resulting boundary
  boundary1 <- fm_segm(mesh1, boundary = TRUE)
  interior1 <- fm_segm(mesh1, boundary = FALSE)
  
  ## Triangulate inner domain
  mesh2 <- fm_rcdt_2d(
    loc = seed_points,
    interior = fm_segm( interior1, nonoverlapping_line_transects$geometry, is.bnd=F),
    boundary = boundary1,
    refine = list(
      min.angle = 27
    ),
    crs = fm_crs(mesh1)
  )
  

  boundary2 <- fm_segm(mesh2, boundary = TRUE)
  interior2 <- fm_segm(mesh2, boundary = FALSE)
  
  mesh3 <- fm_rcdt_2d(
    loc = rbind(seed_points, fm_vertices(mesh2)),
    interior = fm_segm(boundary2, interior2, is.bnd=F),
    refine = list(
      min.angle = 27
    ),
    crs = fm_crs(mesh1)
  )
  
  if (param == 0){
    mesh3
  } else{
    fm_subdivide(mesh3, param)
  }
}

make_spatially_varying_mesh3 <- function(hex_length, subdivisions=0, coarse_edge_len = 2*3*5){
  
  # buffer the transects a bit farther than the halfwidth
  bt <- st_buffer(mexdolphin_sf$samplers, 8*1.25)
  
  # make a fine lattice within the transect segments
  hex_points <- fm_hexagon_lattice(bt, edge_len = hex_length)
  
  # make a coarse lattice over the whole region to attempt to keep regularity
  extra_hex_points <- fm_hexagon_lattice(mexdolphin_sf$ppoly, edge_len = coarse_edge_len)
  
  # combine the two as seed points
  seed_points <- c(hex_points, extra_hex_points)
  
  extension_amount = 60 # seems to work well in practice

  ext <- fm_extensions(
    mexdolphin_sf$ppoly,
    convex = extension_amount,
    crs = st_crs(mexdolphin_sf$ppoly)
  )
  
  # make a relatively simple mesh from the lattices' points
  mm <- fm_mesh_2d(
    loc=seed_points,
    loc.domain = ext[[1]],
    min.angle = 27,
    crs = fm_crs(mexdolphin_sf$mesh)
  )
  
  
  
  if (subdivisions == 0){
    mm
  } else {
    fm_subdivide(mm, n = subdivisions)
  }
}

test_param_tradeoff <- function(
    func_name, params, subtitle, sim_info,
    verbose=T
  ){
  
  mesh_maker <-get(func_name)
  
  file_name <-  paste0(
    format(Sys.time(), "%d-%m-%Y %H-%M"),
    func_name,
    "tradeoff.pdf"
  )
  
  fit_time <- NULL
  error <- NULL
  
  for (param in params){
    cat("\n       ", func_name, "param =", param, "\n")
    
    rm(m, ips, some_model_fit); invisible(gc())
    
    start <- proc.time()
    cat("making mesh and ips \n")
    
    m <- mesh_maker( param )

    approx_dist <- eval_mesh(m)
  
    ips <- fm_int(
      m, samplers = sim_info$buffered_transects
    )
    ips$distance <- dist_to_nearest_line_transect(ips$geometry, sim_info$line_transects)
    
    end <- proc.time()
    print( (end-start)[3])
    
    #plot the error
    if (verbose){
      mpe <- plot_results(approx_dist)
      mtext(paste("param=", param, "# vertices =", nrow(fm_vertices(m)) ))
    } else {
      mpe <- mean( abs(actual_distance-approx_dist)/actual_distance )
    }
    
    
    #fit a simple model with spatial random effect and hn detection function
    start <- proc.time()
    some_model_fit <- fit_model_with_it(m, ips)
    end <- proc.time()
    run_time <- (end-start)[[3]]
    print(run_time)
  
    
    #record stats
    fit_time <- c(fit_time, run_time)
    error <- c(error, mpe)
    
  }
  
  # save the plot
  pdf(file_name)
  
  # i dont see a way of doing this in ggplot, nor do i see a potential benefit
  plot(
    error, 
    fit_time,
    # accross all methods and params, 80 puts them roughly on the same scale
    ylim = range(0, 75, fit_time), 
    xlim = range(0:1, error),
    type = "c", # lines near but not connecting each point
    xlab = "MAPE as %",
    ylab = "fitting time in seconds",
    main = "Trade off between computational time and mesh quality"
  )
  text(
    error, 
    fit_time,
    labels = params,
    col = rainbow(length(error)),
    cex=1.5
  )
  subtitle <- paste("Number indicates", subtitle)
  mtext(subtitle)
  
  # save the plot
  dev.off()
  
  data.frame(fit_time=fit_time, error=error, param=param)

}


if (T){
  
  ### get some data
  set.seed(123)
  sim_info <- simulate_lcgp(
    true_beta0 = -5,
    true_rho = 500, true_sigma_GRF=1
  )
  
  pts <- sim_info$samples_df
  actual_distance <- pts$distance
  
  set.seed(123)
  true_sigmaA <- 1.75; true_sigmaB <- 2; true_gammaA <- 1; true_gammaB <- 4
  true_range <- 500; true_sigma_grf <- 1
  
  sim_info <- simulate_lcgp_dual_obs_HR_thinning(
    true_sigmaA = true_sigmaA, true_sigmaB = true_sigmaB,
    true_gammaA = true_gammaA, true_gammaB = true_gammaB,
    true_beta0 = -4,
    true_rho = true_range, true_sigma_GRF = true_sigma_grf
  )
  
  nonoverlapping <- get_nonoverlapping_samplers()
  
  # V3 5:1
  # V2 0:4
  # V1 1:10*10
  
  different_setups <- list(
    list(func_name="make_spatially_varying_mesh3", params=1, subtitle = "lattice edge lengths within transect segments"),
    list(func_name="make_spatially_varying_mesh2", params=4, subtitle = "further mesh subdivisions"),
    list(func_name="make_spatially_varying_mesh1", params=10*10, subtitle = "number of vertices in each mesh, in 1000s")
  )
  
  setwd("sim_results")
  for (setup in different_setups){
    test_param_tradeoff(
      setup$func_name,
      setup$params,
      setup$subtitle,
      sim_info,
      verbose = T
    )
  }
  setwd("..")
}

# # sanity check all meshes cover the nominal study area
# mesh <- make_spatially_varying_mesh3(60, coarse=60)
# 
# mm <- make_spatially_varying_mesh(10)
# ggplot() + gg(mm) + gg(mesh, edge.color = "red")
# 
# mm <- make_spatially_varying_mesh2(0, get_nonoverlapping_samplers())
# ggplot() + gg(mm) + gg(mesh, edge.color = "red")
# 
# ggplot() + gg(mm) + gg(st_buffer(mexdolphin_sf$samplers, 8), alpha=.1)
# 
# mm <- make_spatially_varying_mesh3(5)
# ggplot() + gg(mm) + gg(mesh, edge.color = "red")



