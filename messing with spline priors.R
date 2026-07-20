rm(list=ls());while(dev.cur()>1){dev.off()};old_par<- par(no.readonly = T, pch=19)

library(INLA)
library(inlabru)
library(sf)
library(fmesher)
library(ggplot2)

dists_grid <- seq(0,8, length.out=1000)

spline_to_detect_func <- function(x){
  # converts the value of the spline at a point to the corresponding value
  # of the detection function at that point
  exp(-x)
}


mesh <- fm_mesh_1d(
  loc = 0:8,
  boundary = c("dirichlet", "free"),
  degree = 2
)
projector <- fm_evaluator(mesh, loc=dists_grid)

ggplot() + geom_fm(data= mesh)

plot_the_variance <- function(
    the_mesh=mesh, distance=dists_grid,
    alpha, rho, sigma
  ){
  
  Q <- fm_matern_precision(the_mesh, alpha, rho, sigma)
  
  basis <- fm_basis(the_mesh, distance)
  
  Sigma <- diag( fm_covariance(Q, basis) )
  
  title = "Marginal variance of the Spline representation"
  subtitle = paste(
    " alpha =",alpha,
   "rho =", rho,
   "sigma = ", sigma
  )
  plot(
    distance, Sigma,
    ylim=0:1,
    type="l",
    main = title
  )
  mtext(subtitle)
}

sample_a_detect_func_prior <- function(alpha, rho, sigma, evaluator=projector, the_mesh=mesh){
  
  prior_weights <- fm_matern_sample(
    the_mesh, 
    alpha=alpha, rho = rho, sigma = sigma,
  
  )
  spline_values <- fm_evaluate(proj=projector, field=prior_weights)
  spline_to_detect_func(spline_values)
  # spline_values
}

# plot(dists_grid, sample_a_detect_func_prior(2,1,1))

sample_detect_funcs_and_plot <- function(
    nreps=100, distances=dists_grid,
    alpha, rho, sigma
  ){
  many_dfs <- lapply(
    1:nreps,
    function(i){
      data.frame(
        x = distances,
        y = sample_a_detect_func_prior(alpha, rho, sigma),
        iter=i
      )
    }
  )
  samples <- do.call(rbind, many_dfs)
  
  main = paste(
    nreps, 
    "realisations of Matern with alpha =",
    alpha,
    "rho =", rho,
    "sigma = ", sigma
  )
  g<- ggplot(samples, aes(x,y, group=iter))  + 
    geom_line(alpha=.4) +
    expand_limits(y=0:1) + ylim(c(0, 2)) +
    labs(title = main)
  print(g)
  #samples
}


set.seed(123)
sample_detect_funcs_and_plot(
  nreps=100,
  alpha=2, rho=3, sigma=.25
)


plot_the_variance(
  alpha=2, rho=2, sigma=.25
)





















