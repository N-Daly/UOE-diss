rm(list=ls());while(dev.cur()>1){dev.off()};old_par<- par(no.readonly = T, pch=19)

library(INLA)
library(inlabru)
library(sf)
library(fmesher)
library(ggplot2)
library(latex2exp)

my_dir <- r"(C:\Users\ND\OneDrive - University of Edinburgh\Dissertation\UOE-diss)"
setwd(my_dir)
source("function Definitions.R")

dists_grid <- seq(0,8, length.out=1000)

mesh <- fm_mesh_1d(
  loc =  c(0,2,4,6,8),
  boundary = c("dirichlet", "free"),
  degree = 2
)
projector <- fm_evaluator(mesh, loc=dists_grid)


plot_the_variance <- function(
    the_mesh=mesh, distance=dists_grid,
    alpha, rho, sigma
  ){
  
  Q <- fm_matern_precision(the_mesh, alpha, rho, sigma)
  
  basis <- fm_basis(the_mesh, distance)
  
  Sigma <- diag( fm_covariance(Q, basis) )
  
  title <- "Marginal variance of the Spline"
  subtitle <- TeX(paste(
    "$\\alpha$ =", alpha,
    "$\\rho$ =", rho,
    "$\\sigma$ = ", sigma
  ))
  
  plot(
    distance, Sigma,
    # ylim=0:1,
    type="l",
    main = title
  )
  mtext(subtitle, cex = 1.5)
}

#not sure what good this is without accounting for the neumann condition
# plot_the_variance(alpha=2,rho=1,sigma=1)


sample_many_detect_priors <- function(
    rho, sigma,
    the_mesh = mesh, proj = projector, distances = dists_grid,
    alpha = 2, nsims = 100
  ){
  
  Q <- fm_matern_precision(the_mesh, alpha, rho, sigma)
  
  dof <- fm_dof(the_mesh)
  
  prior_weights <- fm_sample(
    n=nsims,
    Q=Q,
    constr = list(
      A = matrix(c(1, rep(0, dof-1) ), 1, dof),
      e = matrix(0,1,1)
    )
  )
  
  spline_values <- fm_evaluate(proj=proj, field=prior_weights)
  
  detect_values <- spline_detect_func(spline_values)

  df <- data.frame(
    spline = c(spline_values),
    prob = c(detect_values),
    dist = rep(distances, times = ncol(detect_values)),
    sim = rep(1:nsims, each = nrow(detect_values))
  )
  df
}

rho <- 3; sigma <- 0.75
ss <- sample_many_detect_priors(rho=rho, sigma = sigma, nsims=200)

title <- TeX(paste(
  "Realisations of spline detection function with $\\alpha = 2, \\rho =$",
  rho,
  "and $\\sigma =$",
  sigma
))
(g <- ggplot(ss, aes(dist,prob, group=sim)) + 
  geom_line(alpha=.4) +
  expand_limits(y=0:1) + 
  ylim( 0, min(6, max(ss$prob)) )  +
  labs(title = title) +
  # geom_line(data=data.frame(prob=1, dist = dists_grid, sim=-1), colour= "red")+
    geom_hline(yintercept= 1, linewidth = 2, colour = "red")
)
