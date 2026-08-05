rm(list=ls());while(dev.cur()>1){dev.off()};old_par<- par(no.readonly = T, pch=19)

library(INLA)
library(inlabru)
library(sf)
library(fmesher)
library(ggplot2)

source("function Definitions.R")

dists_grid <- seq(0,8, length.out=1000)

mesh <- fm_mesh_1d(
  loc =  c(0,2,4,6,8),
  boundary = c("dirichlet", "free"),
  degree = 2
)
projector <- fm_evaluator(mesh, loc=dists_grid)


#not sure what good this is without accounting for the neumann condition
plot_the_variance <- function(
    the_mesh=mesh, distance=dists_grid,
    alpha, rho, sigma
  ){
  
  Q <- fm_matern_precision(the_mesh, alpha, rho, sigma)
  
  basis <- fm_basis(the_mesh, distance)
  
  Sigma <- diag( fm_covariance(Q, basis) )
  
  title = "Marginal variance of the Spline"
  subtitle = paste(
    " alpha =",alpha,
   "rho =", rho,
   "sigma = ", sigma
  )
  plot(
    distance, Sigma,
    # ylim=0:1,
    type="l",
    main = title
  )
  mtext(subtitle)
}
plot_the_variance(alpha=2,rho=1,sigma=1)


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
    prob = c(detect_values),
    dist = rep(distances, times = ncol(detect_values)),
    sim = rep(1:nsims, each = nrow(detect_values))
  )
  df
}

rho <- 3; sigma <- 0.5
ss <- sample_many_detect_priors(rho=rho, sigma = sigma, nsims=200)

title <- paste(
  "Realisations of spline detect function with rho =",
  rho,
  "and sigma =",
  sigma
)
(g <- ggplot(ss, aes(dist,prob, group=sim))  + 
  geom_line(alpha=.4) +
  expand_limits(y=0:1) + 
  ylim( 0, min(10, max(ss$prob)) )  +
  labs(title = title)
)






##############
# some data

hn <- function(distance, sigma){  exp(-0.5 * (distance/sigma)^2 )  }

dists <- seq(0, 8, length.out=1000)

set.seed(123)
toy_data <- sample(dists, 100, replace=T)

sigmaA <- 2; sigmaB <- 4

pA <- hn(toy_data, sigmaA); pB <- hn(toy_data, sigmaB)
# detected
ii <- (runif(length(toy_data)) <= pA) | (runif(length(toy_data)) <= pB)
toy_data <- data.frame(distance = toy_data[ii])


############# bru model
dof <- fm_dof(mesh)
detect_spde <- inla.spde2.pcmatern(
  mesh,
  prior.range = c(rho, 0.9), # P(rho < val) = p0
  prior.sigma = c(sigma, 0.2), # P(sigma > val) = p0
  extraconstr = list(
    A=matrix(c(1, rep(0, dof-1) ), 1, dof),
    e = matrix(0, 1, 1)
  )
)

cmp <- ~ spline_effect(main=distance, model=detect_spde) +
  Intercept(1)

form <- distance ~ -spline_effect + Intercept

fit <- lgcp(
  components = cmp,
  formula = form,
  data = toy_data,
  domain = list(distance=fm_mesh_1d( 0:8 )),
  options = list(bru_verbose=1)
)


dp <- predict(
  fit,
  data.frame(distance=dists),
  formula = ~ exp(-spline_effect)
)


plot(
  dists, 
  dp$mean,
  type = "l",
  lwd = 2,
  ylim = range(0,dp$q0.975)
)
lines(dists, dp$q0.975)
lines(dists, dp$q0.025)

# plot(dists, -log(dp$mean))





