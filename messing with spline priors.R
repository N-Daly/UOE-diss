rm(list=ls());while(dev.cur()>1){dev.off()};old_par<- par(no.readonly = T, pch=19)

library(INLA)
library(inlabru)
library(sf)
library(fmesher)
library(ggplot2)
library(patchwork)
library(latex2exp)
library(dplyr)

source("function Definitions.R")

dists_grid <- seq(0,8, length.out=1000)

mesh <- fm_mesh_1d(
  loc =  c(0,2,4,6,8),
  boundary = c("dirichlet", "free"),
  degree = 2
)
projector <- fm_evaluator(mesh, loc=dists_grid)


plot_the_variance <- function(
    rho, sigma,
    alpha = 2,
    the_mesh=mesh, distance=dists_grid
  ){
  
  # the neumann condition forces the first basis to zero
  # so here we take the conditional precision matrix of the other bases (base-ees?)
  # which happtens to be a submatrix of the joint precision matrix
  Q <- fm_matern_precision(the_mesh, rho, sigma, alpha)[-1, -1]
  
  # consider only the remaining bases
  basis <- fm_basis(the_mesh, distance)[, -1]
  
  Sigma <- diag( fm_covariance(Q, basis) )
  
  title <- TeX("Marginal variance of the quadratic B-spline $G(d)$")
  subtitle <- TeX(paste(
    "$\\alpha$ =", alpha,
    "$\\rho$ =", rho,
    "$\\sigma_{u}$ = ", sigma
  ))
  
  scale = 2
  plot(
    distance, Sigma,
    type="l",
    main = title,
    xlab = TeX("distance, $d$, in km"), # idk if this does anything tbh
    ylab = "variance",
    cex.lab = scale,
    cex.main = scale
  )
  mtext(subtitle, cex = 1.5)
}

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

rho <- 50; sigma <- 2; alpha = 2

# plot the marginal variance of the spline with given matern hyperparameters
setwd("figs")
pdf("splinePriorCovar.pdf", width=10, height=10)
plot_the_variance(rho = rho,sigma = sigma)
dev.off()
setwd("..")

# sample from the prior distribution
set.seed(1)
ss <- sample_many_detect_priors(rho=rho, sigma = sigma, nsims=200)

ss <- ss |> group_by(sim) |> dplyr::filter(all( prob <= 1) ) |> ungroup()

length(unique(ss$sim))
subtitle <- TeX(paste(
  "$\\alpha$ =", alpha,
  "$\\rho$ =", rho,
  "$\\sigma_{u}$ = ", sigma
))

spline_title <- TeX("Realisations of a quadratic B-spline $G(d)$ from its prior distribution")
detect_title <- TeX("Realisations of a spline detection function from its prior distribution" )

detect_plot <- ggplot(ss, aes(dist, prob, group=sim)) + 
  geom_line(alpha=.5) +
  expand_limits(y=0:1) + 
  labs(
    title = detect_title,
    subtitle = subtitle
  ) +
  ylab(TeX("$g(d; \\beta)$")) +
  xlab(TeX("distance, $d$, in km")) +
  theme(
    plot.title = element_text(size=rel(2)),
    plot.subtitle = element_text(size=rel(2)),
    axis.title.y  = element_text(size=rel(2)),
    axis.title.x  = element_text(size=rel(2)),
    axis.text.x = element_text(size=rel(1.5)),
    axis.text.y = element_text(size=rel(1.5))
  )

detect_plot

setwd("figs")
ggsave("BsplinePriors.pdf", w =10, h=10)
setwd("..")

spline_plot <- ggplot(ss, aes(dist,spline, group=sim)) + 
  geom_line(alpha=.4) +
  expand_limits(y=0:1) + 
  labs(
    title = spline_title,
    subtitle = subtitle
  ) +
  ylab(TeX("$G(d)$")) +
  xlab(TeX("distance, $d$, in km")) +
  theme(
    plot.title = element_text(size=rel(2)),
    plot.subtitle = element_text(size=rel(2)),
    axis.title.y  = element_text(size=rel(2)),
    axis.title.x  = element_text(size=rel(2)),
    axis.text.x = element_text(size=rel(1.5)),
    axis.text.y = element_text(size=rel(1.5))
  )

spline_plot

setwd("figs")
ggsave("splineDetectPriors.pdf", w =10, h=10)
setwd("..")


print(spline_plot + detect_plot)
