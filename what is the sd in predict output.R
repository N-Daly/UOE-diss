rm(list=ls());while(dev.cur()>1){dev.off()};old_par<- par(no.readonly = T, pch=19);options(digits=3)

library(INLA)
library(inlabru)
library(sf)
library(fmesher)
library(ggplot2)
library("patchwork")

set.seed(1)
n <- 100
xx <- runif(n)
yy <- rpois(n, lambda= exp(.5*xx))

fit <- bru(
  components = ~ beta(x),
  formula = y ~ beta,
  data = data.frame(y=yy, x=xx),
  family = "poisson"
)

fit

pred1 <- predict(
  fit, 
  data.frame(x=1:10),
  formula = ~  beta ,
  n.samples=1000
)
pred2 <- predict(
  fit, 
  data.frame(x=1:10),
  formula = ~  exp(beta) ,
  n.samples=1000
)
head(pred1)
pred1$sd / pred1$x
head(pred2)
pred2$sd / pred2$x

cbind(exp(pred1$sd), pred2$sd)
