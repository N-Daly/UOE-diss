rm(list=ls());while(dev.cur()>1){dev.off()};old_par<- par(no.readonly = T, pch=19);options(digits=3);invisible(gc())

library(INLA)
library(inlabru)
library(sf)
library(fmesher)
library(ggplot2)
library(patchwork)
library(dplyr)
library(ggspatial)


my_dir <- r"(C:\Users\ND\OneDrive - University of Edinburgh\Dissertation\UOE-diss)"
setwd(my_dir)
source("simulate 2D ground truth.R")
source("observer models.R")
source("function Definitions.R")
source("mesh resolution from dist discrepancies.R")

# this is the one used for the model's spatial effects
mesh <- make_spatially_varying_mesh3(60, coarse=60)

tran_segm <- st_buffer(mexdolphin_sf$samplers, 8, endCapStyle = "FLAT")

ggplot() +
  gg(mesh, edge.color = "black") + 
  gg(tran_segm, alpha = .1, colour = "blue") + 
  gg(mexdolphin_sf$ppoly, alpha = 0.2,  color="red") +
  ggspatial::annotation_scale(data=mexdolphin_sf$ppoly) +
  ggtitle("Mesh of the study region") + 
  theme(
    plot.title=element_text(size=rel(2),face="bold"),
    plot.margin=grid::unit(c(0,0,0,0), "mm")
    )

ggsave("modelMesh.pdf", height = 10, width = 20)


# demonstrate mesh construct method 3 -  the two hexagonal lattices
mm <- make_spatially_varying_mesh3(5)

ggplot() + gg(mm) + 
  ggtitle("A mesh composed of fine and coarse hexagonal lattices") + 
  theme(
    plot.title=element_text(size=rel(2),face="bold")
  ) +
  gg(mexdolphin_sf$ppoly, col = "red", alpha= .2)

ggsave("meshConstructionMethod3.pdf", height = 15 ,width = 20)


r <- readRDS("01-08-2026 12-15 amended simulation results.rda")

g <- ggplot(r, aes(fill=model)) +
  expand_limits(y=0) +
  geom_abline(intercept = 0, slope = 0, colour="red", linewidth=2) +
  theme(
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    legend.text = element_text(size=rel(2)),
    legend.title = element_text(size=rel(2)),
    axis.title.y  = element_text(size=rel(2)),
    plot.title = element_text(size=rel(2))#,
    # plot.margin=grid::unit(c(0,0,0,0), "mm")
    ) + 
  scale_fill_discrete(
    name="Model",
    breaks=c("one_obs_HR", "one_obs_spline" , "two_obs_hn", "two_obs_HR"),
    labels=c("HR", "Spline", "Dual HN", "Dual HR")
  )
    
ds_ll <- g + geom_boxplot(aes(y=DS_loglambda)) +
  labs(
    title="Integrated Dawid-Sebastiani score on log lambda",
    y = "Differenced Dawid-Sebastiani score"
  ) 

se_ll <- g + geom_boxplot(aes(y=loglambda_SE)) + 
  labs(
    title="Integrated Squared Error on mean log lambda",
    y = "Differenced Squared Error score"
  )
se_ll <- se_ll + coord_cartesian(ylim=c(NA, 1e06))

ae_l <- g + geom_boxplot(aes(y=lambda_AE)) + 
  labs(
    title = "Integrated Absolute Error on median lambda",
    y = "Differenced Absolute Error score"  
  )

ds_avg_prob <- g + geom_boxplot(aes(y=DS_avg_prob)) + 
  labs(
    title="Dawid Sebastiani score on average detection probability",
    y = "Differenced Dawid-Sebastiani score"
  )

h <- 11; w<- 9
ds_ll
ggsave("results_DS_loglambda.pdf", height = h, width = w)

se_ll
ggsave("results_SE_loglambda.pdf", height = h, width = w)

ae_l
ggsave("results_AE_lambda.pdf", height = h, width = w)

ds_avg_prob
ggsave("results_DS_avg_prob.pdf", height = h, width = w)




