rm(list=ls());while(dev.cur()>1){dev.off()};old_par<- par(no.readonly = T, pch=19);options(digits=3);invisible(gc())

my_dir <- r"(C:\Users\ND\OneDrive - University of Edinburgh\Dissertation\UOE-diss)"
setwd(my_dir)
source("function Definitions.R")
library(reshape2)
library(ggplot2)

dd <- seq(0, 8, length.out=1000)

pa <- hr(dd, 1.75, 1)
pb <- hr(dd, 2, 4)
pany <- 1 - (1-pa)*(1-pb)

# I tried plotting this in ggplot 
# it was too much pain trying to place the legend within the plot

setwd("figs")
pdf("true hr detect.pdf", width = 10, height = 7)
par(cex=1)
plot(
  dd,
  pany,
  ylim = 0:1,
  ylab = "Detection probability",
  xlab = "Perpendicular distance from observers in km",
  lwd = 2,
  type="l",
  main = "True probability of detection within the transect segment"
)
lines(dd, pa, col = "red", lwd = 2)
lines(dd, pb, col = "blue", lwd = 2)

legend(
  "topright",
  legend = c("Observer A and/or B", "Observer A", "Observer B"),
  col = c("black", "red", "blue"),
  lwd = "3"
)
dev.off()
setwd("..")

