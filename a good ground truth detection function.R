source("function Definitions.R")

dd <- seq(0, 8, length.out=1000)

pa <- hr(dd, 1.75, 1)
pb <- hr(dd, 2, 4)
pdual <- 1 - (1-pa)*(1-pb)

plot(
  dd,
  pdual,
  ylim = 0:1,
  ylab = "",
  xlab = "distance",
  lwd = 1
)
lines(dd, pa, col = "red")
lines(dd, pb, col = "blue")
