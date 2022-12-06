## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- fig.show='hold'---------------------------------------------------------
plot(1:10)
plot(10:1)

## ---- echo=FALSE, results='asis'----------------------------------------------
knitr::kable(head(mtcars, 10))

## ---- eval = FALSE------------------------------------------------------------
#  library(roleLite)
#  library(pika)
#  
#  S <- 500
#  aa <- matrix(runif(S^2, 0, 5), nrow = S)
#  diag(aa) <- runif(S, 1, 10)
#  
#  # aa <- matrix(1, nrow = S, ncol = S)
#  
#  x <- rlseries(S, 0.001)
#  
#  boo <- compSim(J = 1000, nu = 0, m = 0.01, alpha = aa,
#                 niter = 100000, niterTimestep = 100,
#                 initMeta = x)
#  
#  # boo <- untbSim(J = 1000, nu = 0, m = 0.1,
#  #                niter = 10000, niterTimestep = 100,
#  #                initMeta = x)
#  
#  doo <- getAbund(boo)
#  
#  
#  plot(doo$tstep, doo$abund, type = 'n', log = 'y')
#  for(i in 1:S) {
#      lines(doo$tstep[doo$spID == i],
#            doo$abund[doo$spID == i],
#            col = viridis::viridis(S)[i])
#  }
#  
#  plot(sad(getAbund(boo, nrow(boo))$abund), ptype = 'rad', log = 'y')

