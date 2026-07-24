###########  detection functions

hn <- function(distance, sigma){  exp(-0.5 * (distance/sigma)^2 )  }
log_hn <-  function(distance, sigma){ -0.5 * (distance/sigma)^2 }

hr <- function(distance, sigma, gamma) {
  # catt("hr arg lengths:", length(distance), length(sigma), length(gamma))
  
  if (any(is.na(sigma))){
    catt("hr : got ", sum(is.na(sigma)), " NAs")
  }
  
  1 - exp( -(distance / sigma)^-gamma )
}
log_hr <- function(distance, sigma, gamma){
  # do (dist/sigma)^-b on log scale
  # frac_part <- -exp( -gamma*log(distance) + gamma*log(sigma) )
  # 
  # log1p( -exp(frac_part) )
  log(hr(distance, sigma, gamma))
}


######## the PDFs for both detection state and distance

log_g_2observer_hn <- function(distance, detected, sigmaA, sigmaB, eps=1e-9){
  
  if ( length(sigmaA) == 1 & length(distance) > 1 ){
    sigmaA <- rep(sigmaA, length(distance))
    sigmaB <- rep(sigmaB, length(distance))  
  }
  
  log_terms <- numeric(length(detected))
  # A alone
  ii <- detected==1 
  log_terms[ii] = log_hn(distance[ii], sigmaA) + log1p(-hn(distance[ii], sigmaB)*(1-eps) )
  # B alone
  ii <- detected==2
  log_terms[ii] = log1p(-hn(distance[ii], sigmaA)*(1-eps) ) + log_hn(distance[ii], sigmaB)
  # A,B
  ii <- detected==3
  log_terms[ii] = log_hn(distance[ii], sigmaA) + log_hn(distance[ii], sigmaB)
  
  log_terms
}

log_g_2observer_hr <- function(distance, detected, sigmaA, sigmaB, gammaA, gammaB, eps=1e-9){
  # 
  # cat(
  #   "log g arg lengths: dist", length(distance),
  #   "detect", length(detected),
  #   "sigmaA, sigmaB", length(sigmaA), length(sigmaB),
  #   "gammaA, gammaB", length(gammaA), length(gammaB),
  #   "\n"
  # )
  # 

  if( length(sigmaA) == 1 & length(distance) > 1){
    sigmaA <- rep(sigmaA, length(distance))
    sigmaB <- rep(sigmaB, length(distance))
  }

  if( length(gammaA) == 1 & length(sigmaA) > 1){
    gammaA <- rep(gammaA, length(sigmaA))
    gammaB <- rep(gammaB, length(sigmaB))
  }
  
  log_terms <- numeric(length(detected))
  # A alone
  ii <- detected==1
  log_terms[ii] = log_hr(distance[ii], sigmaA[ii], gammaA[ii]) + log1p( -hr(distance[ii], sigmaB[ii], gammaB[ii])*(1-eps) ) 
  # B alone
  ii <- detected==2
  log_terms[ii] = log1p( -hr(distance[ii], sigmaA[ii], gammaA[ii])*(1-eps) ) + log_hr(distance[ii], sigmaB[ii], gammaB[ii])
  # Both A,B
  ii <- detected==3
  log_terms[ii] =  log_hr(distance[ii], sigmaA[ii], gammaA[ii]) + hr(distance[ii], sigmaB[ii], gammaB[ii])
  
  log_terms
}

spline_detect_func <- function(spline_effect){ exp(-spline_effect) }
ln_spline_detect_func <- function(spline_effect){ -spline_effect }

####### the probability of detection at a given distance with given parameters

detect_func_2observer_sigma <- function(distance, sigmaA, sigmaB){
  #ln( (1-pA)(1-pB) )
  log_terms <- log1p( -hn(distance, sigmaA) ) + log1p( -hn(distance, sigmaB) )
  # cant think of a better way to do this bit
  #  1-terms  =  1- e^log_terms =  -1*( e^logterms - 1 ) 
  # my assumption is that if the probs are high this will be stable for the small log terms
  # and if the probs are low this will be also fine for the bigger log terms
  -expm1(log_terms)
}

detect_func_2_observer_hr <- function(distance, sigmaA, sigmaB, gammaA, gammaB){
  # g = 1 - (1-pA)*(1-pB)
  #  for the hazard rate pdf
  # pA = 1-exp( frac )
  # so g becomes 1 - exp[ log(1-pA) + log(1-pB) ]
   # = 1 - exp( fracA + fracB )
  
  # do the fraction (dist/sigma)^-b itself on the log scale
  fracA <- -exp( -gammaA*log(distance) + gammaA*log(sigmaA) )
  fracB <- -exp( -gammaB*log(distance) + gammaB*log(sigmaB) )
  
  1 - exp(fracA + fracB)
}

######## Scoring rules

dawid_sebastiani_score <- function(post_pred, true_value){
  E <- post_pred$mean
  # i am not sure if this step is correct
  V <- post_pred$sd^2
  ( (true_value - E)^2 / V)  + log(V)
}


######### distance to transect
dist_to_nearest_line_transect <- function(pts, transects = mexdolphin_sf$samplers){
  
  # this should be in kilometres as the crs is in km units
  apply(st_distance(pts, transects), 1, min)
  
}


catt <- function(...){ cat(..., "\n")}# pet peeve


