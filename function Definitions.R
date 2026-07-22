###########  detection functions

hn <- function(distance, sigma){  exp(-0.5 * (distance/sigma)^2 )  }
log_hn <-  function(distance, sigma){ -0.5 * (distance/sigma)^2 }


######## the PDFs for both detection state and distance

log_g_2observer <- function(distance, detected, sigmaA, sigmaB, eps=1e-6){
  
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


