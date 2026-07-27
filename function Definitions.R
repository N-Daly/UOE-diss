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
  log_terms[ii] = log_hn(distance[ii], sigmaA[ii]) + log1p(-hn(distance[ii], sigmaB[ii])*(1-eps) )
  # B alone
  ii <- detected==2
  log_terms[ii] = log1p(-hn(distance[ii], sigmaA[ii])*(1-eps) ) + log_hn(distance[ii], sigmaB[ii])
  # A,B
  ii <- detected==3
  log_terms[ii] = log_hn(distance[ii], sigmaA[ii]) + log_hn(distance[ii], sigmaB[ii])
  
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
  V <- post_pred$sd^2
  ( (true_value - E)^2 / V)  + log(V)
}

get_scoring_differences <- function(
    models,
    ips,
    true_detect,
    true_loglambda
){
  base_mod <- models[[1]] # should always be the one observer hn
  
  # get the scores of this base model
  base_mod$DS <- dawid_sebastiani_score(base_mod$pred_loglambda, true_loglambda)
  
  base_mod$loglambda_SE <- (base_mod$pred_loglambda$mean - true_loglambda)^2
  
  base_mod$lambda_AE <- abs(base_mod$pred_lambda$median - true_loglambda)
  
  base_mod$detect_AE <- abs(base_mod$pred_detect$median - true_detect)
  
  score_diffs <- NULL
  for (i in 2:length(models)){
    mod <- models[[i]]
    
    #get this model' scores 
    mod$DS <- dawid_sebastiani_score(mod$pred_loglambda, true_loglambda)
    mod$loglambda_SE <- (mod$pred_loglambda$mean - true_loglambda)^2
    mod$lambda_AE <- abs(mod$pred_lambda$median - true_loglambda)
    mod$detect_AE <- abs(mod$pred_detect$median - true_detect)
    
    # now take the difference in scores between this and the base
    # and integrate it over the domain
    mod_diff <- list(
      DS = sum(ips$weight * (base_mod$DS - mod$DS)),
      loglambda_SE = sum(ips$weight * (base_mod$loglambda_SE - mod$loglambda_SE)),
      lambda_AE = sum(ips$weight * (base_mod$lambda_AE - mod$lambda_AE)),
      detect_AE = mean(base_mod$detect_AE - mod$detect_AE)
    )
    
    mod_diff$model <- mod$name
    
    score_diffs <- rbind(score_diffs, mod_diff)
  }
  
  as.data.frame(score_diffs, row.names = F)
}

######### general utitlities
dist_to_nearest_line_transect <- function(pts, transects = mexdolphin_sf$samplers){
  
  # this should be in kilometres as the crs is in km units
  apply(st_distance(pts, transects), 1, min)
  
}


catt <- function(...){ cat(..., "\n")}# pet peeve

# just plot the prediction of detection probabilites
plot_detect_pred <- function(pred_df, main = "", true=true_detect_prob, d= dists){
  plot(
    d,
    pred_df$mean,
    ylim = range(0:1, pred_df$q0.975),
    main = main,
    ylab = "Prob of detection",
    xlab = "distance"
  )
  lines(d, pred_df$q0.975)
  lines(d, pred_df$q0.025)
  
  lines(d, true, col = "red", lwd =2)
}

