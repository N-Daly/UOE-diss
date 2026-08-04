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
  
  # this is the normalising constant to make the terms probabilities
  # i.e what is the prob of detected=i given the object was detected at all?
  log_any_detection <- log1p( -(1-eps)* (1-hn(distance, sigmaA))*(1-hn(distance, sigmaB)) )
  
  log_terms - log_any_detection
}

log_g_2observer_hr <- function(distance, detected, sigmaA, sigmaB, gammaA, gammaB, eps=1e-9){

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
  
  # this is the normalising constant to make the terms probabilities
  # i.e what is the prob of detected=i given the object was detected at all?
  
  # do the fraction (dist/sigma)^-gamma 
  fracA <- -(distance/sigmaA)^-gammaA
  fracB <- -(distance/sigmaB)^-gammaB
  
  log_any_detection <- log1p( -exp(fracA + fracB) )
                                
  log_terms #- log_any_detection
}

spline_detect_func <- function(spline_effect){ exp(-spline_effect) }
ln_spline_detect_func <- function(spline_effect){ -spline_effect }

####### the probability of detection at a given distance with given parameters

detect_func_2observer_hn <- function(distance, sigmaA, sigmaB){
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
  
  # do the fraction (dist/sigma)^-gamma itself on the log scale
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
    true_detect,
    true_loglambda,
    true_avg_prob_detect = mean(true_detect),
    ips_for_pred = fm_int(list(geometry = mexdolphin_sf$mesh), samplers=mexdolphin_sf$ppoly)
){
  base_mod <- models[[1]] # should always be the one observer hn
  
  # get the scores of this base model
  base_mod$DS_loglambda <- dawid_sebastiani_score(base_mod$pred_loglambda, true_loglambda)
  
  base_mod$loglambda_SE <- (base_mod$pred_loglambda$mean - true_loglambda)^2
  
  base_mod$lambda_AE <- abs(base_mod$pred_lambda$median - exp(true_loglambda) )
  
  base_mod$detect_AE <- abs(base_mod$pred_detect$median - true_detect)
  
  base_mod$DS_avg_prob <- dawid_sebastiani_score(base_mod$pred_avg_prob, true_avg_prob_detect)
  
  score_diffs <- NULL
  for (i in 2:length(models)){
    mod <- models[[i]]
    
    #get this model's scores 
    mod$DS_loglambda <- dawid_sebastiani_score(mod$pred_loglambda, true_loglambda)
    mod$loglambda_SE <- (mod$pred_loglambda$mean - true_loglambda)^2
    mod$lambda_AE <- abs(mod$pred_lambda$median - exp(true_loglambda) )
    mod$detect_AE <- abs(mod$pred_detect$median - true_detect)
    mod$DS_avg_prob <- dawid_sebastiani_score(mod$pred_avg_prob, true_avg_prob_detect)
    
    # now take the difference in scores between this and the base
    # and integrate it over the domain
    mod_diff <- list(
      DS_loglambda = sum(ips_for_pred$weight * (mod$DS_loglambda - base_mod$DS_loglambda)),
      loglambda_SE = sum(ips_for_pred$weight * (mod$loglambda_SE - base_mod$loglambda_SE)),
      lambda_AE = sum(ips_for_pred$weight * (mod$lambda_AE - base_mod$lambda_AE)),
      detect_AE = mean(mod$detect_AE - base_mod$detect_AE),
      DS_avg_prob = mod$DS_avg_prob - base_mod$DS_avg_prob
    )
    
    mod_diff$model <- mod$name
    
    score_diffs <- rbind(score_diffs, mod_diff)
  }
  
  res <- as.data.frame(score_diffs)
  
  #  wrangling with classes, for some reason the columns are lists
  j <- which(colnames(res)=="model")
  res[,j] <- as.character(res[,j])
  res[,-j] <- apply(res[,-j], 2, as.numeric)
  
  rownames(res) <- NULL

  res
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

