###########  detection functions

hn <- function(distance, sigma){  exp(-0.5 * (distance/sigma)^2 )  }
log_hn <-  function(distance, sigma){ -0.5 * (distance/sigma)^2 }

hr <- function(distance, sigma, gamma) {
  # catt("hr arg lengths:", length(distance), length(sigma), length(gamma))
  
  if (any(is.na(sigma))){
    catt("hr : got ", sum(is.na(sigma)), " NAs")
  }
  
  1 - exp( -(distance / sigma)^(-gamma) )
}
log_hr <- function(distance, sigma, gamma){
  # do (dist/sigma)^-b on log scale
  # frac_part <- -exp( -gamma*log(distance) + gamma*log(sigma) )
  # 
  # log1p( -exp(frac_part) )
  log(hr(distance, sigma, gamma))
}

log_1mhr <-  function(distance, sigma, gamma, eps=1e-25){
  # ln( 1 - (1- exp(frac) ) = ln(exp(frac)) = frac
  frac <- -(distance/sigma)^(-gamma)
  
  out <- pmax(-5519956, frac)
  
  ii <- !is.finite(frac)
  catt("log1mhr range of output is", range(frac), "clipped to", range(out) )
  if (any(ii)){
    catt("nonfinite results from dists", distance[ii])
  }
  else{
    catt("all outputs finite")
  }
  # incase the true probability of non detection is ~= 0 
  # we output log(epsilon) instead of log(0) for small epsilon > 0
  out 
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

log_g_2observer_hr <- function(distance, detected, sigmaA, sigmaB, gammaA, gammaB, eps=1e-6){

  # catt("range of distances given to detect func", range(distance))
  
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
  log_terms[ii] =  log_hr(distance[ii], sigmaA[ii], gammaA[ii]) + log_hr(distance[ii], sigmaB[ii], gammaB[ii])

  
  # # # A alone
  # ii <- detected==1
  # log_terms[ii] = log_hr(distance[ii], sigmaA[ii], gammaA[ii]) + log_1mhr(distance[ii], sigmaB[ii], gammaB[ii], eps)
  # # B alone
  # ii <- detected==2
  # log_terms[ii] = log_1mhr(distance[ii], sigmaA[ii], gammaA[ii], eps) + log_hr(distance[ii], sigmaB[ii], gammaB[ii])
  # # Both A,B
  # ii <- detected==3
  # log_terms[ii] =  log_hr(distance[ii], sigmaA[ii], gammaA[ii]) + hr(distance[ii], sigmaB[ii], gammaB[ii])

  
  log_terms
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

######### general utilities
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

pretty_print_seconds <- function(secs){
  secs <- round(secs)
  
  hrs <-  secs %/% (60*60)
  secs <- secs %% (60*60)
  
  mins <- secs %/% 60
  secs <- secs %% 60
  
  txt <- ""
  
  if (hrs > 0){
    txt <- paste0(txt, hrs, "hrs", " ")
  }
  if(mins >0){
    txt <- paste0(txt, mins, "mins", " ")
  }
  if(secs >0){
    txt <- paste0(txt, secs, "secs", " ")
  }
  
  txt
}


############ examine LGCP simulation realisation
lets_have_a_look_at_you <- function(
    sim,
    detect_df,
    approx_intensity_range = c(0, 0.225) # so all plots have the same colour scale
){
  
  all_pts <- sim$unthinned_samples_df
  obs_pts <- sim$samples_df
  
  col_states <- c("red", "blue", "purple")
  name_states <- c("Observer A only", "Observer B only", "Observer A and B")

  ##########  plot all points and just observed points
  animal_plot <- ggplot() +
    geom_sf(data = all_pts$geometry, col = "black") +
    labs(title = paste("True abundance =", nrow(all_pts)) ) 
  
  obs_animal_plot <- ggplot() + 
    geom_sf(
      data = obs_pts,
      mapping = aes(colour=as.character(obs_pts$detected)) 
    ) +
    scale_color_manual(
      name = "Type of detection",
      breaks = c("1", "2", "3"),
      values = col_states,
      labels = name_states
    ) +
    labs(title = paste("Total observed = ", nrow(obs_pts)))
  

  ########## plot lambda
  pxls <- fm_pixels(sim_info$the_mesh, mask = sim_info$boundary, dims = c(200, 200))
  
  pxls$loglambda <- fm_evaluate(sim_info$the_mesh, loc= pxls, field=sim_info$log_lambda )
  pxls$lambda <- exp(pxls$loglambda)
  
  intensity_plot <- ggplot() +
    gg(sim_info$boundary, alpha= .1) +
    geom_tile(
      data = pxls,
      aes(geometry = geometry, fill = lambda),
      stat = "sf_coordinates"  
    ) +
    geom_sf(data = sim_info$boundary, alpha = 0.1) +
    labs(title = "True underlying animal density or intensity") +
    scale_fill_continuous(
      name = "Density",
      limits = range(approx_intensity_range, pxls$lambda)
      ) + 
    theme( 
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    ) 
  
  print(intensity_plot/ animal_plot / obs_animal_plot)
  
  ##########  histogram of observed distances, per observer and all animals' distances
  
  
  # the detection functions arent pdfs so theyre only proportional to the true
  # prob of detecting each distance - scaling by a constant for line transects
  rescale_to_overlay <- function(obs_dists, prob){
    
    # the highest point in the histogram,
    hist_max <- max(hist(obs_dists,30, plot=F)$density)
    
    # will now coincide with the highest point in the pdf curve
    factor <- hist_max/max(prob)
    
    factor*prob
  }
  
  detect_df$Amarginal <- rescale_to_overlay(
    obs_pts$distance[obs_pts$detectA],
    detect_df$Amarginal
  )
  detect_df$Bmarginal <- rescale_to_overlay(
    obs_pts$distance[obs_pts$detectB],
    detect_df$Bmarginal
  )
  detect_df$detected1 <- rescale_to_overlay(
    obs_pts$distance[obs_pts$detected==1],
    detect_df$detected1
  )
  detect_df$detected2 <- rescale_to_overlay(
    obs_pts$distance[obs_pts$detected==2],
    detect_df$detected2
  )
  detect_df$detected3 <- rescale_to_overlay(
    obs_pts$distance[obs_pts$detected==3],
    detect_df$detected3
  )
  detect_df$any <- rescale_to_overlay(
    obs_pts$distance,
    detect_df$any
  )
  
  # hist of all animals distances
  
  # make sure all detection plots are in the same frame
  g <- ggplot() + expand_limits(y= c(0, 1) )
  
  all_animals_d <- g + 
    geom_histogram(
      data = all_pts, 
      mapping=aes(distance, y = after_stat(density) ),
      bins=30
    ) 
  
  obs_animals_d <- g + 
    geom_histogram(
      data = obs_pts, 
      mapping=aes(distance, y = after_stat(density) ),
      bins=30
    ) +
    geom_line(
      data = detect_df,
      aes(distance, any),
      linewidth = 2,
      col = "green" # idk what colour is the union of blue and red
    ) +
    labs(
      title = "Observed distances by either observer",
      subtitle = "Theoretical probability of detection overlain"
    )
  
  # hist of those in each detection state and theoretical prob overlaid
  detect_dist_plots <- sapply(
    1:3,
    function(detect_state){
      g + 
        geom_histogram(
          data = obs_pts[obs_pts$detected == detect_state,],
          mapping=aes(distance, y = after_stat(density) ),
          bins = 30
        )+
        geom_line(
          data = detect_df,
          # see https://ggplot2.tidyverse.org/reference/aes_.html
          aes(distance, .data[[paste0("detected", detect_state)]]),
          linewidth = 2,
          col = col_states[detect_state]
        ) +
        labs(
          title = paste("Observed distances of", name_states[detect_state]),
          subtitle = "Theoretical probability of detection overlain"
        )
    }
  )
  
  
  detect_states_plot <- (obs_animals_d  + detect_dist_plots[[1]])/ (detect_dist_plots[[2]] + detect_dist_plots[[3]]) 
  
  
  # hist of distances detected marginally by each observer
  
  g <- ggplot(mapping=aes(distance, y = after_stat(density) )) + 
    expand_limits(y= c(0, 1) )
  
  detect_A_plot <- g + 
    geom_histogram(
      data = obs_pts[obs_pts$detectA,],
      bins = 30,
      closed = "left"
    ) +
    geom_line(
      data = detect_df,
      aes(distance, Amarginal),
      linewidth = 2,
      col = "red"
    ) + 
    labs(title= "Observed distances by Observer A marginally")
  
  detect_B_plot <- g + 
    geom_histogram(
      data = obs_pts[obs_pts$detectB,],
      bins = 30,
      closed = "left"
    ) +
    geom_line(
      data = detect_df,
      aes(distance, Bmarginal),
      linewidth = 2,
      col = "blue"
    ) + 
    labs(title= "Observed distances by Observer B marginally")
  
  print(  detect_states_plot | (detect_A_plot/detect_B_plot) )
}



