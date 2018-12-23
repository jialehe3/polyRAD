RADprocessingPop <- function(object){
  alefq <- list()
  realestimate <- list()
  name1 <- list()
  name2 <- list()
  better_bias <- list()
  if(is.null(object$alleleDepth) == 0){
    object <- list(object)
  }
  for (m in 1:length(object)){
    
    tryCatch({
      real_iterate <- IteratePopStruct(object[[m]])
      alefq[[m]] <- real_iterate$alleleFreq
      name1[[m]] <- real_iterate$locTable
      name2[[m]] <- real_iterate$alleleNucleotides
      better_bias[[m]] <- ((1-real_iterate$normalizedDepthProp)/(1-alefq[[m]]))/(real_iterate$normalizedDepthProp/alefq[[m]])
      
      real_estigeno <- GetProbableGenotypes(real_iterate, omitCommonAllele = FALSE)
      realestimate[[m]] <- real_estigeno
    },error = function(e){cat("ERROR :",conditionMessage(e), "\n")} ) 
  }
  
  alratio <- list()
  tot_depth <- list()
  for (n in 1:length(alefq)){
    alratio[[n]] <- colSums(object[[n]]$alleleDepth)/(colSums(object[[n]]$alleleDepth)+colSums(object[[n]]$antiAlleleDepth))
    tot_depth[[n]] <- colSums(object[[n]]$alleleDepth)+colSums(object[[n]]$antiAlleleDepth)
  }
  
  if(length(alratio) != length(tot_depth)){
    stop("NULL probably exists, dataframe cannot be built")
  }
  
  Ratio_Freq_pop <- data.frame(unlist(alratio),unlist(alefq))
  colnames(Ratio_Freq_pop) <- c("alleleRatio","alleleFrequency")
  
  Ratio_Freq_pop$difference_alefq_alratio <- abs(Ratio_Freq_pop$alleleRatio - Ratio_Freq_pop$alleleFrequency)
  
  
  diff_ratio <- list()
  diff_freq <- list()
  bias_calc <- list()
  for (o in 1:length(object)) {
    diff_ratio[[o]] <- abs(colMeans(object[[o]]$depthRatio, na.rm = TRUE)-alratio[[o]])
    diff_freq[[o]] <- abs(colMeans(object[[o]]$depthRatio, na.rm = TRUE)-alefq[[o]])
    bias_calc[[o]] <- colMeans(object[[o]]$depthRatio, na.rm = TRUE)/alefq[[o]]
  }
  
  Ratio_Freq_pop$diff_ratio <- unlist(diff_ratio)
  Ratio_Freq_pop$diff_freq <- unlist(diff_freq)
  Ratio_Freq_pop$bias <- unlist(bias_calc)
  Ratio_Freq_pop$totaldepth <- unlist(tot_depth)
  Ratio_Freq_pop$betterBias <- unlist(better_bias)
  
  plot(log(object$betterbias),object$diff_freq)
}

