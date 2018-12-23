library(VariantAnnotation)
library(polyRADBetaBinomial)
vcf_list <- list()
bg_list <- list()
setwd("/Users/hejiale/polyRAD")
file_list <- list.files("/Users/hejiale/Desktop/For Jiale/Filtered/")

for(i in 1:11) {
  vcf_list[[i]] <- paste("/Users/hejiale/Desktop/For Jiale/Filtered/", file_list[i], sep = "")
  bg_list[[i]] <- paste("/Users/hejiale/Desktop/For Jiale/Filtered/", file_list[i], ".bgz", sep = "")
}

mysvp <- ScanVcfParam(fixed = "ALT", info = "DP", geno = c("GT", "AD"))

myRAD_rice <- list()
for (j in 1 : length(vcf_list)){
  
    samples <- VariantAnnotation::samples(VariantAnnotation::scanVcfHeader(vcf_list[[j]]))
    tryCatch({
      rice_svp <- ScanVcfParam(fixed = "ALT", info = NA, geno = "AD",
                             samples = samples)
    myRAD_rice[[j]] <- VCF2RADdata(vcf_list[[j]], svparam = rice_svp, yieldSize = NA_integer_)
    },error = function(e){cat("ERROR :",conditionMessage(e), "\n")} )
}

alefq_rice <- list()
realestimate <- list()
name1 <- list()
name2 <- list()
better_bias <- list()

for (m in 1:length(myRAD_rice)){
  tryCatch({
    # can pull allele frequency out of the iterate file: $alleleFreq
    real_iterate <- IteratePopStruct(myRAD_rice[[m]])
    name1[[m]] <- real_iterate$locTable
    name2[[m]] <- real_iterate$alleleNucleotides
    alefq_rice[[m]] <- real_iterate$alleleFreq
    # alleledepth/(alleledepth + antialleledepth) through colsum
    real_estigeno <- GetProbableGenotypes(real_iterate, omitCommonAllele = FALSE)
    realestimate[[m]] <- real_estigeno
    better_bias[[m]] <- ((1-real_iterate$normalizedDepthProp)/(1-alefq_rice[[m]]))/(real_iterate$normalizedDepthProp/alefq_rice[[m]])
    
  },error = function(e){cat("ERROR :",conditionMessage(e), "\n")} )
}

alratio_rice <- list()
tot_depth <- list()
for (n in 1:length(alefq_rice)){
  alratio_rice[[n]] <- colSums(myRAD_rice[[n]]$alleleDepth)/(colSums(myRAD_rice[[n]]$alleleDepth)+colSums(myRAD_rice[[n]]$antiAlleleDepth))
  tot_depth[[n]] <- colSums(myRAD_rice[[n]]$alleleDepth)+colSums(myRAD_rice[[n]]$antiAlleleDepth)
}

Ratio_freq_rice <- data.frame(unlist(alratio_rice),unlist(alefq_rice))
colnames(Ratio_freq_rice) <- c("alleleRatio","alleleFrequency")
Ratio_freq_rice$difference <- abs(Ratio_freq_rice$alleleRatio - Ratio_freq_rice$alleleFrequency)


diff_ratio <- list()
diff_freq <- list()
bias_calc <- list()
for (o in 1:length(myRAD_rice)) {
  diff_ratio[[o]] <- abs(colMeans(myRAD_rice[[o]]$depthRatio, na.rm = TRUE)-alratio_rice[[o]])
  diff_freq[[o]] <- abs(colMeans(myRAD_rice[[o]]$depthRatio, na.rm = TRUE)-alefq_rice[[o]])
  bias_calc[[o]] <- colMeans(myRAD_rice[[o]]$depthRatio, na.rm = TRUE)/alefq_rice[[o]]
}

Ratio_freq_rice$diff_ratio <- unlist(diff_ratio)
Ratio_freq_rice$diff_freq <- unlist(diff_freq)
Ratio_freq_rice$bias <- unlist(bias_calc)
Ratio_freq_rice$totaldepth <- unlist(tot_depth)
Ratio_freq_rice$betterbias <- unlist(better_bias)


# AddNormalizedDepthProp.RADdata <- function(object, ...){
#   # get the total number of reads per individual
#   depthPerInd <- rowSums(object$locDepth)
#   # get proportion of total reads for a taxon belonging to each allele
#   propDepthAl <- sweep(object$alleleDepth, 1, depthPerInd, "/")
#   totAl <- colSums(propDepthAl, na.rm = TRUE)
#   # get proportion of total reads for a taxon belonging to other alleles at locus
#   propDepthAnti <- sweep(object$antiAlleleDepth, 1, depthPerInd, "/")
#   totAnti <- colSums(propDepthAnti, na.rm = TRUE)
#   # get normalized proportion of reads for a locus belonging to an allele
#   object$normalizedDepthProp <- totAl/(totAl + totAnti)
#   
#   return(object)
# }

# better_bias <- list()
# for (i in 1:length(myRAD_rice)){
#   myRAD_rice[[i]] <- AddNormalizedDepthProp.RADdata(myRAD_rice[[i]])
#   better_bias[[i]] <- ((1-myRAD_rice[[i]]$normalizedDepthProp)/(1-alefq_rice[[i]]))/(myRAD_rice[[i]]$normalizedDepthProp/alefq_rice[[i]])
# }

save(Ratio_freq_rice, file = "ratioXfreq_iteratepop_rice_biasincorp.Rdata")

mean(Ratio_freq_rice$betterbias <= 0.67 | Ratio_freq_rice$betterbias >= 1.5)


plot(Ratio_freq_rice$diff_ratio, Ratio_freq_rice$diff_freq)
abline(lm(Ratio_freq_rice$diff_freq~Ratio_freq_rice$diff_ratio), col = "green")

plot(log(Ratio_freq_rice$betterbias),Ratio_freq_rice$diff_freq)

ggplot(Ratio_freq_rice, aes(x = log(betterbias), y = log(totaldepth)))+
  geom_point()+
  geom_density2d()+
  xlim(-6, 6)

# remove.packages("polyRADbias")
# install.packages("polyRADbias", repos = NULL, type = "source")
# detach(package:polyRADbias)
