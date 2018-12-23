library(VariantAnnotation)
library(polyRADBetaBinomial)
myvcf <- "/Users/hejiale/Desktop/For Jiale/170608Msi_PstI_genotypes.vcf"
clengths <- read.csv("/Users/hejiale/polyRAD/datasets/chromosome_lengths.CSV")

clengths$chr <- gsub(clengths$chr, patter = "Chr", replacement = "")
clengths$chr <- gsub(clengths$chr, patter="scaffold",replacement = "SCAFFOLD")
cchrom <- clengths$chr
clen <- clengths$length


myRAD <- list()
count <- 1
for (j in cchrom[1:19]){
  
  k <- match(j, clengths$chr)
  times <- ceiling(clen[k]/10000000)
  for (l in 0:times){
    # extracting headers from myvcf
    samples <- VariantAnnotation::samples(VariantAnnotation::scanVcfHeader(myvcf))
    tryCatch({
      mysvp3 <- ScanVcfParam(fixed = "ALT", info = NA, geno = "AD", 
                             which = GRanges(j, 
                                             IRanges((l*10000000+1), (l+1)*10000000)),
                             samples = samples)
      # get chuncks of svps using the which argument within Scanvcfparam
      # then goto VCF2RADdata
      # then use polyRAD to convert each chunk
      # chombine results in matrix particles
      
      myRAD[[count]] <- VCF2RADdata(myvcf, svparam = mysvp3, yieldSize = NA_integer_)
      count <- count + 1
    },error = function(e){cat("ERROR :",conditionMessage(e), "\n")} )
  }
}







RADprocessing <- function(object){
alefq <- list()
realestimate <- list()
name1 <- list()
name2 <- list()
better_bias <- list()

for (m in 1:length(myRAD)){
  tryCatch({
  # can pull allele frequency out of the iterate file: $alleleFreq
  real_iterate <- IteratePopStruct(myRAD[[m]])
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
  alratio[[n]] <- colSums(myRAD[[n]]$alleleDepth)/(colSums(myRAD[[n]]$alleleDepth)+colSums(myRAD[[n]]$antiAlleleDepth))
  tot_depth[[n]] <- colSums(myRAD[[n]]$alleleDepth)+colSums(myRAD[[n]]$antiAlleleDepth)
}



Ratio_Freq_pop <- data.frame(unlist(alratio),unlist(alefq))
colnames(Ratio_Freq_pop) <- c("alleleRatio","alleleFrequency")

Ratio_Freq_pop$difference_alefq_alratio <- abs(Ratio_Freq_pop$alleleRatio - Ratio_Freq_pop$alleleFrequency)
# compare allele frequency 


diff_ratio <- list()
diff_freq <- list()
bias_calc <- list()
for (o in 1:length(myRAD)) {
  diff_ratio[[o]] <- abs(colMeans(myRAD[[o]]$depthRatio, na.rm = TRUE)-alratio[[o]])
  diff_freq[[o]] <- abs(colMeans(myRAD[[o]]$depthRatio, na.rm = TRUE)-alefq[[o]])
  bias_calc[[o]] <- colMeans(myRAD[[o]]$depthRatio, na.rm = TRUE)/alefq[[o]]
}

# diff_ratio_filtered <- diff_ratio[-88]
Ratio_Freq_pop$diff_ratio <- unlist(diff_ratio)
Ratio_Freq_pop$diff_freq <- unlist(diff_freq)
Ratio_Freq_pop$bias <- unlist(bias_calc)
Ratio_Freq_pop$totaldepth <- unlist(tot_depth)
Ratio_Freq_pop$betterBias <- unlist(better_bias)
}
save(Ratio_Freq_pop, file = "ratioXfreq_pop_biasincorp.Rdata")



# which(Ratio_Freq_pop$diff_freq > 0.1)
# which(Ratio_Freq_pop$diff_ratio > 0.1)
# SNO.1.1$X <- ifelse((SNO.1.1$X.1 %in% lookat$X.1), "lookat", "no")
# ggplot(SNO.1.1)+
#   geom_violin(aes(x = X, y = Depth))+
#   coord_trans(y = "log")


# all values with higher diff_freq values, how often(portion out of the whole data), threshold = 0.1
# separate the targeted value with chromosome and observe positions on chr maybe?
# bias = (colmean(depthratio))/alefreq
length(which(Ratio_Freq_pop$diff_freq >= 0.1))/nrow(Ratio_Freq_pop)
high_diffq <- Ratio_Freq_pop[which(Ratio_Freq_pop$diff_freq >= 0.1),]

pos_in_chr <- list()
for (p in unique(substr(rownames(high_diffq),1, 3))){
  pos_in_chr[[which(unique(substr(rownames(high_diffq),1, 3)) == p)]] <-grep(p, rownames(high_diffq))
}

all_loc <- gsub("S.*._","", gsub("_[ACGT]*$","",rownames(high_diffq)))

high_diffq$pos <- as.integer(all_loc)
high_diffq$chr <- gsub("S","", gsub("_.*","",rownames(high_diffq)))

ggplot(high_diffq)+
  geom_violin(aes(x = chr, y = log(pos)))


# biased if with "bias" value smaller than 0.67 and bigger than 1.5
mean(Ratio_Freq_pop$bias <= 0.67 | Ratio_Freq_pop$bias >= 1.5)
# write a summary
# use sweep to see the allele depth
# sum read depth of each individual, sum of depth of allele/ sum of the depth of the individual, make it into new matrix

rdpmatrix <- matrix(nrow = length(myRAD[[1]]))
# then use colsum_of_the_new_matrix/(allele_depth + anti_allele_depth)

# yd_func <- function(x) {
#   r <- nrow(x$alleleDepth)
#   c <- ncol(x$alleleDepth)
#   allele_matrix <- matrix(NA, r, c)
#   rsum <- rowsum(x$alleleDepth)
#   csum <- colsum(x$alleleDepth)
#   for (q in 1:c) {
#     allele_matrix[, q] <- csum[q]/rsum
#   }
# }
# 
# yd <- sapply(myRAD, function(x) {
#   r <- nrow(x$alleleDepth)
#   c <- ncol(x$alleleDepth)
#   allele_matrix <- matrix(NA, r, c)
#   rsum <- rowSums(x$alleleDepth)
#   csum <- colSums(x$alleleDepth)
#   for (q in 1:c) {
#     allele_matrix[, q] <- csum[q]/rsum
#   }
#   return(allele_matrix)
# })

# newmatrix <- do.call(cbind, yd)

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
# 
# better_bias <- list()
# for (i in 1:length(myRAD)){
#   myRAD[[i]] <- AddNormalizedDepthProp.RADdata(myRAD[[i]])
#   better_bias[[i]] <- ((1-myRAD[[i]]$normalizedDepthProp)/(1-alefq[[i]]))/(myRAD[[i]]$normalizedDepthProp/alefq[[i]])
# }


plot(Ratio_Freq_pop$diff_ratio, Ratio_Freq_pop$diff_freq)
abline(lm(Ratio_Freq_pop$diff_freq~Ratio_Freq_pop$diff_ratio), col = "green")

plot(log(Ratio_Freq_pop$betterBias),Ratio_Freq_pop$diff_freq)

ggplot(Ratio_Freq_pop, aes(x = log(betterBias), y = log(totaldepth)))+
  geom_point()+
  geom_density2d()+
  xlim(-6, 6)

message = FALSE
warning = FALSE
