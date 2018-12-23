library(VariantAnnotation)
myvcf <- "/Users/hejiale/Desktop/For Jiale/170608Msi_PstI_genotypes.vcf"
mybg <- "/Users/hejiale/Desktop/For Jiale/170608Msi_PstI_genotypes.vcf.bgz"
library(polyRAD)
# indexTabix(mybg,"vcf")
hdr <- scanVcfHeader(myvcf)
mysvp <- ScanVcfParam(fixed = "ALT", info = "DP", geno = c("GT", "AD"))
mybg <- gsub ("$", ".bgz", myvcf) 
mytabix <- TabixFile(mybg, yieldSize = 1000)

# start here
open(mytabix)

# make empty list collections
p_collect11 <- list()
chrom_collect <- list()
pos_collect <- list()
MBM_collect <- list()
MAF_collect <- list()
DP_collect <- list()
p_collect13 <- list()
p_collect31 <- list()
AltRef_collect <- list()

count <- 1
# read again for another 1000 lines
while(nrow(vcf_chunk <- readVcf(mytabix, param = mysvp))){
  if (count == 11){
    break
    }else{
# show "."s
missingByMarker <- rowMeans(geno(vcf_chunk)$GT == ".")

# find heterozygotes
nummat <- matrix(NA_integer_, nrow = nrow(geno(vcf_chunk)$GT), ncol = ncol(geno(vcf_chunk)$GT), dimnames = dimnames(geno(vcf_chunk)$GT))
nummat[geno(vcf_chunk)$GT == "0/0"] <- 0L
nummat[geno(vcf_chunk)$GT == "0/1"] <- 1L
nummat[geno(vcf_chunk)$GT == "1/0"] <- 1L
nummat[geno(vcf_chunk)$GT == "1/1"] <- 2L
myMAF <- rowMeans (nummat, na.rm = TRUE)/2

# filter "."
SNPs_to_examine <- which(missingByMarker<= 0.5 & myMAF >= 0.05 & myMAF <=0.95)

# fill variables with NAs to avoid 0s
p_val11 <- rep(NA_real_, length(SNPs_to_examine))
p_val13 <- rep(NA_real_, length(SNPs_to_examine))
p_val31 <- rep(NA_real_, length(SNPs_to_examine))
chrom <- rep(NA_real_, length(SNPs_to_examine))
pos <- rep(NA_real_, length(SNPs_to_examine))
MBM <- rep(NA_real_, length(SNPs_to_examine))
MAF <- rep(NA_real_, length(SNPs_to_examine))
DP <- rep(NA_real_, length(SNPs_to_examine))
AltRef <- rep(NA_real_, length(SNPs_to_examine))

# small loop for each row
for( i in SNPs_to_examine){
  hets <- which(nummat[i,] == 1)
  if(length(hets) ==0){
    p_val11[match(i,SNPs_to_examine)] <- NA
    p_val13[match(i,SNPs_to_examine)] <- NA
    p_val31[match(i,SNPs_to_examine)] <- NA
  } else {
    het_AD <- geno(vcf_chunk)$AD[i,hets]
    ref_tot <- sum(sapply(het_AD, function(x)x[1]))
    alt_tot <- sum(sapply(het_AD, function(x)x[2]))
    p_val11[match(i,SNPs_to_examine)] <- chisq.test(matrix(c(ref_tot,alt_tot),nrow=2, ncol=1))$p.value
    p_val13[match(i,SNPs_to_examine)] <- chisq.test(matrix(c(ref_tot,alt_tot),nrow=2, ncol=1), p = c(0.25,0.75))$p.value
    p_val31[match(i,SNPs_to_examine)] <- chisq.test(matrix(c(ref_tot,alt_tot),nrow=2, ncol=1), p = c(0.75,0.25))$p.value
    AltRef[match(i,SNPs_to_examine)] <- alt_tot/(ref_tot + alt_tot)

# get all DPs
DP_all <- info(vcf_chunk)$DP
DP[match(i,SNPs_to_examine)] <- DP_all[i]

# get all chromosomes no.
chrom_all <- as.character(seqnames(rowRanges(vcf_chunk)))
if (is.na(chrom_all)){
  break}
chrom[match(i,SNPs_to_examine)] <- chrom_all[i]


# get all positions
pos_all <- start(ranges(rowRanges(vcf_chunk)))
pos[match(i,SNPs_to_examine)] <- pos_all[i]

# get all missing by markers
MBM[match(i,SNPs_to_examine)] <- missingByMarker[i]

# get all minor allele frequency
MAF[match(i,SNPs_to_examine)] <- myMAF[i]

  }}


# put all values in separate collections for each big loop
p_collect11[[count]] <- p_val11
p_collect13[[count]] <- p_val13
p_collect31[[count]] <- p_val31
chrom_collect[[count]] <- chrom
pos_collect[[count]] <- pos
MBM_collect[[count]] <- MBM
MAF_collect[[count]] <- MAF
DP_collect[[count]] <- DP
AltRef_collect[[count]] <- AltRef

count <- count+1
}
  }
# end
close(mytabix)

# unlist all
pc11_u <- unlist(p_collect11)
pc13_u <- unlist(p_collect13)
pc31_u <- unlist(p_collect31)
chrom_u <- unlist(chrom_collect)
pos_u <- unlist(pos_collect)
MBM_u <- unlist(MBM_collect)
MAF_u <- unlist(MAF_collect)
DP_u <- unlist(DP_collect)
AltRef_u <- unlist(AltRef_collect)
 # make .csv file 
write.csv(data.frame(chromosome = chrom_u, position = pos_u, MissingByMarker = MBM_u, 
                     MinorAlleleFrequency = MAF_u, Depth = DP_u, pvalue1to1_ref.over.alt = pc11_u,
                     pvalue1to3_ref.over.alt = pc13_u,
                     pvalue3to1_ref.over.alt = pc31_u, alt.over.total = AltRef_u),file= "SNO.1.csv")
a <- sapply(chrom_collect,length)
b <- sapply(pos_collect,length)
c <- sapply(MBM_collect,length)
d <- sapply(DP_collect,length)
e <- sapply(MAF_collect,length)
f <- sapply(p_collect11,length)

err <- which(is.na(match(a,c)))

### polyRAD processing
myvcf <- "/Users/hejiale/Desktop/For Jiale/170608Msi_PstI_genotypes.vcf"
clengths <- read.csv("/Users/hejiale/polyRAD/datasets/chromosome_lengths.CSV")
library(VariantAnnotation)
library(polyRAD)
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

alefq <- list()
realestimate <- list()
name1 <- list()
name2 <- list()
for (m in 1:length(myRAD)){
  # can pull allele frequency out of the iterate file: $alleleFreq
  real_iterate <- IterateHWE(myRAD[[m]])
  name1[[m]] <- real_iterate$locTable
  name2[[m]] <- real_iterate$alleleNucleotides
  ###alefq[[m]] <- real_iterate$alleleFreq
  # alleledepth/(alleledepth + antialleledepth) through colsum
  ###real_estigeno <- GetProbableGenotypes(real_iterate, omitCommonAllele = FALSE)
  ###realestimate[[m]] <- real_estigeno
}

alratio <- list()
tot_depth <- list()
for (n in 1:length(alefq)){
  #alratio[[n]] <- colSums(myRAD[[n]]$alleleDepth)/(colSums(myRAD[[n]]$alleleDepth)+colSums(myRAD[[n]]$antiAlleleDepth))
  tot_depth[[n]] <- colSums(myRAD[[n]]$alleleDepth)+colSums(myRAD[[n]]$antiAlleleDepth)
}
names12 <- data.frame(unlist(name1),unlist(name2))
Ratio_Freq <- data.frame(unlist(alratio),unlist(alefq))
colnames(Ratio_Freq) <- c("alleleRatio","alleleFrequency")
Ratio_Freq$difference <- abs(Ratio_Freq$alleleRatio - Ratio_Freq$alleleFrequency)
# compare allele frequency 
save(Ratio_Freq, file = "ratioXfreq.Rdata")

plot(Ratio_Freq$difference)

# look into loci with more than 0.4 difference
# collection of read depth (for a single read depth at one locus)

sigbias <- Ratio_Freq[which(Ratio_Freq$difference >= 0.4),]
RADloctable <- sapply(myRAD, function(x){x[3]})
uninames <- unique(gsub("_[ACGT]*$","",rownames(sigbias)))
SNO.1.1 <- read.csv("/Users/hejiale/polyRAD/datasets/SNO.1.1.CSV")
SNO.1.1$chromosome <- formatC(SNO.1.1$chromosome, width = 2, flag = "0")
SNO.1.1$uninames <- paste("S",SNO.1.1$chromosome,"_", SNO.1.1$position, sep = "")
loclookat <- na.omit(match(uninames,SNO.1.1$uninames))
lookat <- SNO.1.1[loclookat,]

# compare depth of the lookat with typical depth
# avg of depthRtio in myRAD, if it's close to tot1/(tot1+tot2), then biased, if close to freqency, no real biased.

diff_ratio <- list()
diff_freq <- list()
bias_calc <- list()
for (o in 1:length(myRAD)) {
  #diff_ratio[[o]] <- abs(colMeans(myRAD[[o]]$depthRatio, na.rm = TRUE)-alratio[[o]])
  #diff_freq[[o]] <- abs(colMeans(myRAD[[o]]$depthRatio, na.rm = TRUE)-alefq[[o]])
  bias_calc[[o]] <- colMeans(myRAD[[o]]$depthRatio, na.rm = TRUE)/alefq[[o]]
}

Ratio_Freq$diff_ratio <- unlist(diff_ratio)
Ratio_Freq$diff_freq <- unlist(diff_freq)

which(Ratio_Freq$diff_freq > 0.1)
which(Ratio_Freq$diff_ratio > 0.1)

SNO.1.1$X <- ifelse((SNO.1.1$X.1 %in% lookat$X.1), "lookat", "no")

ggplot(SNO.1.1)+
  geom_violin(aes(x = X, y = Depth))+
  coord_trans(y = "log")

plot(Ratio_Freq$diff_ratio, Ratio_Freq$diff_freq)
abline(lm(Ratio_Freq$diff_freq~Ratio_Freq$diff_ratio), col = "green")

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
plot(high_diffq$diff_freq[pos_in_chr[[1]]], log(as.integer(all_loc[pos_in_chr[[1]]])))
high_diffq$pos <- as.integer(all_loc)
high_diffq$chr <- gsub("S","", gsub("_.*","",rownames(high_diffq)))
ggplot(high_diffq)+
  geom_violin(aes(x = chr, y = log(pos)))

Ratio_Freq$bias <- unlist(bias_calc)
plot(log(Ratio_Freq$bias),Ratio_Freq$diff_freq)

Ratio_Freq$totaldepth <- unlist(tot_depth)
densitybias <- ggplot(Ratio_Freq, aes(x = log(bias), y = log(totaldepth)))+
  geom_point()+
  geom_density2d()+
  xlim(-2, 8)

# biased if with "bias" value smaller than 0.67 and bigger than 1.5
mean(Ratio_Freq$bias <= 0.67 | Ratio_Freq$bias >= 1.5)
# write a summary
# use sweep to see the allele depth
# sum read depth of each individual, sum of depth of allele/ sum of the depth of the individual, make it into new matrix
for(q in length(myRAD)){
  myRAD[[q]]$alleleDepth
}
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

AddNormalizedDepthProp.RADdata <- function(object, ...){
  # get the total number of reads per individual
  depthPerInd <- rowSums(object$locDepth)
  # get proportion of total reads for a taxon belonging to each allele
  propDepthAl <- sweep(object$alleleDepth, 1, depthPerInd, "/")
  totAl <- colSums(propDepthAl, na.rm = TRUE)
  # get proportion of total reads for a taxon belonging to other alleles at locus
  propDepthAnti <- sweep(object$antiAlleleDepth, 1, depthPerInd, "/")
  totAnti <- colSums(propDepthAnti, na.rm = TRUE)
  # get normalized proportion of reads for a locus belonging to an allele
  object$normalizedDepthProp <- totAl/(totAl + totAnti)
  
  return(object)
}
