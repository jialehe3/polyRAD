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
write.csv(data.frame(chromosome = chrom_u, position = pos_u, MissingByMarker = MBM_u, MinorAlleleFrequency = MAF_u, Depth = DP_u, pvalue1to1_ref.over.alt = pc11_u, pvalue1to3_ref.over.alt = pc13_u, pvalue3to1_ref.over.alt = pc31_u, alt.over.total = AltRef_u), file= "SNO.1.csv")

a <- sapply(chrom_collect,length)
b <- sapply(pos_collect,length)
c <- sapply(MBM_collect,length)
d <- sapply(DP_collect,length)
e <- sapply(MAF_collect,length)
f <- sapply(p_collect11,length)

err <- which(is.na(match(a,c)))

### polyRAD processing
myvcf <- "/Users/hejiale/Desktop/For Jiale/170608Msi_PstI_genotypes.vcf"
clengths <- read.csv("/Users/hejiale/chromosome_lengths.CSV")
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
  times <- ceiling(clen[k]/1000000)
  for (l in 0:times){
    # extracting headers from myvcf
    samples <- VariantAnnotation::samples(VariantAnnotation::scanVcfHeader(myvcf))
    tryCatch({
    mysvp3 <- ScanVcfParam(fixed = "ALT", info = NA, geno = "AD", 
                           which = GRanges(j, 
                                           IRanges((l*1000000+1), (l+1)*1000000)),
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

realestimate <- list()
for (m in 1:length(myRAD)){
  # can pull allele frequency out of the iterate file: $alleleFreq
  real_iterate <- IterateHWE(myRAD[[m]])
  # alleledepth/alleledepth + antialleledepth through colsum
  real_estigeno <- GetProbableGenotypes(real_iterate, omitCommonAllele = FALSE)
  realestimate[[m]] <- real_estigeno
}
# compare allele frequency 
# git test 2
