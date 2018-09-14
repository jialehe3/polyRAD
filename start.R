library(VariantAnnotation)
myvcf <- "/Users/hejiale/Desktop/For Jiale/170608Msi_PstI_genotypes.vcf"
mybg <- bgzip(myvcf)
indexTabix(mybg,"vcf")
hdr <- scanVcfHeader(myvcf)
hdr
mysam <- samples(hdr)
mysam [1:20]
geno(hdr)
mysvp <- ScanVcfParam(fixed = "ALT", info = NA, geno = c("GT", "AD"))

mybg <- gsub ("$", ".bgz", myvcf) 
mytabix <- TabixFile(mybg, yieldSize = 1000)

# start here
open(mytabix)

# read again for another 1000 lines
vcf_chunk <- readVcf(mytabix, param = mysvp)

rowRanges (vcf_chunk)
close(mytabix)

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
# show heterozygotes
hets <- which(nummat[SNPs_to_examine[1],] == 1)

het_AD <- geno(vcf_chunk)$AD[SNPs_to_examine[1],hets]
ref_tot <- sum(sapply(het_AD, function(x)x[1]))
alt_tot <- sum(sapply(het_AD, function(x)x[2]))
ref_tot
alt_tot

chisq.test(matrix(c(ref_tot,alt_tot),nrow=2, ncol=1))





while(nrow(vcf_chunk <- readVcf(mytabix,param = mysvp))){}

mytabix <- TabixFile(mybg, yieldSize = NA_integer_)
mysvp <- ScanVcfParam(fixed = "ALT", infor = NA, geno= c ("GT", "AD"), which = GRanges("06", IRanges(1000000,2000000)))
mygenome <- c(`06` = "Chr06")
vcf_chunk_06 <- readVcf(mytabix,genome = mygenome,param = mysvp)
vcf_chunk_06

library(polyRAD)
mysvp <- ScanVcfParam(fixed= "ALT", info = NA, geno = "AD", which = GRanges("06",IRanges(1000000,2000000)), samples = samples(scanVcfHeader(myvcf)))
myRAD <- VCF2RADdata(myvcf,svparam = mysvp, yieldSize= NA_integer_)
myRAD
myRAD$AlleleDepth[1:20,1:4]
myRAD$alleles2loc[1:30]

geno(vcf_chunk)$AD[SNPs_to_examine[1], 1:100]
geno(vcf_chunk)
ref_tot

install.packages("qqman")
install.packages("ggplot2")

