# Generalized function to simulate RAD-seq reads at one locus, given known 
# genotypes.  Can adapt to any ploidy.
# `genotypes` is a vector of numeric genotypes ranging from zero to the ploidy.
# A matrix with reads for each of two alleles is returned.
goldengate <- read.csv("goldengate_data.csv", row.names = 1, stringsAsFactors = FALSE)
conv <- c(0, 1, 1, 2)
names(conv) <- c("A/A", "A/B", "B/A", "B/B")
golden_gate <- matrix(conv[as.matrix(goldengate)], nrow = nrow(goldengate), ncol = ncol(goldengate), dimnames = dimnames(goldengate))
bias_group <- exp(seq(-2, 2, by = 0.25))


SimulateRADReads <- function(genotypes, ploidy = 2, mndepth_shape = 2, 
                             mndepth_scale = 5, inddepth_scale = 10,
                             overdispersion = 9, contamRate = 0.001,
                             bias = 1){
  require(rmutil)
  # randomly select a mean read depth for this locus
  meanDepth <- ceiling(rgamma(1, shape = mndepth_shape, scale = mndepth_scale))
  inddepth_shape <- meanDepth/inddepth_scale
  # randomly generate total locus depth for each taxon
  indDepth <- as.integer(round(rgamma(genotypes, shape = inddepth_shape, scale = inddepth_scale)))
  # reads to output
  outreads <- matrix(integer(length(genotypes)*2), nrow = length(genotypes), ncol = 2)
  # prob. of sampling second allele for each genotype
  alfreq <- mean(genotypes, na.rm = TRUE)/ploidy # frequency of second allele
  alprobs2 <- (0:ploidy)/ploidy * (1 - contamRate) + alfreq * contamRate
  alprobs2 <- alprobs2/(bias * (1 - alprobs2) + alprobs2)

  for(i in 1:length(genotypes)){
    if(is.na(genotypes[i])) next
    reads2 <- rbetabinom(1, indDepth[i], alprobs2[genotypes[i]+1],
                         overdispersion)
    outreads[i,2] <- as.integer(reads2)
    outreads[i,1] <- as.integer(indDepth[i] - reads2)
  }
  
  return(outreads)
}


golden_p <- vector()  
golden_r <- vector()
simu2_collect <- list()
simu2 <- matrix(0L, nrow = nrow(golden_gate), ncol = ncol(golden_gate) *2,
                dimnames = list(rownames(golden_gate), 
                                paste(rep(colnames(golden_gate), each = 2), 0:1, sep = "_")))
for (i in 1:length(bias_group)){
for (j in 1:ncol(golden_gate)){
  if (length(which(golden_gate[,j] ==1)) == 0 ) next
  one <- as.numeric(which(golden_gate[,j] == 1))
  simu <- SimulateRADReads(golden_gate[one,j])
  ####
  simu2[,c(2*j-1,2*j)] <- SimulateRADReads(golden_gate[,j], bias=bias_group[i])
  ad1sum <- sum(simu[,1])
  ad2sum <- sum(simu[,2])
  if (ad1sum == 0 & ad2sum == 0) next
  golden_r[j] <- ad1sum/(ad1sum+ad2sum)
  golden_p[j] <- chisq.test(matrix(c(ad1sum,ad2sum)))$p.value
}
  simu2_collect[[i]] <- simu2
}

golden_align <- read.csv("GoldenGate_Msi7_alignments.csv")
golden_locatable <- data.frame(row.names= colnames(golden_gate),
                               Chr = golden_align$Chr,
                               Pos = golden_align$Pos)
golden_rmse <- vector()
for (k in 1:length(simu2_collect)){
diploid_sim <- RADdata(simu2_collect[[k]], rep(1:ncol(golden_gate), each = 2),
                       golden_locatable, list(2L), 0.001,
                       rep(c("A", "C"), times = ncol(golden_gate)))
golden_iterate <- IterateHWE(diploid_sim)
golden_estigeno <- GetProbableGenotypes(golden_iterate, omitCommonAllele = FALSE)
golden_rmse[k] <- sqrt(mean((golden_gate - golden_estigeno[["genotypes"]])^2, na.rm= TRUE))
}
# changing bias values just like with prob
