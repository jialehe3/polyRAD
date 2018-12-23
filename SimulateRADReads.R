setwd("/Users/hejiale/polyRAD/datasets")
# read csv as string instead of factor
goldengate <- read.csv("goldengate_data.csv", row.names = 1, stringsAsFactors = FALSE)
# making a dictionary
conv <- c(0, 1, 1, 2)
names(conv) <- c("A/A", "A/B", "B/A", "B/B")
# transform dataframe to matrix, apply dictionary and get a vector. transform back to matrix and remodle matrix by nrow, ncol and dimnames.
golden_gate <- matrix(conv[as.matrix(goldengate)], nrow = nrow(goldengate), ncol = ncol(goldengate), dimnames = dimnames(goldengate))

# function to simiulate RAD-seq reads at one locus, given known genotypes ####
# genotypes here formatted as 0, 1, 2.  
# For the GoldenGate data, 0 = A/A, 1 = A/B, and 2 = B/B.
# The genotypes vector should be the genotypes for all individuals at one locus.
SimulateRADReads <- function(genotypes, prob = c(0.5, 0.5),
                             mndepth_shape = 2, mndepth_scale = 5, inddepth_scale = 10){
  # randomly select a mean read depth for this locus 
  meanDepth <- ceiling(rgamma(1, shape = mndepth_shape, scale = mndepth_scale))
  inddepth_shape <- meanDepth/inddepth_scale
  # randomly generate total locus depth for each taxon DP
  indDepth <- as.integer(round(rgamma(genotypes, shape = inddepth_shape, scale = inddepth_scale)))
  # reads to output
  outreads <- matrix(integer(length(genotypes)*2), nrow = length(genotypes), ncol = 2)
  # AD
  for(i in 1:length(genotypes)){
    if(is.na(genotypes[i])) next
    if(genotypes[i] == 0) outreads[i,1] <- indDepth[i]
    if(genotypes[i] == 2) outreads[i,2] <- indDepth[i]
    if(genotypes[i] == 1){
      thesereads <- sample(c(0,1), indDepth[i], replace = TRUE, prob = prob)
      outreads[i,1] <- sum(thesereads == 0)
      outreads[i,2] <- sum(thesereads == 1)
    }
  }
  
  return(outreads)
}

locus <- colnames(golden_gate)
golden_pcollection <- data.frame(row.names = locus)

# loop to get p_values from the sum of all type "1" individuals' AD
golden_p <- vector()  
golden_r <- vector()
simu2_collect <- list()
simu2 <- matrix(0L, nrow = nrow(golden_gate), ncol = ncol(golden_gate) *2,
                dimnames = list(rownames(golden_gate), 
                                paste(rep(colnames(golden_gate), each = 2), 0:1, sep = "_")))
# p_ratio is the bias
p_ratio <- list(c(0.9,0.1), c(0.8,0.2), c(0.7,0.3), c(0.6,0.4), c(0.5,0.5), c(0.4,0.6), c(0.3,0.7), c(0.2,0.8), c(0.1,0.9))
for (l in 1:length(p_ratio)){
for (j in 1:ncol(golden_gate)){
  if (length(which(golden_gate[,j] ==1)) == 0 ) next
  one <- as.numeric(which(golden_gate[,j] == 1))
  simu <- SimulateRADReads(golden_gate[one,j], prob = p_ratio[[l]])
  ####
  simu2[,c(2*j-1,2*j)] <- SimulateRADReads(golden_gate[,j],prob = p_ratio[[l]])
  ad1sum <- sum(simu[,1])
  ad2sum <- sum(simu[,2])
  if (ad1sum == 0 & ad2sum == 0) next
  golden_r[j] <- ad1sum/(ad1sum+ad2sum)
  golden_p[j] <- chisq.test(matrix(c(ad1sum,ad2sum)))$p.value
}
simu2_collect[[l]] <- simu2

# create new data from with locus and p_values only
golden_pcollection <- data.frame(golden_pcollection, data.frame(golden_p))
golden_pcollection <- data.frame(golden_pcollection, data.frame(golden_r))
# change column name
names(golden_pcollection)[names(golden_pcollection) == "golden_p"] <- paste("p_ratio", toString(p_ratio[[l]]))
names(golden_pcollection)[names(golden_pcollection) == "golden_r"] <- paste("alt.over.total", toString(p_ratio[[l]]))
}
# change prob = c() ratio
# do samething as did for previous VCF files

par(mfrow = c(3,3))
for (m in seq(from = 1, to = ncol(golden_pcollection), by = 2)){
  n <- golden_pcollection[, m]
  qqm <- qq(n, main = colnames(golden_pcollection)[m])
}

for (r in seq(from = 2, to = ncol(golden_pcollection), by = 2)){
  s <- golden_pcollection[, r]
  histr <- hist(s, main = colnames(golden_pcollection)[r])
}

####
golden_align <- read.csv("GoldenGate_Msi7_alignments.csv")


golden_locatable <- data.frame(row.names= colnames(golden_gate),
                               Chr = golden_align$Chr,
                               Pos = golden_align$Pos)

library(polyRAD)
golden_rmse <- vector()
for (s in 1:length(simu2_collect)){
diploid_sim <- RADdata(simu2_collect[[s]], rep(1:ncol(golden_gate), each = 2),
                      golden_locatable, list(2L), 0.001,
                      rep(c("A", "C"), times = ncol(golden_gate)))
golden_iterate <- IterateHWE(diploid_sim)
golden_estigeno <- GetProbableGenotypes(golden_iterate, omitCommonAllele = FALSE)
golden_rmse[s] <- sqrt(mean((golden_gate - golden_estigeno[["genotypes"]])^2, na.rm= TRUE))
}

SNO.1.1 <- read.csv("SNO.1.1.csv")
oneone <- which(SNO.1.1$conclusion == '1:1')
SNO.11 <- SNO.1.1[oneone,]

par(mfrow = c(1,2))
qq(SNO.11$pvalue1to1_ref.over.alt, main = "real_pval1to1")
qq(golden_pcollection$ratio..0.5..0.5, main = "goldengate_pval1to1")

par(mfrow = c(1,1))
qq(SNO.1.1$pvalue1to1_ref.over.alt)
