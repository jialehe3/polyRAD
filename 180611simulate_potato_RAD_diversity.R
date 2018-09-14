setwd("/Users/hejiale")
# import data from SolCAP - http://solcap.msu.edu/potato_infinium.shtml
# Potato diversity panel dosage genotype scores

# Import genotypes ####
# read in as character strings
potato_char <- as.matrix(read.csv("potato_5cluster_genotype_scores_3763m_250g.csv",
                        row.names = 1, header = TRUE, stringsAsFactors = FALSE))

# set up matrix to convert to integer
potato_num <- matrix(NA_integer_, nrow = nrow(potato_char), ncol = ncol(potato_char),
                     dimnames = dimnames(potato_char))
# loop through markers to convert to integer
for(i in 1:nrow(potato_num)){
  # find the two alleles, and which one is first
  for(j in 1:ncol(potato_char)){
    if(potato_char[i,j] == ".") next
    nuc1 <- substring(potato_char[i,j], 1, 1)
    nuc2 <- substring(potato_char[i,j], 4, 4)
    if(nuc1 != nuc2) break
  }
  # make sure they are what is expected
  if(!all(c(nuc1, nuc2) %in% c('A', 'C', 'G', 'T'))){
    stop("Unexpected nucleotide")
  }
  
  # set up expected genotypes
  exp_gen <- character(5)
  for(g in 1:5){
    exp_gen[g] <- paste(rep(c(nuc1, nuc2), times = c(5 - g, g - 1)),
                        collapse = "")
  }
  # conversion vector
  conv <- c(0:4, NA_integer_)
  names(conv) <- c(exp_gen, ".")
  
  # convert genotypes and add to matrix
  potato_num[i,] <- conv[potato_char[i,]]
}

# allele frequencies
hist(rowMeans(potato_num, na.rm = TRUE)/4) # very even dist. (probably due to selection of informative SNPs)

# function to simiulate RAD-seq reads at one locus, given known genotypes ####
# genotypes here formatted as 0, 1, 2, 3, 4
SimulateRADReadsTetra <- function(genotypes, ratio = c(0.5, 0.5), mndepth_shape = 2, mndepth_scale = 5, inddepth_scale = 10){
  # randomly select a mean read depth for this locus
  meanDepth <- ceiling(rgamma(1, shape = mndepth_shape, scale = mndepth_scale))
  inddepth_shape <- meanDepth/inddepth_scale
  # randomly generate total locus depth for each taxon
  indDepth <- as.integer(round(rgamma(genotypes, shape = inddepth_shape, scale = inddepth_scale)))
  # reads to output
  outreads <- matrix(integer(length(genotypes)*2), nrow = length(genotypes), ncol = 2)
  
  for(i in 1:length(genotypes)){
    if(is.na(genotypes[i])) next
    if(genotypes[i] == 0) outreads[i,1] <- indDepth[i]
    if(genotypes[i] == 4) outreads[i,2] <- indDepth[i]
    if(genotypes[i] %in% c(1,2,3)){
      alprob <- genotypes[i]/4
      thesereads <- sample(c(0,1), indDepth[i], replace = TRUE, prob = c(1-alprob, alprob)* ratio) 
      outreads[i,1] <- sum(thesereads == 0)
      outreads[i,2] <- sum(thesereads == 1)
    }
  }
  
  return(outreads)
}

# simulate RAD-seq read depth for all 3763 markers ####
potato_sim_reads <- matrix(0L, nrow = ncol(potato_num), ncol = nrow(potato_num) * 2,
                           dimnames = list(colnames(potato_num), 
                                           paste(rep(rownames(potato_num), each = 2), 0:1, sep = "_")))

tetra_ratio <- list(c(0.9,0.1), c(0.8,0.2), c(0.7,0.3), c(0.6,0.4), c(0.5,0.5), c(0.4,0.6), c(0.3,0.7), c(0.2,0.8), c(0.1,0.9))
tetra_ratio2 <- list(c(0.55,0.45), c(0.54,0.46), c(0.53,0.47), c(0.52,0.48), c(0.51, 0.48), c(0.5,0.5), c(0.49,0.51), c(0.48,0.52), c(0.47,0.53), c(0.46,0.54), c(0.45,0.55))
potato_sim_collect <- list()
for (m in 1:length(tetra_ratio2)){
for(i in 1:nrow(potato_num)){
  potato_sim_reads[,c(2*i - 1, 2*i)] <- SimulateRADReadsTetra(potato_num[i,], ratio = tetra_ratio2[[m]])
}
  potato_sim_collect[[m]] <- potato_sim_reads
}

potato_sim_reads[1:10,1:4]
t(potato_num[1:2,1:10])

# get sum of depth at all ratios
sum_depth <- vector()
sum_depth_collect <- list()
for (n in 1:length(tetra_ratio2)){
for (k in 1:ncol(potato_sim_reads)){
 sum_depth[k] <- sum(potato_sim_collect[[n]][,k]) 
}
  sum_depth_collect[[n]] <- sum_depth
}

# get all ad1/(ad1+ad2) and p_values under all ratios
allele_adratio <- vector()
allele_adratio_collect <- list()
allele_pval <- vector()
allele_pval_collect <- list()
for (o in 1:length(sum_depth_collect)){
for (l in 1:(length(sum_depth)/2)){
  sum1 <- sum_depth_collect[[o]][l*2 - 1]
  sum2 <- sum_depth_collect[[o]][l*2]
  allele_adratio[l] <- sum1/(sum1+sum2)
  allele_pval[l] <- chisq.test(matrix(c(sum1,sum2)))$p.value
}
  allele_adratio_collect[[o]] <- allele_adratio 
  allele_pval_collect[[o]] <- allele_pval
}

# get allele frequencies
allele_freq <- vector()
for (j in 1:nrow(potato_num)){
  allele_freq[j] <- sum(potato_num[j,], na.rm= T)/(sum(!is.na(potato_num[j,])) * 4)
}

# combine allele frequency and allele ad ratio 
allele_tetra <- list()
for (r in 1:length(allele_adratio_collect)){
allele_tetra[[r]] <- cbind(allele_freq,allele_adratio_collect[[r]])
}
# filter out duplicates and bind the rest allele frequency and ad ratio
allele_tetra_final <- data.frame(allele_tetra[[1]],row.names = row.names(potato_num))
for (s in 2:length(allele_tetra)){
allele_tetra_final <- cbind(allele_tetra_final, allele_tetra[[s]][,2])
}

# change col names
for (t in 2:ncol(allele_tetra_final)){
  colnames(allele_tetra_final)[t] <- tetra_ratio2[t-1]
}

# plot all 9 scatter plots on one page
par(mfrow = c(3,3))
for (u in 2:ncol(allele_tetra_final)){
  scattertot <- plot(allele_tetra_final$allele_freq, allele_tetra_final[,u], main = tetra_ratio2[u-1], xlab= "frequency", ylab= "ADratio")
}

qq(allele_pval_collect)

# import alignment data ####
potato_align <- read.csv("solcap_snp_tagdigger.csv", stringsAsFactors = FALSE)
str(potato_align)HWE
align_match <- match(rownames(potato_num),
                     potato_align$Old.marker.name)
sum(is.na(align_match))

potato_loctable <- data.frame(row.names = rownames(potato_num),
                              Chr = potato_align$Chr[align_match],
                              Pos = potato_align$SNP_Pos[align_match])

# convert to RADdata
library(polyRAD)
potato_rmse <- vector()
for (v in 1:length(potato_sim_collect)){
potato_sim <- RADdata(potato_sim_collect[[v]], rep(1:nrow(potato_num), each = 2),
                      potato_loctable, list(4L), 0.001,
                      rep(c("A", "C"), times = nrow(potato_num)))

#save(potato_sim, potato_num, file = "180611potato_sim_diversity.RData")

# get estimated genos from simulated read depth
# population frequency (AlleleFreq)  inside potato_iterate
# tot1(depth)/tot1+tot2(depth) = freq1, depth in real data are directly from pipeline, frequency are generate polyRAD by using genotypes generated from ployRAD
potato_iterate <- IterateHWE(potato_sim)
# omitCommonAllele to filter out common alleles and used to look for rare alleles, not helpful while calculating rmse.
potato_estigeno <- GetProbableGenotypes(potato_iterate, omitCommonAllele = FALSE)

# calculate rmse between estimated genos and real genos
potato_rmse[v] <- sqrt(mean((t(potato_num) - potato_estigeno[["genotypes"]])^2, na.rm= TRUE))
}

#123


