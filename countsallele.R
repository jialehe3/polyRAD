
allelecount <- list()
for (i in 1:length(myRAD)){
  readdepth <- myRAD[[i]]$alleleDepth
  alleleloc <- myRAD[[i]]$alleles2loc
  # readdepth[,alleleloc == x] >0 will make a matrix of true of false, in this case it 585*couple thousands filled with True/False, rowSums count on the number of Trues on each row.
  # done by each loci(group of alleles) then each individual.
  # allele2loc is just an index, so every chunck will start from 1
  allelecount[[i]] <- sapply(1:nLoci(myRAD[[i]]),function(x)rowSums(readdepth[,alleleloc == x] >0))
  colnames(allelecount[[i]]) <- GetLoci(myRAD[[i]])
}
# ?"["
# View(allelecount[[1]][343,,drop=FALSE])

mt2individual <- list()
mt2loci <- list()
mt2individualsum <- numeric(nrow(readdepth))

for  (j in 1:length(allelecount)){
  mt2individual[[j]] <- rowSums(allelecount[[j]] > 2 )
  mt2loci[[j]] <- colSums(allelecount[[j]] > 2)
  mt2individualsum <- mt2individualsum + mt2individual[[j]]
}
mt2locisum <- unlist(mt2loci)

individualout <- boxplot(mt2individualsum)$out
lociout <- boxplot(mt2locisum)$out
lociout1 <- boxplot(boxplot(mt2locisum)$out)$out
lociout2 <- boxplot(boxplot(boxplot(mt2locisum)$out)$out)$out

locitemp <- list()
locicount <- list()
for (k in 1:length(myRAD)){
  readdepth <- myRAD[[k]]$alleleDepth
  alleleloc <- myRAD[[k]]$alleles2loc
  locitemp[[k]] <- sapply(1:nLoci(myRAD[[k]]),function(x)rowSums(readdepth[,alleleloc == x]))
  locicount[[k]] <- data.frame(colSums(locitemp[[k]]))
  rownames(locicount[[k]]) <- GetLoci(myRAD[[k]])
}
locicountsum <- do.call(rbind, unname(locicount))
locicountsum <- data.matrix(locicountsum, rownames.force = NA)
names(locicountsum) <- rownames(locicountsum)

threshold50000 <- locicountsum[which(locicountsum > 50000)]
threshold1000mt2 <- mt2locisum[which(mt2locisum > locicountsum/1000)]

thresholdunique <- unique(c(names(threshold1000mt2),names(threshold50000)))
locicountWthresh <- locicountsum[which(is.na(match(names(locicountsum),thresholdunique)) == TRUE)]
mt2locisumWthresh <- mt2locisum[which(is.na(match(names(locicountsum),thresholdunique)) == TRUE)]
lociCWT50000 <- locicountsum[which(is.na(match(names(locicountsum),names(threshold50000))) == TRUE)]
mt2lociCWT50000 <- mt2locisum[which(is.na(match(names(locicountsum),names(threshold50000))) == TRUE)]
lociCWT1k <- locicountsum[which(is.na(match(names(locicountsum),names(threshold1000mt2))) == TRUE)]
mt2lociCWT1k <- mt2locisum[which(is.na(match(names(locicountsum),names(threshold1000mt2))) == TRUE)]


# no threshold used
hist(log(mt2locisum/locicountsum),breaks = 100,main = "origin")
# 50k locisum threshold only
hist(log(mt2lociCWT50000/lociCWT50000),breaks = 100, main = "threshold_50000")
# 1k mt2locisum threshold only
hist(log(mt2lociCWT1k/lociCWT1k),breaks = 100, main = "threshold_1k_mt2")
# both 50k locisum and 1k mt2locisum thresholds apply
hist(log(mt2locisumWthresh/locicountWthresh),breaks = 100, main = "both threshold")


posterprobtest <- IterateHWE(myRAD[[1]])
View(posterprobtest$posteriorProb[[1]][1,,][1:10])

posterprob <- list()
for (m in 1:length(myRAD)){
  tryCatch({
    # can pull allele frequency out of the iterate file: $alleleFreq
    real_iterate <- IterateHWE(myRAD[[m]])
    posterprob[[m]] <- real_iterate$posteriorProb
  },error = function(e){cat("ERROR :",conditionMessage(e), "\n")} )
}

# posterprob[[1]][[1]][1,1:10,1:10] will show colnames and rownames

PPequalzero <- list()
for (n in 1:length(posterprob)) {
  alleleloc <- myRAD[[n]]$alleles2loc
  PPequalzero[[n]] <- sapply(1:nLoci(myRAD[[n]]),function(x)rowSums(posterprob[[n]][[1]][1,,alleleloc == x] < 0.1))
  colnames(PPequalzero[[n]]) <- GetLoci(myRAD[[n]])
}

diff_alc_pp <- list()
count <- 1
for (o in 1:length(allelecount)) {
  for (p in 1:nrow(allelecount[[o]])){
    diff_alc_pp[[count]] <- sum(allelecount[[o]][p,] != PPequalzero[[o]][p,])
    count <- count+1
  }
}

# align posterprob with locicountsum and 