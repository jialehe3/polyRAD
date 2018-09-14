library(polyRAD)
myvcf <- "/Users/hejiale/Desktop/For Jiale/170608Msi_PstI_genotypes.vcf"
clengths <- read.csv("chromosome_lengths.CSV")
clengths$chr <- gsub(clengths$chr, patter = "Chr", replacement = "")
clengths$chr <- gsub(clengths$chr, patter="scaffold",replacement = "SCAFFOLD")
looptimes <- clengths$length/10000000
count <- 1
# seting memories to make loop faster
out <- list()
length(out) <- 200000
lessthan2 <- list()
length(lessthan2) <- 200000

RADpos <- data.frame (
  Chr = character(200000),
  Pos = integer(200000),
  stringsAsFactors = FALSE
)
# big loop
for(j in 1:nrow(clengths)) {
  k=looptimes[j]
  for (l in 0:as.integer(k)){
    # tryCatch errors when nothing pass through the threshold
    out[[count]] <- tryCatch({
      #if (count == 5){
      #  break
      #}else{
      mysvp2 <- ScanVcfParam(fixed = "ALT", info = NA, geno = "AD", which = GRanges(clengths$chr[j], IRanges(l*10000000,10000000+l*10000000-1)), samples = samples(scanVcfHeader(myvcf)))
      
      myRAD <- VCF2RADdata(myvcf, svparam = mysvp2, yieldSize = NA_integer_)
      
      Loci <- nLoci(myRAD)
      
      RADpos[count + (1:nrow(myRAD$locTable))-1,] <-  myRAD$locTable
      
      # small loop, loop by loci
      for(i in 1:Loci){
        my.matrix <- myRAD$alleleDepth[,myRAD$alleles2loc == i]
        # true for less or equal to 2 alleles
        my.target <- rowSums(my.matrix != 0) <= 2
        lessthan2[[count]] <- my.target
        count <- count +1
      }
      
      }, error= function(e){list(c("i"= i,"j"= j,"k" = k,"l" = l), e)})
  }}
#}

# get rid off null values created by the step that preset memories
lessthan2 <- Filter(Negate(is.null), lessthan2)
falsecount <- sapply(lessthan2, function(x) {sum(x == FALSE)})

# get rid off empty characters and empty integers created by the preset step
RAD.chr <- RADpos$Chr[RADpos$Chr !=""]
RAD.pos <- RADpos$Pos[RADpos$Pos !=0]

write.csv(data.frame(Chr = RAD.chr, Pos= RAD.pos, num.of.false = falsecount), file = "polyRAD.csv")

morethan2 <- subset(read.csv("polyRAD.CSV"), select = -c(X))
morethan2$Chr <- as.numeric(morethan2$Chr)
morethan2$Pos <- as.numeric(morethan2$Pos)

manhattan(morethan2, chr = "Chr", bp = "Pos", p = "num.of.false", ylim = c(0,500),logp = FALSE)

# matching SNO.1.2 and polyRAD/morethan2, getting all significant matches
a <- left_join(morethan2, SNO.1.2, by= c("Pos" = "position"))
csv.match <- subset(na.omit(a),select = -c(X))

# matching SNO.11 and morethan2, all 11 matches
b <- left_join(morethan2, SNO.11, by= c("Pos" = "position"))
csv.match11 <- subset(na.omit(b),select = -c(X))

# matching SNO.13 and mt2, all 13 matches
c <- left_join(morethan2, SNO.13, by= c("Pos" = "position"))
csv.match13 <- subset(na.omit(c),select = -c(X))

# matching SNO.31 and mt2, all 31 matches
d <- left_join(morethan2, SNO.31, by= c("Pos" = "position"))
csv.match31 <- subset(na.omit(d),select = -c(X))

# matching SNO.1.1 and mt2, all signicifant, 11, 13 and 31 matches
e <- left_join(morethan2, SNO.1.1, by= c("Pos" = "position"))
csv.matchall <- subset(na.omit(e), select = -c(X))


manhattan(csv.match, chr = "Chr", bp = "Pos", p = "num.of.false", ylim = c(0,500),logp = FALSE)


plot(csv.match13$alt.over.total,log(csv.match13$Depth),col = "#00000050")

ggplot(csv.match, aes(x=log(Depth)))+
  geom_density(fill="blue", alpha=0.5)+
    geom_density(data=csv.match11,fill="purple",aes(x=log(Depth)),alpha=0.5)+
    geom_density(data=csv.match13,fill="yellow",aes(x=log(Depth)),alpha=0.5)+
    geom_density(data=csv.match31,fill="red",aes(x=log(Depth)),alpha=0.5)+
  labs(title= "densityplot", x="Depth")

boxplot(csv.matchall$num.of.false~csv.matchall$conclusion)

ggplot(csv.matchall, aes(x=log(num.of.false+1), fill=conclusion))+
  geom_density(alpha=0.5)+
  labs(title= "densityplot", x="No.false")

# look into values 0 falsecount and diferent CATs
look.into <- which(csv.matchall$num.of.false == 0 & csv.matchall$conclusion == "significant")
look.result <- csv.matchall[look.into,]

look.into11 <- which(csv.match11$num.of.false == 0)
look.result11 <- csv.match11[look.into11,]

look.into13 <- which(csv.match13$num.of.false == 0)
look.result13 <- csv.match13[look.into13,]

look.into31 <- which(csv.match31$num.of.false == 0)
look.result31 <- csv.match31[look.into31,]

# plot all four into one and compare
par(mfrow = c(2,2))
plot(look.result$MissingByMarker,log(look.result$Depth),col = "#00000050")
plot(look.result11$MissingByMarker,log(look.result11$Depth),col = "#00000050")
plot(look.result13$MissingByMarker,log(look.result13$Depth),col = "#00000050")
plot(look.result31$MissingByMarker,log(look.result31$Depth),col = "#00000050")

# MAF around 0.5, MBM and Depth
MAF.look <- which(csv.matchall$MinorAlleleFrequency>0.45 & csv.matchall$MinorAlleleFrequency< 0.55 )
MAF.result <- csv.matchall[MAF.look,]

par(mfrow = c(1,1))
boxplot(csv.matchall$MinorAlleleFrequency ~ csv.matchall$conclusion)
boxplot(MAF.result$MinorAlleleFrequency ~ MAF.result$conclusion)

ggplot(csv.matchall)+
  geom_violin(aes(x = conclusion, y = MinorAlleleFrequency))
  
ggplot(MAF.result)+
  geom_violin(aes(x = conclusion, y = MinorAlleleFrequency))


plot(MAF.result$MinorAlleleFrequency,log(MAF.result$Depth),col = "#00000050")

par(mfrow = c(2,2))
plot(lm(look.result$MissingByMarker ~ log(look.result$Depth)))
