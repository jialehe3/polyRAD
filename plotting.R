SNO.1 <- read.csv("SNO.1.CSV")
SNO.1$chromosome <- as.integer(gsub(SNO.1$chromosome, patter = "SCAFFOLD", replacement = "999"))
dont <- -which(is.na(SNO.1$chromosome))
SNO.1.1 <- SNO.1[dont,]

library(qqman)

# getting values from csv
names(SNO.1.1)
pval11 <- SNO.1.1$pvalue1to1_ref.over.alt
pval13 <- SNO.1.1$pvalue1to3_ref.over.alt
pval31 <- SNO.1.1$pvalue3to1_ref.over.alt
MAFplt <- SNO.1.1$MinorAlleleFrequency
Rdepth <- SNO.1.1$Depth
chromplt <- SNO.1.1$chromosome
posplt <- SNO.1.1$position
ratio<- SNO.1.1$alt.over.total

# get values fit 11
pselect11 <- pval11 >= 0.05
SNO.1.1$match11 = pselect11
# get values not fit 11 but 13
pselect13 <- pval13 >= 0.05
SNO.1.1$match13 = pselect13

# get values not fit 11 but 31
pselect31 <- pval31 >= 0.05
SNO.1.1$match31 = pselect31

# do qqplots with above values
qq(pselect11)
qq(pselect13)
qq(pselect31)

# make a conclusion column
SNO.1.1$conclusion <- ifelse(
  (SNO.1.1$match11 == TRUE), "1:1" ,
  ifelse((SNO.1.1$match13 == TRUE), "1:3" ,
         ifelse((SNO.1.1$match31 == TRUE), "3:1", "significant"
         )))
# getting significant
sig <- which(SNO.1.1$conclusion == 'significant')
SNO.1.2 <- SNO.1.1[sig,]

oneone <- which(SNO.1.1$conclusion == '1:1')
SNO.11 <- SNO.1.1[oneone,]

onethree <- which(SNO.1.1$conclusion == "1:3")
SNO.13 <- SNO.1.1[onethree,]

threeone <- which(SNO.1.1$conclusion == "3:1")
SNO.31 <- SNO.1.1[threeone,]

write.csv(SNO.1.1,file = "SNO.1.1.CSV")
write.csv(SNO.1.2, file = "SNO1.2.sig.CSV")

# plottings
pval11plt <- qq(pval11, main = "p_value 1 to 1")
pval13plt <- qq(pval13, main = "p_value 1 to 3")
pval31plt <- qq(pval31, main = "p_value 3 to 1")
pval11overMAF <- plot(log(pval11),log(MAFplt))
pval11overRD <- plot(log(pval11),log(Rdepth))
pval13overMAF <- plot(log(pval13),log(MAFplt))
pval13overRD <- plot(log(pval13),log(Rdepth))
pval31overMAF <- plot(log(pval31),log(MAFplt))
pval31overRD <- plot(log(pval31),log(Rdepth))


# violin plots
ggplot(SNO.1.1)+
  geom_violin(aes(x = conclusion, y = Depth))+
  coord_trans(y = "log")

# separate plots by 
ggplot(SNO.1.1, aes(x=log(Depth), fill=conclusion))+
  geom_density(alpha=0.5)+
    labs(title= "densityplot", x="Depth")

# make plots color transparent
plot(SNO.1.2$alt.over.total,log(SNO.1.2$Depth),col = "#00000050")
plot(SNO.1.1$alt.over.total,log(SNO.1.1$Depth),col = "#00000020")

# boxplot
boxplot(log(SNO.1.1$Depth)~SNO.1.1$conclusion)

manhattan(SNO.1.1, chr = "chromosome", bp = "position", p = "pvalue1to1_ref.over.alt", ylim = c(0,350))
