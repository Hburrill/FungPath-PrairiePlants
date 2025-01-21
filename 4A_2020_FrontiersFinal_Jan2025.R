# 2020 re-sampling

# load packages ----
library(vegan)
library(stats)
library(emmeans)
library(ggplot2)
library(tidyverse)
library(writexl)
library(ggpubr)
library(devtools)
library(factoextra)
library(ggfortify)
library(dplyr)
library(ecodist)
library(plyr)
if(Sys.getenv("JAVA_HOME")!=""){
  Sys.setenv(JAVA_HOME="")
}
library(rJava)
library(glmulti)
library(lme4)
library(nlme)

# summarySE code source: https://stackoverflow.com/questions/30485154/how-to-write-summaryse-function-in-r
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

# load otus ----
### fungi ----
setwd("C:/Users/haley/OneDrive - University of Kansas/2020data/Fungi1mismatch_output")
metafung <- read.csv("DP4A_meta_2020fungnoKOE.csv", header=T)
metafung$PlantFamily <- factor(metafung$PhyloFam, levels=c('AST','FAB','POA','MIX'))
metafung$PlantFamily <- revalue(metafung$PlantFamily, c("AST"="Asteraceae", "FAB"="Fabaceae", "POA"="Poaceae", "MIX"="Mixture"))
metafung$Precipitation <- factor(metafung$Precip, levels=c('50%','150%'))

metafung_2020 <- read.csv("DP4A_meta_2020fung-w2020.csv", header=T)
# merge in R, pairing using OTU_ID
# first remove line 1 in notepad++
# add _ between OTU and ID
# import table
#otu <- read.table("fungtaxmerge2020_funguild_tab.txt", header=TRUE)
#head(otu)
# # Make row 1 colnames
# rownames(otu) <- otu$OTU_ID
# otu$OTU_ID <- NULL
# head(otu)
# import taxonomy
#funguild <- read.table("fungtaxmerge2020_funguild.taxa.guilds.csv", head=T, sep=",")
#head(funguild)
#funguildOTU <- merge(otu, funguild, by="OTU_ID", all.x=TRUE)
#head(funguildOTU)
# export merged file for funguild
#write.table(funguildOTU, file='funguild_output2020.txt')
# # in excel remove 'confidence' column, upload to funguild, matched & unmatched
allfung <- read.delim("fungaltraits_allfungi_Feb28redo.csv", sep=",", header=T)
head(allfung)

rownames(allfung) <- allfung$OTU_ID
allfung$OTU_ID <- NULL
allfungOTU <- allfung[,-c(90,110,173,193,240:257)]
metafung$div <- diversity(t(allfungOTU))
# remove F115, 135, 198, 218
# create new 'pathogen' column; 0 = not plant pathogen, 1 = yes
# primary_lifestyle
# pathogens: plant_pathogen, unspecified_pathotroph

allfung$ppath2 <- ifelse(allfung$primary_lifestyle=="plant_pathogen" | allfung$Secondary_lifestyle=="plant_pathogen", "1", "0")

## calculate relative abundance of pathogens/sample
d <- colnames(allfung[,sapply(allfung, is.numeric)])
n<-length(d)
columns = c("pprel_ab","pprich","ppdiv") 
outputpath = data.frame(matrix(nrow = n, ncol = length(columns)))
colnames(outputpath) = columns

for(i in 1:n){
  s <- sum(allfung[which(allfung$ppath2=="1"), d[i]])
  df <- as.data.frame(s)
  outputpath$pprich[i] <- df[[1]]}
for(i in 1:n){
  s <- sum(allfung[which(allfung$ppath2=="1"), d[i]])/sum(allfung[,d[i]])
  df <- as.data.frame(s)
  outputpath$pprel_ab[i] <- df[[1]]}
for(i in 1:n){
  s <- diversity(allfung[which(allfung$ppath2=="1"), d[i]])
  df <- as.data.frame(s)
  outputpath$ppdiv[i] <- df[[1]]}

outputpath$Sample <- d # now you have a dataframe with sample ID as one column and relative abundance of rhizobia reads as the other column

# fungpathOTU <- read.delim("fungaltraitsfungpath.csv", sep=",", header=T)
# head(fungpathOTU)
# rownames(fungpathOTU) <- fungpathOTU$OTU_ID
# fungpathOTU$OTU_ID <- NULL

ppath <- allfung[which(allfung$ppath2=="1"),]
ppath <- ppath[,-c(90,110,173,193,240:260)]
ppath <- t(ppath)
ppath_rclr <- transform(x=ppath, "rclr")
ppath_rclr2 <- mutate_all(ppath_rclr, function(x) as.numeric(as.character(x)))
rowSums(ppath_rclr2[,-1])
#ppath_rclr3 <- ppath_rclr2[rowSums(ppath_rclr2[, -1])>0, ]
ppath_ATCH <- vegdist(ppath_rclr2[,-1], method="robust.aitchison")
metafung$pathdiv <- diversity(ppath)
outputpath <- outputpath[-c(90,110,173,193),]
metafung$pathrel_ab <- outputpath$pprel_ab
#meta_ppath <- merge(metafung, output, by="Sample", sort=F)


# create new 'saprotroph' column; 0 = not sap, 1 = yes
# saprotrophs: litter_saprotroph, soil_saprotroph, wood_saprotroph, unspecified_saprotroph

allfung$sap <- ifelse(allfung$primary_lifestyle=="litter_saprotroph"| allfung$primary_lifestyle=="soil_saprotroph" | allfung$primary_lifestyle=="wood_saprotroph" , "1", "0")

## calculate relative abundance of pathogens/sample
e <- colnames(allfung[,sapply(allfung, is.numeric)])
o<-length(e)
columns = c("saprel_ab","saprich","sapdiv") 
outputsap = data.frame(matrix(nrow = o, ncol = length(columns)))
colnames(outputsap) = columns

for(i in 1:n){
  s <- sum(allfung[which(allfung$sap=="1"), d[i]])
  df <- as.data.frame(s)
  outputsap$saprich[i] <- df[[1]]}
for(i in 1:n){
  s <- sum(allfung[which(allfung$sap=="1"), d[i]])/sum(allfung[,d[i]])
  df <- as.data.frame(s)
  outputsap$saprel_ab[i] <- df[[1]]}
for(i in 1:n){
  s <- diversity(allfung[which(allfung$sap=="1"), d[i]])
  df <- as.data.frame(s)
  outputsap$sapdiv[i] <- df[[1]]}

outputsap$Sample <- e # now you have a dataframe with sample ID as one column and relative abundance of rhizobia reads as the other column

# fungsapOTU <- read.delim("fungaltraitsfungsap.csv", sep=",", header=T)
# head(fungsapOTU)
# rownames(fungsapOTU) <- fungsapOTU$OTU_ID
# fungsapOTU$OTU_ID <- NULL

sapsub <- allfung[which(allfung$sap=="1"),]
sapsub <- sapsub[,-c(90,110,173,193,240:259)]
sapsub <- t(sapsub)
fsap_rclr <- transform(x=sapsub, "rclr")
fsap_rclr2 <- mutate_all(fsap_rclr, function(x) as.numeric(as.character(x)))
fsap_ATCH <- vegdist(fsap_rclr2[,-1], method="robust.aitchison")
metafung$sapdiv <- diversity(sapsub)
outputsap <- outputsap[-c(90,110,173,193),]
metafung$saprel_ab <- outputsap$saprel_ab

#meta_fsap <- merge(metafung, outputsap, by="Sample", sort=F)

#meta_fsap$Sap2pathRich <- meta_fsap$saprich/meta_ppath$pprich

# sap <- read.table("fungsaps_Dim4A2020noKOE.csv", sep=",", header=T)
# head(sap)
# rownames(sap) <- sap$OTU_ID
# sap$OTU_ID <- NULL
# sap <- t(sap)
# # add sp rich
# 
# # perm realized density as predictors (yr 2) 
# # then do yr 3 and see how the pathogen load differs
# # for cortin and chafas
# 
# pathfung <- read.table("fungpath4A2020_noKOE.csv", sep=",", header=T)
# rownames(pathfung) <- pathfung$OTU_ID
# pathfung$OTU_ID <- NULL
# pathfung <- t(pathfung)
# 
# Fungdiv <- diversity(sap)
# metafung$sapdiv <- Fungdiv
# Fungpathdiv <- diversity(pathfung)
# metafung$Fungpathdiv <- Fungpathdiv

# from vegdist() vignette:
# Aitchison (1986) distance is equivalent to Euclidean distance between CLR-transformed samples ("clr") and deals with positive compositional data. Robust Aitchison distance by Martino et al. (2019) uses robust CLR ("rlcr"), making it applicable to non-negative data including zeroes (unlike the standard Aitchison).
# sapf_rclr <- transform(sap, "rclr")
# sapfun_ATCH <- vegdist(sapf_rclr, method="robust.aitchison")
# 
# # fungal pathogens
# Pfun_rclr <- transform(pathfung, "rclr")
# Pfun_ATCH <- vegdist(Pfun_rclr, method="robust.aitchison")

### OOMYCETES ----
setwd("C:/Users/haley/OneDrive - University of Kansas/2020data/Oomy1mismatch_output")
oomy <- read.table("asv_Oomy2020-filt4AnoKOE.txt", header=T)
head(oomy)
rownames(oomy) <- oomy$OTU_ID
oomy$OTU_ID <- NULL
oomy <- t(oomy)

oomy_rclr <- transform(oomy, "rclr")
oomy_ATCH <- vegdist(oomy_rclr, method="robust.aitchison")

metaoo <- read.csv("DP4A_meta_2020oomynoKOE.csv", header=T)
metaoo$PlantFamily <- factor(metaoo$PhyloFam, levels=c('AST','FAB','POA','MIX'))
metaoo$PlantFamily <- revalue(metaoo$PlantFamily, c("AST"="Asteraceae", "FAB"="Fabaceae", "POA"="Poaceae", "MIX"="Mixture"))
metaoo$Precipitation <- factor(metaoo$Precip, levels=c('50%','150%'))

metaoo_2020 <- read.csv("DP4A_meta_2020oomy-w2020.csv", header=T)
oomydiv <- diversity(oomy)
metaoo$oomydiv <- oomydiv
metaoo_2020$oomydiv <- oomydiv

### BACTERIA ----
setwd("C:/Users/haley/OneDrive - University of Kansas/2020data/Bact1mismatch_output")
# bact <- read.table("otu_tableDP2020-bact_4AnoKOE.txt", header=T)
# head(bact)
# rownames(bact) <- bact$OTU_ID
# bact$OTU_ID <- NULL
# bact <- t(bact)
# bacts_rclr <- transform(x=bact, "rclr")
# bacts_ATCH <- vegdist(bacts_rclr[-1], method="robust.aitchison")

metab <- read.csv("DP4A_meta_2020bactnoKOE.csv", header=T)
metab$PlantFamily <- factor(metab$PhyloFam, levels=c('AST','FAB','POA','MIX'))
metab$PlantFamily <- revalue(metab$PlantFamily, c("AST"="Asteraceae", "FAB"="Fabaceae", "POA"="Poaceae", "MIX"="Mixture"))
metab$Precipitation <- factor(metab$Precip, levels=c('50%','150%'))

metab_2020 <- read.csv("DP4A_meta_2020bact-w2020.csv", header=T)

## merge taxonomy and OTU tables
otus <- read.delim("otu_tableDP2020-bact_4AnoKOE.txt", header=TRUE)
tax <- read.delim("DP2020bact-taxonomy.tsv", header=TRUE)

# bact div
bactOTU <- otus
rownames(bactOTU) <- bactOTU$OTU_ID
bactOTU$OTU_ID <- NULL
bactOTU <- t(bactOTU)
bactdiv <- diversity(bactOTU)
metab$bacdiv <- bactdiv
metab_2020$bacdiv <- bactdiv

# merge data frames by OTU_ID/Feature.ID
all <- merge(otus, tax, by.x = "OTU_ID", by.y = "Feature.ID")

## assign N-fixer status based on known genera
# split Taxon string into multiple columns
library(reshape2)
newColNames <- c("D0", "D1", "D2", "D3", "D4", "D5", "D6")
newCols <- colsplit(all$Taxon, ";", newColNames)
all_tax <- cbind(all, newCols)
all_tax$Taxon <- NULL
all_tax <- all_tax[,-239] # remove 'confidence' column

# create new 'nfixer' column; 0 = not rhizobia, 1 = rhizobia
all_tax$nfixer <- ifelse(all_tax$D5=="D_5__Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium"| 
                           all_tax$D5=="D_5__Bradyrhizobium" | all_tax$D5=="D_5__Ensifer" | 
                           all_tax$D5=="D_5__Burkholderia-Caballeronia-Paraburkholderia" | all_tax$D5=="D_5__Mesorhizobium", "1", "0")


## calculate relative abundance of N-fixers/sample
d <- colnames(all_tax[,sapply(all_tax, is.numeric)])
n<-length(d)
columns = c("rel_ab","richness","diversity") 
output = data.frame(matrix(nrow = n, ncol = length(columns)))
colnames(output) = columns

for(i in 1:n){
  s <- sum(all_tax[which(all_tax$nfixer=="1"), d[i]])
  df <- as.data.frame(s)
  output$richness[i] <- df[[1]]}
for(i in 1:n){
  s <- sum(all_tax[which(all_tax$nfixer=="1"), d[i]])/sum(all_tax[,d[i]])
  df <- as.data.frame(s)
  output$rel_ab[i] <- df[[1]]}
for(i in 1:n){
  s <- diversity(all_tax[which(all_tax$nfixer=="1"), d[i]])
  df <- as.data.frame(s)
  output$diversity[i] <- df[[1]]}

output$Sample <- d # now you have a dataframe with sample ID as one column and relative abundance of rhizobia reads as the other column

rhizobia <- all_tax[which(all_tax$nfixer=="1"),]
rownames(rhizobia) <- rhizobia$OTU_ID
rhizobia$OTU_ID <- NULL
rhizobia <- t(rhizobia[,-c(1,234:241)])
rhizobia_rclr <- transform(x=rhizobia, "rclr")
str(rhizobia_rclr)
rhizobia_rclr <- mutate_all(rhizobia_rclr, function(x) as.numeric(as.character(x)))
rhiz_rclr2 <- rhizobia_rclr[-c(233:236),-1]
rhizobia_ATCH <- vegdist(rhizobia_rclr[-c(233:236),-1], method="robust.aitchison")
metarhiz <- merge(metab, output[-234,], by="Sample", sort=F)
metarhiz$Precipitation <- factor(metarhiz$Precip, levels=c('50%','150%'))

bact <- all_tax[which(all_tax$nfixer=="0"),]
head(bact)
rownames(bact) <- bact$OTU_ID
bact$OTU_ID <- NULL
bact <- t(bact[,-c(234:241)])
bact_rclr <- transform(x=bact, "rclr")
bact_rclr <- mutate_all(bact_rclr, function(x) as.numeric(as.character(x)))
bacts_ATCH <- vegdist(bact_rclr[-1], method="robust.aitchison")

### AMF ----
setwd("C:/Users/haley/OneDrive - University of Kansas/2020data/AMF1mismatch_output")
AMF <- read.table("2020AMF-otu97table_clean_4A-noKOE.txt", header=T)
head(AMF)
rownames(AMF) <- AMF$OTU_ID
AMF$OTU_ID <- NULL
AMF <- t(AMF)
AMF_rclr <- transform(x=AMF, "rclr")
AMF_ATCH <- vegdist(AMF_rclr[-1], method="robust.aitchison")

metamf <- read.csv("DP4A_meta_2020amfnoKOE.csv", header=T)
metamf$PlantFamily <- factor(metamf$PhyloFam, levels=c('AST','FAB','POA','MIX'))
metamf$PlantFamily <- revalue(metamf$PlantFamily, c("AST"="Asteraceae", "FAB"="Fabaceae", "POA"="Poaceae", "MIX"="Mixture"))
metamf$Precipitation <- factor(metamf$Precip, levels=c('50%','150%'))

metamf_2020 <- read.csv("DP4A_meta_2020amf-w2020.csv", header=T)
AMFdiv <- diversity(AMF)
metamf$AMFdiv <- AMFdiv
metamf_2020$amfdiv <- AMFdiv


# diversity models ----
setwd("C:/Users/haley/OneDrive - University of Kansas/2020data")

## fungi ----
allfundiv <- aov(div ~ Block + PhyloFam*R_SpRich*Precip +  ANDGER + ELYCAN + BOUGRA + PANVIR + AMOCAN + DALCAN + DALPUR + DESILL + DESCAN + CHAFAS + LIAPYC + CORTIN + ECHPAL + EUPALT + SILINT + HELMOL, data=metafung)
anova(allfundiv)
fungdiv_out <- as.data.frame(anova(allfundiv))
write_xlsx(fungdiv_out, "fungdiv_out.xlsx")

sapdiv <- aov(sapdiv ~ Block + PhyloFam*R_SpRich*Precip +  ANDGER + ELYCAN + BOUGRA + PANVIR + AMOCAN + DALCAN + DALPUR + DESILL + DESCAN + CHAFAS + LIAPYC + CORTIN + ECHPAL + EUPALT + SILINT + HELMOL, data=metafung)
anova(sapdiv)
sapdiv_out <- as.data.frame(anova(sapdiv))
write_xlsx(sapdiv_out, "sapdiv_out.xlsx")

res <- glmulti(sapdiv ~ Block + PhyloFam*R_SpRich*Precip +  ANDGER + ELYCAN + BOUGRA + PANVIR + AMOCAN + DALCAN + DALPUR + DESILL + DESCAN + CHAFAS + LIAPYC + CORTIN + ECHPAL + EUPALT + SILINT + HELMOL, data=metafung, level=1, fitfunction=glm, crit="aicc") # 1|Subblock*Precip

# best fit:
# sapdiv ~ 1 + Precip + ELYCAN + BOUGRA + DALPUR

sapdivBF <- lmer(sapdiv ~ 1 + Precip + ELYCAN + BOUGRA + DALPUR + (1|Subblock:Precip), data=metafung)
sapdivcoeff <- as.data.frame(coef(summary(sapdivBF)))
sapdivcoeff$p <- 2 * (1 - pnorm(abs(sapdivcoeff$`t value`)))
write_xlsx(sapdivcoeff, "sapdiv_bestfit.xlsx")

sapRA <- aov(saprel_ab ~ Block + PhyloFam*R_SpRich*Precip + SCHSCO + ANDGER + ELYCAN + BOUGRA + PANVIR + AMOCAN + DALCAN + DALPUR + DESILL + DESCAN + CHAFAS + LIAPYC + CORTIN + ECHPAL + EUPALT + SILINT + HELMOL, data=metafung)
anova(sapRA)
sapRA_out <- as.data.frame(anova(sapRA))
write_xlsx(sapRA_out, "sapRA_out.xlsx")

RAres <- glmulti(saprel_ab ~ Block + PhyloFam*R_SpRich*Precip +  ANDGER + ELYCAN + BOUGRA + PANVIR + AMOCAN + DALCAN + DALPUR + DESILL + DESCAN + CHAFAS + LIAPYC + CORTIN + ECHPAL + EUPALT + SILINT + HELMOL, data=metafung, level=1, fitfunction=glm, crit="aicc")
# best model:
# saprel_ab ~ 1 + Block + R_SpRich + Precip + ANDGER + ELYCAN + DALPUR + DESILL + CORTIN

sapRAbf <- lmer(saprel_ab ~ 1 + Block + R_SpRich + Precip + ANDGER + ELYCAN + DALPUR + DESILL + CORTIN + (1|Subblock:Precip), data=metafung)
sapRAcoeff <- as.data.frame(coef(summary(sapRAbf)))
sapRAcoeff$p <- 2 * (1 - pnorm(abs(sapRAcoeff$`t value`)))
write_xlsx(sapRAcoeff, "sapRA_bestfit.xlsx")

ggplot(metafung, aes(x=Precip,y=log(saprel_ab), color=Precip))+ geom_boxplot(fill=NA) +geom_point(alpha=0.3,aes(color=Precip),position=position_jitter())+labs(title="Differences in fungal saprobe diversity between precipitation treatments")

metafung$Precip <- factor(metafung$Precip, levels=c('50%','150%'))
### sapRA ----
emm_options(rg.limit = 2097152)
ref.123<-emmeans(lm(saprel_ab ~ Block + PlantFamily*Precipitation, data=metafung),~PlantFamily*Precipitation)
ref.table.123<-as.data.frame(ref.123)
p123<-ggplot(ref.table.123, aes(x=PlantFamily, y=emmean, color=PlantFamily)) + 
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), width=0.5,size=1, position=position_dodge(1)) + facet_wrap(~Precipitation) +
  labs(title="Fungal saprotrophs", x="Plant Family Composition", y = "Mean relative abundance") + theme(axis.text.x = element_text(angle = 60, vjust = 0.5), legend.position = "none")
p123
p124<-ggplot(ref.table.123, aes(x=Precipitation, y=log(emmean), color=Precipitation)) +
  geom_boxplot() + 
  labs(title="Saprotroph relative abundance by precipitation treatment", x="Precipitation treatments", y = "Mean relative abundance") + theme(axis.text.x = element_text(angle = 60, vjust = 0.5), legend.position = "none")
p124


meansapRA <- summarySE(metafung, measurevar="saprel_ab", groupvars=c("Precip","PlantFamily"))

RA1 <- ggplot(meansapRA, aes(x=Precip,y=saprel_ab, color=Precip))+ geom_boxplot(fill=NA) +geom_point(alpha=0.3,aes(color=Precip),position=position_jitter())+labs(title="Fungal saprotroph relative abundance between precipitation treatments", x="Precipitation Treatment", y = "Mean Relative Abundance")

# vs realized
fundivR <- aov(metafung$sapdiv ~ Block + PhyloFam*R_SpRich*Precip + SCHSCO_bio + ANDGER_bio + ELYCAN_bio + BOUGRA_bio + PANVIR_bio + AMOCAN_bio + DALCAN_bio + DALPUR_bio + DESILL_bio + DESCAN_bio + CHAFAS_bio + LIAPYC_bio + CORTIN_bio + ECHPAL_bio + EUPALT_bio + SILINT_bio + HELMOL_bio, data=metafung_2020)
summary(fundivR)
cor.test(sapdiv$residuals, fundivR$residuals)
# p-value < 2.2e-16, cor 0.9701726

# pathogens
funpathdiv <- aov(pathdiv ~ Block + PhyloFam*R_SpRich*Precip + SCHSCO + ANDGER + ELYCAN + BOUGRA + PANVIR + AMOCAN + DALCAN + DALPUR + DESILL + DESCAN + CHAFAS + LIAPYC + CORTIN + ECHPAL + EUPALT + SILINT + HELMOL, data=metafung)
summary(funpathdiv)
funpathdiv_out <- as.data.frame(anova(funpathdiv))
write_xlsx(funpathdiv_out, "funpathdiv_out.xlsx")

res2 <- glmulti(pathdiv ~ Block + PhyloFam*R_SpRich*Precip +  ANDGER + ELYCAN + BOUGRA + PANVIR + AMOCAN + DALCAN + DALPUR + DESILL + DESCAN + CHAFAS + LIAPYC + CORTIN + ECHPAL + EUPALT + SILINT + HELMOL, data=metafung, level=1, fitfunction=glm, crit="aicc")
# best model:
# pathdiv~1+Block+PhyloFam+BOUGRA+PANVIR

funpathdivBF <- lmer(pathdiv ~ 1+Block+PhyloFam+BOUGRA+PANVIR + (1|Subblock:Precip), data=metafung)
funpathdivcoeff <- as.data.frame(coef(summary(funpathdivBF)))
funpathdivcoeff$p <- 2 * (1 - pnorm(abs(funpathdivcoeff$`t value`)))
write_xlsx(funpathdivcoeff, "funpathdiv_bestfit.xlsx")

ggplot(metafung, aes(x=PhyloFam,y=pathdiv, color=PhyloFam))+ geom_boxplot(fill=NA) +geom_point(alpha=0.3,aes(color=PhyloFam),position=position_jitter())+labs(title="Fungal pathogen diversity among plant family treatments", x="Plant Family Composition", y = "Diversity (H')")

funpathRA <- aov(pathrel_ab ~ Block + PhyloFam*R_SpRich*Precip + SCHSCO + ANDGER + ELYCAN + BOUGRA + PANVIR + AMOCAN + DALCAN + DALPUR + DESILL + DESCAN + CHAFAS + LIAPYC + CORTIN + ECHPAL + EUPALT + SILINT + HELMOL, data=metafung)
summary(funpathRA)
funpathRA_out <- as.data.frame(anova(funpathRA))
write_xlsx(funpathRA_out, "funpathRA_out.xlsx")

RAres2 <- glmulti(pathrel_ab ~ Block + PhyloFam*R_SpRich*Precip +  ANDGER + ELYCAN + BOUGRA + PANVIR + AMOCAN + DALCAN + DALPUR + DESILL + DESCAN + CHAFAS + LIAPYC + CORTIN + ECHPAL + EUPALT + SILINT + HELMOL, data=metafung, level=1, fitfunction=glm, crit="aicc")
# best model:
# pathrel_ab ~ 1 + PhyloFam + BOUGRA

funpathRAbf <- aov(pathrel_ab ~ 1 + PhyloFam + BOUGRA + Subblock*Precip, data=metafung)
summary(funpathRAbf)
funpathRA_out <- as.data.frame(anova(funpathRAbf))
write_xlsx(funpathRA_out, "funpathRA_bf.xlsx")

funpathRAbf <- lmer(pathrel_ab ~ 1 + PhyloFam + BOUGRA + (1|Subblock:Precip), data=metafung)
funpathRAcoeff <- as.data.frame(coef(summary(funpathRAbf)))
funpathRAcoeff$p <- 2 * (1 - pnorm(abs(funpathRAcoeff$`t value`)))
write_xlsx(funpathRAcoeff, "funpathRA_bf.xlsx")

meanpathRA <- summarySE(metafung, measurevar="pathrel_ab", groupvars=c("Precip","PhyloFam"))

RA2 <- ggplot(meanpathRA, aes(x=PhyloFam,y=pathrel_ab, color=PhyloFam))+ geom_boxplot(fill=NA) +geom_point(alpha=0.3,aes(color=PhyloFam),position=position_jitter())+labs(title="Fungal pathogen relative abundance between plant family treatments", x="Plant Family Composition", y = "Mean Relative Abundance")

# meanpathdiv <- summarySE(metafung, measurevar="pathdiv", groupvars=c("Precip","PlantFamily"))

emm_options(rg.limit = 2097152)
### pathRA ----
ref.122<-emmeans(lm(pathrel_ab ~ Block + PlantFamily*R_SpRich*Precipitation, data=metafung),~PlantFamily*R_SpRich*Precipitation)
ref.table.122<-as.data.frame(ref.122)
p122<-ggplot(ref.table.122, aes(x=PlantFamily, y=emmean, color=PlantFamily)) + geom_point(position=position_dodge(1), size =4) + geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), width=0.5,size=1, position=position_dodge(1)) + facet_wrap(~Precipitation) + labs(title="Fungal plant pathogens", x="Plant Family Composition", y = "Mean relative abundance") + theme(axis.text.x = element_text(angle = 60, vjust = 0.5), legend.position = "none")
p122

p125<-ggplot(ref.table.122, aes(x=Precip, y=log(emmean), color=Precip)) + 
  geom_boxplot() + 
  labs(title="Pathogen relative abundance by precipitation treatment", x="Precipitation treatments", y = "Mean relative abundance")
p125

funpathdivR <- aov(metafung$pathdiv ~ Block + PhyloFam*R_SpRich*Precip + SCHSCO_bio + ANDGER_bio + ELYCAN_bio + BOUGRA_bio + PANVIR_bio + AMOCAN_bio + DALCAN_bio + DALPUR_bio + DESILL_bio + DESCAN_bio + CHAFAS_bio + LIAPYC_bio + CORTIN_bio + ECHPAL_bio + EUPALT_bio + SILINT_bio + HELMOL_bio, data=metafung_2020)
summary(funpathdivR)
cor.test(funpathdiv$residuals, funpathdivR$residuals)
# p-value < 2.2e-16, cor 0.9762391

#sap2pathmod <- aov(Sap2pathRich ~ Block + PhyloFam*R_SpRich*Precip + SCHSCO + ANDGER + ELYCAN + BOUGRA + PANVIR + AMOCAN + DALCAN + DALPUR + DESILL + DESCAN + CHAFAS + LIAPYC + CORTIN + ECHPAL + EUPALT + SILINT + HELMOL, data=meta_fsap)
#summary(sap2pathmod)


## oomy ----
oomydiv <- aov(oomydiv ~ Block + PhyloFam*R_SpRich*Precip + SCHSCO + ANDGER + ELYCAN + BOUGRA + PANVIR + AMOCAN + DALCAN + DALPUR + DESILL + DESCAN + CHAFAS + LIAPYC + CORTIN + ECHPAL + EUPALT + SILINT + HELMOL, data=metaoo)
summary(oomydiv) # precip and phylo
oomydiv_out <- as.data.frame(anova(oomydiv))
write_xlsx(oomydiv_out, "oomydiv_out.xlsx")

res3 <- glmulti(oomydiv ~ Block + PhyloFam*R_SpRich*Precip +  ANDGER + ELYCAN + BOUGRA + PANVIR + AMOCAN + DALCAN + DALPUR + DESILL + DESCAN + CHAFAS + LIAPYC + CORTIN + ECHPAL + EUPALT + SILINT + HELMOL, data=metaoo, level=1, fitfunction=glm, crit="aicc")
# best model:
# oomydiv~1+PhyloFam+Precip+ANDGER

oomydivBF <- lmer(oomydiv ~ 1+PhyloFam+Precip+ANDGER + (1|Subblock:Precip), data=metaoo)
oomydivcoeff <- as.data.frame(coef(summary(oomydivBF)))
oomydivcoeff$p <- 2 * (1 - pnorm(abs(oomydivcoeff$`t value`)))
write_xlsx(oomydivcoeff, "oomydiv_bestfit.xlsx")

metaoo$Precip <- factor(metaoo$Precip, levels=c('50%','150%'))
ggplot(metaoo, aes(x=Precip,y=oomydiv, color=Precip))+ geom_boxplot(fill=NA) +geom_point(alpha=0.3,aes(color=Precip),position=position_jitter())+labs(title="Differences in oomycete diversity by precip")

ggplot(metaoo, aes(x=PhyloFam,y=oomydiv, color=Precip))+ geom_boxplot(fill=NA) +geom_point(alpha=0.3,aes(color=Precip),position=position_jitter())+labs(title="Differences in oomycete diversity by family mixes")

ggplot(metaoo, aes(x=PhyloFam,y=oomydiv, color=PhyloFam))+ geom_boxplot(fill=NA) +geom_point(alpha=0.3,aes(color=PhyloFam),position=position_jitter())+labs(title="Differences in oomycete diversity by family mixes")

metaoo50 <- metaoo[which(metaoo$Precip=="drought"),]
ggplot(metaoo50, aes(x=R_SpRich,y=oomydiv, color=PhyloFam))+ geom_smooth(method = 'glm', level=0.95, formula = y ~ x, se=T, aes(fill=PhyloFam)) +geom_point(alpha=0.3,aes(color=PhyloFam),position=position_jitter())+labs(title="Oomycete diversity response to Sp Rich and Family, 50% precip")

metaoo150 <- metaoo[(metaoo$Precip=="overwater"),]
ggplot(metaoo150, aes(x=R_SpRich,y=oomydiv, color=PhyloFam))+ geom_smooth(method = 'glm', level=0.95, formula = y ~ x, se=T, aes(fill=PhyloFam)) +geom_point(alpha=0.3,aes(color=PhyloFam),position=position_jitter())+labs(title="Oomycete diversity response to Sp Rich and Family, 150% precip")

# vs realized
oomydivR <- aov(oomydiv ~ Block + PhyloFam*R_SpRich*Precip + SCHSCO_bio + ANDGER_bio + ELYCAN_bio + BOUGRA_bio + PANVIR_bio + AMOCAN_bio + DALCAN_bio + DALPUR_bio + DESILL_bio + DESCAN_bio + CHAFAS_bio + LIAPYC_bio + CORTIN_bio + ECHPAL_bio + EUPALT_bio + SILINT_bio + HELMOL_bio, data=metaoo_2020)
summary(oomydivR)
cor.test(oomydiv$residuals, oomydivR$residuals)

## bact ----
# rhizobia
rhizRA <- aov(log(rel_ab) ~ Block + PhyloFam*R_SpRich*Precip + SCHSCO + ANDGER + ELYCAN + BOUGRA + PANVIR + AMOCAN + DALCAN + DALPUR + DESILL + DESCAN + CHAFAS + LIAPYC + CORTIN + ECHPAL + EUPALT + SILINT + HELMOL, data=metarhiz)
summary(rhizRA)
rhizRA_out <- as.data.frame(anova(rhizRA))
write_xlsx(rhizRA_out, "rhizRA_out.xlsx")

RAres3 <- glmulti(rel_ab ~ Block + PhyloFam*R_SpRich*Precip +  ANDGER + ELYCAN + BOUGRA + PANVIR + AMOCAN + DALCAN + DALPUR + DESILL + DESCAN + CHAFAS + LIAPYC + CORTIN + ECHPAL + EUPALT + SILINT + HELMOL, data=metarhiz, level=1, fitfunction=glm, crit="aicc")
# best model:
# rel_ab ~ 1+Block+PhyloFam+ANDGER+DALPUR

rhizRAbf <- lmer(rel_ab ~ 1+Block+PhyloFam+ANDGER+DALPUR + (1|Subblock:Precip), data=metarhiz)
rhizRAcoeff <- as.data.frame(coef(summary(rhizRAbf)))
rhizRAcoeff$p <- 2 * (1 - pnorm(abs(rhizRAcoeff$`t value`)))
write_xlsx(rhizRAcoeff, "rhizRA_out.xlsx")

meanrhizRA <- summarySE(metarhiz, measurevar="rel_ab", groupvars=c("Precip","PhyloFam"))
RA3 <- ggplot(meanrhizRA, aes(x=PhyloFam,y=rel_ab, color=PhyloFam))+ geom_boxplot(fill=NA) +geom_point(alpha=0.3,aes(color=PhyloFam),position=position_jitter())+labs(title="Rhizobia bacteria relative abundance between plant family treatments", x="Plant Family Composition", y = "Mean Relative Abundance")

emm_options(rg.limit = 2097152)
### RhizRA ----
ref.10<-emmeans(lm(rel_ab ~ Block + PlantFamily*R_SpRich*Precipitation, data=metarhiz),~PlantFamily*R_SpRich*Precipitation)
ref.table.10<-as.data.frame(ref.10)
p10<-ggplot(ref.table.10, aes(x=PlantFamily, y=emmean, color=PlantFamily)) + geom_point(position=position_dodge(1), size =4) + geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), width=0.5,size=1, position=position_dodge(1)) + facet_wrap(~Precipitation) + labs(title="Rhizobia bacteria", x="Plant Family Composition", y = "Mean relative abundance") + theme(axis.text.x = element_text(angle = 60, vjust = 0.5), legend.position = "none")
p10
ggsave("rhizRA_phylo.png", p10, scale = 1, width = 10, height = 8, units ="in", dpi = 300, limitsize = TRUE)

rhizdiv <- aov(diversity ~ Block + PhyloFam*R_SpRich*Precip + SCHSCO + ANDGER + ELYCAN + BOUGRA + PANVIR + AMOCAN + DALCAN + DALPUR + DESILL + DESCAN + CHAFAS + LIAPYC + CORTIN + ECHPAL + EUPALT + SILINT + HELMOL, data=metarhiz)
summary(rhizdiv)
write_xlsx(as.data.frame(anova(rhizdiv)), "rhizdiv_out.xlsx")

res4 <- glmulti(diversity ~ Block + PhyloFam*R_SpRich*Precip +  ANDGER + ELYCAN + BOUGRA + PANVIR + AMOCAN + DALCAN + DALPUR + DESILL + DESCAN + CHAFAS + LIAPYC + CORTIN + ECHPAL + EUPALT + SILINT + HELMOL, data=metarhiz, level=1, fitfunction=glm, crit="aicc")
# best model:
# diversity~1+Block+Precip+BOUGRA+DALPUR

rhizdivBF <- lmer(rel_ab ~ 1+Block+Precip+BOUGRA+DALPUR + (1|Subblock:Precip), data=metarhiz)
rhizdivcoeff <- as.data.frame(coef(summary(rhizdivBF)))
rhizdivcoeff$p <- 2 * (1 - pnorm(abs(rhizdivcoeff$`t value`)))
write_xlsx(rhizdivcoeff, "rhizdiv_bestfit.xlsx")

metarhiz$Precip <- factor(metarhiz$Precip, levels=c('50%','150%'))
ggplot(metarhiz, aes(x=Precip,y=diversity, color=Precip))+ geom_boxplot(fill=NA) +geom_point(alpha=0.3,aes(color=Precip),position=position_jitter())+labs(title="Rhizobia species diversity between precipitation treatments")

bactdiv <- aov(bactdiv ~ Block + PhyloFam*R_SpRich*Precip + SCHSCO + ANDGER + ELYCAN + BOUGRA + PANVIR + AMOCAN + DALCAN + DALPUR + DESILL + DESCAN + CHAFAS + LIAPYC + CORTIN + ECHPAL + EUPALT + SILINT + HELMOL, data=metab)
summary(bactdiv)
bactdiv_out <- as.data.frame(anova(bactdiv))
write_xlsx(bactdiv_out, "bactdiv_out.xlsx")

res5 <- glmulti(bactdiv ~ Block + PhyloFam*R_SpRich*Precip +  ANDGER + ELYCAN + BOUGRA + PANVIR + AMOCAN + DALCAN + DALPUR + DESILL + DESCAN + CHAFAS + LIAPYC + CORTIN + ECHPAL + EUPALT + SILINT + HELMOL, data=metab, level=1, fitfunction=glm, crit="aicc")
# best model:
# bactdiv~1+Block+Precip+AMOCAN+DALCAN

bactdivBF <- lmer(bactdiv ~ 1+Block+Precip+AMOCAN+DALCAN + (1|Subblock:Precip), data=metab)
bactdivcoeff <- as.data.frame(coef(summary(bactdivBF)))
bactdivcoeff$p <- 2 * (1 - pnorm(abs(bactdivcoeff$`t value`)))
write_xlsx(bactdivcoeff, "bactdiv_bestfit.xlsx")

metab$Precip <- factor(metab$Precip, levels=c('50%','150%'))
ggplot(metab, aes(x=Precip,y=bactdiv, color=Precip))+ geom_boxplot(fill=NA) +geom_point(alpha=0.3,aes(color=Precip),position=position_jitter())+labs(title="Differences in bacteria diversity by precip")

metab50 <- metab[which(metab$Precip=="drought"),]
ggplot(metab50, aes(x=R_SpRich,y=bactdiv, color=PhyloFam))+ geom_smooth(method = 'glm', level=0.95, formula = y ~ x, se=T, aes(fill=PhyloFam)) +geom_point(alpha=0.3,aes(color=PhyloFam),position=position_jitter())+labs(title="Bacteria diversity response to Sp Rich and Family, 50% precip")

metab150 <- metab[(metab$Precip=="overwater"),]
ggplot(metab150, aes(x=R_SpRich,y=bactdiv, color=PhyloFam))+ geom_smooth(method = 'glm', level=0.95, formula = y ~ x, se=T, aes(fill=PhyloFam)) +geom_point(alpha=0.3,aes(color=PhyloFam),position=position_jitter())+labs(title="Bacteria diversity response to Sp Rich and Family, 150% precip")

metarhiz2020 <- merge(output[-234,], metab_2020, by="Sample", sort=F)
rhizdivR <- aov(diversity ~ Block + PhyloFam*R_SpRich*Precip + SCHSCO_bio + ANDGER_bio + ELYCAN_bio + BOUGRA_bio + PANVIR_bio + AMOCAN_bio + DALCAN_bio + DALPUR_bio + DESILL_bio + DESCAN_bio + CHAFAS_bio + LIAPYC_bio + CORTIN_bio + ECHPAL_bio + EUPALT_bio + SILINT_bio + HELMOL_bio, data=metarhiz2020)
summary(rhizdivR)
cor.test(rhizdiv$residuals, rhizdivR$residuals)

bactdivR <- aov(bactdiv ~ Block + PhyloFam*R_SpRich*Precip + SCHSCO_bio + ANDGER_bio + ELYCAN_bio + BOUGRA_bio + PANVIR_bio + AMOCAN_bio + DALCAN_bio + DALPUR_bio + DESILL_bio + DESCAN_bio + CHAFAS_bio + LIAPYC_bio + CORTIN_bio + ECHPAL_bio + EUPALT_bio + SILINT_bio + HELMOL_bio, data=metab_2020)
summary(bactdivR)
cor.test(bactdiv$residuals, bactdivR$residuals)

## amf ----
amfdiv <- aov(AMFdiv ~ Block + PhyloFam*R_SpRich*Precip + SCHSCO + ANDGER + ELYCAN + BOUGRA + PANVIR + AMOCAN + DALCAN + DALPUR + DESILL + DESCAN + CHAFAS + LIAPYC + CORTIN + ECHPAL + EUPALT + SILINT + HELMOL, data=metamf)
summary(amfdiv)
amfdiv_out <- as.data.frame(anova(amfdiv))
write_xlsx(amfdiv_out, "amfdiv_out.xlsx")

res6 <- glmulti(AMFdiv ~ Block + PhyloFam*R_SpRich*Precip +  ANDGER + ELYCAN + BOUGRA + PANVIR + AMOCAN + DALCAN + DALPUR + DESILL + DESCAN + CHAFAS + LIAPYC + CORTIN + ECHPAL + EUPALT + SILINT + HELMOL, data=metamf, level=1, fitfunction=glm, crit="aicc")
# best model:
# AMFdiv~1+Block+Precip+ELYCAN+PANVIR

amfdivBF <- lmer(AMFdiv ~ 1+Block+Precip+AMOCAN+DALCAN + (1|Subblock:Precip), data=metamf)
amfdivcoeff <- as.data.frame(coef(summary(amfdivBF)))
amfdivcoeff$p <- 2 * (1 - pnorm(abs(amfdivcoeff$`t value`)))
write_xlsx(amfdivcoeff, "amfdiv_bestfit.xlsx")


metamf$Precip <- factor(metamf$Precip, levels=c('50%','150%'))
ggplot(metamf, aes(x=Precip,y=AMFdiv, color=Precip))+ geom_boxplot(fill=NA) +geom_point(alpha=0.3,aes(color=Precip),position=position_jitter())+labs(title="Differences in AMF diversity by precip")

amfdivR <- aov(AMFdiv ~ Block + PhyloFam*R_SpRich*Precip + SCHSCO_bio + ANDGER_bio + ELYCAN_bio + BOUGRA_bio + PANVIR_bio + AMOCAN_bio + DALCAN_bio + DALPUR_bio + DESILL_bio + DESCAN_bio + CHAFAS_bio + LIAPYC_bio + CORTIN_bio + ECHPAL_bio + EUPALT_bio + SILINT_bio + HELMOL_bio, data=metamf_2020)
summary(amfdivR)
cor.test(amfdiv$residuals, amfdivR$residuals)

## diversity figure ----
# # code used to save diversity response
# write_xlsx(metafung,"MetaSoilFungwDiv.xlsx")
# write_xlsx(metaoo,"MetaSoilOomywDiv.xlsx")
# write_xlsx(metab,"MetaSoilBactwDiv.xlsx")
# write_xlsx(metaas,"MetaSoilAMFwDiv.xlsx")

#fungcoeff <- as.data.frame(summary(fundiv)$coefficients)
#write_xlsx(fungcoeff,"fungcoeff.xlsx")
#Rfungcoeff <- as.data.frame(summary(fundivR)$coefficients)
#write_xlsx(Rfungcoeff,"Rfungcoeff.xlsx")

meansapdivprecip <- summarySE(metafung, measurevar="sapdiv", groupvars="Precip")
p1<-ggplot(meansapdivprecip, aes(x=Precip, y=sapdiv, color=Precip)) +  geom_point(position=position_dodge(1), size =4) + geom_errorbar(aes(ymin=sapdiv-se,ymax=sapdiv+se), width=0.5,size=1, position=position_dodge(1)) + labs(title="Fungal saprobe diversity response to precipitation", x="Precipitation treatment (% ambient)", y = "Shannon-Index of Fungal Saprotrophs") + scale_color_manual(name="Precip", values=c("red","blue")) + theme(legend.position = "none")
p1

#fungpathcoeff <- as.data.frame(summary(funpathdiv)$coefficients)
#write_xlsx(fungpathcoeff,"fungpathcoeff.xlsx")
#Rfungpathcoeff <- as.data.frame(summary(funpathdivR)$coefficients)
#write_xlsx(Rfungpathcoeff,"Rfungpathcoeff.xlsx")
meanpathdiv <- summarySE(metafung, measurevar="pathdiv", groupvars="PlantFamily")
pathdiv2 <- ggplot(meanpathdiv, aes(x=PlantFamily,y=pathdiv, color=PlantFamily))  +  geom_point(position=position_dodge(1), size =4) + geom_errorbar(aes(ymin=pathdiv-se,ymax=pathdiv+se), width=0.5,size=1, position=position_dodge(1))+labs(title="Fungal pathogen diversity between plant family treatments", x="Plant Family Composition", y = "Shannon-Index of Fungal Pathogens")  + theme(axis.text.x = element_text(angle = 60, vjust = 0.5), plot.title = element_text(size=13), legend.position = "none")
pathdiv2


#oomycoeff <- as.data.frame(summary(oomydiv)$coefficients)
#write_xlsx(oomycoeff,"oomycoeff.xlsx")
#Roomycoeff <- as.data.frame(summary(oomydivR)$coefficients)
#write_xlsx(Roomycoeff,"Roomycoeff.xlsx")
ref.3<-emmeans(lm(oomydiv ~ Block + PlantFamily*R_SpRich*Precip, data=metaoo),~PlantFamily*R_SpRich*Precip)
ref.table.3<-as.data.frame(ref.3)
p3<-ggplot(ref.table.3, aes(x=Precip, y=emmean, color=PlantFamily)) + 
  geom_point(position=position_dodge(1), size =4) + 
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), width=0.5,size=1, position=position_dodge(1)) +
  labs(title="Oomycete diversity response to precipitation", x="Precipitation treatment (% ambient)", y = "Shannon-Index of Oomycetes")   + theme(plot.title = element_text(size=13), legend.position = "bottom")
p3

ref.8<-emmeans(lm(oomydiv ~ Block + PlantFamily*R_SpRich, data=metaoo),~PlantFamily*R_SpRich)
ref.table.8<-as.data.frame(ref.8)
p8<-ggplot(ref.table.8, aes(x=PlantFamily, y=emmean, color=PlantFamily)) +  geom_point(position=position_dodge(1), size =4) + geom_errorbar(aes(ymin=lower.CL,ymax=upper.CL), width=0.5,size=1, position=position_dodge(1)) +
  labs(title="Oomycete Diversity response", x="Plant family composition", y = "Shannon-Index of Oomycetes")
p8

p6 <- ggplot(ref.table.3, aes(x=Precip, y=emmean, color=Precip)) +  geom_point() + geom_errorbar(aes(ymin=lower.CL,ymax=upper.CL)) + labs(title="Oomycete diversity response to precipitation", x="Precipitation treatment", y = "Shannon-Index of Oomycetes") + scale_color_manual(name="Precip", values=c("red","blue")) + theme(legend.position = "none") + facet_wrap(~factor(PlantFamily, levels=c('Asteraceae','Fabaceae','Poaceae','Mixture')), ncol=2, nrow=2)
p6

#bactdivcoeff <- as.data.frame(summary(bactdiv)$coefficients)
#write_xlsx(bactdivcoeff,"bactcoeff.xlsx")
#Rbactdivcoeff <- as.data.frame(summary(bactdivR)$coefficients)
#write_xlsx(Rbactdivcoeff,"Rbactcoeff.xlsx")
meanrhizdiv <- summarySE(metarhiz, measurevar="diversity", groupvars="Precip")
p4<-ggplot(meanrhizdiv, aes(x=Precip, y=diversity, color=Precip)) +  geom_point(position=position_dodge(1), size =4) + geom_errorbar(aes(ymin=diversity-se,ymax=diversity+se), width=0.5,size=1, position=position_dodge(1)) + 
  labs(title="Rhizobia diversity response to precipitation", x="Precipitation treatment (% ambient)", y = "Shannon-Index of Rhizobia Bacteria") +
  scale_color_manual(name="Precip", values=c("red","blue")) + theme(legend.position = "none")
p4

meanbactdiv <- summarySE(metab, measurevar="bacdiv", groupvars="Precip")
p5<-ggplot(meanbactdiv, aes(x=Precip, y=bacdiv, color=Precip))+  geom_point(position=position_dodge(1), size =4) + geom_errorbar(aes(ymin=bacdiv-se,ymax=bacdiv+se), width=0.5,size=1, position=position_dodge(1)) + 
  labs(title="Bacteria diversity response to precipitation", x="Precipitation treatment (% ambient)", y = "Shannon-Index of Bacteria") +
  scale_color_manual(name="Precip", values=c("red","blue")) + theme(legend.position = "none")
p5

p7<- ggplot(metab, aes(x=R_SpRich, y=bactdiv, color=Precip)) +  geom_point() + geom_smooth(method="lm", se=F) + labs(title="Bacteria diversity response to precipitation", x="Precipitation treatment", y = "Shannon-Index of Bacteria") + scale_color_manual(name="Precip", values=c("red","blue")) + facet_wrap(~factor(PhyloFam, levels=c('AST','FAB','POA','MIX')), ncol=2, nrow=2)
p7

#amfcoeff <- as.data.frame(summary(amfdiv)$coefficients)
#write_xlsx(amfcoeff,"amfcoeff.xlsx")
#Ramfcoeff <- as.data.frame(summary(amfdivR)$coefficients)
#write_xlsx(Ramfcoeff,"Ramfcoeff.xlsx")
meanamfdiv <- summarySE(metamf, measurevar="AMFdiv", groupvars="Precip")
p9<-ggplot(meanamfdiv, aes(x=Precip, y=AMFdiv, color=Precip)) +  geom_point(position=position_dodge(1), size =4) + geom_errorbar(aes(ymin=AMFdiv-se,ymax=AMFdiv+se), width=0.5,size=1, position=position_dodge(1)) + 
  labs(title="AMF diversity response to precipitation", x="Precipitation treatment (% ambient)", y = "Shannon-Index of AMF") +
  scale_color_manual(name="Precip", values=c("red","blue")) + theme(legend.position = "none")
p9

## save figs 1 & 2 ----
# combine into 1 fig
theme_set(theme_pubr())
setwd("C:/Users/haley/OneDrive - University of Kansas/2020data")
combdivs <- ggarrange(p1, pathdiv2, p4, p5, p3, p9, labels = c("A", "B","C","D","E","F"), ncol = 2, nrow = 3)
ggsave(combdivs, filename="combdivs_wrhizOct29.jpg", scale = 1, width = 15, height = 20, units ="in", dpi = 300, limitsize = TRUE)

combdiv_int <- ggarrange(p6, p7, p8, labels = c("A", "B", "C"), ncol = 2, nrow = 2)
ggsave(combdiv_int, filename="combdiv_int.png", scale = 1, width = 20, height = 16, units ="in", dpi = 300, limitsize = TRUE)

# RA fig
theme_set(theme_pubr())
combRA <- ggarrange(p123, p122, p10, labels = c("A", "B","C"), ncol = 1, nrow = 3)
ggsave(combRA, filename="combRA-Nov19.jpg", scale = 1, width = 10, height = 15, units ="in", dpi = 300, limitsize = TRUE)

# to save individually
# ggsave(p1, filename="fungdiv.png", path = "C:/Users/hburrill/Desktop/KU/Bever lab/Dimensions/4A/2020seqreads/figures", scale = 1, width = 4, height = 3, units ="in", dpi = 300, limitsize = TRUE)

# permanova ----
set.seed(1)
sap_th <- adonis2(fsap_ATCH ~ Block + PhyloFam*R_SpRich*Precip +  ANDGER + ELYCAN + BOUGRA + PANVIR + AMOCAN + DALCAN + DALPUR + DESILL + DESCAN + CHAFAS + LIAPYC + CORTIN + ECHPAL + EUPALT + SILINT + HELMOL, data=metafung)
sap_th
FSperm <- as.data.frame(sap_th)
write_xlsx(FSperm, "FSperm.xlsx")

mod0 <- rda(sapsub ~ 1, metafung)  # Model with intercept only
mod1 <- rda(sapsub ~ Block + PhyloFam*R_SpRich*Precip +  ANDGER + ELYCAN + BOUGRA + PANVIR + AMOCAN + DALCAN + DALPUR + DESILL + DESCAN + CHAFAS + LIAPYC + CORTIN + ECHPAL + EUPALT + SILINT + HELMOL, metafung)  # Model with all explanatory variables

## With scope present, the default direction is "both"
ordistep(mod0, scope = formula(mod1), perm.max = 200)

## Example without scope. Default direction is "backward"
ordistep(mod1, perm.max = 200)

sap_fin <- adonis2(fsap_ATCH ~ R_SpRich + Precip + DALCAN + EUPALT, data=metafung, by = "terms")
sap_fin
FSperm <- as.data.frame(sap_fin)
write_xlsx(FSperm, "FSperm_ordistep.xlsx")

sap_cov <- adonis2(fsap_ATCH ~ Block + PhyloFam*R_SpRich*Precip +  ANDGER_cov + ELYCAN_cov + BOUGRA_cov + PANVIR_cov + AMOCAN_cov + DALCAN_cov + DALPUR_cov + DESILL_cov + DESCAN_cov + CHAFAS_cov + LIAPYC_cov + CORTIN_cov + ECHPAL_cov + EUPALT_cov + SILINT_cov + HELMOL_cov, data=metafung_2020)
sap_cov

cor.test(sap_th$R2, sap_cov$R2) # 0.999 without SCHSCO

# betadisp
FungalSaprotrophDispersion <- betadisper(fsap_ATCH, metafung$PhyloFam)
anova(FungalSaprotrophDispersion) # ns precip # phylofam 0.06
plot(FungalSaprotrophDispersion, hull=FALSE, ellipse=TRUE) ##sd ellipse

# PERMANOVA pathogens
funpath_th <- adonis2(ppath_ATCH ~ Block + PhyloFam*R_SpRich*Precip +  ANDGER + ELYCAN + BOUGRA + PANVIR + AMOCAN + DALCAN + DALPUR + DESILL + DESCAN + CHAFAS + LIAPYC + CORTIN + ECHPAL + EUPALT + SILINT + HELMOL, data=metafung)
funpath_th
FPperm <- as.data.frame(funpath_th)
write_xlsx(FPperm, "FPperm.xlsx")

# mod0 <- rda(ppath ~ 1, metafung)  # Model with intercept only
mod1 <- rda(ppath ~ Block + PhyloFam*R_SpRich*Precip +  ANDGER + ELYCAN + BOUGRA + PANVIR + AMOCAN + DALCAN + DALPUR + DESILL + DESCAN + CHAFAS + LIAPYC + CORTIN + ECHPAL + EUPALT + SILINT + HELMOL, metafung)  # Model with all explanatory variables

## Example without scope. Default direction is "backward"
ordistep(mod1, perm.max = 200)

path_fin <- adonis2(ppath_ATCH ~ Block + PhyloFam + Precip + PhyloFam:Precip, data=metafung, by = "terms")
path_fin
FPperm <- as.data.frame(path_fin)
write_xlsx(FPperm, "FPperm_ordistep.xlsx")

funpath_cov <- adonis2(ppath_ATCH ~ Block + PhyloFam*R_SpRich*Precip +  ANDGER_cov + ELYCAN_cov + BOUGRA_cov + PANVIR_cov + AMOCAN_cov + DALCAN_cov + DALPUR_cov + DESILL_cov + DESCAN_cov + CHAFAS_cov + LIAPYC_cov + CORTIN_cov + ECHPAL_cov + EUPALT_cov + SILINT_cov + HELMOL_cov, data=metafung_2020)
funpath_cov

cor.test(funpath_th$R2, funpath_cov$R2) # 0.999

# betadisp
pathdisp <- betadisper(ppath_ATCH, metafung$Precip)
anova(pathdisp) # precip ns
plot(pathdisp, hull=FALSE, ellipse=TRUE) ##sd ellipse

# oomycetes
oomy_th <- adonis2(oomy_ATCH ~ Block + PhyloFam*R_SpRich*Precip + ANDGER + ELYCAN + BOUGRA + PANVIR + AMOCAN + DALCAN + DALPUR + DESILL + DESCAN + CHAFAS + LIAPYC + CORTIN + ECHPAL + EUPALT + SILINT + HELMOL, data=metaoo)
oomy_th
OOperm <- as.data.frame(oomy_th)
write_xlsx(OOperm, "OOperm.xlsx")

mod1 <- rda(oomy ~ Block + PhyloFam*R_SpRich*Precip +  ANDGER + ELYCAN + BOUGRA + PANVIR + AMOCAN + DALCAN + DALPUR + DESILL + DESCAN + CHAFAS + LIAPYC + CORTIN + ECHPAL + EUPALT + SILINT + HELMOL, metaoo)  # Model with all explanatory variables

## Example without scope. Default direction is "backward"
ordistep(mod1, perm.max = 200)

oomy_fin <- adonis2(oomy_ATCH ~ Block + PhyloFam + Precip +  PhyloFam:Precip, data=metaoo, by = "terms")
oomy_fin
oomyperm <- as.data.frame(oomy_fin)
write_xlsx(oomyperm, "oomyperm_ordistep.xlsx")

oomy_cov <- adonis2(oomy_ATCH ~ Block + PhyloFam*R_SpRich*Precip + ANDGER_cov + ELYCAN_cov + BOUGRA_cov + PANVIR_cov + AMOCAN_cov + DALCAN_cov + DALPUR_cov + DESILL_cov + DESCAN_cov + CHAFAS_cov + LIAPYC_cov + CORTIN_cov + ECHPAL_cov + EUPALT_cov + SILINT_cov + HELMOL_cov, data=metaoo_2020)
oomy_cov

cor.test(oomy_th$R2, oomy_cov$R2) # 0.999

OomyceteDispersion <- betadisper(oomy_ATCH, metaoo$Precip)
anova(OomyceteDispersion) # precip betadisp 0.067
plot(OomyceteDispersion, hull=FALSE, ellipse=TRUE) ##sd ellipse

# bacteria
rhiz_th <- adonis2(rhizobia_ATCH ~ Block + PhyloFam*R_SpRich*Precip + ANDGER + ELYCAN + BOUGRA + PANVIR + AMOCAN + DALCAN + DALPUR + DESILL + DESCAN + CHAFAS + LIAPYC + CORTIN + ECHPAL + EUPALT + SILINT + HELMOL, data=metarhiz[-c(1:2),])
rhiz_th
write_xlsx(as.data.frame(rhiz_th), "rhizperm.xlsx")

mod1 <- rda(rhizobia ~ Block + PhyloFam*R_SpRich*Precip +  ANDGER + ELYCAN + BOUGRA + PANVIR + AMOCAN + DALCAN + DALPUR + DESILL + DESCAN + CHAFAS + LIAPYC + CORTIN + ECHPAL + EUPALT + SILINT + HELMOL, metarhiz[-1,])  # Model with all explanatory variables

## Example without scope. Default direction is "backward"
ordistep(mod1, perm.max = 200)

rhiz_fin <- adonis2(rhizobia_ATCH ~ Block + Precip + AMOCAN, data=metarhiz[-1,], by = "terms")
rhiz_fin
rhizperm <- as.data.frame(rhiz_fin)
write_xlsx(rhizperm, "rhizperm_ordistep.xlsx")

rhiz_cov <- adonis2(rhizobia_ATCH ~ Block + PhyloFam*R_SpRich*Precip + ANDGER_cov + ELYCAN_cov + BOUGRA_cov + PANVIR_cov + AMOCAN_cov + DALCAN_cov + DALPUR_cov + DESILL_cov + DESCAN_cov + CHAFAS_cov + LIAPYC_cov + CORTIN_cov + ECHPAL_cov + EUPALT_cov + SILINT_cov + HELMOL_cov, data=metarhiz2020[-c(1:2),])
rhiz_cov
cor.test(rhiz_th$R2, rhiz_cov$R2) # 0.999
RhizobiaDispersion <- betadisper(rhizobia_ATCH, metarhiz[-c(1),]$Precip)
anova(RhizobiaDispersion) # .0398
plot(RhizobiaDispersion, hull=FALSE, ellipse=TRUE) ##sd ellipse

bact_th <- adonis2(bacts_ATCH ~ Block + PhyloFam*R_SpRich*Precip + ANDGER + ELYCAN + BOUGRA + PANVIR + AMOCAN + DALCAN + DALPUR + DESILL + DESCAN + CHAFAS + LIAPYC + CORTIN + ECHPAL + EUPALT + SILINT + HELMOL, data=metab)
bact_th
bacperm <- as.data.frame(bact_th)
write_xlsx(bacperm, "bacperm.xlsx")

mod1 <- rda(bact ~ Block + PhyloFam*R_SpRich*Precip +  ANDGER + ELYCAN + BOUGRA + PANVIR + AMOCAN + DALCAN + DALPUR + DESILL + DESCAN + CHAFAS + LIAPYC + CORTIN + ECHPAL + EUPALT + SILINT + HELMOL, metab)  # Model with all explanatory variables

## Example without scope. Default direction is "backward"
ordistep(mod1, perm.max = 200)

bact_fin <- adonis2(bacts_ATCH ~ Block + PhyloFam + Precip + BOUGRA + PANVIR + LIAPYC, data=metab, by = "terms")
bact_fin
bactperm <- as.data.frame(bact_fin)
write_xlsx(bactperm, "bactperm_ordistep.xlsx")

bact_cov <- adonis2(bacts_ATCH ~ Block + PhyloFam*R_SpRich*Precip + ANDGER_cov + ELYCAN_cov + BOUGRA_cov + PANVIR_cov + AMOCAN_cov + DALCAN_cov + DALPUR_cov + DESILL_cov + DESCAN_cov + CHAFAS_cov + LIAPYC_cov + CORTIN_cov + ECHPAL_cov + EUPALT_cov + SILINT_cov + HELMOL_cov, data=metab_2020)
bact_cov

cor.test(bact_th$R2, bact_cov$R2) # 0.999

bacdisp <- betadisper(bacts_ATCH, metab$Precip)
anova(bacdisp) # ns
bacdisp2 <- betadisper(bacts_ATCH, metab$PhyloFam)
anova(bacdisp2) # ns

amf_th <- adonis2(AMF_ATCH ~ Block + PhyloFam*R_SpRich*Precip + ANDGER + ELYCAN + BOUGRA + PANVIR + AMOCAN + DALCAN + DALPUR + DESILL + DESCAN + CHAFAS + LIAPYC + CORTIN + ECHPAL + EUPALT + SILINT + HELMOL, data=metamf)
amf_th
AMFperm <- as.data.frame(amf_th)
write_xlsx(AMFperm, "AMFperm.xlsx")

mod1 <- rda(AMF ~ Block + PhyloFam + R_SpRich + Precip + DALCAN + DESILL + LIAPYC + CORTIN, metamf)  # Model with all explanatory variables

mod0 <- rda(AMF ~ 1, metamf)  # Model with intercept only

## With scope present, the default direction is "both"
ordistep(mod0, scope = formula(mod1), perm.max = 200)

## Example without scope. Default direction is "backward"
ordistep(mod1, perm.max = 200)

amf_fin <- adonis2(AMF_ATCH ~ Block + PhyloFam + R_SpRich + Precip + DALCAN + DESILL + LIAPYC + CORTIN, data=metamf, by = "terms")
amf_fin
amfperm <- as.data.frame(amf_fin)
write_xlsx(amfperm, "amfperm_ordistep.xlsx")

amf_cov <- adonis2(AMF_ATCH ~ Block + PhyloFam*R_SpRich*Precip + ANDGER_cov + ELYCAN_cov + BOUGRA_cov + PANVIR_cov + AMOCAN_cov + DALCAN_cov + DALPUR_cov + DESILL_cov + DESCAN_cov + CHAFAS_cov + LIAPYC_cov + CORTIN_cov + ECHPAL_cov + EUPALT_cov + SILINT_cov + HELMOL_cov, data=metamf_2020)
amf_cov

cor.test(amf_th$R2, amf_cov$R2) # 0.999

AMFungalDispersion <- betadisper(AMF_ATCH, metamf$Precip)
anova(AMFungalDispersion) # precip 0.04, PhyloFam ns
plot(AMFungalDispersion, hull=FALSE, ellipse=TRUE)

## perm tables ----
setwd("C:/Users/haley/OneDrive - University of Kansas/2020data/figures")

sap_th %>%
  kbl(caption = "Saprotroph fungi Theoretical Response") %>%
  kable_classic(full_width = F, html_font = "Cambria") %>%
  column_spec(6, bold=T) %>%
  save_kable("fungsapTH.png")
sap_cov %>%
  kbl(caption = "Saprotroph fungi Realized Response") %>%
  kable_classic(full_width = F, html_font = "Cambria") %>%
  column_spec(6, bold=T) %>%
  save_kable("fungsapR.png")

funpath_th %>%
  kbl(caption = "Fungal Pathogen Theoretical Response") %>%
  kable_classic(full_width = F, html_font = "Cambria") %>%
  column_spec(6, bold=T) %>%
  save_kable("figures/fungpathTH.png")
funpath_cov %>%
  kbl(caption = "Fungal Pathogen Realized Response") %>%
  kable_classic(full_width = F, html_font = "Cambria") %>%
  column_spec(6, bold=T) %>%
  save_kable("figures/fungpathR.png")

oomy_th %>%
  kbl(caption = "Oomycete Theoretical Response") %>%
  kable_classic(full_width = F, html_font = "Cambria") %>%
  column_spec(6, bold=T) %>%
  save_kable("figures/oomyTH.png")
oomy_cov %>%
  kbl(caption = "Oomycete Realized Response") %>%
  kable_classic(full_width = F, html_font = "Cambria") %>%
  column_spec(6, bold=T) %>%
  save_kable("figures/oomyR.png")

bact_th %>%
  kbl(caption = "Bacteria Theoretical Response") %>%
  kable_classic(full_width = F, html_font = "Cambria") %>%
  column_spec(6, bold=T) %>%
  save_kable("figures/bactTH.png")
bact_cov %>%
  kbl(caption = "Bacteria Realized Response") %>%
  kable_classic(full_width = F, html_font = "Cambria") %>%
  column_spec(6, bold=T) %>%
  save_kable("figures/bactR.png")

amf_th %>%
  kbl(caption = "AMF Theoretical Response") %>%
  kable_classic(full_width = F, html_font = "Cambria") %>%
  column_spec(6, bold=T) %>%
  save_kable("figures/amfTH.png")
amf_cov %>%
  kbl(caption = "AMF Realized Response") %>%
  kable_classic(full_width = F, html_font = "Cambria") %>%
  column_spec(6, bold=T) %>%
  save_kable("figures/amfR.png")

# PCA 2020 only ----
## fungi ----
fungpath.pca <- prcomp(na.omit(ppath_ATCH), center = TRUE)
summary(fungpath.pca)
summary(aov(fungpath.pca$x[,6] ~ Block + PhyloFam*R_SpRich*Precip +  ANDGER + ELYCAN + BOUGRA + PANVIR + AMOCAN + DALCAN + DALPUR + DESILL + DESCAN + CHAFAS + LIAPYC + CORTIN + ECHPAL + EUPALT + SILINT + HELMOL, data=metafung)) # phylo 1, 3 # precip 2, 4 # sp rich 3, 5
#pf1 <- fungpath.pca$scores[,1]
#metafung$pf1 <- pf1

metafung$pf1 <- fungpath.pca$x[,1]
metafung$pf2 <- fungpath.pca$x[,2]
metafung$pf3 <- fungpath.pca$x[,3]
metafung$pf4 <- fungpath.pca$x[,4]
metafung$pf6 <- fungpath.pca$x[,6]

funpath_phylo <- ggplot(data=metafung, aes(x=pf1, y=pf3, color=PlantFamily)) + geom_point() + stat_ellipse(geom="polygon", aes(fill = PlantFamily), alpha = 0.05, show.legend = T, level = 0.95) + labs(title="Fungal pathogen composition among plant family mixtures",x="Dim 1 (71.9%)", y = "Dim 3 (2.3%)") + theme(plot.title = element_text(size=13), legend.position = "none")
funpath_phylo
metafung$Precip <- factor(metafung$Precip, level=c("50%","150%"))
funpath_precip <- ggplot(data=metafung, aes(x=pf4, y=pf6, color=Precip)) + geom_point() + stat_ellipse(geom="polygon", aes(fill=Precip), alpha = 0.05, show.legend = T, level = 0.95) + labs(title="Fungal plant pathogens",x="PC4 (1.9%)", y = "PC6 (1.3%)") + theme(plot.title = element_text(size=13), legend.position = "none")
funpath_precip

fungsap.pca <- prcomp(fsap_ATCH, center = TRUE)
summary(fungsap.pca)
summary(aov(fungsap.pca$x[,4] ~ Block + PhyloFam*R_SpRich*Precip +  ANDGER + ELYCAN + BOUGRA + PANVIR + AMOCAN + DALCAN + DALPUR + DESILL + DESCAN + CHAFAS + LIAPYC + CORTIN + ECHPAL + EUPALT + SILINT + HELMOL, data=metafung)) # phylo 1, 2 # precip 2, 3
#sf1 <- fungsap.pca$scores[,1] # 55.9%
#metafung$sf1 <- sf1

metafung$sf1 <- fungsap.pca$x[,1]
metafung$sf2 <- fungsap.pca$x[,2]
metafung$sf3 <- fungsap.pca$x[,3]
metafung$sf4 <- fungsap.pca$x[,4]

funsap_phylo <- ggplot(data=metafung, aes(x=sf1, y=sf2, color=PhyloFam)) + geom_point() + stat_ellipse(geom="polygon", aes(fill = PhyloFam), alpha = 0.05, show.legend = T, level = 0.95) +labs(title="Fungal saprotrophs",x="PC1 (59.5%)", y = "PC2 (4.6%)") + theme(plot.title = element_text(size=13), legend.position = "none")
funsap_phylo
funsap_precip <- ggplot(data=metafung, aes(x=sf3, y=sf4, color=Precip)) + geom_point() + stat_ellipse(geom="polygon", aes(fill = Precip), alpha = 0.05, show.legend = T, level = 0.95) +labs(title="Fungal saprotrophs",x="PC3 (3.4%)", y = "PC4 (1.7%)") + theme(plot.title = element_text(size=13), legend.position = "none")
funsap_precip

# metafung$PlntSpRich <- as.factor(metafung$R_SpRich)
# funsap_sprich <- ggplot(data=metafung, aes(x=sf3, y=sf4, color=PlntSpRich)) + geom_point() + stat_ellipse(geom="polygon", aes(fill = PlntSpRich), alpha = 0.05, show.legend = T, level = 0.95) + labs(title="Fungal saprotroph composition among plant species richness",x="Dim 3 (4.2%)", y = "Dim 4 (2.0%)") + theme(plot.title = element_text(size=13), legend.position = "top")
# funsap_sprich

# theme_set(theme_pubr())
# fungpcas <- ggarrange(pca1, pca2, labels = c("A", "B"), ncol = 2, nrow = 1)
# ggsave(fungpcas, filename="fungpcas.png", path = "C:/Users/hburrill/Desktop/KU/Bever lab/Dimensions/4A/2020seqreads/figures", scale = 1, width = 8, height = 3, units ="in", dpi = 300, limitsize = TRUE)

## oomy ----
oomy.pca<- prcomp(oomy_ATCH, center=T)
summary(oomy.pca)
summary(aov(oomy.pca$x[,2] ~ Block + PhyloFam*R_SpRich*Precip +  ANDGER + ELYCAN + BOUGRA + PANVIR + AMOCAN + DALCAN + DALPUR + DESILL + DESCAN + CHAFAS + LIAPYC + CORTIN + ECHPAL + EUPALT + SILINT + HELMOL, data=metaoo)) # precip 1, 2 # precip & phylo axes 7 & 8
metaoo$o1 <- oomy.pca$x[,1]
metaoo$o2 <- oomy.pca$x[,2]
metaoo$o7 <- oomy.pca$x[,7]
metaoo$o8 <- oomy.pca$x[,8]

metaoo$Precip <- factor(metaoo$Precip, level=c("50%","150%"))
oomy_precip <- ggplot(data=metaoo, aes(x=o1, y=o2, color=Precip)) + geom_point() + stat_ellipse(geom="polygon", aes(fill = Precip), alpha = 0.05, show.legend = T, level = 0.95) +labs(title="Oomycetes",x="PC1 (66.8%)", y = "PC2 (3.7%)") + theme(plot.title = element_text(size=13), legend.position = "none")
oomy_precip

precip.leg <- ggplot(data=metaoo, aes(x=o1, y=o2, color=Precip)) + geom_point() + stat_ellipse(geom="polygon", aes(fill = Precip), alpha = 0.05, show.legend = T, level = 0.95) +labs(title="Oomycetes",x="PC1 (66.8%)", y = "PC2 (3.7%)") + theme(plot.title = element_text(size=13), legend.position = "bottom") + scale_color_manual(values=c('red','blue'))
ggsave(precip.leg, filename="precipleg.jpg", scale = 1, width = 8, height = 10, units ="in", dpi = 300, limitsize = TRUE)

oomy_phylo <- ggplot(data=metaoo, aes(x=o1, y=o2, color=PhyloFam)) + geom_point() + stat_ellipse(geom="polygon", aes(fill = PhyloFam), alpha = 0.05, show.legend = T, level = 0.95) +labs(title="Oomycete composition between plant family treatments",x="Dim 1 (66.8%)", y = "Dim 2 (3.7%)") + theme(plot.title = element_text(size=13), legend.position = "none")
oomy_phylo

oomy_int50 <- ggplot(data=metaoo[(metaoo$Precip=="50%"),], aes(x=o7, y=o8, color=PhyloFam)) + geom_point() + stat_ellipse(geom="polygon", aes(fill = PhyloFam), alpha = 0.05, show.legend = T, level = 0.95) +labs(title="Precipitation 50% ambient",x="Dim 7 (1.3%)", y = "Dim 8 (1.1%)") + theme(plot.title = element_text(size=13), legend.position = "top")
oomy_int50

oomy_int150 <- ggplot(data=metaoo[(metaoo$Precip=="150%"),], aes(x=o7, y=o8, color=PhyloFam)) + geom_point() + stat_ellipse(geom="polygon", aes(fill = PhyloFam), alpha = 0.05, show.legend = T, level = 0.95) +labs(title="Precipitation 150% ambient",x="Dim 7 (1.3%)", y = "Dim 8 (1.1%)") + theme(plot.title = element_text(size=13), legend.position = "top")
oomy_int150

theme_set(theme_pubr())
oomyints_pc <- ggarrange(oomy_int50, oomy_int150, labels = c("A", "B"), ncol = 2, nrow = 1)
annotate_figure(oomyints_pc, top=text_grob("Oomycete composition between plant family and precipitation treatments", size=10))
ggsave(annotate_figure(oomyints_pc, top=text_grob("Oomycete composition between plant family and precipitation treatments", size=15)), filename="oomyints_pc.png", scale = 1, width = 8, height = 3, units ="in", dpi = 300, limitsize = TRUE)

## bact ----
bac.pca <- prcomp(bacts_ATCH, center=T)
summary(bac.pca)
head(summary(aov(bac.pca$x[,5] ~ Block + PhyloFam*R_SpRich*Precip +  ANDGER + ELYCAN + BOUGRA + PANVIR + AMOCAN + DALCAN + DALPUR + DESILL + DESCAN + CHAFAS + LIAPYC + CORTIN + ECHPAL + EUPALT + SILINT + HELMOL, data=metab))) # phylo 4, 5 # precip 1, 4
metab$b1 <- bac.pca$x[,1]
metab$b4 <- bac.pca$x[,4]

bac_phylo <- ggplot(data=metab, aes(x=b1, y=b4, color=PlantFamily)) + geom_point() + stat_ellipse(geom="polygon", aes(fill = PlantFamily), alpha = 0.05, show.legend = T, level = 0.95) +labs(title="Bacteria composition",x="PC1 (57.3%)", y = "PC4 (1.0%)") + theme(plot.title = element_text(size=13), legend.position = "none")
bac_phylo

metab$Precip <- factor(metab$Precip, level=c("50%","150%"))
bac_precip <- ggplot(data=metab, aes(x=b1, y=b4, color=Precip)) + geom_point() + stat_ellipse(geom="polygon", aes(fill = Precip), alpha = 0.05, show.legend = T, level = 0.95) +labs(title="Bacteria composition",x="PC1 (57.3%)", y = "PC4 (1.0%)")+ theme(plot.title = element_text(size=13), legend.position = "none")
bac_precip

rhiz.pca <- prcomp(rhizobia_ATCH, center=T)
summary(rhiz.pca)
head(summary(aov(rhiz.pca$x[,5] ~ Block + PhyloFam*R_SpRich*Precip +  ANDGER + ELYCAN + BOUGRA + PANVIR + AMOCAN + DALCAN + DALPUR + DESILL + DESCAN + CHAFAS + LIAPYC + CORTIN + ECHPAL + EUPALT + SILINT + HELMOL, data=metarhiz[-1,]))) # phylo 4, 5 # precip 1, 4
metarhiz_pca <- metarhiz[-1,]
metarhiz_pca$r1 <- rhiz.pca$x[,1]
metarhiz_pca$r5 <- rhiz.pca$x[,5]

metarhiz_pca$Precip <- factor(metarhiz_pca$Precip, level=c("50%","150%"))
rhiz_precip <- ggplot(data=metarhiz_pca, aes(x=r1, y=r5, color=Precip)) + geom_point() + stat_ellipse(geom="polygon", aes(fill = Precip), alpha = 0.05, show.legend = T, level = 0.95) +labs(title="Rhizobia composition",x="PC1 (63.1%)", y = "PC5 (2.1%)")+ theme(plot.title = element_text(size=13), legend.position = "none")
rhiz_precip

# theme_set(theme_pubr())
# bacpcas <- ggarrange(pca6, pca7, labels = c("G", "H"), ncol = 2, nrow = 1)
# ggsave(bacpcas, filename="bacpcas.png", path = "C:/Users/hburrill/Desktop/KU/Bever lab/Dimensions/4A/2020seqreads/figures", scale = 1, width = 8, height = 3, units ="in", dpi = 300, limitsize = TRUE)

## amf ----
amf.pca <- prcomp(AMF_ATCH, center=T)
summary(amf.pca)
summary(aov(amf.pca$x[,7] ~ Block + PhyloFam*R_SpRich*Precip +  ANDGER + ELYCAN + BOUGRA + PANVIR + AMOCAN + DALCAN + DALPUR + DESILL + DESCAN + CHAFAS + LIAPYC + CORTIN + ECHPAL + EUPALT + SILINT + HELMOL, data=metamf)) # phylo 4, 7 # precip 1, 2

metamf$pc1 <- amf.pca$x[,1]
metamf$pc2 <- amf.pca$x[,2]
metamf$pc4 <- amf.pca$x[,4]
metamf$pc7 <- amf.pca$x[,7]

amf_phylo <- ggplot(data=metamf, aes(x=pc4, y=pc7, color=PlantFamily)) + geom_point() + stat_ellipse(geom="polygon", aes(fill = PlantFamily), alpha = 0.05, show.legend = T, level = 0.95) +labs(title="AM Fungi",x="PC4 (1.9%)", y = "PC7 (1.4%)")+ theme(plot.title = element_text(size=13), legend.position = "bottom")
amf_phylo

metamf$Precip <- factor(metamf$Precip, level=c("50%","150%"))
amf_precip <- ggplot(data=metamf, aes(x=pc1, y=pc2, color=Precip)) + geom_point() + stat_ellipse(geom="polygon", aes(fill = Precip), alpha = 0.05, show.legend = T, level = 0.95) +labs(title="AM Fungi",x="PC1 (32.6%)", y = "PC2 (4.6%)")+ theme(plot.title = element_text(size=13), legend.position = "bottom")
amf_precip

# amfpcas <- ggarrange(pca8, pca9, labels = c("I", "J"), ncol = 2, nrow = 1)
# ggsave(amfpcas, filename="amfpcas.png", path = "C:/Users/hburrill/Desktop/KU/Bever lab/Dimensions/4A/2020seqreads/figures", scale = 1, width = 8, height = 3, units ="in", dpi = 300, limitsize = TRUE)

### combine fig ----
theme_set(theme_pubr())
combpcas1 <- ggarrange(funsap_precip, funsap_phylo, funpath_precip, funpath_precip, oomy_precip, oomy_precip, bac_precip, bac_phylo, rhiz_precip, rhiz_precip, amf_precip, amf_phylo, labels = c("A", "B","C", "", "D", "", "E","F","G", "","H","I"), ncol = 2, nrow = 6)
ggsave(combpcas1, filename="combpcas_aug21-2.jpg", scale = 1, width = 15, height = 20, units ="in", dpi = 300, limitsize = TRUE)

combpcas1 <- ggarrange(funpath_precip, funsap_phylo, funsap_precip, funsap_phylo, oomy_precip, oomy_phylo, bac_precip, bac_phylo, amf_precip, amf_phylo, labels = c("A", "","B","C","D","","E","F","G","H"), ncol = 2, nrow = 4)


# PCA 2018-2020 ----
setwd("C:/Users/haley/OneDrive - University of Kansas/2020data/2018-20_mergepcoa")

## fungi ----
# fungi merge for funguild
#setwd("C:/Users/haley/OneDrive - University of Kansas/2020data/2018-20_bioinform/final-otu-tables")
#otu <- read.table("otu_table-wtax.txt", header=TRUE)
#head(otu)
# funguild
#guild <- read.table("otu_table-wtax.taxa.guilds.txt", header=T)
#head(guild)
#funguilds <- merge(otu, guild, by="OTU_ID", all.x=TRUE)
#head(funguilds)
#write.table(funguilds, file="4AFUNGuild_2018-20.txt")

setwd("C:/Users/haley/OneDrive - University of Kansas/2020data/2018-20_mergepcoa")
# fung saps
metafun <- read.csv("DP2018-20fungsapmeta.csv", header=T)
head(metafun)
fungsap1820 <- read.table("otu_DP2018-20_fungsapredo.csv", sep=",", header=T)
head(fungsap1820)
rownames(fungsap1820) <- fungsap1820$OTU_ID
fungsap1820$OTU_ID <- NULL
fungsap1820 <- t(fungsap1820)
fungsap1820_rclr <- transform(fungsap1820, "rclr")
fungsap1820_ATCH <- vegdist(fungsap1820_rclr, method="robust.aitchison")

fungsap.pca <- prcomp(fungsap1820_ATCH, center=T) # phylo 1, 2 # precip 2, 3
metafun$pc1 <- fungsap.pca$x[,1]
metafun$pc2 <- fungsap.pca$x[,2]
metafun$pc3 <- fungsap.pca$x[,3]
metafun$PlantFamily <- metafun$PhyloFam
DP2018FS <- metafun[which(metafun$library=='DP2018'),]

fs1 <- ggplot(data=DP2018FS, aes(x=pc1, y=pc2, color=PlantFamily, shape=PlantFamily)) + geom_point() + stat_ellipse(geom="polygon", aes(fill = PlantFamily), alpha = 0.001, show.legend = T, level = 0.95) +labs(title="Fungal saprobe composition 2018",x="PC1", y = "PC2") + theme(plot.title = element_text(size=13), legend.position = "none")
fs1

DP2020FS <- metafun[(metafun$library=='DP2020-150perc') | (metafun$library=='DP2020-50perc'),]

fs2 <- fs1 + geom_point(data=DP2020FS, color="grey", aes(x=pc1, y=pc2))
fs2

fs3 <- ggplot(data=DP2020FS, aes(x=pc1, y=pc2, color=PlantFamily, shape=PlantFamily)) + geom_point() + stat_ellipse(geom="polygon", aes(fill = PlantFamily), alpha = 0.001, show.legend = T, level = 0.95) +labs(title="Fungal saprobe composition 2020",x="PC1", y = "PC2") + theme(plot.title = element_text(size=13), legend.position = "none")
fs4 <- fs3 + geom_point(data=DP2018FS, color="grey", aes(x=pc1, y=pc2))

DP2020FS50 <- metafun[(metafun$library=='DP2020-50perc'),]
DP2020FS150 <- metafun[(metafun$library=='DP2020-150perc'),]

fs5 <- ggplot(data=DP2020FS50, aes(x=pc1, y=pc2, color="red")) + geom_point() + stat_ellipse(geom="polygon", alpha = 0.001, show.legend = T, level = 0.95) +labs(title="Fungal saprobe composition 2020",x="PC1", y = "PC2") + theme(plot.title = element_text(size=13), legend.position = "none")
fs6 <- fs5 + geom_point(data=DP2020FS150, color="blue", aes(x=pc1, y=pc2)) + stat_ellipse(data=DP2020FS150, geom="polygon", alpha = 0.001, show.legend = T, level = 0.95, color="blue") + geom_point(data=DP2018FS, color="grey", aes(x=pc1, y=pc2))
fs6

# fung path
metafungpath <- read.csv("DP2018-20fungpathmeta.csv", header=T)
head(metafungpath)
fungpath1820 <- read.csv("otu_DP2018-20_fungpathredo.csv", sep=",", header=T)
head(fungpath1820)
rownames(fungpath1820) <- fungpath1820$OTU_ID
fungpath1820$OTU_ID <- NULL
fungpath1820 <- t(fungpath1820)
fungpath1820_2 <- fungpath1820[,-c(448:959)]
fungpath1820_rclr <- transform(fungpath1820_2, "rclr")
fungpath1820_ATCH <- vegdist(fungpath1820_rclr, method="robust.aitchison")

fungpath.pca <- prcomp(fungpath1820_ATCH, center=T) # phylo 1, 3 # precip 2, 4 # sp rich 2
metafun$pcp1 <- fungpath.pca$x[,1]
metafun$pcp2 <- fungpath.pca$x[,2]
metafun$pcp3 <- fungpath.pca$x[,3]
metafun$pcp4 <- fungpath.pca$x[,4]
DP2018FP <- metafun[which(metafun$library=='DP2018'),]

fp1 <- ggplot(data=DP2018FP, aes(x=pcp1, y=pcp2, color=PlantFamily, shape=PlantFamily)) + geom_point() + stat_ellipse(geom="polygon", aes(fill = PlantFamily), alpha = 0.001, show.legend = T, level = 0.95) +labs(title="Fungal pathogen composition 2018",x="PC1", y = "PC2") + theme(plot.title = element_text(size=13), legend.position = "none")
fp1

DP2020FP <- metafun[(metafun$library=='DP2020-150perc') | (metafun$library=='DP2020-50perc'),]

fp2 <- fp1 + geom_point(data=DP2020FP, color="grey", aes(x=pcp1, y=pcp2))
fp2

fp3 <- ggplot(data=DP2020FP, aes(x=pcp1, y=pcp2, color=PlantFamily, shape=PlantFamily)) + geom_point() + stat_ellipse(geom="polygon", aes(fill = PlantFamily), alpha = 0.001, show.legend = T, level = 0.95) +labs(title="Fungal pathogen composition 2020",x="PC1", y = "PC2") + theme(plot.title = element_text(size=13), legend.position = "none")
fp4 <- fp3 + geom_point(data=DP2020FP, color="grey", aes(x=pcp1, y=pcp2))
fp4

DP2020FP50 <- metafun[(metafun$library=='DP2020-50perc'),]
DP2020FP150 <- metafun[(metafun$library=='DP2020-150perc'),]

fp5 <- ggplot(data=DP2020FP50, aes(x=pcp1, y=pcp2, color="red")) + geom_point() + stat_ellipse(geom="polygon", alpha = 0.001, show.legend = T, level = 0.95) +labs(title="Fungal pathogen composition 2020",x="PC1", y = "PC2") + theme(plot.title = element_text(size=13), legend.position = "none")
fp6 <- fp5 + geom_point(data=DP2020FP150, color="blue", aes(x=pcp1, y=pcp2)) + stat_ellipse(data=DP2020FP150, geom="polygon", alpha = 0.001, show.legend = T, level = 0.95, color="blue") + geom_point(data=DP2018FP, color="grey", aes(x=pcp1, y=pcp2))
fp6

## oomy ----
metaoomy <- read.csv("oomymeta2018-20.csv", header=T)
head(metaoomy)
oomy1820 <- read.table("otu_oomy2018-20-nokoe.csv", sep=",", header=T)
head(oomy1820)
rownames(oomy1820) <- oomy1820$OTU_ID
oomy1820$OTU_ID <- NULL
oomy1820 <- t(oomy1820)
oomy1820_rclr <- transform(oomy1820, "rclr")
oomy1820_ATCH <- vegdist(oomy1820_rclr, method="robust.aitchison")

oomy.pca <- prcomp(oomy1820_ATCH, center=T) # precip 1, 2
metaoomy$pc1 <- oomy.pca$x[,1]
metaoomy$pc2 <- oomy.pca$x[,2]
metaoomy$PlantFamily <- metaoomy$PhyloFam
DP2018oomy <- metaoomy[which(metaoomy$library=='DP2018'),]

o1 <- ggplot(data=DP2018oomy, aes(x=pc1, y=pc2, color=PlantFamily, shape=PlantFamily)) + geom_point() + stat_ellipse(geom="polygon", aes(fill = library), alpha = 0.001, show.legend = T, level = 0.95) +labs(title="Oomycete composition 2018",x="PC1", y = "PC2") + theme(plot.title = element_text(size=13), legend.position = "none")
o1

DP2020oomy <- metaoomy[(metaoomy$library=='DP2020-150perc') | (metaoomy$library=='DP2020-50perc'),]

o2 <- o1 + geom_point(data=DP2020oomy, color="grey", aes(x=pc1, y=pc2))

o3 <- ggplot(data=DP2020oomy, aes(x=pc1, y=pc2, color=PlantFamily, shape=PlantFamily)) + geom_point() + stat_ellipse(geom="polygon", aes(fill = PlantFamily), alpha = 0.001, show.legend = T, level = 0.95) +labs(title="Oomycete composition 2020",x="PC1", y = "PC2") + theme(plot.title = element_text(size=13), legend.position = "none")
o4 <- o3 + geom_point(data=DP2018oomy, color="grey", aes(x=pc1, y=pc2))

DP2020oomy50 <- metaoomy[(metaoomy$library=='DP2020-50perc'),]
DP2020oomy150 <- metaoomy[(metaoomy$library=='DP2020-150perc'),]

o5 <- ggplot(data=DP2020oomy50, aes(x=pc1, y=pc2, color="red")) + geom_point() + stat_ellipse(geom="polygon", alpha = 0.001, show.legend = T, level = 0.95) +labs(title="Oomycete composition 2020",x="PC1", y = "PC2") + theme(plot.title = element_text(size=13), legend.position = "none")
o6 <- o5 + geom_point(data=DP2020oomy150, color="blue", aes(x=pc1, y=pc2)) + stat_ellipse(data=DP2020oomy150, geom="polygon", alpha = 0.001, show.legend = T, level = 0.95, color="blue") + geom_point(data=DP2018oomy, color="grey", aes(x=pc1, y=pc2))
o6

## bacteria ----
metabact <- read.csv("bactmeta2018-20.csv", header=T)
head(metabact)

## merge taxonomy and OTU tables
otus <- read.table("DP2018-20bact-grouped-nokoe.csv", sep=",", header=T)

# merge data frames by OTU_ID/Feature.ID
all <- merge(otus, tax, by.x = "OTU_ID", by.y = "Feature.ID")

## assign N-fixer status based on known genera
# split Taxon string into multiple columns
library(reshape2)
newColNames <- c("D0", "D1", "D2", "D3", "D4", "D5", "D6")
newCols <- colsplit(all$Taxon, ";", newColNames)
all_tax <- cbind(all, newCols)
all_tax$Taxon <- NULL
all_tax <- all_tax[,-239] # remove 'confidence' column

# create new 'nfixer' column; 0 = not rhizobia, 1 = rhizobia
all_tax$nfixer <- ifelse(all_tax$D5=="D_5__Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium"| 
                           all_tax$D5=="D_5__Bradyrhizobium" | all_tax$D5=="D_5__Ensifer" | 
                           all_tax$D5=="D_5__Burkholderia-Caballeronia-Paraburkholderia" | all_tax$D5=="D_5__Mesorhizobium", "1", "0")


## calculate relative abundance of N-fixers/sample
d <- colnames(all_tax[,sapply(all_tax, is.numeric)])
n<-length(d)
columns = c("rel_ab","richness","diversity") 
output = data.frame(matrix(nrow = n, ncol = length(columns)))
colnames(output) = columns

for(i in 1:n){
  s <- sum(all_tax[which(all_tax$nfixer=="1"), d[i]])
  df <- as.data.frame(s)
  output$richness[i] <- df[[1]]}
for(i in 1:n){
  s <- sum(all_tax[which(all_tax$nfixer=="1"), d[i]])/sum(all_tax[,d[i]])
  df <- as.data.frame(s)
  output$rel_ab[i] <- df[[1]]}
for(i in 1:n){
  s <- diversity(all_tax[which(all_tax$nfixer=="1"), d[i]])
  df <- as.data.frame(s)
  output$diversity[i] <- df[[1]]}

output$Sample <- d # now you have a dataframe with sample ID as one column and relative abundance of rhizobia reads as the other column

rhizobia1820 <- all_tax[which(all_tax$nfixer=="1"),]
rownames(rhizobia1820) <- rhizobia1820$OTU_ID
rhizobia1820$OTU_ID <- NULL
rhizobia1820 <- t(rhizobia1820)
rhizobia1820_rclr <- transform(x=rhizobia1820[-c(353:361),], "rclr")
str(rhizobia1820_rclr)
rhizobia1820_rclr <- mutate_all(rhizobia1820_rclr[,-1], function(x) as.integer(as.character(x)))
rhizobia1820_rclrno0 <- rhizobia1820_rclr[rowSums(rhizobia1820_rclr[])>0,]
rhizobia1820_ATCH <- vegdist(rhizobia1820_rclrno0, method="robust.aitchison")
showna <- as.data.frame(lapply(lapply(rhizobia1820_rclr, is.na), table))
min(rowSums(rhizobia1820_rclr)) # df1[rowSums(df1[])>0,]
rhiz.pca <- prcomp(rhizobia1820_ATCH, center=T) # phylo 4, 5 # precip 1, 4

rhizmeta1820 <- metabact[-c(2,60,192,238),]
rhizmeta1820$pc1 <- rhiz.pca$x[,1]
rhizmeta1820$pc2 <- rhiz.pca$x[,2]
rhizmeta1820$pc4 <- rhiz.pca$x[,4]
rhizmeta1820$pc5 <- rhiz.pca$x[,5]
summary(aov(rhiz.pca$scores[,5] ~ Block + PhyloFam*PlntDiv*library +  ANDGER + ELYCAN + BOUGRA + PANVIR + AMOCAN + DALCAN + DALPUR + DESILL + DESCAN + CHAFAS + LIAPYC + CORTIN + ECHPAL + EUPALT + SILINT + HELMOL, data=rhizmeta1820))
rhizmeta1820$PlantFamily <- rhizmeta1820$PhyloFam
DP2018rhiz <- rhizmeta1820[which(rhizmeta1820$library=='DP2018'),]

r1 <- ggplot(data=DP2018rhiz, aes(x=pc1, y=pc2, color=PlantFamily, shape=PlantFamily)) + geom_point() + stat_ellipse(geom="polygon", aes(fill = PlantFamily), alpha = 0.001, show.legend = T, level = 0.95) +labs(title="Rhizobia composition 2018",x="PC1", y = "PC2") + theme(plot.title = element_text(size=13), legend.position="none")

DP2020rhiz <- rhizmeta1820[(rhizmeta1820$library=='DP2020-150perc') | (rhizmeta1820$library=='DP2020-50perc'),]
r2 <- r1 + geom_point(data=DP2020rhiz, color="grey", aes(x=pc1, y=pc2))

r3 <- ggplot(data=DP2020rhiz, aes(x=pc1, y=pc2, color=PlantFamily, shape=PlantFamily)) + geom_point() + stat_ellipse(geom="polygon", aes(fill = PlantFamily), alpha = 0.001, show.legend = T, level = 0.95) +labs(title="Rhizobia composition 2020",x="PC1", y = "PC2") + theme(plot.title = element_text(size=13), legend.position="none")
r4 <- r3 + geom_point(data=DP2018rhiz, color="grey", aes(x=pc1, y=pc2))

DP2020rhiz50 <- rhizmeta1820[(rhizmeta1820$library=='DP2020-50perc'),]
DP2020rhiz150 <- rhizmeta1820[(rhizmeta1820$library=='DP2020-150perc'),]

r5 <- ggplot(data=DP2020rhiz50, aes(x=pc1, y=pc2, color="red")) + geom_point() + stat_ellipse(geom="polygon", alpha = 0.001, show.legend = T, level = 0.95) +labs(title="Rhizobia composition 2020",x="PC1", y = "PC2") + theme(plot.title = element_text(size=13), legend.position = "none")
r6 <- r5 + geom_point(data=DP2020rhiz150, color="blue", aes(x=pc1, y=pc2)) + stat_ellipse(data=DP2020rhiz150, geom="polygon", alpha = 0.001, show.legend = T, level = 0.95, color="blue") + geom_point(data=DP2018rhiz, color="grey", aes(x=pc1, y=pc2))
r6

bact1820 <- all_tax[which(all_tax$nfixer=="0"),]
head(bact1820)
rownames(bact1820) <- bact1820$OTU_ID
bact1820$OTU_ID <- NULL
bact1820 <- t(bact1820)
bact1820_rclr <- transform(x=bact1820[-c(353:361),], "rclr")
bact1820_rclr <- mutate_all(bact1820_rclr, function(x) as.numeric(as.character(x)))
bact1820_ATCH <- vegdist(bact1820_rclr[-1], method="robust.aitchison")

bact1820.pca <- prcomp(bact1820_ATCH, center=T) # phylo 4, 5 # precip 1, 4
metabact1820 <- metabact[-238,]
metabact1820$pc1 <- bact1820.pca$x[,1]
metabact1820$pc2 <- bact1820.pca$x[,2]
metabact1820$pc4 <- bact1820.pca$x[,4]
metabact1820$pc5 <- bact1820.pca$x[,5]
metabact1820$PlantFamily <- metabact1820$PhyloFam
DP2018bac <- metabact1820[which(metabact1820$library=='DP2018'),]

b1 <- ggplot(data=DP2018bac, aes(x=pc1, y=pc2, color=PlantFamily, shape=PlantFamily)) + geom_point() + stat_ellipse(geom="polygon", aes(fill = PlantFamily), alpha = 0.001, show.legend = T, level = 0.95) +labs(title="Bacteria composition 2018",x="PC1", y = "PC2") + theme(plot.title = element_text(size=13), legend.position="none")

DP2020bac <- metabact1820[(metabact1820$library=='DP2020-150perc') | (metabact1820$library=='DP2020-50perc'),]
b2 <- b1 + geom_point(data=DP2020bac, color="grey", aes(x=pc1, y=pc2))

b3 <- ggplot(data=DP2020bac, aes(x=pc1, y=pc2, color=PlantFamily, shape=PlantFamily)) + geom_point() + stat_ellipse(geom="polygon", aes(fill = PlantFamily), alpha = 0.001, show.legend = T, level = 0.95) +labs(title="Bacteria composition 2020",x="PC1", y = "PC2") + theme(plot.title = element_text(size=13), legend.position="none")
b4 <- b3 + geom_point(data=DP2018bac, color="grey", aes(x=pc1, y=pc2))

DP2020bac50 <- metabact1820[(metabact1820$library=='DP2020-50perc'),]
DP2020bac150 <- metabact1820[(metabact1820$library=='DP2020-150perc'),]

b5 <- ggplot(data=DP2020bac50, aes(x=pc1, y=pc2, color="red")) + geom_point() + stat_ellipse(geom="polygon", alpha = 0.001, show.legend = T, level = 0.95) +labs(title="Bacteria composition 2020",x="PC1", y = "PC2") + theme(plot.title = element_text(size=13), legend.position = "none")
b6 <- b5 + geom_point(data=DP2020bac150, color="blue", aes(x=pc1, y=pc2)) + stat_ellipse(data=DP2020bac150, geom="polygon", alpha = 0.001, show.legend = T, level = 0.95, color="blue") + geom_point(data=DP2018bac, color="grey", aes(x=pc1, y=pc2))
b6

## AMF ----
metaAM <- read.csv("AMFmeta2018-20.csv", header=T)
head(metaAM)
#metaAM[,2:320]
amf1820 <- read.table("Dp2018-20AMFotus.csv", sep=",", header=T)
head(amf1820)
rownames(amf1820) <- amf1820$OTU_ID
amf1820$OTU_ID <- NULL
amf1820 <- t(amf1820)
amf1820_rclr <- transform(amf1820, "rclr")
amf1820_ATCH <- vegdist(amf1820_rclr, method="robust.aitchison")

amf.pca <- prcomp(amf1820_ATCH, center=T) # phylo 4, 7 # precip 1, 2
metaAM$pc1 <- amf.pca$x[,1]
metaAM$pc2 <- amf.pca$x[,2]
metaAM$pc4 <- amf.pca$x[,4]
metaAM$pc7 <- amf.pca$x[,7]
metaAM$PlantFamily <- metaAM$PhyloFam
metaAM$PlantFamily <- revalue(metaAM$PlantFamily, c("AST"="Asteraceae", "FAB"="Fabaceae", "POA"="Poaceae", "MIX"="Mixture"))
DP2018amf <- metaAM[which(metaAM$library=='DP2018'),]
DP2020amf <- metaAM[(metaAM$library=='DP2020-150perc') | (metaAM$library=='DP2020-50perc'),]

a1 <- ggplot(data=DP2018amf, aes(x=pc1, y=pc2, color=PlantFamily, shape=PlantFamily)) + geom_point() + stat_ellipse(geom="polygon", aes(fill = PlantFamily), alpha = 0.001, show.legend = T, level = 0.95) +labs(title="AMF composition 2018",x="PC1", y = "PC2") + theme(plot.title = element_text(size=13), legend.position = "bottom")
a1
a2 <- a1 + geom_point(data=DP2020amf, color="grey", aes(x=pc1, y=pc2))
a2

a3 <- ggplot(data=DP2020amf, aes(x=pc1, y=pc2, color=PlantFamily, shape=PlantFamily)) + geom_point() + stat_ellipse(geom="polygon", aes(fill = PlantFamily), alpha = 0.001, show.legend = T, level = 0.95) +labs(title="AMF composition 2020",x="PC1", y = "PC2") + theme(plot.title = element_text(size=13), legend.position = "bottom")
a4 <- a3 + geom_point(data=DP2018amf, color="grey", aes(x=pc1, y=pc2))

DP2020amf50 <- metaAM[(metaAM$library=='DP2020-50perc'),]
DP2020amf150 <- metaAM[(metaAM$library=='DP2020-150perc'),]

a5 <- ggplot(data=DP2020amf50, aes(x=pc1, y=pc2, color="red")) + geom_point() + stat_ellipse(geom="polygon", alpha = 0.001, show.legend = T, level = 0.95) +labs(title="AMF composition 2020",x="PC1", y = "PC2") + theme(plot.title = element_text(size=13), legend.position = "bottom")
a6 <- a5 + geom_point(data=DP2020amf150, color="blue", aes(x=pc1, y=pc2)) + stat_ellipse(data=DP2020amf150, geom="polygon", alpha = 0.001, show.legend = F, level = 0.95, color="blue") + geom_point(data=DP2018amf, color="grey", aes(x=pc1, y=pc2))
a6

### combine fig ----
theme_set(theme_pubr())
combpcas_fin <- ggarrange(fs2, fs4, fs6, fp2, fp4, fp6, b2, b4, b6, r2, r4, r6, o2, o4, o6, a2, a4, a6, labels = c("A", "","","B","","","C","","","D", "", "", "E", "", "", "F", "", ""), ncol = 3, nrow = 6)
ggsave(combpcas_fin, filename="combpcas_jan21.jpg", scale = 1, width = 20, height = 26, units ="in", dpi = 300, limitsize = TRUE)