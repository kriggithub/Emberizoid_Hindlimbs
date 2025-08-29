#setwd()

###############################################################
### 1. Read in full dataset and explore morphological variation
all <- read.csv("AllBlackbirdData.csv")

hist(all$L.Tarsus) 
hist(all$W.Tarsus) # 16116 P. decumanus largest sizes for multiple measurements
hist(all$D.Tarsus)
hist(all$L.Midtoe) # 11930 C. cela 0 midtoe length
hist(all$W.Midtoe.P1)
hist(all$L.Midtoe.Claw)
hist(all$L.Hallux)
hist(all$W.Hallux.P1)
hist(all$L.Hallux.Claw)

bad <- c(11930)

all <- all[all$SpecimenID%in%bad==F,] # remove specimens w/ incorrect measurements for now

rownames(all)<-all$SpecimenID

############################################
### 2. Perform PCA on hind limb measurements
morpho <- log(all[,5:13]) # extract and ln-transform measurements

#View(morpho)

rownames(morpho) <- all$SpecimenID
morpho <- morpho[complete.cases(morpho),] # remove specimens w/ any missing measurements (removed 1)

all <- all[complete.cases(all[,5:13]),]

pca <- prcomp(morpho) # perform PCA on covariance matrix

pct.var <- pca$sdev^2 / sum(pca$sdev^2) # extract % variation for each axis
# PC1 = 77.3%, PC2 = 8.9%, PC3 = 5.6%

eigvec <- pca$rotation # extract eigenvector matrix
#view(eigvec)
plot(eigvec[,1:2],pch=19,xlab="PC1 Loading",ylab="PC2 Loading",cex=1.5,
     xlim=c(-0.55,0.55),ylim=c(-0.6,0.8))
abline(h=0,lty=2,lwd=0.5)
abline(v=0,lty=2,lwd=0.5)
text(eigvec[,1:2],rownames(eigvec),pos=2)

plot(eigvec[,3:2],pch=19,xlab="PC3 Loading",ylab="PC2 Loading",cex=1.5,
     xlim=c(-0.75,0.6),ylim=c(-0.6,0.8))
abline(h=0,lty=2,lwd=0.5)
abline(v=0,lty=2,lwd=0.5)
text(eigvec[,3:2],rownames(eigvec),pos=2)

#PC1 - all positively loaded; a size axis
#PC2 - widths - negative, lengths - positive; a width vs length axis
#      width of hallux, tarsus, midtoe and depth of midtoe negative;
#      length of midtoe claw, tarsus, midtoe, hallux claw, hallux positive
#PC3 - another length of leg vs width and depth of tarsus/toes
#      width of midtoe, tarsus, length of hallux, depth of tarsus negative
#      length of tarsus, midtoe, hallux claw, width of hallux, and length of midtoe claw positive

#write.csv(rbind(eigvec,pct.var),file="PCA_Eigenvectorblackbird.csv",row.names = T)

scores <- data.frame(pca$x) # extract PC scores
plot(scores$PC1,scores$PC2)
plot(scores$PC3,scores$PC2)

all.pc <- cbind(all,scores[,1:3]) # combine PC1-3 with all other data

#write.csv(all.pc,file="AllBlackbirdData_PCA.csv",row.names = F)


#################################################################################
### 3. Read in phylogenetic tree and ensure that dataset and tree match perfectly
require(ape)
require(phytools)

tree <- read.nexus("~/RegisRProjects/ImfeldResearch/Blackbird_hindlimbs/Icteridae_SpeciesTree_MCC.tre") # read in phylogeny

missing <- tree$tip.label[tree$tip.label%in%all.pc$Species==F] # ID species in tree but not in dataset

tree <- drop.tip(tree,missing) # drop missing species
all.pc <- all.pc[all.pc$Species%in%tree$tip.label==T,] # drop specimens not in tree (8)


#############################################################
### 4. Aggregate species-level data for phylogenetic analyses (for blackbirds have no scratching, but locomotive behavior)
require(dplyr)

all.sp.morph <- all.pc %>% group_by(Species) %>% 
  summarize_at(vars("Mass..g.","L.Tarsus","W.Tarsus","D.Tarsus","L.Midtoe",
                    "W.Midtoe.P1","L.Midtoe.Claw","L.Hallux","W.Hallux.P1",
                    "L.Hallux.Claw","PC1","PC2","PC3"),mean) # aggregate means for quantitative variables
all.sp.eco <- all.pc %>% group_by(Species) %>% 
  summarize_at(vars("Habitat","Foraging.Score","Locomotion","Singing.Score",
                    "Nesting.Score", "Parasite"),unique) # extract categorical values for category variables

all.sp <- merge(all.sp.morph,all.sp.eco,by="Species") # merge trait and ecological data together

all.sp$Habitat <- as.factor(all.sp$Habitat)
all.sp$Foraging.Score <- as.factor(all.sp$Foraging.Score)
all.sp$Locomotion <- as.factor(all.sp$Locomotion)
all.sp$Singing.Score <- as.factor(all.sp$Singing.Score)
all.sp$Nesting.Score <- as.factor(all.sp$Nesting.Score)
all.sp$Parasite <- as.factor(all.sp$Parasite)


rownames(all.sp) <- all.sp$Species


###########################################################
### 5. Calculate phylogenetic signal for measurements & PCs
phylo.sig <- function(column,dataset,tree){
  this.trait <- dataset[,column] # extract the trait 
  names(this.trait) <- dataset$Species # assign species as names
  lambda <- phylosig(tree,this.trait,method="lambda",test=T) # calculate Pagel's lambda
  k <- phylosig(tree,this.trait,method="K",test=T) # calculate Blomberg's K
  results <- c(colnames(dataset)[column],lambda$lambda,lambda$P,k$K,k$P) # assemble results
  names(results) <- c("Trait","Lambda","P.Lambda","K","P.K") # assign names to results
  results # print results
}


ps.results <- matrix(nrow=12,ncol=5) # set up phylogenetic signal data matrix
colnames(ps.results) <- c("Trait","Lambda","P.Lambda","K","P.K")

for(i in 3:14){ # calculate phylogenetic signal for all continuous variables
  ps.results[i-2,] <- phylo.sig(i,all.sp,tree)
}

ps.results<-as.data.frame(ps.results) # reorganize as data frame
ps.results[,2:5] <- apply(ps.results[2:5],2,as.numeric) # make continuous data numeric

#write.csv(ps.results,file="Phylogenetic_Signal_Results_blackbirds.csv",row.names=F)


#################################
### 6. Fit pANOVA and PGLM models
require(phylolm)
require(geomorph)
require(caper)

PC1 <- all.sp$PC1
PC2 <- all.sp$PC2
PC3 <- all.sp$PC3
names(PC1) <- all.sp$Species
names(PC2) <- all.sp$Species
names(PC3) <- all.sp$Species


# Habitat (removed n=1 shrubland)
all.sp.noshrub <- subset(all.sp, Habitat != "Shrubland")
all.sp.noshrub$Habitat <- droplevels(all.sp.noshrub$Habitat)
all.sp.noshrub$Habitat <- as.factor(all.sp.noshrub$Habitat)
rownames(all.sp.noshrub) <- all.sp.noshrub$Species

PC1noshrub <- all.sp.noshrub$PC1
PC2noshrub <- all.sp.noshrub$PC2
PC3noshrub <- all.sp.noshrub$PC3
names(PC1noshrub) <- all.sp.noshrub$Species
names(PC2noshrub) <- all.sp.noshrub$Species
names(PC3noshrub) <- all.sp.noshrub$Species

Habitatnoshrub <- all.sp.noshrub$Habitat
names(Habitatnoshrub) <- all.sp.noshrub$Species
treenoshrub <- drop.tip(tree, "Icterus_parisorum")


hab1noshrub <- phylANOVA(treenoshrub, Habitatnoshrub, PC1noshrub)
#NONE significant
hab2noshrub <- phylANOVA(treenoshrub, Habitatnoshrub, PC2noshrub)
#Grassland vs. Forest, Wetland, Woodland significant
hab3noshrub <- phylANOVA(treenoshrub, Habitatnoshrub, PC3noshrub)
#NONE Significant

hab.cols.fig <- c("forestgreen","grey66","palegreen2", "wheat", "dodgerblue", "burlywood4")
boxplot(PC2noshrub~Habitatnoshrub,all.sp.noshrub, col = hab.cols.fig, ylab="PC2 Score",cex.lab=1.25)




# Foraging
Foraging <- as.numeric(all.sp$Foraging.Score)
names(Foraging) <- all.sp$Species
f.cd <- comparative.data(tree,all.sp[,c(1,3:14,16)],names.col="Species")

for1 <- pgls(PC1~Foraging.Score,f.cd)
summary(for1)
# Not significant
for2 <- pgls(PC2~Foraging.Score,f.cd)
summary(for2)
# Not significant
for3 <- pgls(PC3~Foraging.Score,f.cd)
summary(for3)
# Borderline significant
plot(PC3~Foraging.Score,all.sp,pch=19)



# Locomotion
Locomotion <- all.sp$Locomotion
names(Locomotion) <- all.sp$Species

loco1 <- phylANOVA(tree,Locomotion,PC1)
#NONE significant
loco2 <- phylANOVA(tree,Locomotion,PC2)
#NONE significant
loco3 <- phylANOVA(tree,Locomotion,PC3)
#NONE significant




# Singing
s.cd <- comparative.data(tree,all.sp[,c(1,3:14,18)],names.col="Species")
sing1 <- pgls(PC1~Singing.Score,s.cd)
summary(sing1)
# Not significant
sing2 <- pgls(PC2~Singing.Score,s.cd)
summary(sing2)
# Not significant
sing3 <- pgls(PC3~Singing.Score,s.cd)
summary(sing3)
# Not significant




# Nesting
Nesting <- all.sp$Nesting.Score
names(Nesting) <- all.sp$Species
n.cd <- comparative.data(tree,all.sp[,c(1,3:14,19)],names.col="Species")

nest1 <- pgls(PC1~Nesting.Score,n.cd)
summary(nest1)
# Not significant
nest2 <- pgls(PC2~Nesting.Score,n.cd)
summary(nest2)
# Not significant
nest3 <- pgls(PC3~Nesting.Score,n.cd)
summary(nest3)
# Not significant




##################################
### 7. Fit non-phylogenetic models
# Habitat
hab.mod1 <- lm(PC1~Habitat,all.pc)
summary(hab.mod1)
TukeyHSD(aov(PC1~Habitat,all.pc))
#Woodland vs Generalist, Woodland vs. Grassland, Woodland vs. Wetland significant

hab.mod2 <- lm(PC2~Habitat,all.pc)
summary(hab.mod2)
TukeyHSD(aov(PC2~Habitat,all.pc))
#Generalist vs Forest, Wetland vs. Forest, Woodland vs Generalist, all grassland comparisons

hab.mod3 <- lm(PC3~Habitat,all.pc)
summary(hab.mod3)
TukeyHSD(aov(PC3~Habitat,all.pc))
#ALL wetland comparisons significant



# Foraging
for.mod1 <- lm(PC1~Foraging.Score,all.pc)
summary(for.mod1)
#Significant

for.mod2 <- lm(PC2~Foraging.Score,all.pc)
summary(for.mod2)
#Significant

for.mod3 <- lm(PC3~Foraging.Score,all.pc)
summary(for.mod3)
#Not significant

plot(PC1~Foraging.Score,all.pc)
plot(PC2~Foraging.Score,all.pc)




# Locomotion
loco.mod1 <- lm(PC1~Locomotion,all.pc)
summary(loco.mod1)
TukeyHSD(aov(PC1~Locomotion,all.pc))
#Significant

loco.mod2 <- lm(PC2~Locomotion,all.pc)
summary(loco.mod2)
TukeyHSD(aov(PC2~Locomotion,all.pc))
#Significant

loco.mod3 <- lm(PC3~Locomotion,all.pc)
summary(loco.mod3)
TukeyHSD(aov(PC3~Locomotion,all.pc))
#Not Significant




# Singing
sing.mod1 <- lm(PC1~Singing.Score,all.pc)
summary(sing.mod1)
#Not significant

sing.mod2 <- lm(PC2~Singing.Score,all.pc)
summary(sing.mod2)
#Significant

sing.mod3 <- lm(PC3~Singing.Score,all.pc)
summary(sing.mod3)
# Not Significant
boxplot(PC2~Singing.Score,all.pc)




# Nesting
nest.mod1 <- lm(PC1~Nesting.Score,all.pc)
summary(nest.mod1)
#Significant

nest.mod2 <- lm(PC2~Nesting.Score,all.pc)
summary(nest.mod2)
#Significant

nest.mod3 <- lm(PC3~Nesting.Score,all.pc)
summary(nest.mod3)
#Not significant
boxplot(PC1~Nesting.Score,all.pc)
boxplot(PC2~Nesting.Score,all.pc)




##################################
### 8. MANOVA and phyMANOVA (categorical variables)
require(geiger)
manova.morph <-as.matrix(all.sp.morph[,3:11]) # create named matrix of response variables
rownames(manova.morph) <- all.sp.morph$Species # assign species names to row names
manova.morph <- log(manova.morph)


#Habitat
Habitat <- all.sp$Habitat # create named vector for predictor
names(Habitat) <- all.sp$Species
### Pull out p-value from phyMANOVA simulations, and compute average across 10 iterations
hab.man.pvaluesafter <- numeric(10) #initialize p-value vector for after phylogeny
for (i in 1:10){
  hab.man <- aov.phylo(manova.morph ~ Habitat, phy = tree, nsim=1000) # run MANOVA with phylMANOVA simulations
  
  hab.man.pvaluesafter[i] <- attributes(hab.man)$summary$`Pr(>F) given phy`[1] # yoink p value to store
}
mean(hab.man.pvalues) # 0.2444555
sd(hab.man.pvalues) # 0.01649646
#Significant until given phylogeny




#Locomotion
Locomotion <- all.sp$Locomotion
names(Locomotion) <- all.sp$Species

loco.man.pvalues <- numeric(10) #initialize p-value vector
for (i in 1:10){
  loco.man <- aov.phylo(manova.morph ~ Locomotion, phy = tree, nsim=1000)
  
  loco.man.pvalues[i] <- attributes(loco.man)$summary$`Pr(>F) given phy`[1] 
}
mean(loco.man.pvalues) # 0.798002
sd(loco.man.pvalues) # 0.01345915
#Significant until given phylogeny




#Parasite
Parasite <- all.sp$Parasite
names(Parasite) <- all.sp$Species
para.man.pvalues <- numeric(10) #initialize p-value vector
for (i in 1:10){
  para.man <- aov.phylo(manova.morph ~ Parasite, phy = tree, nsim=1000)
  
  para.man.pvalues[i] <- attributes(para.man)$summary$`Pr(>F) given phy`[1] 
}
mean(para.man.pvalues) # 0.9964036
sd(para.man.pvalues) # 0.002009071
#Not Significant





##################################
### 9. Multivariate regression with and without accounting for phylogeny (continuous predictors)
### Multivariate regressions for continuous predictors (Foraging, Singing, Nesting):
require(car) # car contains Anova() command to parse out multivariate regressions

# Foraging
Foraging <- as.numeric(all.sp$Foraging.Score)
names(Foraging) <- all.sp$Species

for.mult <- lm(manova.morph~Foraging) # fit the multivariate regression
summary(Anova(for.mult)) # examine summary of multivariate regression
# Significant 

require(mvMORPH)
Foraging <- Foraging[!is.na(Foraging)] # need to remove NAs from this object for mvMORPH
for.morph <- manova.morph[names(Foraging),] # make morpho match ecological data
for.tree <- drop.tip(tree,tree$tip.label[tree$tip.label%in%names(Foraging)==F]) # same for tree
mult.list <- list(For.Morpho = for.morph[for.tree$tip.label,1:9], 
                  Foraging=Foraging[for.tree$tip.label]) # assemble list with morpho & ecological data sorted to match order of tip labels in tree
for.pmult <- mvgls(For.Morpho~Foraging, data = mult.list, tree = for.tree, model="BM") # fit phylogenetic multivariate GLS
for.pmult.null <- mvgls(For.Morpho~1, data = mult.list, tree = for.tree, model="BM") # fit null model for comparison
GIC(for.pmult)$GIC - GIC(for.pmult.null)$GIC
#58.0852


# Singing
Singing <- as.numeric(all.sp$Singing)
names(Singing) <- all.sp$Species

sing.mult <- lm(manova.morph~Singing) # fit the multivariate regression
summary(Anova(sing.mult)) # examine summary of multivariate regression
# Not Significant

Singing <- Singing[!is.na(Singing)] # need to remove NAs from this object for mvMORPH
sing.morph <- manova.morph[names(Singing),] # make morpho match ecological data
sing.tree <- drop.tip(tree,tree$tip.label[tree$tip.label%in%names(Singing)==F]) # same for tree
mult.list <- list(Sing.Morpho = sing.morph[sing.tree$tip.label,1:9], 
                  Singing=Singing[sing.tree$tip.label]) # assemble list with morpho & ecological data sorted to match order of tip labels in tree
sing.pmult <- mvgls(Sing.Morpho~Singing, data = mult.list, tree = sing.tree, model="BM") # fit phylogenetic multivariate GLS
sing.pmult.null <- mvgls(Sing.Morpho~1, data = mult.list, tree = sing.tree, model="BM") # fit null model for comparison
GIC(sing.pmult)$GIC - GIC(sing.pmult.null)$GIC
# 50.84088




# Nesting
Nesting <- as.numeric(all.sp$Nesting.Score)
names(Nesting) <- all.sp$Species

nest.mult <- lm(manova.morph~Nesting) # fit the multivariate regression
summary(Anova(nest.mult)) # examine summary of multivariate regression
#Significant 

Nesting <- Nesting[!is.na(Nesting)] # need to remove NAs from this object for mvMORPH
nest.morph <- manova.morph[names(Nesting),] # make morpho match ecological data
nest.tree <- drop.tip(tree,tree$tip.label[tree$tip.label%in%names(Nesting)==F]) # same for tree
mult.list <- list(nest.Morpho = nest.morph[nest.tree$tip.label,1:9], 
                  Nesting=Nesting[nest.tree$tip.label]) # assemble list with morpho & ecological data sorted to match order of tip labels in tree
nest.pmult <- mvgls(nest.Morpho~Nesting, data = mult.list, tree = nest.tree, model="BM") # fit phylogenetic multivariate GLS
nest.pmult.null <- mvgls(nest.Morpho~1, data = mult.list, tree = nest.tree, model="BM") # fit null model for comparison
GIC(nest.pmult)$GIC - GIC(nest.pmult.null)$GIC
# 64.52815


################################
### 9. Ancestral character estimation

# Habitat
hab <- as.factor(all.sp$Habitat) # pull out habitat as a factor variable
names(hab) <- all.sp$Species # assign names to it
hab <- na.omit(hab) # remove any species with missing habitat data
hab.tree <- drop.tip(tree,tree$tip.label[tree$tip.label%in%names(hab)==F]) # then prune tree to have them match

hab.simmap <- make.simmap(hab.tree,hab,model="ER",nsim=1000, pi = "estimated") # perform ACE with equal rates model
hab.simmap2 <- make.simmap(hab.tree,hab,model="ARD",nsim=1000, pi="estimated") # same, but with all rates different model
hab.simmap3 <- make.simmap(hab.tree,hab,model="SYM",nsim=1000, pi="estimated") # lastly, symmetrical rates model
AIC(hab.simmap[[1]], hab.simmap2[[1]], hab.simmap3[[1]]) # compare AIC scores for 3 fitted models
# ER model best supported, dAIC = 5.443

hab.cols <- c("forestgreen","grey66","palegreen2", "wheat", "dodgerblue", "burlywood4")
# assign colors to habitats - you may need to change this because you've got more habitats for the blackbirds
names(hab.cols) <- c("Forest", "Generalist", "Grassland","Shrubland","Wetland", "Woodland")

# Be sure to change hab.simmap3 below to whatever your best model is
hab.ace <- summary(hab.simmap,plot=F) # generate consensus ACE from the 1000 simulations of the best model
pdf(file="Habitat_Tree_Fig.pdf",height=11,width=8.5,useDingbats = F) # set up PDF to save image
plot(hab.ace,cex=c(0.5,0.3),colors=hab.cols) # plot the ACE on the phylogeny
legend(x = 12.5, y = 18, legend = names(hab.cols), fill = hab.cols, cex = 0.8, bty = "n")
dev.off() #export and save the PDF


# Locomotion
loco <- as.factor(all.sp$Locomotion) # pull out habitat as a factor variable
names(loco) <- all.sp$Species # assign names to it
loco <- na.omit(loco) # remove any species with missing habitat data
loco.tree <- drop.tip(tree,tree$tip.label[tree$tip.label%in%names(loco)==F]) # then prune tree to have them match

loco.simmap <- make.simmap(loco.tree,loco,model="ER",nsim=1000, pi = "estimated") # perform ACE with equal rates model
loco.simmap2 <- make.simmap(loco.tree,loco,model="ARD",nsim=1000, pi="estimated") # same, but with all rates different model
loco.simmap3 <- make.simmap(loco.tree,loco,model="SYM",nsim=1000, pi="estimated") # lastly, symmetrical rates model
AIC(loco.simmap[[1]], loco.simmap2[[1]], loco.simmap3[[1]]) # compare AIC scores for 3 fitted models
# ER & SYM model best supported, dAIC = 0

loco.cols.sim<-c("gray30","gray75")
# assign colors to habitats - you may need to change this because you've got more habitats for the blackbirds
names(loco.cols.sim) <- c("Hop", "Walk")

# Be sure to change hab.simmap3 below to whatever your best model is
loco.ace <- summary(loco.simmap,plot=F) # generate consensus ACE from the 1000 simulations of the best model
pdf(file="Locomotion_Tree_Fig.pdf",height=11,width=8.5,useDingbats = F) # set up PDF to save image
plot(loco.ace,cex=c(0.5,0.3),colors=loco.cols.sim) # plot the ACE on the phylogeny
legend(x = 12.5, y = 18, legend = names(loco.cols.sim), fill = loco.cols.sim, cex = 0.8, bty = "n")
dev.off() #export and save the PDF










################################
### 10. Generate relevant figures (One Size Graph Using PC1 violin plots, one shape graph using PC2 and PC3)
require(vioplot)
require(ggplot2)
require(phylolm)
require(geomorph)
require(caper)
require(plyr)


find_hull <- function(df) df[chull(df$PC2, df$PC3), ]

# Habitat
hab.cols.fig <- c("forestgreen","grey66","palegreen2", "wheat", "dodgerblue", "burlywood4")
pdf(file="HabitatSize.pdf",width=8,height=8,useDingbats = F)
vioplot(PC1~Habitat,all.pc,col=hab.cols.fig,border=F,lineCol=NA,
        rectCol=NA,plotCentre="line",
        ylab="Size Score (PC1)",cex.lab=1.5)
legend("topright",legend=c("p = 0.917"),bty='n')
dev.off()

hullshab <- ddply(all.sp, "Habitat", find_hull)


fig1hab <- ggplot(data=all.sp, aes(PC2, PC3, color = Habitat, fill = Habitat)) + 
  geom_point() +
  scale_color_manual(values = hab.cols.fig) + 
  scale_fill_manual(values = hab.cols.fig) + 
  labs(x = "PC2", y = "PC3", color = "Habitat", fill = "Habitat", title = "Blackbird Feet Shape (Habitat)") +
  geom_polygon(data=hullshab, alpha=.4) +
  annotate(geom = "text", x=-0.17, y = -0.17, label = "Shrubland", size = 3) +
  annotate(geom = "text", x = 0.4, y = 0.32, label = "PC2: p = 0.010", size = 4)+
  annotate(geom = "text", x = 0.4, y = 0.3, label = "PC3: p = 0.632", size = 4)+
  theme_minimal()


ggsave(filename = "HabitatShape.pdf", plot = fig1hab, width = 8, height = 8)


# Foraging
scale.cols<-c("black","gray30","gray50","gray75","gray95")
pdf(file="ForagingSize.pdf",width=8,height=8,useDingbats = F)
vioplot(PC1~Foraging.Score,all.pc,col=scale.cols,border=F,lineCol=NA,
        rectCol=NA,plotCentre="line",
        ylab="Size Score (PC1)",xlab="Foraging Score",cex.lab=1.5)
legend("topright",legend=c("p = 0.118"),bty='n')
dev.off()

all.sp.for <- subset(all.sp, !is.na(Foraging.Score))
hullsfor <- ddply(all.sp.for, "Foraging.Score", find_hull)


fig1for <- ggplot(data=all.sp.for, aes(PC2, PC3, color = Foraging.Score, fill = Foraging.Score)) + 
  geom_point() +
  scale_color_manual(values = scale.cols) + 
  scale_fill_manual(values = scale.cols) + 
  labs(x = "PC2", y = "PC3", color = "Foraging Score", fill = "Foraging Score", title = "Blackbird Feet Shape (Foraging)") +
  geom_polygon(data=hullsfor, alpha=.4) +
  annotate(geom = "text", x = 0.4, y = 0.32, label = "PC2: p = 0.230", size = 4) +
  annotate(geom = "text", x = 0.4, y = 0.3, label = "PC3: p = 0.040", size = 4) +
  theme_minimal()


ggsave(filename = "ForagingShape.pdf", plot = fig1for, width = 8, height = 8)



# Locomotion
loco.cols<-c("gray30","gray75")
all.sp.loco <- subset(all.sp, !is.na(Locomotion))
pdf(file="LocomotionSize.pdf",width=8,height=8,useDingbats = F)
vioplot(PC1~Locomotion,all.pc,col=loco.cols,border=F,lineCol=NA,
        rectCol=NA,plotCentre="line",
        ylab="Size Score (PC1)",xlab="Locomotion",cex.lab=1.5)
legend("topright",legend=c("p = 0.185"),bty='n')
dev.off()


hullsloco <- ddply(all.sp.loco, "Locomotion", find_hull)


fig1loco <- ggplot(data=all.sp.loco, aes(PC2, PC3, color = Locomotion, fill = Locomotion)) + 
  geom_point() +
  scale_color_manual(values = loco.cols) + 
  scale_fill_manual(values = loco.cols) + 
  labs(x = "PC2", y = "PC3", color = "Locomotion", fill = "Locomotion", title = "Blackbird Feet Shape (Locomotion)") +
  geom_polygon(data=hullsloco, alpha=.4) +
  annotate(geom = "text", x = 0.4, y = 0.32, label = "PC2: p = 0.302", size = 4) +
  annotate(geom = "text", x = 0.4, y = 0.3, label = "PC3: p = 0.613", size = 4) +
  theme_minimal()


ggsave(filename = "LocomotionShape.pdf", plot = fig1loco, width = 8, height = 8)









# Singing
pdf(file="SingingSize.pdf",width=8,height=8,useDingbats = F)
sing.cols <- c("black","gray50","gray75")
vioplot(PC1~Singing.Score,all.pc,col=sing.cols,border=F,lineCol=NA,
        rectCol=NA,plotCentre="line",
        ylab="Size Score (PC1)",xlab="Singing Score",cex.lab=1.5)
legend("topright",legend=c("p = 0.694"),bty='n')
dev.off()

all.sp.sing <- subset(all.sp, !is.na(Singing.Score))
hullssing <- ddply(all.sp.sing, "Singing.Score", find_hull)


fig1sing <- ggplot(data=all.sp.sing, aes(PC2, PC3, color = Singing.Score, fill = Singing.Score)) + 
  geom_point() +
  scale_color_manual(values = sing.cols) + 
  scale_fill_manual(values = sing.cols) + 
  labs(x = "PC2", y = "PC3", color = "Singing Score", fill = "Singing Score", title = "Blackbird Feet Shape (Singing)") +
  geom_polygon(data=hullssing, alpha=.4) +
  annotate(geom = "text", x = 0.4, y = 0.32, label = "PC2: p = 0.131", size = 4) +
  annotate(geom = "text", x = 0.4, y = 0.3, label = "PC3: p = 0.423", size = 4) +
  theme_minimal()


ggsave(filename = "SingingShape.pdf", plot = fig1sing, width = 8, height = 8)







# Nesting
scale.cols2<-c("black","gray30","gray50","gray75","gray85")
pdf(file="NestingSize.pdf",height=8.5,width=8.5,useDingbats = F)
vioplot(PC1~Nesting.Score,all.pc,col=scale.cols2,border=F,lineCol=NA,
        rectCol=NA,plotCentre="line",
        ylab="Size Score (PC1)",xlab="Nesting Score",cex.lab=1.5)
legend("topright",legend=c("p = 0.560"),bty='n')
dev.off()


all.sp.nest <- subset(all.sp, !is.na(Nesting.Score))
hullsnest <- ddply(all.sp.nest, "Nesting.Score", find_hull)


fig1nest <- ggplot(data=all.sp.nest, aes(PC2, PC3, color = Nesting.Score, fill = Nesting.Score)) + 
  geom_point() +
  scale_color_manual(values = scale.cols2) + 
  scale_fill_manual(values = scale.cols2) + 
  labs(x = "PC2", y = "PC3", color = "Nesting Score", fill = "Nesting Score", title = "Blackbird Feet Shape (Nesting)") +
  geom_polygon(data=hullsnest, alpha=.4) +
  annotate(geom = "text", x = 0.4, y = 0.32, label = "PC2: p = 0.124", size = 4) +
  annotate(geom = "text", x = 0.4, y = 0.3, label = "PC3: p = 0.567", size = 4) +
  theme_minimal()


ggsave(filename = "NestingShape.pdf", plot = fig1nest, width = 8, height = 8)









