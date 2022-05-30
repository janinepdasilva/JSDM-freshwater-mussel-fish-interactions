# JSDM-freshwater-mussel-fish-interactions
#Data and code for the article "Happily ever after: co-occurrence patterns between freshwater mussels and their fish hosts uncovered by Joint Species Distribution Models"


################################################################################
#################### Joint Species Distribution Models #########################
#################### Freshwater mussels and fish hosts #########################
################################################################################


library(Hmsc)
library(dplyr)

setwd("D:/R_Folders/JSDMs_hmsc_HyRIV")


#Read data
data <- read.csv("douro_comm_env_HyRIVERS.csv", stringsAsFactors = TRUE, header = TRUE)

#Organize data
#Construct XData dataframe (environmental variables) that will be used as fixed effects
XData <- data.frame(dplyr::select(data, tmp_dc_cmx, ele_mt_cav, pre_mm_cyr, ORD_STRA, crp_pc_use, S_CaCO3_av))
head(XData)

#Construct the YData dataframe (species data)
YData <- data.frame(dplyr::select(data, Mmar, Aana, Plit, Udel, Aarc, Agal, Aolig, Ccal, Cpal, Lboc, Pdur, Salb, Scar, Stru))

#Analyse species richness and prevalence
library(ggplot2)
S = as.data.frame(rowSums(YData > 0))
max(S)
colnames(S) <- c("Species richness")
freq_S <- as.data.frame(table(S))
freq_S

P = as.data.frame(colMeans(YData > 0))

# x- and y-coordinates of sampling sites. We store these as the xy-matrix to be able to fit a spatial model

xy <- as.matrix(cbind(data$X, data$Y))

colnames(xy) <- c("x-coordinate", "y-coordinate")
rownames(xy)=data$HYRIV_ID
head(xy)
par(mfrow=c(1,1))
plot(xy, asp=1) # show the map (NB., equal aspect ratio in the map)


# The data that we will utilize for the Hmsc model consist of the community data matrix Y,
# the environmental data frame XData,
# the coordinates of the sampling locations xy,
# but could also include the species trait data frame TrData,
# and the phylogenetic tree PhyloTree.

# To define a spatial random effect at the level of the route, we need to include the river ID in the studyDesign

studyDesign <- data.frame(as.factor(data$HYRIV_ID))
colnames(studyDesign) <- c("HYRIV_ID")
rL <- HmscRandomLevel(sData = xy)

#IMPORTANT! Row names of xy correspond to the units of the studyDesign. 
#This is how HMSC understands how the units of the random effect (rows of xy) correspond to the sampling units

nChains = 3
nParallel = 3 # optional setting of nParallel
thin = 10 # Run 1, 10, 100 to achieve model convergence
samples = 1500*thin
transient = 500*thin

#Run 2 models: Model FULL includes both the environmental covariates and the spatial random effect of sampling site; 
#Model space includes the spatial random effect but no environmental covariates. No species traits and phylogenetic relationships.

XFormula = ~ tmp_dc_cmx + ele_mt_cav + pre_mm_cyr + ORD_STRA + crp_pc_use + S_CaCO3_av

m.FULL = Hmsc(Y = YData, XData = XData, XFormula = XFormula,
              distr = "probit",
              studyDesign = studyDesign,
              ranLevels = list(HYRIV_ID=rL))


m.SPACE = Hmsc(Y = YData, XData = XData, XFormula = ~1,
               distr = "probit",
               studyDesign = studyDesign,
               ranLevels = list(HYRIV_ID= rL))

m.FULL
m.SPACE

#Combine the three models

models <- list(m.FULL, m.SPACE)
for (i in 1:2) {
  models[[i]] = sampleMcmc(models[[i]], thin = thin, samples = samples, transient = transient, nChains = nChains)
  
}


#Examine MCMC convergence

mpost <- convertToCodaObject(models[[1]], spNamesNumbers = c(T,F),
                             covNamesNumbers = c(T,F))
psrf.beta <- gelman.diag(mpost$Beta, multivariate = FALSE)$psrf
tmp <- mpost$Omega[[1]]
z <- dim(tmp[[1]]) [2]

for (i in 1:length(tmp)) {
  tmp[[i]] = tmp[[i]]
}

psrf.omega <- gelman.diag(tmp, multivariate = FALSE)$psrf
par(mfrow=c(1,2))

hist(psrf.beta, xlab = expression(paste("Potential reduction factor (parameter ", beta, ")")))
hist(psrf.omega, xlab = expression(paste("Potential reduction factor (parameter ", Omega, ")")))



#########################################################################################

#EVALUATING MODEL FIT AND COMPARING MODELS


partition = createPartition(models[[1]], nfolds = 2, column = "HYRIV_ID")
MF = list()
MFCV = list()
for (i in 1:2){
preds = computePredictedValues(models[[i]])
MF[[i]] = evaluateModelFit(hM = models[[i]], predY = preds)
preds = computePredictedValues(models[[i]],
 partition = partition)
 MFCV[[i]] = evaluateModelFit(hM = models[[i]], predY = preds)
}


WAIC = unlist(lapply(models, FUN = computeWAIC))
WAIC = as.data.frame(WAIC)

#########################################################################################


#Exploring Parameter estimates

VP = list()
for (i in 1:2){
  VP[[i]] = computeVariancePartitioning(models[[i]])
}


#Plot Variance partitioning with Hmsc

plotVariancePartitioning(models[[1]], VP[[1]], args.legend = list(x = "topright", bty = "n", inset=c(-0.15, 0)))


#Plot heatmap

postBeta <- getPostEstimate(models[[1]], parName = "Beta")
plotBeta(models[[1]], post=postBeta, supportLevel = 0.95, param = "Sign", plotTree = FALSE, spNamesNumbers = c(TRUE, FALSE), 
         covNamesNumbers = c(TRUE, FALSE), colors = colorRampPalette(c("grey","white","black")), colorLevels = 3)


# species associations revealed by the random effects with the corrplot function
par(mfrow = c(1,2))
#Model FULL - residual associations that account for the responses of species to environmental conditions
library(corrplot)
OmegaCor = computeAssociations(models[[1]])
supportLevel = 0.85
toPlot = ((OmegaCor[[1]]$support>supportLevel)
          + (OmegaCor[[1]]$support<(1-supportLevel))>0)*OmegaCor[[1]]$mean




corrplot(toPlot, method="color", type = "lower", tl.srt = 45, col=c("grey", "white", "black"), title = "Residual co-occurrences", 
         tl.cex = 0.8, tl.col = "black", font = 3, mar=c(0,0,1,0))


#Model SPATIAL raw associations

OmegaCor = computeAssociations(models[[2]])
supportLevel = 0.85
toPlot = ((OmegaCor[[1]]$support>supportLevel)
          + (OmegaCor[[1]]$support<(1-supportLevel))>0)*OmegaCor[[1]]$mean

corrplot(toPlot, method="color", type = "lower", tl.srt = 45, col=c("grey", "white", "black"), title = "Raw co-occurrences", 
         tl.cex = 0.8, tl.col = "black", font = 3, mar=c(0,0,1,0))
