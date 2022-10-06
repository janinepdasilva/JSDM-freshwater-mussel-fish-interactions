################################################################################
#################### Joint Species Distribution Models #########################
#################### Freshwater mussels and fish hosts #########################
################################################################################
rm(list = ls())
library(Hmsc)
library(corrplot)
library(tidyverse)
#setwd("D:/R_Folders/JSDMs_hmsc_HyRIV")


#Read data
data <- read.csv("../douro_comm_env_HyRIVERS.csv", stringsAsFactors = TRUE, header = TRUE)

#Organize data
#Construct XData dataframe (environmental variables) that will be used as fixed effects
XData <- data.frame(dplyr::select(data, tmp_dc_cmx, ele_mt_cav, pre_mm_cyr, ORD_STRA, crp_pc_use, S_CaCO3_av))
head(XData)

#Construct the YData dataframe (species data)
YData <- data.frame(dplyr::select(data, Mmar, Aana, Plit, Udel, Aarc, Agal, Aolig, Ccal, Cpal, Lboc, Pdur, Salb, Scar, Stru))

#Analyse species richness and prevalence
S = as.data.frame(rowSums(YData > 0))
max(S)
colnames(S) <- c("Species richness")
freq_S <- as.data.frame(table(S))
freq_S

P = as.data.frame(colMeans(YData > 0))

# x- and y-coordinates of sampling sites. We store these as the xy-matrix to be able to fit a spatial model
xy <- as.matrix(cbind(data$NEAR_X, data$NEAR_Y))

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
m_FULL = Hmsc(Y = YData, XData = XData, XFormula = XFormula,
              distr = "probit",
              studyDesign = studyDesign,
              ranLevels = list(HYRIV_ID=rL))


m_SPACE = Hmsc(Y = YData, XData = XData, XFormula = ~1,
               distr = "probit",
               studyDesign = studyDesign,
               ranLevels = list(HYRIV_ID= rL))

m_ENV = Hmsc(Y = YData, XData = XData, XFormula = XFormula, 
             distr = "probit")


m_FULL
m_SPACE
m_ENV


################################################################################
#FIT MODELS

#Combine the three models
# models_description <- list(m_FULL, m_SPACE, m_ENV)
models_fitted <- list()
# for (i in 1:length(models)) {
#   models_fitted[[i]] = sampleMcmc(models_description[[i]], thin = thin, samples = samples, transient = transient, nChains = nChains, nParallel = 3)
#   }
models_fitted[[1]] = sampleMcmc(m_ENV, thin = thin, samples = samples, transient = transient, nChains = nChains, nParallel = 3)


#Examine MCMC convergence
mpost <- convertToCodaObject(models_fitted[[1]], spNamesNumbers = c(T,F), covNamesNumbers = c(T,F))
psrf_beta <- gelman.diag(mpost$Beta, multivariate = FALSE)$psrf
tmp <- mpost$Omega[[1]]
z <- dim(tmp[[1]]) [2]

for (i in 1:length(tmp)) {
  tmp[[i]] = tmp[[i]]
  }

psrf_omega <- gelman.diag(tmp, multivariate = FALSE)$psrf
par(mfrow=c(1,2))

hist(psrf_beta, xlab = expression(paste("Potential reduction factor (parameter ", beta, ")")))
hist(psrf_omega, xlab = expression(paste("Potential reduction factor (parameter ", Omega, ")")))



#########################################################################################
#EVALUATING MODEL FIT AND COMPARING MODELS
partition = createPartition(models_fitted[[1]], nfolds = 2, column = "HYRIV_ID")
MF = list()
MFCV = list()
for (i in 1:2){
   preds = computePredictedValues(models_fitted[[i]])
   MF[[i]] = evaluateModelFit(hM = models_fitted[[i]], predY = preds)
   preds = computePredictedValues(models_fitted[[i]], partition = partition)
   MFCV[[i]] = evaluateModelFit(hM = models_fitted[[i]], predY = preds)
   }


WAIC = unlist(lapply(models_fitted, FUN = computeWAIC))
WAIC = as.data.frame(WAIC)


################################################################################
#Calculate between species correlation due to environmental covariates. Note this is done based on the posterior distribution of the species-specific regression coefficients, as opposed to being based on the prior distribution [which then be proportion to V]

calc_Xcor <- function(codasamples, model, est = "mean", prob = 0.95) {
   num_spp <- ncol(model$Y)
   all_Beta_samples <- do.call(rbind.data.frame, codasamples$Beta)
   num_mcmc <- nrow(all_Beta_samples)
   Xmatrix <- model$XScaled
   
   enviro_cor_mat <- enviro_cor_mat_cilower <- enviro_cor_mat_ciupper <- enviro_cov_mat <- sig_enviro_cor_mat <- enviro_cor_support_mat <- matrix(0, num_spp, num_spp)
   rownames(enviro_cor_mat) <- rownames(enviro_cor_mat_cilower) <- rownames(enviro_cor_mat_ciupper) <- rownames(enviro_cov_mat) <- rownames(sig_enviro_cor_mat) <- rownames(enviro_cor_support_mat) <- colnames(model$Y)
  colnames(enviro_cor_mat) <- colnames(enviro_cor_mat_cilower) <- colnames(enviro_cor_mat_ciupper) <- colnames(enviro_cov_mat) <- colnames(sig_enviro_cor_mat) <- colnames(enviro_cor_support_mat) <- colnames(model$Y)
   all_cors <- all_support <- array(NA, dim = c(num_spp, num_spp, num_mcmc))
   
   message("Calculating (posterior) between species correlation matrix due to measured covariates...")
   pb <- txtProgressBar(initial = 0, max = num_mcmc, style =3) 

   
   for(k0 in 1:num_mcmc) {
      cw_Beta <- matrix(as.vector(unlist(all_Beta_samples[k0,])), nrow = num_spp, byrow = TRUE)
      rownames(cw_Beta) <- colnames(model$Y)
      colnames(cw_Beta) <- colnames(model$XScaled)
      
      all_cors[,,k0] <- cor(tcrossprod(as.matrix(Xmatrix), cw_Beta))
      all_support[,,k0] <- all_cors[,,k0] > 0
      setTxtProgressBar(pb, k0)
      } 
   close(pb)
    
   for(j in 1:num_spp) {
      for(j2 in 1:num_spp) { ## Average/Median over the MCMC samples
         if(est == "median") {
            enviro_cov_mat[j, j2] <- median(all_cors[j, j2, ])
            enviro_cor_mat[j, j2] <- median(all_cors[j, j2, ])
            }
         if(est == "mean") {
            enviro_cov_mat[j, j2] <- mean(all_cors[j, j2, ])
            enviro_cor_mat[j, j2] <- mean(all_cors[j, j2, ])
            }
         
         sig_enviro_cor_mat[j, j2] <- enviro_cor_mat[j, j2]
         get.hpd.cors <- HPDinterval(as.mcmc(all_cors[j, j2,]), prob = prob)
         enviro_cor_mat_cilower[j, j2] <- get.hpd.cors[1]
         enviro_cor_mat_ciupper[j, j2] <- get.hpd.cors[2]
         if(0 > get.hpd.cors[1] & 0 < get.hpd.cors[2]) {
            sig_enviro_cor_mat[j, j2] <- 0
            }
         
         enviro_cor_support_mat[j, j2] <- mean(all_support[j, j2, ])
         }
   }

  return(list(cor = enviro_cor_mat, 
              cor_lower = enviro_cor_mat_cilower, 
              cor_upper = enviro_cor_mat_ciupper, 
              sig_cor = sig_enviro_cor_mat,
              support = enviro_cor_support_mat))
   }   


Xcorrelations_ENV <- calc_Xcor(codasamples = convertToCodaObject(models_fitted[[1]], spNamesNumbers = c(T,F), covNamesNumbers = c(T,F)),
                               model = m_ENV)


supportLevel = 0.95
toPlot = ((Xcorrelations_ENV$support > supportLevel)
          + (Xcorrelations_ENV$support < (1-supportLevel))>0) * Xcorrelations_ENV$cor


#########################################################################################
#Exploring Parameter estimates
VP = list()
for (i in 1:length(models_fitted)) {
  VP[[i]] = computeVariancePartitioning(models_fitted[[i]])
  }
   

#Plot Variance partitioning with Hmsc
plotVariancePartitioning(models_fitted[[1]], VP[[1]], args.legend = list(x = "topright", bty = "n", inset=c(-0.15, 0)))


#Plot heatmap
postBeta <- getPostEstimate(models_fitted[[1]], parName = "Beta")
plotBeta(models_fitted[[1]], post=postBeta, supportLevel = 0.95, param = "Sign", plotTree = FALSE, spNamesNumbers = c(TRUE, FALSE), 
         covNamesNumbers = c(TRUE, FALSE), colors = colorRampPalette(c("grey","white","black")), colorLevels = 3)


# species associations revealed by the random effects with the corrplot function
par(mfrow = c(1,2))
#Model FULL - residual associations that account for the responses of species to environmental conditions
OmegaCor = computeAssociations(models_fitted[[1]])
supportLevel = 0.85
toPlot = ((OmegaCor[[1]]$support>supportLevel)
          + (OmegaCor[[1]]$support<(1-supportLevel))>0)*OmegaCor[[1]]$mean


corrplot(toPlot, method="color", type = "lower", tl.srt = 45, col=c("grey", "white", "black"), title = "Residual co-occurrences", 
         tl.cex = 0.8, tl.col = "black", font = 3, mar=c(0,0,1,0))


#Model SPATIAL raw associations
OmegaCor = computeAssociations(models_fitted[[2]])
supportLevel = 0.85
toPlot = ((OmegaCor[[1]]$support>supportLevel)
          + (OmegaCor[[1]]$support<(1-supportLevel))>0)*OmegaCor[[1]]$mean


corrplot(toPlot, method="color", type = "lower", tl.srt = 45, col=c("grey", "white", "black"), title = "Raw co-occurrences", 
         tl.cex = 0.8, tl.col = "black", font = 3, mar=c(0,0,1,0))
