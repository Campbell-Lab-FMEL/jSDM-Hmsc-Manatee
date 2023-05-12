setwd("C:/Users/amely/OneDrive - University of Florida/B_Campbell_Lab/E_HMSC_Manatee_paper_20162020/hmsc")

library(Hmsc)
library(ggplot2)
library(tidyverse)
library(vioplot)

set.seed(1)


#############################################
#### S2 Model convergence and fit ###########
#############################################

## Load the model 
load("mod_pland2500m_chains3_samples250_thin_1000.Rdata")


#### Checking MCMC convergence diagnostics 

## extract posterior distribution from model object
post = convertToCodaObject(mod, spNamesNumbers = c(T,F), covNamesNumbers = c(T,F))

# Beta parameters (species-environment)
# Gamma parameters (trait-environment)
par(mfrow=c(1,2))

# effective sample size
hist(effectiveSize(post$Beta), main="ess(beta)") 
hist(effectiveSize(post$Gamma), main="ess(gamma)")

# potential scale reduction factor
summary(gelman.diag(post$Beta, multivariate = FALSE)$psrf)
# point mean: 1.0031, upper CI mean: 1.0138 
summary(gelman.diag(post$Gamma, multivariate = FALSE)$psrf)
# point mean: 1.0022, upper CI mean: 1.0113 

vioplot(gelman.diag(post$Beta, multivariate = FALSE)$psrf,
        main="psrf(beta)") 
vioplot(gelman.diag(post$Gamma,multivariate = FALSE)$psrf,
        main="psrf(gamma)")

#############################################
#### S3 Inspecting model fit ################
#############################################

# explanatory power - MF
load("mod_pland2500m_chains3_samples250_thin1000_MFexpl.Rdata")
# predictive power - MFCV
load("mod_pland2500m_chains3_samples250_thin1000_MFpred.Rdata")


## model fit: explanatory power
#plot discriminatory power Tjur R2
#plot discriminatory power AUC
plot(colSums(mod$Y, na.rm = T)/mod$ny, MF$AUC, 
     main=paste("AUC . mean = ", round(mean(MF$AUC,na.rm = TRUE),2), sep=""), 
     xlab = "Prevalence")
plot(colSums(mod$Y, na.rm = T)/mod$ny, MF$TjurR2,
     main=paste("TjurR2 . Mean = ", round(mean(MF$TjurR2,na.rm = TRUE),3), sep=""), 
     xlab = "Prevalence")
par(mfrow=c(1,1))


# mean RSME, AUC, and Tjur R2
tmp = c(mean(MF$AUC,na.rm = TRUE), mean(MF$TjurR2, na.rm = TRUE), 
        mean(MFCV$AUC, na.rm = TRUE), mean(MFCV$TjurR2, na.rm = TRUE))
names(tmp)=c("AUC", "TjurR2", "AUC (CV)", "TjurR2 (CV)")
tmp

# plot predictive over explanatory power for each species
par(mfrow=c(1,2))

## AUC
plot(MF$AUC,MFCV$AUC,xlim=c(0,1),ylim=c(0,1),
     xlab = "explanatory power",
     ylab = "predictive power",
     main = paste0("AUC\n",
                   "mean(MF) = ",as.character(round(mean(MF$AUC, na.rm=TRUE), digits = 2)),
                   ", mean(MFCV) = ",as.character(round(mean(MFCV$AUC, na.rm=TRUE), digits = 2))),
     abline(0,1, col="darkblue"))
  abline(v=0.5,  col="grey")
  abline(h=0.5,  col="grey")

## TjurR2
plot(MF$TjurR2,MFCV$TjurR2,xlim=c(0,1),ylim=c(0,1),
     xlab = "explanatory power",
     ylab = "predictive power",
     main = paste0("Tjur R2\n",
                   "mean(MF) = ",as.character(round(mean(MF$TjurR2, na.rm=TRUE), digits = 3)),
                   ", mean(MFCV) = ",as.character(round(mean(MFCV$TjurR2, na.rm=TRUE), digits = 3))),
     abline(0,1, col="darkblue"))

par(mfrow=c(1,1))

####
prev <- colSums(mod$Y, na.rm = T)/mod$ny
AUC.ex <- MF$AUC 
AUC.pred <- MFCV$AUC
TR2.ex <- MF$TjurR2
TR2.pred <- MFCV$TjurR2

spp_fit <- as.data.frame(cbind(prev, AUC.ex, AUC.pred, TR2.ex, TR2.pred))
write.csv(spp_fit, "results/pland2500m_Tr2_chains3_samples250_thin_1000_spp_fit.csv")

#############################################
#### S4 Variance partitioning ###############
#############################################


##### proportion of explained variance ####
#### Partitioning of the variance explained (as captured by Tjur R2)
load("mod_pland2500m_chains3_samples250_thin_1000.Rdata")

par(mar=c(4,4,4,4))

## show variance partitioning between all random and fixed effects
VP = computeVariancePartitioning(mod, na.ignore = T)
VP.t = VP # to calculate prop of total variance explained


VP$vals <- VP$vals[,-c(11,21)] # removes An. albimanus and Cx. nigripalpus 
mycols = rainbow(nrow(VP$vals))
plotVariancePartitioning(mod , VP = VP,
                         cols = mycols,
                         las = 2, # las = 2: y-axis labels horizontal, x-axis labels vertical
                         main = "Variance Partitioning \nProportion of Explained Variance (Tjur R2)") 

vpd <- VP$vals # proportion of explained variance for covariate and species
write.csv(vpd, "var_par_expl_manatee2016_20.csv") 


## show variance partitioning with landscape variables grouped:
mod$covNames # check all covariates 
group_LC=c(1,1,1,1,1,1) # Group partitioning based on included covariates, grouping includes intercept
groupnames_LC = c("land cover")

VP_LC = computeVariancePartitioning(mod, group = group_LC, groupnames = groupnames_LC, na.ignore = T) 

VP_LC$vals <- VP_LC$vals[,-c(11,21)] # removes An. albimanus and Cx. nigripalpus 
plotVariancePartitioning(mod, VP = VP_LC,
                         cols = mycols,
                         las = 2, 
                         main = "Variance Partitioning \nProportion of Explained Variance (Tjur R2)")

vpd_LC <- VP_LC$vals # proportion of explained variance for covariate and species
write.csv(vpd_LC, "var_par_expl_LC_manatee2016_20.csv") 

######

##### proportion of total variance ####
## partitioning the total variance (as captured by Tjur R2)

# calculate total explained variance 
for(k in 1:mod$ns){
  VP.t$vals[,k] = MF$TjurR2[k]*VP.t$vals[,k]
  }

vpd.t <- VP.t$vals %>% rbind(MF$TjurR2) # can use TjurR2 to order by explained variance 
row.names(vpd.t)[8] <- "TjurR2"
vpd.t <- vpd.t[,-c(11,21)] # removes An. albimanus and Cx. nigripalpus 

write.csv(vpd.t, "var_par_tot_manatee2016_20.csv")




#############################################
#### S5 pred over land cover gradients ######
#############################################

## total: fixes non-focal covariates to their most likely value, given the value of focal variable (linear relationship)

    # Gradient.t = constructGradient(mod, focalVariable = "pl_1")
    # predY = predict(mod, Gradient=Gradient.t, expected = TRUE)

# marginal: fixes non-focal covariates to the mean of their respective values
  # Problem: with land cover percentages, total land cover accummulates to more than 100% 

    # Gradient.m = constructGradient(mod, focalVariable = "pl_1", non.focalVariables = 1)
    # predY2 = predict(mod, Gradient=Gradient.m, expected = TRUE) 

  
## Load the model 
load("mod_pland2500m_chains3_samples250_thin_1000.Rdata")

drkgrey <- rgb(0.5,0.5,0.5, alpha=0.5) # data points color: black, transparency: 50%

par(mfrow=c(2,1))

# water  ###########
# cicol <- confidence intervall. rgb as percent values
# RGB: 47,84,151

Gradient.m = constructGradient(mod, focalVariable = "pl_1", non.focalVariables = 1)
predY2 = predict(mod, Gradient=Gradient.m, expected = TRUE)

Gradient.t = constructGradient(mod, focalVariable = "pl_1")
predY = predict(mod, Gradient=Gradient.t, expected = TRUE)

### visualization
# species richness
plotGradient(mod, Gradient.m, pred=predY2, measure="S", xlabel="Water [%]", 
             cicol= rgb(0.18,0.33,0.59, alpha = 0.5), pointcol = drkgrey,
             showData = TRUE, 
             showPosteriorSupport = T, 
             main = "Species presence-absence model: summed response (marginal effect)")
plotGradient(mod, Gradient.t, pred=predY, 
             measure="S", 
             xlabel="Water [%]",
             cicol= rgb(0.18,0.33,0.59, alpha = 0.5),
             pointcol = drkgrey,
             showData = TRUE, 
             showPosteriorSupport = T, 
             main = "Species presence-absence model: summed response (total effect)")

# CWMT WNV
plotGradient(mod, Gradient.m, pred=predY2, 
             measure="T", index = 2, # plots WNV
             xlabel="Water [%]",
             cicol= rgb(0.18,0.33,0.59, alpha = 0.5),
             pointcol = drkgrey,
             showData = TRUE, 
             showPosteriorSupport = T, 
             main = "Species presence-absence model: summed response (marginal effect)")
plotGradient(mod, Gradient.t, pred=predY, 
             measure="T", index = 2, # plots WNV
             xlabel="Water [%]",
             cicol= rgb(0.18,0.33,0.59, alpha = 0.5),
             pointcol = drkgrey,
             showData = TRUE, 
             showPosteriorSupport = T, 
             main = "Species presence-absence model: summed response (total effect)")


# developed ###########
# RGB: 202,32,38

Gradient.m = constructGradient(mod, focalVariable = "pl_2", non.focalVariables = 1)
predY2 = predict(mod, Gradient=Gradient.m, expected = TRUE)

Gradient.t = constructGradient(mod, focalVariable = "pl_2")
predY = predict(mod, Gradient=Gradient.t, expected = TRUE)

## test: reduce predictions along gradient
#EpredY = Reduce("+", predY2)/length(predY2)

### visualization
# species richness
plotGradient(mod, Gradient.m, pred=predY2, 
             measure="S", 
             xlabel="Developed [%]",
             cicol= rgb(0.79,0.13,0.15, alpha = 0.5),
             pointcol = drkgrey,
             showData = TRUE, 
             showPosteriorSupport = T, 
             main = "Species presence-absence model: summed response (marginal effect)")
plotGradient(mod, Gradient.t, pred=predY, 
             measure="S", 
             xlabel="Developed [%]",
             cicol= rgb(0.79,0.13,0.15, alpha = 0.5),
             pointcol = drkgrey,
             showData = TRUE, 
             showPosteriorSupport = T, 
             main = "Species presence-absence model: summed response (total effect)")

# CWMT WNV
plotGradient(mod, Gradient.m, pred=predY2, 
             measure="T", index=2, # plots WNV
             xlabel="Developed [%]",
             cicol= rgb(0.79,0.13,0.15, alpha = 0.5),
             pointcol = drkgrey,
             showData = TRUE, 
             showPosteriorSupport = T, 
             main = "Species presence-absence model: summed response (marginal effect)")
plotGradient(mod, Gradient.t, pred=predY, 
             measure="T", index = 2, # plots WNV
             xlabel="Developed [%]",
             cicol= rgb(0.79,0.13,0.15, alpha = 0.5),
             pointcol = drkgrey,
             showData = TRUE, 
             showPosteriorSupport = T, 
             main = "Species presence-absence model: summed response (total effect)")


# cropland  ###########
# RGB: 234,165,38

Gradient.m = constructGradient(mod, focalVariable = "pl_13", non.focalVariables = 1)
predY2 = predict(mod, Gradient=Gradient.m, expected = TRUE)

Gradient.t = constructGradient(mod, focalVariable = "pl_13")
predY = predict(mod, Gradient=Gradient.t, expected = TRUE)

### visualization
# species richness
plotGradient(mod, Gradient.m, pred=predY2, 
             measure="S", 
             xlabel="Cropland [%]",
             cicol= rgb(0.92,0.65,0.15, alpha = 0.5),
             pointcol = drkgrey,
             showData = TRUE, 
             showPosteriorSupport = T, 
             main = "Species presence-absence model: summed response (marginal effect)")
plotGradient(mod, Gradient.t, pred=predY, 
             measure="S", 
             xlabel="Cropland [%]",
             cicol= rgb(0.92,0.65,0.15, alpha = 0.5),
             pointcol = drkgrey,
             showData = TRUE, 
             showPosteriorSupport = T, 
             main = "Species presence-absence model: summed response (total effect)")

# CWMT WNV
plotGradient(mod, Gradient.m, pred=predY2, 
             measure="T", index=2, # plots WNV
             xlabel="Cropland [%]",
             cicol= rgb(0.92,0.65,0.15, alpha = 0.5),
             pointcol = drkgrey,
             showData = TRUE, 
             showPosteriorSupport = T, 
             main = "Species presence-absence model: summed response (marginal effect)")
plotGradient(mod, Gradient.t, pred=predY, 
             measure="T", index = 2, # plots WNV
             xlabel="Cropland [%]",
             cicol= rgb(0.92,0.65,0.15, alpha = 0.5),
             pointcol = drkgrey,
             showData = TRUE, 
             showPosteriorSupport = T, 
             main = "Species presence-absence model: summed response (total effect)")


# herb_wetland   ###########
# RGB: 101,179,214

Gradient.m = constructGradient(mod, focalVariable = "pl_15", non.focalVariables = 1)
predY2 = predict(mod, Gradient=Gradient.m, expected = TRUE)

Gradient.t = constructGradient(mod,focalVariable = "pl_15")
predY = predict(mod, Gradient=Gradient.t, expected = TRUE)

### visualization
# species richness
plotGradient(mod, Gradient.m, pred=predY2, 
             measure="S", 
             xlabel="Herbaceous Wetland [%]",
             cicol= rgb(0.4,0.7,0.84, alpha = 0.5),
             pointcol = drkgrey,
             showData = TRUE, 
             showPosteriorSupport = T, 
             main = "Species presence-absence model: summed response (marginal effect)")
plotGradient(mod, Gradient.t, pred=predY, 
             measure="S", 
             xlabel="Herbaceous Wetland [%]",
             cicol= rgb(0.4,0.7,0.84, alpha = 0.5),
             pointcol = drkgrey,
             showData = TRUE, 
             showPosteriorSupport = T, 
             main = "Species presence-absence model: summed response (total effect)")

# CWMT WNV
plotGradient(mod, Gradient.m, pred=predY2, 
             measure="T", index=2, # plots WNV
             xlabel="Herbaceous Wetland [%]",
             cicol= rgb(0.4,0.7,0.84, alpha = 0.5),
             pointcol = drkgrey,
             showData = TRUE, 
             showPosteriorSupport = T, 
             main = "Species presence-absence model: summed response (marginal effect)")
plotGradient(mod, Gradient.t, pred=predY, 
             measure="T", index = 2, # plots WNV
             xlabel="Herbaceous Wetland [%]",
             cicol= rgb(0.4,0.7,0.84, alpha = 0.5),
             pointcol = drkgrey,
             showData = TRUE, 
             showPosteriorSupport = T, 
             main = "Species presence-absence model: summed response (total effect)")


# woody_wetland   ###########
# RGB: 98,146,63

Gradient.m = constructGradient(mod, focalVariable = "pl_16", non.focalVariables = 1)
predY2 = predict(mod, Gradient=Gradient.m, expected = TRUE)

Gradient.t = constructGradient(mod, focalVariable = "pl_16")
predY = predict(mod, Gradient = Gradient.t, expected = TRUE)

### visualization
# species richness
plotGradient(mod, Gradient.m, pred=predY2, 
             measure="S", 
             xlabel="Woody Wetland [%]",
             cicol= rgb(0.38,0.57,0.25, alpha = 0.5),
             pointcol = drkgrey,
             showData = TRUE, 
             showPosteriorSupport = T, 
             main = "Species presence-absence model: summed response (marginal effect)")
plotGradient(mod, Gradient.t, pred=predY, 
             measure="S", 
             xlabel="Woody Wetland [%]",
             cicol= rgb(0.38,0.57,0.25, alpha = 0.5),
             pointcol = drkgrey,
             showData = TRUE, 
             showPosteriorSupport = T, 
             main = "Species presence-absence model: summed response (total effect)")

# CWMT WNV
plotGradient(mod, Gradient.m, pred=predY2, 
             measure="T", index=2, # plots WNV
             xlabel="Woody Wetland [%]",
             cicol= rgb(0.38,0.57,0.25, alpha = 0.5),
             pointcol = drkgrey,
             showData = TRUE, 
             showPosteriorSupport = T, 
             main = "Species presence-absence model: summed response (marginal effect)")
plotGradient(mod, Gradient.t, pred=predY, 
             measure="T", index = 2, # plots WNV
             xlabel="Woody Wetland [%]",
             cicol= rgb(0.38,0.57,0.25, alpha = 0.5),
             pointcol = drkgrey,
             showData = TRUE, 
             showPosteriorSupport = T, 
             main = "Species presence-absence model: summed response (total effect)")



#############################################
#### S6 spatial predictions #################
#############################################

### load model 
hmsc.mod <- "mod_pland2500m_chains3_samples250_thin_1000.Rdata"

load(hmsc.mod)

mod.name <-  hmsc.mod %>% 
  str_remove(pattern = "mod") %>% 
  str_remove(pattern = ".Rdata") 

### load grid
grid2 = read.csv("pland_1000m_grid.csv")

# set up output folder
out_path <- "results/xy_pred"
if (!dir.exists(out_path)) dir.create(out_path)


# change here for resp. output name:
#ex.out <- "T"  # expected = T -- predicts occurrence probabilities (e.g 0.2, )
ex.out <- "F"  # expected = F -- predicts occurence (0 or 1)
########################


### make predictions to smaller subsets of grid to be more efficient ####

# find adequate sub-units by identifying run time for subset prediction:
# startTime <- Sys.time()
# ## function() # 
# endTime <- Sys.time()
# print(endTime - startTime)

n_sub <- 7  # total number of subsets

# point.range <- nrow(grid2)/n_sub
# 2275 / 7 = 325
# 325 points: 38.20671 secs 


#### loop1: predict communities at prediction points ##############

### Loop preparation: count variables

# loop progress counter
loop.index <- 0

# update prediction points
point.range <- nrow(grid2)/n_sub

# start values for point range:
n1 <- 1
n2 <- point.range

# empty dfs to write into:
grid.c = data.frame() # control pred point selection
pred.all = data.frame() # gather pres and grid info into one df

### Run loop #####
for(i in 1:n_sub) {
  
  if (i == 1) startTime <- Sys.time()
  
  # move counter along
  loop.index <- loop.index + 1
  
  ## select grid subset
  grid = grid2[n1:n2, ]
  
  # gather everything into one df part1 - grid :
  g.temp1 = grid
  g.temp1$loop.n <- loop.index
  grid.c = rbind(grid.c, g.temp1)
  
  xy.grid = as.matrix(cbind(grid$x, grid$y))
  
  XData.grid = data.frame(pl_1 = grid$pland_1, pl_2 = grid$pland_2, pl_13 = grid$pland_13, pl_15 = grid$pland_15, pl_16 = grid$pland_16, stringsAsFactors = T)
  
  ### prepare gradient for prediction
  Gradient = prepareGradient(mod, XDataNew = XData.grid, sDataNew = list(ID=xy.grid))
  
  # Make prediction and reduce
  nParallel=2
  predY = predict(mod, Gradient=Gradient, expected = T, nParallel=nParallel)
  # expected = T -- predicts occurrence probabilities (e.g 0.2, )
  # expected = F -- predicts occurence (0 or 1)
  
  EpredY=Reduce("+",predY)/length(predY) # EpredY is a matrix of posterior mean occurrence probabilities
  
  # save subset prediction object:
  EpredY.name <- paste(out_path, "/EpredY_sub_",n1,"_", n2, "_exp_", ex.out, mod.name, ".rds", sep = "")
  saveRDS(EpredY,  file = EpredY.name )
  
  # populate gather df
  g.temp1 <- cbind(g.temp1, EpredY)
  pred.all = rbind(pred.all, g.temp1)
  
  
  # progress message:
  message <- paste("\n successfully looped through subset ", loop.index, " of ", n_sub,
                   "\n predicted communities to ", n2, " out of ", nrow(grid2), " prediction points", sep = "")
  cat(message)
  
  
  # update point ranges:
  n1 <- n1 + point.range
  n2 <- n2 + point.range
  
  
  if (i == n_sub){
    # save full dataset:
    all.name <- paste(out_path, "/PredY_xy_cov_all_exp_", ex.out,  mod.name, ".rds", sep = "")
    saveRDS(pred.all,  file = all.name)
    
    all.name2 <- paste(out_path, "/PredY_xy_all_exp_", ex.out, mod.name, ".rds", sep = "")
    saveRDS(pred.all[c(1:2, 9:53)],  file = all.name2)
    
    endTime <- Sys.time()
    T_tot <- (endTime - startTime)
    cat("\n Duration of spatial prediction loop: ", round(T_tot, 2), " minutes", sep = "")
  }
  
}



# # check selected pred. points::
# library(ggplot2)
# library(viridis)
# 
# ggplot(data = pred.all, aes(x=x, y=y, fill=loop.n)) +
#   geom_raster() +
#   scale_fill_viridis()

########################


#### Prediction Maps ##########################################################################################

## load combined spatial dataset prediction data for model:
pred.all <- readRDS(list.files(path = out_path, pattern = str_c("xy_all_exp_", ex.out, mod.name, ".rds$"), full.names = TRUE))

### spatial information
xy = pred.all[,1:2]
### predicted species values
preddat <- as.matrix(pred.all[, 3:47])

### calculate community-based values:
S = rowSums(preddat) # species richness
CWM = (preddat%*%mod$Tr)/matrix(rep(S, mod$nt),ncol=mod$nt) # Community-weighted mean prop of WNV vector spp

# write df for data viz:
mapData = data.frame(xy,S,CWM[, 2:3], stringsAsFactors=TRUE)


#### Visualizations:

##### basic maps (no smoothing) ####
### map out predicted species richness 
library(viridis)

ggplot(data = mapData, aes(x=x, y=y)) + 
  geom_raster(aes(fill = S)) + 
  scale_fill_viridis(option = "G", 
                     end=0.97, # end = 0.9x tone down bright turquoise top value a little
                     na.value = "transparent") + 
  xlab("Longitude") +
  ylab("Latitude") +
  coord_equal() +
  theme_bw()


### map out predicted CWMTs:
library(wesanderson)
#define wes anderson color palette
pal <- wes_palette("Zissou1", 100, type = "continuous")

# CWMT range:
# 0.3046071 - 0.5135398

ggplot(data = mapData, aes(x=x, y=y)) +
  geom_raster(aes(fill = WNV)) +
  #  ggtitle("Predicted proportion of WNV vector species") + 
  xlab("Longitude") +
  ylab("Latitude") +
  scale_fill_gradientn(colours = pal, limits = c(0.28,0.53), 
                       #name="WNV competent \nspecies in %"
  ) +
  coord_equal() +
  theme_bw()


