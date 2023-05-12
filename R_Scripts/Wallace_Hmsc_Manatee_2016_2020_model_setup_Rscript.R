

setwd("/home/abauer/hmsc_manatee")


library(Hmsc)
set.seed(1)
library(raster)
library(maptools)
library(parallel)
library(corrplot)

########################
# read in the data

dat<-read.csv("manatee_mosq_by_site_pland_2016to20.csv", stringsAsFactors = F)

dim(dat)
names(dat)

### environmental features
XData = data.frame(ID = dat$ID_final, pl_1 = dat$pland_1, pl_2 = dat$pland_2, pl_13 = dat$pland_13, pl_15 = dat$pland_15, pl_16 = dat$pland_16, year = dat$year)

XData$ID <- as.factor(XData$ID)
XData$year <- as.factor(XData$year)


### species matrix ###
#truncated to presence/absence:
Y = as.matrix(dat[,3:47]) >0 
Y = apply(Y,MARGIN = 2,FUN = as.numeric) # some post-processing functions assume Y has num 1/0 values


### spatial data ###
coords <- read.csv("manatee_coords_all_IDs.csv", stringsAsFactors = F)

xy = as.matrix(cbind(coords$Long_final,coords$Lat_final))
rownames(xy)=coords$ID_final
colnames(xy)=c("x-coordinate","y-coordinate")


### trait matrix ###
alltraits = read.csv("Manatee_2016_2020_traits.csv", stringsAsFactors = T)

Tr = as.data.frame(alltraits)
rownames(Tr) <- Tr$Species


#######################

### define spatial random effect by including sData
rL = HmscRandomLevel(sData=xy)


### study design ###
# define random effect at the level of ID and year

studyDesign = data.frame(ID = XData$ID, year = XData$year)

ID = HmscRandomLevel(units = studyDesign$ID)
year = HmscRandomLevel(units = studyDesign$year)

ranlevels = list(ID = rL, year = year) 
ranlevels

# set XFormula

#used VIF (Zuur et al. 2009) to check no multicollinearity of variables:
XFormula = ~ pl_1 + pl_2 + pl_13 + pl_15 + pl_16

# set Trait formula
TrFormula = ~ WNV


####################
# set up the model

simul <- Hmsc(Y=Y, XData = XData, XFormula = XFormula, TrData = Tr, TrFormula = TrFormula, studyDesign = studyDesign, ranLevels=ranlevels, distr = "probit")  
# probit distribution: presence/absence data

simul # explores model object

# RUN ON CLUSTER ONLY 
thin = 1000
samples = 250 
nChains = 3
transient = 50*thin


mod = sampleMcmc(simul, samples = samples, transient = transient, thin = thin, nChains = nChains, nParallel = nChains, alignPost = TRUE)

save(mod,file="mod_pland2500m_chains3_samples1000_thin_1000.Rdata") # run on cluster

