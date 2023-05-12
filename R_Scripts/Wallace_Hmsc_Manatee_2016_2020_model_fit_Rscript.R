setwd("/home/abauer/hmsc_manatee")

library(Hmsc)
set.seed(1)

## Load the model 
load("mod_pland2500m_chains3_samples250_thin_1000.Rdata")


#### evauluate model fit ####

#explanatory power
preds = computePredictedValues(mod)
MF = evaluateModelFit(hM=mod, predY=preds)

save(MF, file="mod_pland2500m_chains3_samples250_thin_1000_MFexpl.Rdata")

#predictive power
partition = createPartition(mod, nfolds = 2)
preds2 = computePredictedValues(mod, partition = partition)
MFCV = evaluateModelFit(hM = mod, predY = preds2)

save(MFCV, file="mod_pland2500m_chains3_samples250_thin_1000_MFpred.Rdata")
