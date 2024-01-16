##===========================================================================================================================================================================================================================================================================#
# Simulation code for observational, Linear propensity, use indicator as observed Y, KM ipcw for censoring KM, pseudo reg. as rf, estimate propensity using RF, conditional survival prob: RSF; kernel shap with 100 test point and 100 randomly sampled trained for bg data #
# Author: Na Bo                                                                                                                                                                                                                                                              #
#============================================================================================================================================================================================================================================================================#
### scenario 1
#setwd("your working directory containing simulation functions")
library(survival)
library(randomForestSRC)
library(grf)
library(tidyverse)
library(randomForest)
library(pec)
library(glmnet)
library(kernelshap)

# load functions
source("meta_learner_simulation_censor_survival_prob_kernelshap_ipcw.R")
source("weibull_true.R")

# load simulation data
load(paste0("your working directory containing simulated data","/balanced_setting1_traindata1.RData")) # change data for different settings
load(paste0("your working directory containing simulated data","/observational_10p1000n10000test_balanced_scenario1.RData"))
load(paste0("your working directory containing simulated data","/observational_10p1000n100test_kernelshap_balanced_scenario1.RData"))

set.seed(100871+1+15213) 
tt<-12

# run meta-learners
rsf_X<-rsf_HTE_X(data=dat,testdat=testdat,testdat_kernelshap=testdat_shap,kernelshap_bg=testdat_shap[,which(names(testdat_shap$data)%in%c("V1", "V2","V3","V4", "V5","V6", "V7", "V8","V9", "V10","Treatment"))],ntrees=1000,time.interest=tt,xtype=2,k.folds=10,IPCW.method="KM",IPCW=T,pseudo_reg_separate_trt=T,pseudo_reg="weighted_RF",est_pi=1, propensity_method="regression_forest",propensity_method_X=c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"), propensity=0.5)
rsf_M <- rsf_HTE_M(data=dat,testdat=testdat,testdat_kernelshap=testdat_shap,kernelshap_bg=testdat_shap[,which(names(testdat_shap)%in%c("V1", "V2","V3","V4", "V5","V6", "V7", "V8","V9", "V10","Treatment"))],ntrees=1000,time.interest=tt,IPCW=T,IPCW.method="KM",k.folds=10,impute_type=2,est_pi=1,propensity_method="regression_forest",propensity_method_X=c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"), propensity=0.5,pseudo_reg="weighted_RF")
rsf_DR <- rsf_HTE_DR_Kennedy(data=dat,testdat=testdat,testdat_kernelshap=testdat_shap,kernelshap_bg=testdat_shap[,which(names(testdat_shap)%in%c("V1", "V2","V3","V4", "V5","V6", "V7", "V8","V9", "V10","Treatment"))],ntrees=1000,time.interest=tt,impute_type=2,IPCW=T,IPCW.method = "KM",k.folds=10,est_pi=1,propensity_method="regression_forest",propensity_method_X=c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"),propensity=0.5,pseudo_reg="weighted_RF")
rsf_R <- rsf_HTE_R(data=dat,testdat=testdat,testdat_kernelshap=testdat_shap,kernelshap_bg=testdat_shap[,which(names(testdat_shap)%in%c("V1", "V2","V3","V4", "V5","V6", "V7", "V8","V9", "V10","Treatment"))],ntrees=1000,time.interest=tt,impute_type=2,IPCW=T,k.folds=10,IPCW.method="KM",est_pi=1,propensity_method="regression_forest",propensity_method_X=c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"),propensity=0.5,pseudo_reg="weighted_RF")
rsf_D <- rsf_HTE_D(data=dat,testdat=testdat,testdat_kernelshap=testdat_shap,kernelshap_bg=testdat_shap[,which(names(testdat_shap)%in%c("V1", "V2","V3","V4", "V5","V6", "V7", "V8","V9", "V10","Treatment"))],ntrees=1000,time.interest=tt,impute_type=2,IPCW=T,IPCW.method ="KM",k.folds=10,EA=F,est_pi=1,propensity_method="regression_forest",propensity_method_X=c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"),propensity=0.5,pseudo_reg="weighted_RF")
rsf_D_EA <- rsf_HTE_D(data=dat,testdat=testdat,testdat_kernelshap=testdat_shap,kernelshap_bg=testdat_shap[,which(names(testdat_shap)%in%c("V1", "V2","V3","V4", "V5","V6", "V7", "V8","V9", "V10","Treatment"))],ntrees=1000,time.interest=tt,impute_type=2,IPCW=T,IPCW.method ="KM",k.folds=10,EA=T,est_pi=1,propensity_method="regression_forest",propensity_method_X=c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"),propensity=0.5,pseudo_reg="weighted_RF")

# run Weibull true model
weib.true<-weibull_true(data=dat,testdat=testdat,time.interest=tt,scenario=1)

# run CSF
cs.forest.prob <- causal_survival_forest(dat[,which(names(dat)%in%c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"))], dat$Time, dat$Treatment, dat$Event,target = "survival.probability",horizon = tt,num.trees = 500,mtry = sqrt(15),min.node.size = 15,num.threads=1)
csf.prob = predict(cs.forest.prob, testdat[,which(names(testdat)%in%c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"))])
csa_diff<-data.frame(diff=csf.prob$predictions)
shap_csf <- kernelshap(cs.forest.prob, as.matrix(testdat_shap[,which(names(testdat_shap)%in%c("V1", "V2","V3","V4", "V5","V6", "V7", "V8","V9", "V10"))]), bg_X = as.matrix(testdat_shap[,which(names(testdat_shap)%in%c("V1", "V2","V3","V4", "V5","V6", "V7", "V8","V9", "V10","Treatment"))]))
shapvip_csf <- shap_csf$S

# combine results
diff<-data.frame(rsf_X_diff=rsf_X$diff,rsf_M_diff=rsf_M$diff,rsf_DR_diff=rsf_DR$diff,rsf_R_diff=rsf_R$diff,rsf_D_diff=rsf_D$diff,rsf_D_EA_diff=rsf_D_EA$diff,csf_diff=csa_diff$diff,weib_diff=weib.true$diff,true_diff=mydata$true.diff.survprob)
colnames(diff) <- c("rsf_X_diff","rsf_M_diff","rsf_DR_diff","rsf_R_diff","rsf_D_diff","rsf_D_EA_diff","csf_diff","weib_diff","true_diff")
shapvip <- as.data.frame(cbind(rsf_X$shapvip,rsf_M$shapvip,rsf_DR$shapvip,rsf_R$shapvip,rsf_D$shapvip,rsf_D_EA$shapvip,shapvip_csf))
colnames(shapvip) <- c(paste0("X_V",seq(1:10)),paste0("M_V",seq(1:10)),paste0("DR_V",seq(1:10)),paste0("R_V",seq(1:10)),paste0("D_V",seq(1:10)),paste0("DEA_V",seq(1:10)),paste0("CSF_V",seq(1:10)))

rm(list=ls())

