#--------------------------------------------------- X-learner --------------------------------------------------------#
rsf_HTE_X <- function(data=dat,testdat=mydata$data,testdat_kernelshap=mydata_kernelshap$data,kernelshap_bg=train_kernelshap_bg,ntrees=1000,time.interest=12,xtype=2,k.folds=10,IPCW.method="KM",IPCW=T,pseudo_reg_separate_trt=T,pseudo_reg="weighted_RF",est_pi=1,propensity_method="correct_logistic",propensity_method_X=c("V1","V5"), propensity=0.5){
  data<-data[complete.cases(data),]
  data1 <- data[which(data$Treatment==1),-c(which(names(data)=="Treatment"), which(names(data)=="sim"))]
  data0 <- data[which(data$Treatment==0),-c(which(names(data)=="Treatment"), which(names(data)=="sim"))]
  
  #Step1: Fit 2 RFSRC models
  rfsrc_data1 <-
    grf::survival_forest(X=as.matrix(data1[,which(names(data1) %in% c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"))]),Y=as.vector(data1$Time),D=as.vector(data1$Event),prediction.type="Nelson-Aalen",num.trees = 500,num.threads = 1) # OOB estimates
  
  rfsrc_data0 <-
    grf::survival_forest(X=as.matrix(data0[,which(names(data0) %in% c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"))]),Y=as.vector(data0$Time),D=as.vector(data0$Event),prediction.type="Nelson-Aalen",num.trees = 500,num.threads = 1) # OOB estimates
  
  if (xtype==0){ # use rfsrc_data1 and rfsrc_data0 to impute observed survival probability
    #Step 2: calculate D1 and D2
    # D1=mu1(1)-mu0(1)
    predict_rfsrc_median01 <- predict(rfsrc_data0, newdata = data1[,which(names(data1)%in%c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"))],failure.times=rep(time.interest,dim(data1)[1]), prediction.times="time")$predictions # rep(time.interest,dim(data1)[1])
    predict_rfsrc_median11 <- predict(rfsrc_data1, failure.times=rep(time.interest,dim(data1)[1]), prediction.times="time")$predictions
    data1$d1<- predict_rfsrc_median11-predict_rfsrc_median01
    
    # D0=mu1(0)-mu0(0)
    predict_rfsrc_median00 <- predict(rfsrc_data0, failure.times=rep(time.interest,dim(data0)[1]), prediction.times="time")$predictions
    predict_rfsrc_median10 <- predict(rfsrc_data1, newdata=data0[,which(names(data0)%in%c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"))], failure.times=rep(time.interest,dim(data0)[1]), prediction.times="time")$predictions
    data0$d0<- predict_rfsrc_median10-predict_rfsrc_median00
    
  }
  
  if (xtype==2){ # use IPCW on complete cases
    # KM IPCW: 10 fold cv
    U <- pmin(data$Time, time.interest)                         # truncated follow-up time by t0
    if (IPCW.method == "KM") { # KM as IPCW method
      fold.id <- sample(rep(seq(k.folds), length = nrow(data)))
      C.hat <- rep(NA, length(fold.id))
      for (z in 1:k.folds) {
        c.fit <- survival::survfit(survival::Surv(data$Time[!fold.id == z], 1 - data$Event[!fold.id == z]) ~ 1)
        kmc <- summary(c.fit, times = U[fold.id == z])
        C.hat[fold.id == z] <- kmc$surv[match(U[fold.id == z], kmc$time)]
      }
    }
    
    # 
    data0_surv<-ifelse(data0$Time>=time.interest,1,ifelse(data0$Time<time.interest & data0$Event==1,0,NA))
    data1_surv<-ifelse(data1$Time>=time.interest,1,ifelse(data1$Time<time.interest & data1$Event==1,0,NA))
    
    idx.trt <- which(data$Treatment==1)
    C.hat.trt <- C.hat[idx.trt]
    C.hat.trt.cen <- C.hat.trt[!is.na(data1_surv)]
    
    idx.ctrl <- which(data$Treatment==0)
    C.hat.ctrl <- C.hat[idx.ctrl]
    C.hat.ctrl.cen <- C.hat.ctrl[!is.na(data0_surv)]
    
    predict_rfsrc_median01 <-  predict(rfsrc_data0, newdata = data1[,which(names(data1)%in%c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"))], failure.times=rep(time.interest,dim(data1)[1]),prediction.times="time")$predictions
    data1$d1<- data1_surv-predict_rfsrc_median01
    
    # D0=mu1(0)-mu0(0)
    predict_rfsrc_median10 <-  predict(rfsrc_data1, newdata=data0[,which(names(data0)%in%c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"))], failure.times=rep(time.interest,dim(data0)[1]), prediction.times="time")$predictions
    data0$d0<- predict_rfsrc_median10-data0_surv
    
  }
  
  
  #Step 3: model d0 and d1
  if (IPCW==T){ # xtype must equal to 2.
    if (pseudo_reg_separate_trt==T){
      if (pseudo_reg=="weighted_lasso"){
        modeld0 <- cv.glmnet(y=data0[!is.na(data0$d0),which(names(data0)%in%c("d0"))], x=as.matrix(data0[!is.na(data0$d0),which(names(data0)%in%c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"))]),weights = (1/C.hat.ctrl.cen),alpha=1,nlambda = 20) #family = gaussian(), SL.library = SL_methods
        modeld1 <- cv.glmnet(y=data1[!is.na(data1$d1),which(names(data1)%in%c("d1"))], x=as.matrix(data1[!is.na(data1$d1),which(names(data1)%in%c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"))]),weights = (1/C.hat.trt.cen),alpha=1,nlambda = 20) #family = gaussian(), SL.library = SL_methods
        
      }else if (pseudo_reg=="weighted_RF"){
        modeld0 <- grf::regression_forest(Y=data0[!is.na(data0$d0),which(names(data0)%in%c("d0"))], X=data0[!is.na(data0$d0),which(names(data0)%in%c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"))],sample.weights = (1/C.hat.ctrl.cen),num.trees=ntrees,mtry=sqrt(10),min.node.siz=15,num.threads=1) 
        modeld1 <- grf::regression_forest(Y=data1[!is.na(data1$d1),which(names(data1)%in%c("d1"))], X=data1[!is.na(data1$d1),which(names(data1)%in%c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"))],sample.weights = (1/C.hat.trt.cen),num.trees=ntrees,mtry=sqrt(10),min.node.siz=15,num.threads=1) 
        
      }
      
    }else if (pseudo_reg_separate_trt==F){
      data0$d <- data0$d0
      data1$d <- data1$d1
      data_speudo_reg <- rbind(data0[,which(names(data0)%in%c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10","d"))],data1[,which(names(data1)%in%c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10","d"))])
      C.hat.cen <- rbind(C.hat.ctrl.cen, C.hat.trt.cen)
      
      if (pseudo_reg=="weighted_lasso"){
        modeld <- cv.glmnet(y=data_speudo_reg[!is.na(data_speudo_reg$d),which(names(data_speudo_reg)%in%c("d"))], x=as.matrix(data_speudo_reg[!is.na(data_speudo_reg$d),which(names(data_speudo_reg)%in%c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"))]),weights = (1/C.hat.cen),alpha=1,nlambda = 20) #family = gaussian(), SL.library = SL_methods
        
      }else if (pseudo_reg=="weighted_RF"){
        modeld <- regression_forest(Y=data_speudo_reg[!is.na(data_speudo_reg$d),which(names(data_speudo_reg)%in%c("d"))], X=data_speudo_reg[!is.na(data_speudo_reg$d),which(names(data_speudo_reg)%in%c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"))],sample.weights = (1/C.hat.cen),num.trees=ntrees,mtry=sqrt(10),min.node.siz=15,num.threads=1) 

      }
      
    }
    
  }else if (IPCW==F){ # xtype must not equal to 2.
    if (pseudo_reg_separate_trt==T){
      if (pseudo_reg=="weighted_lasso"){
        modeld0 <- cv.glmnet(y=data0[!is.na(data0$d0),which(names(data0)%in%c("d0"))], x=as.matrix(data0[!is.na(data0$d0),which(names(data0)%in%c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"))]),alpha=1,nlambda = 20) #family = gaussian(), SL.library = SL_methods
        modeld1 <- cv.glmnet(y=data1[!is.na(data1$d1),which(names(data1)%in%c("d1"))], x=as.matrix(data1[!is.na(data1$d1),which(names(data1)%in%c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"))]),alpha=1,nlambda = 20) #family = gaussian(), SL.library = SL_methods
        
      }else if (pseudo_reg=="weighted_RF"){
        modeld0 <- grf::regression_forest(Y=data0[!is.na(data0$d0),which(names(data0)%in%c("d0"))], X=data0[!is.na(data0$d0),which(names(data0)%in%c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"))],num.trees=ntrees,mtry=sqrt(10),min.node.siz=15,num.threads=1) 
        modeld1 <- grf::regression_forest(Y=data1[!is.na(data1$d1),which(names(data1)%in%c("d1"))], X=data1[!is.na(data1$d1),which(names(data1)%in%c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"))],num.trees=ntrees,mtry=sqrt(10),min.node.siz=15,num.threads=1) 
        
      }
      
    }else if (pseudo_reg_separate_trt==F){
      data0$d <- data0$d0
      data1$d <- data1$d1
      data_speudo_reg <- rbind(data0[,which(names(data0)%in%c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10","d"))],data1[,which(names(data1)%in%c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10","d"))])
      
      if (pseudo_reg=="weighted_lasso"){
        modeld <- cv.glmnet(y=data_speudo_reg[!is.na(data_speudo_reg$d),which(names(data_speudo_reg)%in%c("d"))], x=as.matrix(data_speudo_reg[!is.na(data_speudo_reg$d),which(names(data_speudo_reg)%in%c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"))]),alpha=1,nlambda = 20) #family = gaussian(), SL.library = SL_methods
        
      }else if (pseudo_reg=="weighted_RF"){
        modeld <- grf::regression_forest(Y=data_speudo_reg[!is.na(data_speudo_reg$d),which(names(data_speudo_reg)%in%c("d"))], X=data_speudo_reg[!is.na(data_speudo_reg$d),which(names(data_speudo_reg)%in%c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"))],num.trees=ntrees,mtry=sqrt(10),min.node.siz=15,num.threads=1) 
        
      }
      
    }
    
  }
  
  
  ### predict d0 and d1 or d on test data 
  if (pseudo_reg_separate_trt==T) {
    
    if (pseudo_reg=="weighted_lasso"){
      predictiond0 <- predict(modeld0, as.matrix(testdat[,which(names(testdat)%in%c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"))]), type = "response") # This predict() function is for glmnet
      predictiond1 <- predict(modeld1, as.matrix(testdat[,which(names(testdat)%in%c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"))]), type = "response")
    }else if(pseudo_reg=="weighted_RF"){
      predictiond0 <- predict(modeld0, testdat[,which(names(testdat)%in%c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"))], estimate.variance = F)$predictions
      predictiond1 <- predict(modeld1, testdat[,which(names(testdat)%in%c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"))], estimate.variance = F)$predictions
    }
    
    if (est_pi==1 ){
  
      if (propensity_method=="correct_logistic"){ # not cv version
        rf <- glm(Treatment~., family = "binomial"(link = "logit"),data = data[,which(names(data)%in%c(propensity_method_X,"Treatment"))])
        pred <- predict(rf,newdata=testdat[,which(names(testdat)%in%c(propensity_method_X))],type='response')
      }else if(propensity_method=="regression_forest"){ # throw in all X's in random forest; OOB prediction on training data
        rf<-grf::probability_forest(X=as.matrix(data[,which(names(data)%in%c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"))]),Y=as.factor(as.vector(data$Treatment)),num.trees = ntrees,num.threads = 1)
        pred<-predict(rf,testdat[,which(names(testdat)%in%c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"))])$predictions[,2]
      }
      
    }else if (est_pi==0 ){
      pred <- rep(propensity,dim(testdat)[1])
    }
    
    X.diff<-pred*predictiond0+(1-pred)*predictiond1
    
  }else if (pseudo_reg_separate_trt==F){
    if (pseudo_reg=="weighted_lasso"){
      predictiond <- predict(modeld, as.matrix(testdat[,which(names(testdat)%in%c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"))]), type = "response") # This predict() function is for glmnet
    }else if(pseudo_reg=="weighted_RF"){
      predictiond <- predict(modeld, testdat[,which(names(testdat)%in%c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"))], estimate.variance = F)$predictions
    }
    
    X.diff <- predictiond
  }  
  
  ## kernelshap
  if (pseudo_reg_separate_trt==T) {
  system.time(
    shap_d0 <- kernelshap(modeld0, as.matrix(testdat_kernelshap[which(testdat_kernelshap$Treatment==0),which(names(testdat_kernelshap)%in%c("V1", "V2","V3","V4", "V5","V6", "V7", "V8","V9", "V10"))]), bg_X = as.matrix(kernelshap_bg[kernelshap_bg$Treatment==0,which(names(kernelshap_bg)%in%c("V1", "V2","V3","V4", "V5","V6", "V7", "V8","V9", "V10"))]))
  )
  shapvip_d0 <- shap_d0$S
  
  system.time(
    shap_d1 <- kernelshap(modeld1, as.matrix(testdat_kernelshap[which(testdat_kernelshap$Treatment==1),which(names(testdat_kernelshap)%in%c("V1", "V2","V3","V4", "V5","V6", "V7", "V8","V9", "V10"))]), bg_X = as.matrix(kernelshap_bg[kernelshap_bg$Treatment==1,which(names(kernelshap_bg)%in%c("V1", "V2","V3","V4", "V5","V6", "V7", "V8","V9", "V10"))]))
  )
  shapvip_d1 <- shap_d1$S 
  
  }
  
  shapvip_d <- rbind(shapvip_d0,shapvip_d1)
  
  return(list(diff=X.diff, propensity=pred, pred.s1=predictiond1, pred.s0=predictiond0,shapvip_d0 =shapvip_d0,shapvip_d1 =shapvip_d1,shapvip=shapvip_d)) 
  
}

#------------------------------------------------ M-learner ---------------------------------------------------#
rsf_HTE_M <- function(data=dat,testdat=mydata$data,testdat_kernelshap=mydata_kernelshap$data,kernelshap_bg=train_kernelshap_bg,ntrees=1000,time.interest=12,IPCW=T,IPCW.method="KM",k.folds=10,impute_type=2,est_pi=0, propensity=0.5,propensity_method="correct_logistic",propensity_method_X=c("V1","V5"),pseudo_reg="weighted_RF"){ 
  data<-data[complete.cases(data),]

  #### fit propensity model of treatment assignment using out of bag prediction
  ### random forest
  if (est_pi==1){
    
    if (propensity_method=="correct_logistic"){ # not cv version
      rf <- glm(Treatment~., family = "binomial"(link = "logit"),data = data[,which(names(data)%in%c(propensity_method_X,"Treatment"))])
      data$pred.pi <- predict(rf,type='response')
    }else if(propensity_method=="regression_forest"){ # throw in all X's in random forest; OOB prediction on training data
      rf<-grf::probability_forest(X=as.matrix(data[,which(names(data)%in%c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"))]),Y=as.factor(as.vector(data$Treatment)),num.trees = ntrees,num.threads = 1)
      data$pred.pi<-predict(rf)$predictions[,2] 
    }
    
  }else if (est_pi==0){
    data$pred.pi <- rep(propensity,dim(data)[1])
  }
  
  #### modified outcome
  if (impute_type==2){ # > time.interest, 1; < time.interest, 0; ipcw on complete cases
    # KM IPCW: 10 fold cv
    U <- pmin(data$Time, time.interest)                         # truncated follow-up time by t0
    if (IPCW.method == "KM") { # KM as IPCW method
      fold.id <- sample(rep(seq(k.folds), length = nrow(data)))
      C.hat <- rep(NA, length(fold.id))
      for (z in 1:k.folds) {
        c.fit <- survival::survfit(survival::Surv(data$Time[!fold.id == z], 1 - data$Event[!fold.id == z]) ~ 1)
        kmc <- summary(c.fit, times = U[fold.id == z])
        C.hat[fold.id == z] <- kmc$surv[match(U[fold.id == z], kmc$time)]
      }
    }
    
    data_surv<-ifelse(data$Time>=time.interest,1,ifelse(data$Time<time.interest & data$Event==1,0,NA))
    C.hat.cen <- C.hat[!is.na(data_surv)]
    
    if (stable.propensity==T){
      num.pred <- sum(data$Treatment)/length(data$Treatment)
      Y_M <- data_surv * (data$Treatment*(num.pred/data$pred.pi) - (1-data$Treatment)*((1-num.pred)/(1-data$pred.pi))) 
    }else if (stable.propensity==F){
      if (propensity_method !="EB_cal" ){
        Y_M <- data_surv * (data$Treatment/data$pred.pi - (1-data$Treatment)/(1-data$pred.pi)) 
      }else{
        Y_M <- data_surv * (data$Treatment*wi - (1-data$Treatment)*(1-wi)) 
      }
      
    }
    
  }
  
  if (IPCW==T){
    if (pseudo_reg=="weighted_lasso"){
      modeld <- cv.glmnet(y=Y_M[!is.na(Y_M)], x = as.matrix(data[!is.na(Y_M),which(names(data)%in%c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"))]), weights = (1/C.hat.cen),alpha=1,nlambda = 20)
    }else if (pseudo_reg=="weighted_RF"){
      modeld <- grf::regression_forest(Y=Y_M[!is.na(Y_M)], X=data[!is.na(Y_M),which(names(data)%in%c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"))],sample.weights = (1/C.hat.cen),num.trees=ntrees,mtry=sqrt(10),min.node.siz=15,num.threads=1) 
    }
  }
  
  ##### predict on test data
  if(pseudo_reg=="weighted_lasso"){
    predictiond<-predict(modeld,as.matrix(testdat[,which(names(testdat)%in%c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"))]),type="response")
  }else if (pseudo_reg=="weighted_RF"){
    predictiond <- predict(modeld, testdat[,which(names(testdat)%in%c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"))], estimate.variance = F)$predictions
  }
  
  #### kernelshap
  system.time(
    shap_d <- kernelshap(modeld, as.matrix(testdat_kernelshap[,which(names(testdat_kernelshap)%in%c("V1", "V2","V3","V4", "V5","V6", "V7", "V8","V9", "V10"))]), bg_X = as.matrix(kernelshap_bg[,which(names(kernelshap_bg)%in%c("V1", "V2","V3","V4", "V5","V6", "V7", "V8","V9", "V10"))]))
  )
  shapvip_d <- shap_d$S

  return(list(diff=predictiond,propensity=data$pred.pi,shapvip=shapvip_d )) 
  
}


#------------------------------------------------ DR-learner ---------------------------------------------------#
rsf_HTE_DR_Kennedy <- function(data=dat,testdat=mydata$data,testdat_kernelshap=mydata_kernelshap$data,kernelshap_bg=train_kernelshap_bg,ntrees=1000,time.interest=15,impute_type=2,IPCW=T,IPCW.method = "KM",k.folds=10,est_pi=1,propensity_method="correct_logistic",propensity_method_X=c("V1","V5"),propensity=0.5,pseudo_reg="weighted_RF"){
  data<-data[complete.cases(data),]
  data1 <- data[which(data$Treatment==1),-c(which(names(data)=="Treatment"), which(names(data)=="sim"))]
  data0 <- data[which(data$Treatment==0),-c(which(names(data)=="Treatment"), which(names(data)=="sim"))]
  
  #Step1: Fit 2 RFSRC models; OOB estimation
    rfsrc_data1 <-
      grf::survival_forest(X=as.matrix(data1[,which(names(data1) %in% c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"))]),Y=as.vector(data1$Time),D=as.vector(data1$Event),prediction.type="Nelson-Aalen",num.trees = 500,num.threads = 1) # OOB estimates
    
    rfsrc_data0 <-
      grf::survival_forest(X=as.matrix(data0[,which(names(data0) %in% c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"))]),Y=as.vector(data0$Time),D=as.vector(data0$Event),prediction.type="Nelson-Aalen",num.trees = 500,num.threads = 1) # OOB estimates


  #Step 2: construct pseudo outcome
  if (est_pi==1){
    
    if (propensity_method=="correct_logistic"){ # not cv version
      rf <- glm(Treatment~., family = "binomial"(link = "logit"),data = data[,which(names(data)%in%c(propensity_method_X,"Treatment"))])
      data$pi.hat <- predict(rf,type='response')
    }else if(propensity_method=="regression_forest"){ # throw in all X's in random forest; OOB prediction on training data
      rf<-grf::probability_forest(X=as.matrix(data[,which(names(data)%in%c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"))]),Y=as.factor(as.vector(data$Treatment)),num.trees = ntrees,num.threads = 1)
      data$pi.hat<-predict(rf)$predictions[,2]
    }
    
  }else if (est_pi==0){
    data$pi.hat <- rep(propensity,dim(data)[1])
  }
  
  data$predict_rfsrc_median0 <- predict(rfsrc_data0, newdata=data[,which(names(data)%in%c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"))], failure.times=rep(time.interest,dim(data)[1]), prediction.times="time")$predictions 
  data$predict_rfsrc_median1 <- predict(rfsrc_data1, newdata=data[,which(names(data)%in%c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"))], failure.times=rep(time.interest,dim(data)[1]), prediction.times="time")$predictions 
  
  data$I<-ifelse(data$Time>=time.interest,1,ifelse(data$Time<time.interest & data$Event==1,0,NA))
  
  data$pseudo <-(data$Treatment*wi+(1-data$Treatment)*(-1*wi))*(data$I-data$Treatment*data$predict_rfsrc_median1-(1-data$Treatment)*data$predict_rfsrc_median0) + data$predict_rfsrc_median1-data$predict_rfsrc_median0 
  
  # step3: regress pseudo outcome on covariates
  if (impute_type==2){
    # KM IPCW: 10 fold cv
    U <- pmin(data$Time, time.interest)                         # truncated follow-up time by t0
    if (IPCW.method == "KM") { # KM as IPCW method
      fold.id <- sample(rep(seq(k.folds), length = nrow(data)))
      C.hat <- rep(NA, length(fold.id))
      for (z in 1:k.folds) {
        c.fit <- survival::survfit(survival::Surv(data$Time[!fold.id == z], 1 - data$Event[!fold.id == z]) ~ 1)
        kmc <- summary(c.fit, times = U[fold.id == z])
        C.hat[fold.id == z] <- kmc$surv[match(U[fold.id == z], kmc$time)]
      }
    }
    
    C.hat.cen <- C.hat[!is.na(data$I)]
  }
  
  if (IPCW==T){
    if (pseudo_reg=="weighted_lasso"){
      modeld <- cv.glmnet(y=data$pseudo[!is.na(data$pseudo)], x = as.matrix(data[!is.na(data$pseudo),which(names(data)%in%c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"))]), weights = (1/C.hat.cen),alpha=1,nlambda=20)
    }else if (pseudo_reg=="weighted_RF"){
      modeld <- grf::regression_forest(Y=data$pseudo[!is.na(data$pseudo)], X=data[!is.na(data$pseudo),which(names(data)%in%c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"))],sample.weights = (1/C.hat.cen),num.trees=ntrees,mtry=sqrt(10),min.node.siz=15,num.threads=1) 
    }
  }
  
  # step4: predict on test data
  if(pseudo_reg=="weighted_lasso"){
    predictiond<-predict(modeld,as.matrix(testdat[,which(names(testdat)%in%c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"))]),type="response")
  }else if (pseudo_reg=="weighted_RF"){
    predictiond <- predict(modeld, testdat[,which(names(testdat)%in%c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"))], estimate.variance = F)$predictions
  }
  
  #### kernelshap
  system.time(
    shap_d <- kernelshap(modeld, as.matrix(testdat_kernelshap[,which(names(testdat_kernelshap)%in%c("V1", "V2","V3","V4", "V5","V6", "V7", "V8","V9", "V10"))]), bg_X = as.matrix(kernelshap_bg[,which(names(kernelshap_bg)%in%c("V1", "V2","V3","V4", "V5","V6", "V7", "V8","V9", "V10"))]))
  )
  shapvip_d <- shap_d$S
  
  return(list(diff=predictiond,propensity=data$pihat,shapvip=shapvip_d)) 
}

#-------------------------------------------- R-learner -----------------------------------------------#
rsf_HTE_R <- function(data=dat,testdat=mydata$data,testdat_kernelshap=mydata_kernelshap$data,kernelshap_bg=train_kernelshap_bg,ntrees=1000,time.interest=15,impute_type=2,IPCW=T,k.folds=10,IPCW.method="KM",est_pi=1,propensity_method="correct_logistic",propensity_method_X=c("V1","V5"),propensity=0.5,pseudo_reg="weighted_RF"){
  data<-data[complete.cases(data),]
  data1 <- data[which(data$Treatment==1),-c(which(names(data)=="Treatment"), which(names(data)=="sim"))]
  data0 <- data[which(data$Treatment==0),-c(which(names(data)=="Treatment"), which(names(data)=="sim"))]
  
  #### fit propensity model of treatment assignment using out of bag prediction on training data
  ### random forest
  if (est_pi==1){
    
    if (propensity_method=="correct_logistic"){ # not cv version
      rf <- glm(Treatment~., family = "binomial"(link = "logit"),data = data[,which(names(data)%in%c(propensity_method_X,"Treatment"))])
      data$pred.pi <- predict(rf,type='response')
    }else if(propensity_method=="regression_forest"){ # throw in all X's in random forest; OOB prediction on training data
      rf<-grf::probability_forest(X=as.matrix(data[,which(names(data)%in%c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"))]),Y=as.factor(as.vector(data$Treatment)),num.trees = ntrees,num.threads = 1)
      data$pred.pi<-predict(rf)$predictions[,2]
    }
    
  }else if(est_pi==0){
    data$pred.pi <- rep(propensity, dim(data)[1])
  }

  #### Fit 2 RFSRC models
  rfsrc_data1 <-
    grf::survival_forest(X=as.matrix(data1[,which(names(data1) %in% c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"))]),Y=as.vector(data1$Time),D=as.vector(data1$Event),prediction.type="Nelson-Aalen",num.trees = 500,num.threads = 1) # OOB estimates
  
  rfsrc_data0 <-
    grf::survival_forest(X=as.matrix(data0[,which(names(data0) %in% c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"))]),Y=as.vector(data0$Time),D=as.vector(data0$Event),prediction.type="Nelson-Aalen",num.trees = 500,num.threads = 1) # OOB estimates
  
  #### construct m
  data$s0_hat <- predict(rfsrc_data0, newdata=data[,which(names(data)%in%c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"))], failure.times=rep(time.interest,dim(data)[1]), prediction.times="time")$predictions # summary of s0_hat on training data: 0.6 to 0.8
  data$s1_hat <- predict(rfsrc_data1, newdata=data[,which(names(data)%in%c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"))], failure.times=rep(time.interest,dim(data)[1]), prediction.times="time")$predictions # summary of s0_hat on training data: 0.5 to 0.9
  data$m_hat <- data$pred.pi*data$s1_hat+(1-data$pred.pi)*data$s0_hat

  #### create pseudo outcome
  data$I<-ifelse(data$Time>=time.interest,1,ifelse(data$Time<time.interest & data$Event==1,0,NA))
  data$pseudo <- ((data$I) - data$m_hat) / (data$Treatment - data$pred.pi)
  
  #### create censoring distribution
  if (impute_type==2){
    # KM IPCW: 10 fold cv
    U <- pmin(data$Time, time.interest)                         # truncated follow-up time by t0
    if (IPCW.method == "KM") { # KM as IPCW method
      fold.id <- sample(rep(seq(k.folds), length = nrow(data)))
      C.hat <- rep(NA, length(fold.id))
      for (z in 1:k.folds) {
        c.fit <- survival::survfit(survival::Surv(data$Time[!fold.id == z], 1 - data$Event[!fold.id == z]) ~ 1)
        kmc <- summary(c.fit, times = U[fold.id == z])
        C.hat[fold.id == z] <- kmc$surv[match(U[fold.id == z], kmc$time)]
      }
    }

  }

  #### weight
  weight <- (1/C.hat)*((data$Treatment - data$pred.pi)^2)
  weight.cen <-  weight[!is.na(data$I)]
  
  #### solve objective function
  if (IPCW==T){
    if (pseudo_reg=="weighted_lasso"){
      modeld <- cv.glmnet(y=data$pseudo[!is.na(data$pseudo)], x = as.matrix(data[!is.na(data$pseudo),which(names(data)%in%c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"))]), weights = weight.cen,alpha=1,nlambda=20)
    }else if (pseudo_reg=="weighted_RF"){
      modeld <- grf::regression_forest(Y=data$pseudo[!is.na(data$pseudo)], X=data[!is.na(data$pseudo),which(names(data)%in%c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"))],sample.weights = weight.cen,num.trees=ntrees,mtry=sqrt(10),min.node.siz=15,num.threads=1) 
    }
  }
  
  ### predict on test data
  if(pseudo_reg=="weighted_lasso"){
    predictiond<-predict(modeld,as.matrix(testdat[,which(names(testdat)%in%c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"))]),type="response")
  }else if (pseudo_reg=="weighted_RF"){
    predictiond <- predict(modeld, testdat[,which(names(testdat)%in%c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"))], estimate.variance = F)$predictions
  }
  
  #### kernelshap
  system.time(
    shap_d <- kernelshap(modeld, as.matrix(testdat_kernelshap[,which(names(testdat_kernelshap)%in%c("V1", "V2","V3","V4", "V5","V6", "V7", "V8","V9", "V10"))]), bg_X = as.matrix(kernelshap_bg[,which(names(kernelshap_bg)%in%c("V1", "V2","V3","V4", "V5","V6", "V7", "V8","V9", "V10"))]))
  )
  shapvip_d <- shap_d$S
  
  return(list(diff=predictiond,propensity=data$pred.pi,shapvip=shapvip_d)) 
}



#------------------------------------------------------------ D-learner and DEA-learner --------------------------------------------------------------#
rsf_HTE_D <- function(data=dat,testdat=mydata$data,testdat_kernelshap=mydata_kernelshap$data,kernelshap_bg=train_kernelshap_bg,ntrees=1000,time.interest=15,impute_type=2,IPCW=T,IPCW.method ="KM",k.folds=10,EA=F,est_pi=1,propensity_method="correct_logistic",propensity_method_X=c("V1","V5"),propensity=0.5,pseudo_reg="weighted_RF"){
  data<-data[complete.cases(data),]
  
  #### fit propensity model of treatment assignment using out of bag prediction on training data
  ### random forest
  if (est_pi==1){
    
    if (propensity_method=="correct_logistic"){ # not cv version
      rf <- glm(Treatment~., family = "binomial"(link = "logit"),data = data[,which(names(data)%in%c(propensity_method_X,"Treatment"))])
      data$pred.pi <- predict(rf,type='response')
    }else if(propensity_method=="regression_forest"){ # throw in all X's in random forest; OOB prediction on training data
      rf<-grf::probability_forest(X=as.matrix(data[,which(names(data)%in%c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"))]),Y=as.factor(as.vector(data$Treatment)),num.trees = ntrees,num.threads = 1)
      data$pred.pi<-predict(rf)$predictions[,2]
    }
    
    }else if(est_pi==0){
    data$pred.pi <- rep(propensity,dim(data)[1])
  }

  #### create m_hat
  rfsrc_data <-
    grf::survival_forest(X=as.matrix(data[,which(names(data) %in% c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"))]),Y=as.vector(data$Time),D=as.vector(data$Event),prediction.type="Nelson-Aalen",num.trees = 500,num.threads = 1) # OOB estimates
  data$m_hat <- predict(rfsrc_data,failure.times=rep(time.interest,dim(data)[1]),prediction.times="time")$predictions 
  
  #### create pseudo outcome
  data$I<-ifelse(data$Time>=time.interest,1,ifelse(data$Time<time.interest & data$Event==1,0,NA))
  if (EA==F){
    data$pseudo <- 2*(2*data$Treatment-1)*data$I 
  }else if (EA==T){
    data$pseudo <- 2*(2*data$Treatment-1)*(data$I-data$m_hat)
  }
  
  #### create censoring distribution
  if (impute_type==2){
    # KM IPCW: 10 fold cv
    U <- pmin(data$Time, time.interest)                         # truncated follow-up time by t0
    if (IPCW.method == "KM") { # KM as IPCW method
      fold.id <- sample(rep(seq(k.folds), length = nrow(data)))
      C.hat <- rep(NA, length(fold.id))
      for (z in 1:k.folds) {
        c.fit <- survival::survfit(survival::Surv(data$Time[!fold.id == z], 1 - data$Event[!fold.id == z]) ~ 1)
        kmc <- summary(c.fit, times = U[fold.id == z])
        C.hat[fold.id == z] <- kmc$surv[match(U[fold.id == z], kmc$time)]
      }
    }
    
  }
  
  #### weight
  weight <- (1/C.hat)*(2*data$Treatment-1)*(data$Treatment*(1/4)*wi+(1-data$Treatment)*(-1/4)*wi)
  weight.cen <- weight[!is.na(data$I)]

  #### solve objective function
  if (IPCW==T){
    if (pseudo_reg=="weighted_lasso"){
      modeld <- cv.glmnet(y=data$pseudo[!is.na(data$pseudo)], x = as.matrix(data[!is.na(data$pseudo),which(names(data)%in%c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"))]), weights = weight.cen,alpha=1,nlambda=20)
    }else if (pseudo_reg=="weighted_RF"){
      modeld <- grf::regression_forest(Y=data$pseudo[!is.na(data$pseudo)], X=data[!is.na(data$pseudo),which(names(data)%in%c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"))],sample.weights = weight.cen,num.trees=ntrees,mtry=sqrt(10),min.node.siz=15,num.threads=1) 
    }
    
  }

  if(pseudo_reg=="weighted_lasso"){
    predictiond<-predict(modeld,as.matrix(testdat[,which(names(testdat)%in%c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"))]),type="response")
  }else if (pseudo_reg=="weighted_RF"){
    predictiond <- predict(modeld, testdat[,which(names(testdat)%in%c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10"))], estimate.variance = F)$predictions
  }
  
  #### kernelshap
  system.time(
    shap_d <- kernelshap(modeld, as.matrix(testdat_kernelshap[,which(names(testdat_kernelshap)%in%c("V1", "V2","V3","V4", "V5","V6", "V7", "V8","V9", "V10"))]), bg_X = as.matrix(kernelshap_bg[,which(names(kernelshap_bg)%in%c("V1", "V2","V3","V4", "V5","V6", "V7", "V8","V9", "V10"))]))
  )
  shapvip_d <- shap_d$S
  
  return(list(diff=predictiond, propensity=data$pred.pi,shapvip=shapvip_d)) 
}



