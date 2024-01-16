weibull_true <- function(data=dat,testdat=testdat,time.interest=NULL,scenario=1){
  data<-data[complete.cases(data),]
  data1 <- data[which(data$Treatment==1),-c(which(names(data)=="Treatment"))]
  data0 <- data[which(data$Treatment==0),-c(which(names(data)=="Treatment"))]
  
  # Fit Weibull model
  if (scenario==1){
    weib0 <- survreg(Surv(data0$Time, data0$Event)~data0$V1+data0$V3+data0$V6+data0$V7+data0$V10,dist="weibull")
    sigma0<-summary(weib0)$scale
    mu0<-summary(weib0)$coef[1]
    lambda0<-exp(mu0)
    eta0<-1/sigma0
    beta0<-summary(weib0)$coef[-1]*(-1)/sigma0
    
    weib1 <- survreg(Surv(data1$Time, data1$Event)~data1$V1+data1$V2+data1$V3+data1$V5+data1$V6+data1$V7+data1$V8+data1$V10,dist="weibull")
    sigma1<-summary(weib1)$scale
    mu1<-summary(weib1)$coef[1]
    lambda1<-exp(mu1)
    eta1<-1/sigma1
    beta1<-summary(weib1)$coef[-1]*(-1)/sigma1
    # Predict using counterfactual data
    time.interest1<-time.interest0<-time.interest
    predict_weib1 <- exp(-exp(as.matrix(testdat[,c(match(c("V1","V2","V3","V5","V6","V7","V8","V10"),names(testdat)))])%*%beta1)*(time.interest/lambda1)^eta1)
    predict_weib0 <- exp(-exp(as.matrix(testdat[,c(match(c("V1","V3","V6","V7","V10"),names(testdat)))])%*%beta0)*(time.interest/lambda0)^eta0)
  }
  
  if (scenario==2){
    data0$expv1<-exp(data0$V1)
    data0$sqv2 <- (data0$V2)^2
    data0$sinv3 <- sin(data0$V3)
    data0$sinv6 <- sin(data0$V6)
    data0$v1v8 <- data0$V1*data0$V8
    weib0 <- survreg(Surv(data0$Time, data0$Event)~data0$expv1+data0$sqv2+data0$sinv3+data0$V5+data0$sinv6+data0$V7+data0$v1v8+data0$V10,dist="weibull")
    sigma0<-summary(weib0)$scale
    mu0<-summary(weib0)$coef[1]
    lambda0<-exp(mu0)
    eta0<-1/sigma0
    beta0<-summary(weib0)$coef[-1]*(-1)/sigma0
    
    data1$expv1<-exp(data1$V1)
    data1$sqv2 <- (data1$V2)^2
    data1$sinv3 <- sin(data1$V3)
    data1$sinv6 <- sin(data1$V6)
    data1$v1v8 <- data1$V1*data1$V8
    weib1 <- survreg(Surv(data1$Time, data1$Event)~data1$expv1+data1$sqv2+data1$sinv3+data1$V5+data1$sinv6+data1$V7+data1$v1v8+data1$V10+data1$V2+data1$V8,dist="weibull")
    sigma1<-summary(weib1)$scale
    mu1<-summary(weib1)$coef[1]
    lambda1<-exp(mu1)
    eta1<-1/sigma1
    beta1<-summary(weib1)$coef[-1]*(-1)/sigma1
    # Predict using counterfactual data
    time.interest1<-time.interest0<-time.interest
    testdat$expv1<-exp(testdat$V1)
    testdat$sqv2 <- (testdat$V2)^2
    testdat$sinv3 <- sin(testdat$V3)
    testdat$sinv6 <- sin(testdat$V6)
    testdat$v1v8 <- testdat$V1*testdat$V8
    predict_weib1 <- exp(-exp(as.matrix(testdat[,c(match(c("expv1","sqv2","sinv3","V5","sinv6","V7","v1v8","V10","V2","V8"),names(testdat)))])%*%beta1)*(time.interest/lambda1)^eta1)
    predict_weib0 <- exp(-exp(as.matrix(testdat[,c(match(c("expv1","sqv2","sinv3","V5","sinv6","V7","v1v8","V10"),names(testdat)))])%*%beta0)*(time.interest/lambda0)^eta0)
  }
  
  if (scenario==3){
    data0$expv1<-exp(data0$V1)
    data0$sqv2<-data0$V2^2
    data0$sinv3 <- sin(data0$V3)
    data0$sinv6 <- sin(data0$V6)
    data0$v1v8 <- data0$V1*data0$V8
    weib0 <- survreg(Surv(data0$Time, data0$Event)~data0$expv1+data0$sqv2+data0$sinv3+data0$V5+data0$sinv6+data0$V7+data0$v1v8+data0$V10,dist="weibull")
    sigma0<-summary(weib0)$scale
    mu0<-summary(weib0)$coef[1]
    lambda0<-exp(mu0)
    eta0<-1/sigma0
    beta0<-summary(weib0)$coef[-1]*(-1)/sigma0
    
    data1$expv1<-exp(data1$V1)
    data1$sqv2<-data1$V2^2
    data1$sinv3 <- sin(data1$V3)
    data1$sinv6 <- sin(data1$V6)
    data1$v1v8 <- data1$V1*data1$V8
    data1$sqv1 <- (data1$V1)^2
    data1$v2v3 <- data1$V2*data1$V3
    data1$expv5 <- exp(data1$V5)
    data1$v3v8 <- data1$V3*data1$V8
    weib1 <- survreg(Surv(data1$Time, data1$Event)~data1$expv1+data1$sqv2+data1$sinv3+data1$V5+data1$sinv6+data1$V7+data1$v1v8+data1$V10+data1$sqv1+data1$v2v3+data1$expv5+data1$v3v8,dist="weibull")
    sigma1<-summary(weib1)$scale
    mu1<-summary(weib1)$coef[1]
    lambda1<-exp(mu1)
    eta1<-1/sigma1
    beta1<-summary(weib1)$coef[-1]*(-1)/sigma1
    # Predict using counterfactual data
    time.interest1<-time.interest0<-time.interest
    testdat$expv1<-exp(testdat$V1)
    testdat$sqv2<-testdat$V2^2
    testdat$sinv3 <- sin(testdat$V3)
    testdat$sinv6 <- sin(testdat$V6)
    testdat$v1v8 <- testdat$V1*testdat$V8
    testdat$sqv1 <- (testdat$V1)^2
    testdat$v2v3 <- testdat$V2*testdat$V3
    testdat$expv5 <- exp(testdat$V5)
    testdat$v3v8 <- testdat$V3*testdat$V8
    predict_weib1 <- exp(-exp(as.matrix(testdat[,c(match(c("expv1","sqv2","sinv3","V5","sinv6","V7","v1v8","V10","sqv1","v2v3","expv5","v3v8"),names(testdat)))])%*%beta1)*(time.interest/lambda1)^eta1)
    predict_weib0 <- exp(-exp(as.matrix(testdat[,c(match(c("expv1","sqv2","sinv3","V5","sinv6","V7","v1v8","V10"),names(testdat)))])%*%beta0)*(time.interest/lambda0)^eta0)
  }
  
  
  weib.diff <- predict_weib1-predict_weib0
  weib.ratio <- predict_weib1/predict_weib0
  return(list(diff=weib.diff,ratio=weib.ratio,pred0=predict_weib0,pred1=predict_weib1)) 
}
