#cutoff<-c(80,90,100,110,120,130,140,150,160,1000)

library(ggplot2)
library(dplyr)
library(data.table)
library(survival)
library(DiagrammeR)
library(regmedint)

#Set Up
p<<-1.4
#outcome regression AFT
theta0<<-7.5
theta1<<-0.12
theta2<<-0.015#should be negative
theta3<<--0.023
theta4<<--0.02

theta0_ni<<-7.7
theta1_ni<<-0.59
theta2_ni<<-0.0032 #should be negative
cutoff<-c(90,110,130,150,1000)



#invastigate if the adjustment works or not using simulation



n<-50
# Simulation:

# Method: 

rd11<-NULL
set.seed(1)
for(j in 1:n){
  
  tryCatch({
    print(paste("1",j))
    e1.te2<-NULL
    e1.te<-NULL
    e1.pnde<-NULL
    e1.tnde<-NULL
    e1.pnie<-NULL
    e1.tnie<-NULL
    e1.pm<-NULL
    e1.beta0<-NULL
    e1.beta1<-NULL
    e1.beta2<-NULL
    e1.sigma<-NULL
    e1.theta1<-NULL
    e1.theta2<-NULL
    e1.theta3<-NULL
    e1.theta4<-NULL
    
    
    e2.te2<-NULL
    e2.te<-NULL
    e2.pnde<-NULL
    e2.tnde<-NULL
    e2.pnie<-NULL
    e2.tnie<-NULL
    e2.pm<-NULL
    e2.beta0<-NULL
    e2.beta1<-NULL
    e2.beta2<-NULL
    e2.sigma<-NULL
    e2.theta1<-NULL
    e2.theta2<-NULL
    e2.theta3<-NULL
    e2.theta4<-NULL
    
    e3.te2<-NULL
    e3.te<-NULL
    e3.pnde<-NULL
    e3.tnde<-NULL
    e3.pnie<-NULL
    e3.tnie<-NULL
    e3.pm<-NULL
    e3.beta0<-NULL
    e3.beta1<-NULL
    e3.beta2<-NULL
    e3.sigma<-NULL
    e3.theta1<-NULL
    e3.theta2<-NULL
    e3.theta3<-NULL
    e3.theta4<-NULL
    
    cutoff<-c(90,110,130,150,1000)
    for(i in 1:5){
      data<-simulation.function.v42(20000,cutoff = cutoff[i])
      data$ldlc_change<-data$ldlc12 - data$ldlcb
      data$ldlc_change_m<-data$ldlc12m - data$ldlcm
      fit<-lm(data,formula = ldlc_change_m~ldlc_change+ldlcb+drug)
      data$ldlc_change_m_pred<-predict(fit,data)
      #data$cen<-rep(1,20000)
      #data$d<-rep(0,20000)
      #data$logldlcb<-log(data$ldlcb)
      z<-1/(1+exp(-theta0+theta1*data$drug-
                        theta2*data$ldlc_change_m-
                        theta3*data$ldlc_change_m*data$drug-
                        theta4*data$ldlcb))
      pr = 1/(1+exp(-z))         # pass through an inv-logit function
      data$y1 = rbinom(20000,1,pr/10)      # bernoulli response variable(rear outcome)
      
      
      c_median<-summary(data$ldlcb)[3]
      c_mean<-summary(data$ldlcb)[4]
      fit1<-regmedint(data = data,
                      ## Variables
                      yvar = "y1",
                      avar = "drug",
                      mvar = "ldlc_change",
                      cvar = c("ldlcb"),
                      #eventvar = "cen",
                      ## Values at which effects are evaluated
                      a0 = 0,
                      a1 = 1,
                      m_cde = 0,
                      c_cond = 120,
                      ## Model types
                      mreg = "linear",
                      #yreg = "survCox",
                      yreg = "logistic",
                      ## Additional specification
                      interaction = TRUE,
                      
                      casecontrol = FALSE)
      
      summ1<-summary(fit1)
      
      te.fit<-glm(y1~drug+ldlcb,data = data,family = "binomial")
      
      e1.te2[i]<-te.fit$coefficients[2]
      e1.te[i]<-summ1$summary_myreg[6]
      e1.pnde[i]<-summ1$summary_myreg[2]
      e1.tnie[i]<-summ1$summary_myreg[3]
      e1.tnde[i]<-summ1$summary_myreg[4]
      e1.pnie[i]<-summ1$summary_myreg[5]
      e1.pm[i]<-summ1$summary_myreg[7]
      e1.beta0[i]<-fit1$mreg_fit$coefficients[1]
      e1.beta1[i]<-fit1$mreg_fit$coefficients[2]
      e1.beta2[i]<-fit1$mreg_fit$coefficients[3]
      #e1.sigma[i]<-summ1$summary_mreg_fit$sigma
      e1.theta1[i]<-fit1$yreg_fit$coefficients[2]
      e1.theta2[i]<-fit1$yreg_fit$coefficients[3]
      e1.theta3[i]<-fit1$yreg_fit$coefficients[5]
      e1.theta4[i]<-fit1$yreg_fit$coefficients[4]
      
      fit2<-regmedint(data = data,
                      ## Variables
                      yvar = "y1",
                      avar = "drug",
                      mvar = "ldlc_change_m",
                      cvar = c("ldlcb"),
                      #eventvar = "cen",
                      ## Values at which effects are evaluated
                      a0 = 0,
                      a1 = 1,
                      m_cde = 0,
                      c_cond = 120,
                      ## Model types
                      mreg = "linear",
                      #yreg = "survCox",
                      yreg = "logistic",
                      ## Additional specification
                      interaction = TRUE,
                      
                      casecontrol = FALSE)
      
      summ2<-summary(fit2)
      
      e2.te2[i]<-te.fit$coefficients[2]
      e2.te[i]<-summ2$summary_myreg[6]
      e2.pnde[i]<-summ2$summary_myreg[2]
      e2.tnie[i]<-summ2$summary_myreg[3]
      e2.tnde[i]<-summ2$summary_myreg[4]
      e2.pnie[i]<-summ2$summary_myreg[5]
      e2.pm[i]<-summ2$summary_myreg[7]
      e2.beta0[i]<-fit2$mreg_fit$coefficients[1]
      e2.beta1[i]<-fit2$mreg_fit$coefficients[2]
      e2.beta2[i]<-fit2$mreg_fit$coefficients[3]
      #e2.sigma[i]<-summ2$summary_mreg_fit$sigma
      e2.theta1[i]<-fit2$yreg_fit$coefficients[2]
      e2.theta2[i]<-fit2$yreg_fit$coefficients[3]
      e2.theta3[i]<-fit2$yreg_fit$coefficients[5]
      e2.theta4[i]<-fit2$yreg_fit$coefficients[4]
      
      fit3<-regmedint(data = data,
                      ## Variables
                      yvar = "y1",
                      avar = "drug",
                      mvar = "ldlc_change_m_pred",
                      cvar = c("ldlcb"),
                      #eventvar = "cen",
                      ## Values at which effects are evaluated
                      a0 = 0,
                      a1 = 1,
                      m_cde = 0,
                      c_cond = 120,
                      ## Model types
                      mreg = "linear",
                      #yreg = "survCox",
                      yreg = "logistic",
                      ## Additional specification
                      interaction = TRUE,
                      
                      casecontrol = FALSE)
      
      summ3<-summary(fit3)
      
      
      #te.fit<-survreg(Surv(y1,cen)~drug+ldlcb,data = data,dist = "weibull")
      # out.fit<-survreg(Surv(y1,cen)~drug+ldlc_change.m+drug*ldlc_change.m, dist="weibull",data = data)
      # med.fit<-lm(data = data,ldlc_change ~ drug)
      
      
      e3.te2[i]<-te.fit$coefficients[2]
      e3.te[i]<-summ3$summary_myreg[6]
      e3.pnde[i]<-summ3$summary_myreg[2]
      e3.tnie[i]<-summ3$summary_myreg[3]
      e3.tnde[i]<-summ3$summary_myreg[4]
      e3.pnie[i]<-summ3$summary_myreg[5]
      e3.pm[i]<-summ3$summary_myreg[7]
      e3.beta0[i]<-fit3$mreg_fit$coefficients[1]
      e3.beta1[i]<-fit3$mreg_fit$coefficients[2]
      e3.beta2[i]<-fit3$mreg_fit$coefficients[3]
     # e3.sigma[i]<-summ3$summary_mreg_fit$sigma
      e3.theta1[i]<-fit3$yreg_fit$coefficients[2]
      e3.theta2[i]<-fit3$yreg_fit$coefficients[3]
      e3.theta3[i]<-fit3$yreg_fit$coefficients[5]
      e3.theta4[i]<-fit3$yreg_fit$coefficients[4]
      
      
      
      
    }     
    
    results<-data.frame(cutoff,
                        e1.te2,e1.te,e1.pnde,e1.tnie,e1.tnde,e1.pnie,e1.pm,
                        e1.beta0,e1.beta1,e1.beta2,
                        #e1.sigma,
                        e1.theta1,e1.theta2,e1.theta3,e1.theta4,
                        
                        e2.te2,e2.te,e2.pnde,e2.tnie,e2.tnde,e2.pnie,e2.pm,
                        e2.beta0,e2.beta1,e2.beta2,
                        #e2.sigma,
                        e2.theta1,e2.theta2,e2.theta3,e2.theta4,
                        
                        e3.te2,e3.te,e3.pnde,e3.tnie,e3.tnde,e3.pnie,e3.pm,
                        e3.beta0,e3.beta1,e3.beta2,
                        #e3.sigma,
                        e3.theta1,e3.theta2,e3.theta3,e3.theta4,
                        
                        
                        interaction = 1)
    results$iter<-j
    rd11<-rbind(rd11,results)
    
    
  }, error=function(e){} )
  
}
rd11$method<-"Adjust_RGTM_H1"#

#e1 biased
#e2 gound truth
#e3 adjusted


n<-1000
# Simulation:

# Method: 

rd10<-NULL
set.seed(1)
for(j in 1:n){
  
  tryCatch({
    print(paste("2",j))
    e1.te2<-NULL
    e1.te<-NULL
    e1.pnde<-NULL
    e1.tnde<-NULL
    e1.pnie<-NULL
    e1.tnie<-NULL
    e1.pm<-NULL
    e1.beta0<-NULL
    e1.beta1<-NULL
    e1.beta2<-NULL
    e1.sigma<-NULL
    e1.theta1<-NULL
    e1.theta2<-NULL
    e1.theta3<-NULL
    e1.theta4<-NULL
    
    
    e2.te2<-NULL
    e2.te<-NULL
    e2.pnde<-NULL
    e2.tnde<-NULL
    e2.pnie<-NULL
    e2.tnie<-NULL
    e2.pm<-NULL
    e2.beta0<-NULL
    e2.beta1<-NULL
    e2.beta2<-NULL
    e2.sigma<-NULL
    e2.theta1<-NULL
    e2.theta2<-NULL
    e2.theta3<-NULL
    e2.theta4<-NULL
    
    e3.te2<-NULL
    e3.te<-NULL
    e3.pnde<-NULL
    e3.tnde<-NULL
    e3.pnie<-NULL
    e3.tnie<-NULL
    e3.pm<-NULL
    e3.beta0<-NULL
    e3.beta1<-NULL
    e3.beta2<-NULL
    e3.sigma<-NULL
    e3.theta1<-NULL
    e3.theta2<-NULL
    e3.theta3<-NULL
    e3.theta4<-NULL
    
    cutoff<-c(90,110,130,150,1000)
    for(i in 1:5){
      data<-simulation.function.v4(20000,cutoff = cutoff[i])
      data$ldlc_change<-data$ldlc12 - data$ldlcb
      data$ldlc_change_m<-data$ldlc12m - data$ldlcm
      fit<-lm(data,formula = ldlc_change_m~ldlc_change+ldlcb+drug)
      data$ldlc_change_m_pred<-predict(fit,data)
      data$cen<-rep(1,20000)
      #data$d<-rep(0,20000)
      #data$logldlcb<-log(data$ldlcb)
      data$y1<-rweibull(20000,shape = 1/p,scale = exp(theta0+theta1*data$drug+
                                                        theta2*data$ldlc_change_m+
                                                        theta3*data$ldlc_change_m*data$drug+
                                                        theta4*data$ldlcb))
      
      
      
      c_median<-summary(data$ldlcb)[3]
      c_mean<-summary(data$ldlcb)[4]
      fit1<-regmedint(data = data,
                      ## Variables
                      yvar = "y1",
                      avar = "drug",
                      mvar = "ldlc_change",
                      cvar = c("ldlcb"),
                      eventvar = "cen",
                      ## Values at which effects are evaluated
                      a0 = 0,
                      a1 = 1,
                      m_cde = 0,
                      c_cond = 120,
                      ## Model types
                      mreg = "linear",
                      #yreg = "survCox",
                      yreg = "survAFT_weibull",
                      ## Additional specification
                      interaction = TRUE,
                      
                      casecontrol = FALSE)
      
      summ1<-summary(fit1)
      
      te.fit<-survreg(Surv(y1,cen)~drug+ldlcb,data = data,dist = "weibull")
      
      e1.te2[i]<-te.fit$coefficients[2]
      e1.te[i]<-summ1$summary_myreg[6]
      e1.pnde[i]<-summ1$summary_myreg[2]
      e1.tnie[i]<-summ1$summary_myreg[3]
      e1.tnde[i]<-summ1$summary_myreg[4]
      e1.pnie[i]<-summ1$summary_myreg[5]
      e1.pm[i]<-summ1$summary_myreg[7]
      e1.beta0[i]<-fit1$mreg_fit$coefficients[1]
      e1.beta1[i]<-fit1$mreg_fit$coefficients[2]
      e1.beta2[i]<-fit1$mreg_fit$coefficients[3]
      e1.sigma[i]<-summ1$summary_mreg_fit$sigma
      e1.theta1[i]<-fit1$yreg_fit$coefficients[2]
      e1.theta2[i]<-fit1$yreg_fit$coefficients[3]
      e1.theta3[i]<-fit1$yreg_fit$coefficients[5]
      e1.theta4[i]<-fit1$yreg_fit$coefficients[4]
      
      fit2<-regmedint(data = data,
                      ## Variables
                      yvar = "y1",
                      avar = "drug",
                      mvar = "ldlc_change_m",
                      cvar = c("ldlcb"),
                      eventvar = "cen",
                      ## Values at which effects are evaluated
                      a0 = 0,
                      a1 = 1,
                      m_cde = 0,
                      c_cond = 120,
                      ## Model types
                      mreg = "linear",
                      #yreg = "survCox",
                      yreg = "survAFT_weibull",
                      ## Additional specification
                      interaction = TRUE,
                      
                      casecontrol = FALSE)
      
      summ2<-summary(fit2)
      
      e2.te2[i]<-te.fit$coefficients[2]
      e2.te[i]<-summ2$summary_myreg[6]
      e2.pnde[i]<-summ2$summary_myreg[2]
      e2.tnie[i]<-summ2$summary_myreg[3]
      e2.tnde[i]<-summ2$summary_myreg[4]
      e2.pnie[i]<-summ2$summary_myreg[5]
      e2.pm[i]<-summ2$summary_myreg[7]
      e2.beta0[i]<-fit2$mreg_fit$coefficients[1]
      e2.beta1[i]<-fit2$mreg_fit$coefficients[2]
      e2.beta2[i]<-fit2$mreg_fit$coefficients[3]
      e2.sigma[i]<-summ2$summary_mreg_fit$sigma
      e2.theta1[i]<-fit2$yreg_fit$coefficients[2]
      e2.theta2[i]<-fit2$yreg_fit$coefficients[3]
      e2.theta3[i]<-fit2$yreg_fit$coefficients[5]
      e2.theta4[i]<-fit2$yreg_fit$coefficients[4]
      
      fit3<-regmedint(data = data,
                      ## Variables
                      yvar = "y1",
                      avar = "drug",
                      mvar = "ldlc_change_m_pred",
                      cvar = c("ldlcb"),
                      eventvar = "cen",
                      ## Values at which effects are evaluated
                      a0 = 0,
                      a1 = 1,
                      m_cde = 0,
                      c_cond = 120,
                      ## Model types
                      mreg = "linear",
                      #yreg = "survCox",
                      yreg = "survAFT_weibull",
                      ## Additional specification
                      interaction = TRUE,
                      
                      casecontrol = FALSE)
      
      summ3<-summary(fit3)
      
      
      #te.fit<-survreg(Surv(y1,cen)~drug+ldlcb,data = data,dist = "weibull")
      # out.fit<-survreg(Surv(y1,cen)~drug+ldlc_change.m+drug*ldlc_change.m, dist="weibull",data = data)
      # med.fit<-lm(data = data,ldlc_change ~ drug)
      
      
      e3.te2[i]<-te.fit$coefficients[2]
      e3.te[i]<-summ3$summary_myreg[6]
      e3.pnde[i]<-summ3$summary_myreg[2]
      e3.tnie[i]<-summ3$summary_myreg[3]
      e3.tnde[i]<-summ3$summary_myreg[4]
      e3.pnie[i]<-summ3$summary_myreg[5]
      e3.pm[i]<-summ3$summary_myreg[7]
      e3.beta0[i]<-fit3$mreg_fit$coefficients[1]
      e3.beta1[i]<-fit3$mreg_fit$coefficients[2]
      e3.beta2[i]<-fit3$mreg_fit$coefficients[3]
      e3.sigma[i]<-summ3$summary_mreg_fit$sigma
      e3.theta1[i]<-fit3$yreg_fit$coefficients[2]
      e3.theta2[i]<-fit3$yreg_fit$coefficients[3]
      e3.theta3[i]<-fit3$yreg_fit$coefficients[5]
      e3.theta4[i]<-fit3$yreg_fit$coefficients[4]
      
      
    }      
    
    
    
    results<-data.frame(cutoff,
                        e1.te2,e1.te,e1.pnde,e1.tnie,e1.tnde,e1.pnie,e1.pm,
                        e1.beta0,e1.beta1,e1.beta2,e1.sigma,
                        e1.theta1,e1.theta2,e1.theta3,e1.theta4,
                        
                        e2.te2,e2.te,e2.pnde,e2.tnie,e2.tnde,e2.pnie,e2.pm,
                        e2.beta0,e2.beta1,e2.beta2,e2.sigma,
                        e2.theta1,e2.theta2,e2.theta3,e2.theta4,
                        
                        e3.te2,e3.te,e3.pnde,e3.tnie,e3.tnde,e3.pnie,e3.pm,
                        e3.beta0,e3.beta1,e3.beta2,e3.sigma,
                        e3.theta1,e3.theta2,e3.theta3,e3.theta4,
                        
                        
                        interaction = 1)
    results$iter<-j
    rd10<-rbind(rd10,results)
    
  }, error=function(e){} )
  
}
rd10$method<-"Adjust_RGTM_H0"#

setwd("/Users/sichenghao/Documents/GitHub/Mediation-RGTM/")
#fwrite(rd10,file = "rd10.csv")
fwrite(rd11,file = "rd11_binary_5000.csv")


