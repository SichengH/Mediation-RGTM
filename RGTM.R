





#cutoff<-c(80,90,100,110,120,130,140,150,160,1000)
#for(i in 1:10){
library(ggplot2)
library(dplyr)
library(data.table)
library(survival)
library(DiagrammeR)
library(regmedint)

p<<-1.4
#outcome regression AFT
theta0<<-7.5
theta1<<-0.12
theta2<<-0.015#should be negative
theta3<<--0.023
theta4<<--0.02

#outcome regression AFT, no interaction term

theta0_ni<<-7.7
theta1_ni<<-0.59
theta2_ni<<-0.0032 #should be negative
#
beta1<<--47

# #Cox parameter
# theta1<--0.15
# theta2<--0.022
# theta3<-0.032321
# 
# theta1_ni<--0.82
# theta2_ni<--0.005

n<-200
################# PART ONE ########################

#method1

rd11<-NULL
rd10<-NULL
set.seed(1)
for(j in 1:n){
  
  tryCatch({
    print(paste("1",j))
    e.te2<-NULL
    e.te<-NULL
    e.pnde<-NULL
    e.tnde<-NULL
    e.pnie<-NULL
    e.tnie<-NULL
    e.pm<-NULL
    e.beta0<-NULL
    e.beta1<-NULL
    e.beta2<-NULL
    e.sigma<-NULL
    e.theta1<-NULL
    e.theta2<-NULL
    e.theta3<-NULL
    e.theta4<-NULL
    
    e0.te2<-NULL
    e0.te<-NULL
    e0.pnde<-NULL
    e0.tnde<-NULL
    e0.pnie<-NULL
    e0.tnie<-NULL
    e0.pm<-NULL
    e0.beta0<-NULL
    e0.beta1<-NULL
    e0.beta2<-NULL
    e0.sigma<-NULL
    e0.theta1<-NULL
    e0.theta2<-NULL
    #e0.theta3<-NULL
    e0.theta4<-NULL
    
    cutoff<-c(80,90,100,110,120,130,140,150,160,1000)
    for(i in 1:10){
      data<-simulation.function.v4(20000,cutoff = cutoff[i])
      data$ldlc_change<-data$ldlc12 - data$ldlcb
      data$ldlc_change_m<-data$ldlc12m - data$ldlcm
      data$cen<-rep(1,20000)
      #data$d<-rep(0,20000)
      #data$logldlcb<-log(data$ldlcb)
      data$y1<-rweibull(20000,shape = 1/p,scale = exp(theta0+theta1*data$drug+
                                                        theta2*data$ldlc_change_m+
                                                        theta3*data$ldlc_change_m*data$drug+
                                                        theta4*data$ldlcb))
      
      data$y0<-rweibull(20000,shape = 1/p,scale = exp(theta0_ni+theta1_ni*data$drug+
                                                        theta2_ni*data$ldlc_change.new+
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
      
      summ<-summary(fit1)
      
      
      
      te.fit<-survreg(Surv(y1,cen)~drug+ldlcb,data = data,dist = "weibull")
      # out.fit<-survreg(Surv(y1,cen)~drug+ldlc_change.m+drug*ldlc_change.m, dist="weibull",data = data)
      # med.fit<-lm(data = data,ldlc_change ~ drug)
      
      
      e.te2[i]<-te.fit$coefficients[2]
      e.te[i]<-summ$summary_myreg[6]
      e.pnde[i]<-summ$summary_myreg[2]
      e.tnie[i]<-summ$summary_myreg[3]
      e.tnde[i]<-summ$summary_myreg[4]
      e.pnie[i]<-summ$summary_myreg[5]
      e.pm[i]<-summ$summary_myreg[7]
      e.beta0[i]<-fit1$mreg_fit$coefficients[1]
      e.beta1[i]<-fit1$mreg_fit$coefficients[2]
      e.beta2[i]<-fit1$mreg_fit$coefficients[3]
      e.sigma[i]<-summ$summary_mreg_fit$sigma
      e.theta1[i]<-fit1$yreg_fit$coefficients[2]
      e.theta2[i]<-fit1$yreg_fit$coefficients[3]
      e.theta3[i]<-fit1$yreg_fit$coefficients[5]
      e.theta4[i]<-fit1$yreg_fit$coefficients[4]
      
      
      #No interaction
      fit0<-regmedint(data = data,
                      ## Variables
                      yvar = "y0",
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
                      yreg = "survCox",
                      ## Additional specification
                      interaction = FALSE,
                      casecontrol = FALSE)
      
      summ<-summary(fit0)
      
      
      te.fit<-survreg(Surv(y0,cen)~drug+ldlcb,data = data,dist = "weibull")
      #out.fit<-survreg(Surv(y)~drug+ldlc_change+drug*ldlc_change, dist="weibull",data = data)
      #med.fit<-lmw(data = data,ldlc_change ~ drug)
      
      
      e0.te2[i]<-te.fit$coefficients[1]
      e0.te[i]<-summ$summary_myreg[6]
      e0.pnde[i]<-summ$summary_myreg[2]
      e0.tnie[i]<-summ$summary_myreg[3]
      e0.tnde[i]<-summ$summary_myreg[4]
      e0.pnie[i]<-summ$summary_myreg[5]
      e0.pm[i]<-summ$summary_myreg[7]
      e0.beta0[i]<-fit0$mreg_fit$coefficients[1]
      e0.beta1[i]<-fit0$mreg_fit$coefficients[2]
      e0.beta2[i]<-fit0$mreg_fit$coefficients[3]
      e0.sigma[i]<-summ$summary_mreg_fit$sigma
      e0.theta1[i]<-fit0$yreg_fit$coefficients[1]
      e0.theta2[i]<-fit0$yreg_fit$coefficients[2]
      e0.theta4[i]<-fit0$yreg_fit$coefficients[3]
      
    }
    results<-data.frame(cutoff,e.te2,e.te,e.pnde,e.tnie,e.tnde,e.pnie,e.pm,
                        e.beta0,e.beta1,e.beta2,e.sigma,
                        e.theta1,e.theta2,e.theta3,e.theta4,interaction = 1)
    results$iter<-j
    rd11<-rbind(rd11,results)
    
    results0<-data.frame(cutoff,e0.te2,e0.te,e0.pnde,e0.tnie,e0.tnde,e0.pnie,e0.pm,
                         e0.beta0,e0.beta1,e0.beta2,e0.sigma,
                         e0.theta1,e0.theta2,e0.theta4,interaction = 0)
    results0$iter<-j
    rd10<-rbind(rd10,results0)
    
  }, error=function(e){})
  
}
rd11$method<-"H0I1"
rd10$method<-"H0I0"

rd41<-NULL
rd40<-NULL
set.seed(1)
for(j in 1:n){
  
  tryCatch({
    print(paste("2",j))
    e.te2<-NULL
    e.te<-NULL
    e.pnde<-NULL
    e.tnde<-NULL
    e.pnie<-NULL
    e.tnie<-NULL
    e.pm<-NULL
    e.beta0<-NULL
    e.beta1<-NULL
    e.beta2<-NULL
    e.sigma<-NULL
    e.theta1<-NULL
    e.theta2<-NULL
    e.theta3<-NULL
    e.theta4<-NULL
    
    e0.te2<-NULL
    e0.te<-NULL
    e0.pnde<-NULL
    e0.tnde<-NULL
    e0.pnie<-NULL
    e0.tnie<-NULL
    e0.pm<-NULL
    e0.beta0<-NULL
    e0.beta1<-NULL
    e0.beta2<-NULL
    e0.sigma<-NULL
    e0.theta1<-NULL
    e0.theta2<-NULL
    #e0.theta3<-NULL
    e0.theta4<-NULL
    
    cutoff<-c(80,90,100,110,120,130,140,150,160,1000)
    for(i in 1:10){
      data<-simulation.function21.new(20000,cutoff = cutoff[i])
      data$ldlc_change.new<-data$ldlc12 - data$ldlcm
      data$cen<-rep(1,20000)
      data$d<-rep(0,20000)
      data$logldlcb<-log(data$ldlcb)
      data$y1<-rweibull(20000,shape = 1/p,scale = exp(theta0+theta1*data$drug+
                                                        theta2*data$ldlc_change.new+
                                                        theta3*data$ldlc_change.new*data$drug+
                                                        theta4*data$ldlcb))
      
      data$y0<-rweibull(20000,shape = 1/p,scale = exp(theta0_ni+theta1_ni*data$drug+
                                                        theta2_ni*data$ldlc_change.new+
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
      
      summ<-summary(fit1)
      
      
      
      te.fit<-survreg(Surv(y1,cen)~drug+ldlcb,data = data,dist = "weibull")
      # out.fit<-survreg(Surv(y1,cen)~drug+ldlc_change.m+drug*ldlc_change.m, dist="weibull",data = data)
      # med.fit<-lm(data = data,ldlc_change ~ drug)
      
      
      e.te2[i]<-te.fit$coefficients[2]
      e.te[i]<-summ$summary_myreg[6]
      e.pnde[i]<-summ$summary_myreg[2]
      e.tnie[i]<-summ$summary_myreg[3]
      e.tnde[i]<-summ$summary_myreg[4]
      e.pnie[i]<-summ$summary_myreg[5]
      e.pm[i]<-summ$summary_myreg[7]
      e.beta0[i]<-fit1$mreg_fit$coefficients[1]
      e.beta1[i]<-fit1$mreg_fit$coefficients[2]
      e.beta2[i]<-fit1$mreg_fit$coefficients[3]
      e.sigma[i]<-summ$summary_mreg_fit$sigma
      e.theta1[i]<-fit1$yreg_fit$coefficients[2]
      e.theta2[i]<-fit1$yreg_fit$coefficients[3]
      e.theta3[i]<-fit1$yreg_fit$coefficients[5]
      e.theta4[i]<-fit1$yreg_fit$coefficients[4]
      
      
      #No interaction
      fit0<-regmedint(data = data,
                      ## Variables
                      yvar = "y0",
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
                      yreg = "survCox",
                      ## Additional specification
                      interaction = FALSE,
                      casecontrol = FALSE)
      
      summ<-summary(fit0)
      
      
      te.fit<-survreg(Surv(y0,cen)~drug+ldlcb,data = data,dist = "weibull")
      #out.fit<-survreg(Surv(y)~drug+ldlc_change+drug*ldlc_change, dist="weibull",data = data)
      #med.fit<-lmw(data = data,ldlc_change ~ drug)
      
      
      e0.te2[i]<-te.fit$coefficients[1]
      e0.te[i]<-summ$summary_myreg[6]
      e0.pnde[i]<-summ$summary_myreg[2]
      e0.tnie[i]<-summ$summary_myreg[3]
      e0.tnde[i]<-summ$summary_myreg[4]
      e0.pnie[i]<-summ$summary_myreg[5]
      e0.pm[i]<-summ$summary_myreg[7]
      e0.beta0[i]<-fit0$mreg_fit$coefficients[1]
      e0.beta1[i]<-fit0$mreg_fit$coefficients[2]
      e0.beta2[i]<-fit0$mreg_fit$coefficients[3]
      e0.sigma[i]<-summ$summary_mreg_fit$sigma
      e0.theta1[i]<-fit0$yreg_fit$coefficients[1]
      e0.theta2[i]<-fit0$yreg_fit$coefficients[2]
      e0.theta4[i]<-fit0$yreg_fit$coefficients[3]
      
    }
    results<-data.frame(cutoff,e.te2,e.te,e.pnde,e.tnie,e.tnde,e.pnie,e.pm,
                        e.beta0,e.beta1,e.beta2,e.sigma,
                        e.theta1,e.theta2,e.theta3,e.theta4,interaction = 1)
    results$iter<-j
    rd41<-rbind(rd41,results)
    
    results0<-data.frame(cutoff,e0.te2,e0.te,e0.pnde,e0.tnie,e0.tnde,e0.pnie,e0.pm,
                         e0.beta0,e0.beta1,e0.beta2,e0.sigma,
                         e0.theta1,e0.theta2,e0.theta4,interaction = 0)
    results0$iter<-j
    rd40<-rbind(rd40,results0)
    
  }, error=function(e){})
  
}
rd41$method<-"H1I1"
rd40$method<-"H1I0"
rdI1<-rbind(rd11,rd41)
rdI0<-rbind(rd10,rd40)
setwd("/Users/sichenghao/Desktop/Thesis_Simu/RGTM/")

write.csv(rdI1,file = "n1000-1_new2.csv")
write.csv(rdI0,file = "n1000-0_new2.csv")


# ggplot(data = data)+
#   geom_point(aes(x = ldlcb,y = ldlc_change,color = as.factor(drug)))+
#   geom_smooth(method = lm,aes(x = ldlcb,y = ldlc_change,color = as.factor(drug)))




