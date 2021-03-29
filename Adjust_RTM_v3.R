
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
theta1<<-0.0  #direct effect, positive is protective
theta2<<--0.015  #effect * ldlc_change; negative is protective
theta3<<-0 # ldlc_change * drug assume no effect
theta4<<--0.02 # cov

# theta0_ni<<-7.7
# theta1_ni<<-0.59
# theta2_ni<<-0.0032 #should be negative
cutoff<-c(90,110,130,150,1000)

#
data<-simulation.function.v4(20000,cutoff = cutoff[3])
data$ldlc_change<-data$ldlc12 - data$ldlcb

data$ldlc_change_m<-data$ldlc12m - data$ldlcm

#fit<-lm(data,formula = ldlc_change_m~ldlc_change+ldlcb+drug)
#plot(fit$residuals)

#data$ldlc_change_m_pred<-predict(fit,data)
#plot(data$ldlc_change_m,data$ldlc_change_m_pred)#mean prediction works

#using predicted average as mediator


data$cen<-rep(1,20000)

data$y1<-rweibull(20000,shape = 1/p,scale = exp(theta0+theta1*data$drug+
                                                  theta2*data$ldlc_change_m+
                                                  theta3*data$ldlc_change_m*data$drug+
                                                  theta4*data$ldlcb))


#explore
ggplot(data) + geom_point((aes(x = ldlcb,y = ldlc_change,col = drug)))#negative

ggplot(data) + geom_point((aes(x = ldlcb,y = ldlc_change_m,col = drug)))#no relation 

ldlcb.to.outcome<-survreg(Surv(y1,cen)~ldlcb,data = data,dist = "weibull")
summary(ldlcb.to.outcome)

block.all<-survreg(Surv(y1,cen)~drug+ldlc_change+ldlcb,data = data,dist = "weibull")
summary(block.all)

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
summ
#invastigate if the adjustment works or not using simulation



n<-1000
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
      data<-simulation.function.v43(20000,cutoff = cutoff[i])
      data$ldlc_change<-data$ldlc121 - data$ldlcb
      data$ldlc_change2<-data$ldlc12 - data$ldlcb
      data$ldlc_change_m<-data$ldlc12m - data$ldlcm
      #fit<-lm(data,formula = ldlc_change_m~ldlc_change+ldlcb+drug)
      # data$ldlc_change_m_pred<-predict(fit,data)
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
                      mvar = "ldlc_change2",
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
      
      data0<-data%>%filter(drug == 0)
      
      
      m.fit<-lm(data = data0,ldlc_change~ldlcb)
      data$ldlc_change_adj<-data$ldlc_change - 
        data$ldlcb * m.fit$coefficients[2] - 
        rep(m.fit$coefficients[1],nrow(data))
      #ggplot(data,aes(x = ldlcb,y = ldlc_change_adj,col = drug))+
      #  geom_point(alpha = 0.1)+
      #  geom_smooth(method= "lm")
      
      
      fit3<-regmedint(data = data,
                      ## Variables
                      yvar = "y1",
                      avar = "drug",
                      mvar = "ldlc_change_adj",
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
    rd11<-rbind(rd11,results)
    
    
  }, error=function(e){} )
  
}
rd11$method<-"1_0_adj"#

#e1 Interaction
#e2 No interaction
#e3 adjusted with interaction


setwd("/Users/sh/Documents/GitHub/Mediation-RGTM/")

fwrite(rd11,file = "v3-1.csv")

#v2-0: messed up adj
#v2-1: updated adj, correctly, positive theta2(harmful)
#v2-2: CATE beta1(negative theta2)
#v2-3: f4, 
#v2-4: f42, CATE beta1

#v3-1: adjust vs multiple measurements 



