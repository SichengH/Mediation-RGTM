simulation.function.v5<-function(sample.size = 20000,cutoff = 1000){
  n<-sample.size/2
  osp = 2
  n0<-0
  n1<-0
  data0<-NULL
  data1<-NULL
  #control arm
  for(i in 1:100){
    #print(i)
    if(n0<n){
      ldlcm=rnorm(n = n*osp,mean = 116,sd = 14.2)#mean distribution of LDLC
      d0<-data.frame(ldlcm,ldlcb=rnorm(n = n*osp,mean = ldlcm,sd = 16))
      d0<-d0%>%filter(ldlcb<cutoff)
      data0<-rbind(data0,d0)
      n0<-nrow(data0)
    } else {
      break
    }
  }
  
  data0<-data0[1:n,]
  #data0$ldlc12m<-data0$ldlcm+rnorm(n,0,10)#12m has higher variacne then ldlcm due to possible life style change after entering the trial
  data0$ldlc12<-rnorm(n = n,mean = data0$ldlcm,sd = 16)
  #treatment arm
  for(i in 1:100){
    #print(i)
    if(n1<n){
      ldlcm=rnorm(n = n*osp,mean = 116,sd = 14.2)
      d1<-data.frame(ldlcm,ldlcb=rnorm(n = n*osp,mean = ldlcm,sd = 16))
      d1<-d1%>%filter(ldlcb<cutoff)
      data1<-rbind(data1,d1)
      n1<-nrow(data1)
    } else {
      break
    }
  }
  data1<-data1[1:n,]
  data1$ldlc12<-data1$ldlcm - 47 + rnorm(n,0,16)
  data0$drug = 0
  data1$drug = 1
  data = rbind(data0,data1)
  data<-data%>%mutate(ldlc_change = ldlc12 - ldlcb)
  return(data)
}


simulation.function.v42<-function(sample.size = 20000,cutoff = 1000){
  n<-sample.size/2
  osp = 2
  n0<-0
  n1<-0
  data0<-NULL
  data1<-NULL
  #control arm
  for(i in 1:100){
    #print(i)
    if(n0<n){
      ldlcm=rnorm(n = n*osp,mean = 116,sd = 14.2)#mean distribution of LDLC
      d0<-data.frame(ldlcm,ldlcb=rnorm(n = n*osp,mean = ldlcm,sd = 16))
      d0<-d0%>%filter(ldlcb<cutoff)
      data0<-rbind(data0,d0)
      n0<-nrow(data0)
    } else {
      break
    }
  }
  
  data0<-data0[1:n,]
  data0$ldlc12m<-data0$ldlcm+rnorm(n,0,15)#12m has higher variacne then ldlcm due to possible life style change after entering the trial
  data0$ldlc12<-data0$ldlc12m+rnorm(n,0,16)#add vairance of rgtm
  #treatment arm
  for(i in 1:100){
    #print(i)
    if(n1<n){
      ldlcm=rnorm(n = n*osp,mean = 116,sd = 14.2)
      d1<-data.frame(ldlcm,ldlcb=rnorm(n = n*osp,mean = ldlcm,sd = 16))
      d1<-d1%>%filter(ldlcb<cutoff)
      data1<-rbind(data1,d1)
      n1<-nrow(data1)
    } else {
      break
    }
  }
  data1<-data1[1:n,]
  data1$ldlc12m<-data1$ldlcm - 0.7*data1$ldlcb +30 + rnorm(n,0,15)
  data1$ldlc12<-data1$ldlc12m+rnorm(n,0,sd = 16)
  data0$drug = 0
  data1$drug = 1
  data = rbind(data0,data1)
  #data<-data%>%mutate(ldlc_change = ldlc12 - ldlcb,ldlc_change.m = ldlc12m - ldlcb)
  return(data)
}
simulation.function.v4<-function(sample.size = 20000,cutoff = 1000){
  n<-sample.size/2
  osp = 2
  n0<-0
  n1<-0
  data0<-NULL
  data1<-NULL
  #control arm
  for(i in 1:100){
    #print(i)
    if(n0<n){
      ldlcm=rnorm(n = n*osp,mean = 116,sd = 14.2)#mean distribution of LDLC
      d0<-data.frame(ldlcm,ldlcb=rnorm(n = n*osp,mean = ldlcm,sd = 16))
      d0<-d0%>%filter(ldlcb<cutoff)
      data0<-rbind(data0,d0)
      n0<-nrow(data0)
    } else {
      break
    }
  }
  
  data0<-data0[1:n,]
  data0$ldlc12m<-data0$ldlcm+rnorm(n,0,15)#12m has higher variance then ldlcm due to possible life style change after entering the trial
  data0$ldlc12<-data0$ldlc12m+rnorm(n,0,16)#add variance of rgtm
  #treatment arm
  for(i in 1:100){
    #print(i)
    if(n1<n){
      ldlcm=rnorm(n = n*osp,mean = 116,sd = 14.2)
      d1<-data.frame(ldlcm,ldlcb=rnorm(n = n*osp,mean = ldlcm,sd = 16))
      d1<-d1%>%filter(ldlcb<cutoff)
      data1<-rbind(data1,d1)
      n1<-nrow(data1)
    } else {
      break
    }
  }
  data1<-data1[1:n,]
  data1$ldlc12m<-data1$ldlcm - 47 + rnorm(n,0,15)
  data1$ldlc12<-data1$ldlc12m+rnorm(n,0,sd = 16)
  data0$drug = 0
  data1$drug = 1
  data = rbind(data0,data1)
  #data<-data%>%mutate(ldlc_change = ldlc12 - ldlcb,ldlc_change.m = ldlc12m - ldlcb)
  return(data)
}




#First number is what mediator distribution
#1 as normal homogenious, 2 as log-normal homogenious
#3 as normal heteogenious and 4 as log-normal heterogenious

simulation.function11v3.multi<-function(sample.size = 20000,cutoff = 1000){
  n<-sample.size/2
  osp = 2
  n0<-0
  n1<-0
  data0<-NULL
  data1<-NULL
  #control arm
  for(i in 1:100){
    #print(i)
    if(n0<n){
      ldlcm=rnorm(n = n*osp,mean = 116,sd = 14.2)#mean distribution of LDLC
      d0<-data.frame(ldlcm,ldlcb1=rnorm(n = n*osp,mean = ldlcm,sd = 6),
                     ldlcb2=rnorm(n = n*osp,mean = ldlcm,sd = 6),
                     ldlcb3=rnorm(n = n*osp,mean = ldlcm,sd = 6),
                     ldlc12m = ldlcm,
                     ldlc12m1=rnorm(n = n*osp,mean = ldlcm,sd = 6),
                     ldlc12m2=rnorm(n = n*osp,mean = ldlcm,sd = 6),
                     ldlc12m3=rnorm(n = n*osp,mean = ldlcm,sd = 6))#dataframe for control arm
      d0<-d0%>%mutate(ldlcb =(ldlcb1+ldlcb2+ldlcb3)/3)%>%filter(ldlcb<cutoff)
      data0<-rbind(data0,d0)
      n0<-nrow(data0)
    } else {
      break
    }
  }
  
  data0<-data0[1:n,]
  data0$ldlc12m<-data0$ldlcm+rnorm(n,0,6)
  
  #treatment arm
  for(i in 1:100){
    #print(i)
    if(n1<n){
      ldlcm=rnorm(n = n*osp,mean = 116,sd = 14.2)
      d1<-data.frame(ldlcm,ldlcb1=rnorm(n = n*osp,mean = ldlcm,sd = 6),
                     ldlcb2=rnorm(n = n*osp,mean = ldlcm,sd = 6),
                     ldlcb3=rnorm(n = n*osp,mean = ldlcm,sd = 6))
      d1<-d1%>%mutate(ldlcb =(ldlcb1+ldlcb2+ldlcb3)/3)%>%filter(ldlcb<cutoff)
      data1<-rbind(data1,d1)
      n1<-nrow(data1)
    } else {
      break
    }
  }
  data1<-data1[1:n,]
  data1$ldlc12m<-data1$ldlcm - 47+rnorm(n,0,6)
  data1$ldlc12m1<-data1$ldlc12m+rnorm(n,0,sd = 6)
  data1$ldlc12m2<-data1$ldlc12m+rnorm(n,0,sd = 6)
  data1$ldlc12m3<-data1$ldlc12m+rnorm(n,0,sd = 6)
  
  data0$drug = 0
  data1$drug = 1
  
  data = rbind(data0,data1)
  #data$ldlcb<-mean(data$ldlcb1,data$ldlcb2,data$ldlcb3)
  data$ldlc12<-(data$ldlc12m1+data$ldlc12m2+data$ldlc12m3)/3
  data<-data%>%mutate(ldlc_change = ldlc12 - ldlcb,ldlc_change.m = ldlc12m - ldlcm)
  return(data)
}



simulation.function11v2<-function(sample.size = 20000,cutoff = 1000){
  n<-sample.size/2
  osp = 2
  n0<-0
  n1<-0
  data0<-NULL
  data1<-NULL
  #control arm
  for(i in 1:100){
    #print(i)
    if(n0<n){
      ldlcm=rnorm(n = n*osp,mean = 116,sd = 14.2)#mean distribution of LDLC
      d0<-data.frame(ldlcm,ldlcb=rnorm(n = n*osp,mean = ldlcm,sd = 16),
                     ldlc12=rnorm(n = n*osp,mean = ldlcm,sd = 16))#dataframe for control arm
      d0<-d0%>%filter(ldlcb<cutoff)
      data0<-rbind(data0,d0)
      n0<-nrow(data0)
    } else {
      break
    }
  }
  
  data0<-data0[1:n,]
  data0$ldlc12m<-data0$ldlcm+rnorm(n,0,1)
  
  #treatment arm
  for(i in 1:100){
    #print(i)
    if(n1<n){
      ldlcm=rnorm(n = n*osp,mean = 116,sd = 14.2)
      d1<-data.frame(ldlcm,ldlcb=rnorm(n = n*osp,mean = ldlcm,sd = 16))
      d1<-d1%>%filter(ldlcb<cutoff)
      data1<-rbind(data1,d1)
      n1<-nrow(data1)
    } else {
      break
    }
  }
  data1<-data1[1:n,]
  data1$ldlc12m<-data1$ldlcm - 47+rnorm(n,0,1)
  data1$ldlc12<-data1$ldlc12m+rnorm(n,0,sd = 16)
  data0$drug = 0
  data1$drug = 1
  data = rbind(data0,data1)
  data<-data%>%mutate(ldlc_change = ldlc12 - ldlcb,ldlc_change.m = ldlc12m - ldlcm)
  return(data)
}


simulation.function21v2<-function(sample.size = 20000,cutoff = 1000){
  n<-sample.size/2
  osp = 2
  n0<-0
  n1<-0
  data0<-NULL
  data1<-NULL
  #control arm
  for(i in 1:100){
    #print(i)
    if(n0<n){
      ldlcm=rnorm(n = n*osp,mean = 116,sd = 14.2)#mean distribution of LDLC
      d0<-data.frame(ldlcm,ldlcb=rnorm(n = n*osp,mean = ldlcm,sd = 16),
                     ldlc12=rnorm(n = n*osp,mean = ldlcm,sd = 16))#dataframe for control arm
      d0<-d0%>%filter(ldlcb<cutoff)
      data0<-rbind(data0,d0)
      n0<-nrow(data0)
    } else {
      break
    }
  }
  
  data0<-data0[1:n,]
  data0$ldlc12m<-data0$ldlcm+rnorm(n,0,1)
  
  #treatment arm
  for(i in 1:100){
    #print(i)
    if(n1<n){
      ldlcm=rnorm(n = n*osp,mean = 116,sd = 14.2)
      d1<-data.frame(ldlcm,ldlcb=rnorm(n = n*osp,mean = ldlcm,sd = 16))
      d1<-d1%>%filter(ldlcb<cutoff)
      data1<-rbind(data1,d1)
      n1<-nrow(data1)
    } else {
      break
    }
  }
  data1<-data1[1:n,]
  data1$ldlc12m<-data1$ldlcm - 0.7*data1$ldlcm +30+rnorm(n,0,1)
  data1$ldlc12<-data1$ldlc12m+rnorm(n,0,sd = 16)
  data0$drug = 0
  data1$drug = 1
  data = rbind(data0,data1)
  data<-data%>%mutate(ldlc_change = ldlc12 - ldlcb,ldlc_change.m = ldlc12m - ldlcm)
  return(data)
}


############ Normal Distribution (Distribution-1)#############



simulation.function11.new<-function(sample.size = 20000,cutoff = 1000){
  n<-sample.size/2
  osp = 2
  n0<-0
  n1<-0
  data0<-NULL
  data1<-NULL
  #control arm
  for(i in 1:100){
    #print(i)
    if(n0<n){
      ldlcm=rnorm(n = n*osp,mean = 116,sd = 14.2)#mean distribution of LDLC
      d0<-data.frame(ldlcm,ldlcb=rnorm(n = n*osp,mean = ldlcm,sd = 16),
                     ldlc12=rnorm(n = n*osp,mean = ldlcm,sd = 16))#dataframe for control arm
      d0<-d0%>%filter(ldlcb<cutoff)
      data0<-rbind(data0,d0)
      n0<-nrow(data0)
    } else {
      break
    }
  }

  data0<-data0[1:n,]
  data0$ldlc12m<-data0$ldlcm
  
 #treatment arm
   for(i in 1:100){
    #print(i)
    if(n1<n){
      ldlcm=rnorm(n = n*osp,mean = 116,sd = 14.2)
      d1<-data.frame(ldlcm,ldlcb=rnorm(n = n*osp,mean = ldlcm,sd = 16))
      d1<-d1%>%filter(ldlcb<cutoff)
      data1<-rbind(data1,d1)
      n1<-nrow(data1)
    } else {
      break
    }
  }
  data1<-data1[1:n,]
  data1$ldlc12m<-data1$ldlcm - 47
  data1$ldlc12<-data1$ldlc12m+rnorm(n,0,sd = 16)
  data0$drug = 0
  data1$drug = 1
  data = rbind(data0,data1)
  data<-data%>%mutate(ldlc_change = ldlc12 - ldlcb,ldlc_change.m = ldlc12m - ldlcm)
  return(data)
}

#function12: reduce measurement variation 16—>6

simulation.function12.new<-function(sample.size = 20000,cutoff = 1000){
  n<-sample.size/2
  osp = 2
  n0<-0
  n1<-0
  data0<-NULL
  data1<-NULL
  #control arm
  for(i in 1:100){
    #print(i)
    if(n0<n){
      ldlcm=rnorm(n = n*osp,mean = 116,sd = 14.2)#mean distribution of LDLC
      d0<-data.frame(ldlcm,ldlcb=rnorm(n = n*osp,mean = ldlcm,sd = 6),
                     ldlc12=rnorm(n = n*osp,mean = ldlcm,sd = 6))#dataframe for control arm
      d0<-d0%>%filter(ldlcb<cutoff)
      data0<-rbind(data0,d0)
      n0<-nrow(data0)
    } else {
      break
    }
  }
  
  data0<-data0[1:n,]
  data0$ldlc12m<-data0$ldlcm
  
  #treatment arm
  for(i in 1:100){
    #print(i)
    if(n1<n){
      ldlcm=rnorm(n = n*osp,mean = 116,sd = 14.2)
      d1<-data.frame(ldlcm,ldlcb=rnorm(n = n*osp,mean = ldlcm,sd = 6))
      d1<-d1%>%filter(ldlcb<cutoff)
      data1<-rbind(data1,d1)
      n1<-nrow(data1)
    } else {
      break
    }
  }
  data1<-data1[1:n,]
  data1$ldlc12m<-data1$ldlcm -47 
  data1$ldlc12<-data1$ldlc12m+rnorm(n,0,sd = 6)
  data0$drug = 0
  data1$drug = 1
  data = rbind(data0,data1)
  data<-data%>%mutate(ldlc_change = ldlc12 - ldlcb)
  return(data)
}

#function13: reduced individual variation 14.2->5

simulation.function13.new<-function(sample.size = 20000,cutoff = 1000){
  n<-sample.size/2
  osp = 2
  n0<-0
  n1<-0
  data0<-NULL
  data1<-NULL
  #control arm
  for(i in 1:100){
    #print(i)
    if(n0<n){
      ldlcm=rnorm(n = n*osp,mean = 116,sd = 5)#mean distribution of LDLC
      d0<-data.frame(ldlcm,ldlcb=rnorm(n = n*osp,mean = ldlcm,sd = 16),
                     ldlc12=rnorm(n = n*osp,mean = ldlcm,sd = 16))#dataframe for control arm
      d0<-d0%>%filter(ldlcb<cutoff)
      data0<-rbind(data0,d0)
      n0<-nrow(data0)
    } else {
      break
    }
  }
  
  data0<-data0[1:n,]
  data0$ldlc12m<-data0$ldlcm
  
  #treatment arm
  for(i in 1:100){
    #print(i)
    if(n1<n){
      ldlcm=rnorm(n = n*osp,mean = 116,sd = 5)
      d1<-data.frame(ldlcm,ldlcb=rnorm(n = n*osp,mean = ldlcm,sd = 16))
      d1<-d1%>%filter(ldlcb<cutoff)
      data1<-rbind(data1,d1)
      n1<-nrow(data1)
    } else {
      break
    }
  }
  data1<-data1[1:n,]
  data1$ldlc12m<-data1$ldlcm -47
  data1$ldlc12<-data1$ldlc12m+rnorm(n,0,sd = 16)
  data0$drug = 0
  data1$drug = 1
  data = rbind(data0,data1)
  data<-data%>%mutate(ldlc_change = ldlc12 - ldlcb)
  return(data)
}

########## Normal distribution with Heterogeity(Distribution#2) ############

simulation.function21.new<-function(sample.size = 20000,cutoff = 1000){
  n<-sample.size/2
  osp = 2
  n0<-0
  n1<-0
  data0<-NULL
  data1<-NULL
  #control arm
  for(i in 1:100){
    #print(i)
    if(n0<n){
      ldlcm=rnorm(n = n*osp,mean = 116,sd = 14.2)#mean distribution of LDLC
      d0<-data.frame(ldlcm,ldlcb=rnorm(n = n*osp,mean = ldlcm,sd = 16),
                     ldlc12=rnorm(n = n*osp,mean = ldlcm,sd = 16))#dataframe for control arm
      d0<-d0%>%filter(ldlcb<cutoff)
      data0<-rbind(data0,d0)
      n0<-nrow(data0)
    } else {
      break
    }
  }
  
  data0<-data0[1:n,]
  data0$ldlc12m<-data0$ldlcm
  
  #treatment arm
  for(i in 1:100){
    #print(i)
    if(n1<n){
      ldlcm=rnorm(n = n*osp,mean = 116,sd = 14.2)
      d1<-data.frame(ldlcm,ldlcb=rnorm(n = n*osp,mean = ldlcm,sd = 16))
      d1<-d1%>%filter(ldlcb<cutoff)
      data1<-rbind(data1,d1)
      n1<-nrow(data1)
    } else {
      break
    }
  }
  data1<-data1[1:n,]
  data1$ldlc12m<-data1$ldlcm - 0.7*data1$ldlcm +30
  data1$ldlc12<-data1$ldlc12m+rnorm(n,0,sd = 16)
  data0$drug = 0
  data1$drug = 1
  data = rbind(data0,data1)
  data<-data%>%mutate(ldlc_change = ldlc12 - ldlcb,ldlc_change.m = ldlc12m - ldlcm)
  return(data)
}

#16—>6
simulation.function22.new<-function(sample.size = 20000,cutoff = 1000){
  n<-sample.size/2
  osp = 2
  n0<-0
  n1<-0
  data0<-NULL
  data1<-NULL
  #control arm
  for(i in 1:100){
    #print(i)
    if(n0<n){
      ldlcm=rnorm(n = n*osp,mean = 116,sd = 14.2)#mean distribution of LDLC
      d0<-data.frame(ldlcm,ldlcb=rnorm(n = n*osp,mean = ldlcm,sd = 6),
                     ldlc12=rnorm(n = n*osp,mean = ldlcm,sd = 6))#dataframe for control arm
      d0<-d0%>%filter(ldlcb<cutoff)
      data0<-rbind(data0,d0)
      n0<-nrow(data0)
    } else {
      break
    }
  }
  
  data0<-data0[1:n,]
  data0$ldlc12m<-data0$ldlcm
  
  #treatment arm
  for(i in 1:100){
    #print(i)
    if(n1<n){
      ldlcm=rnorm(n = n*osp,mean = 116,sd = 14.2)
      d1<-data.frame(ldlcm,ldlcb=rnorm(n = n*osp,mean = ldlcm,sd = 6))
      d1<-d1%>%filter(ldlcb<cutoff)
      data1<-rbind(data1,d1)
      n1<-nrow(data1)
    } else {
      break
    }
  }
  data1<-data1[1:n,]
  data1$ldlc12m<-data1$ldlcm - 0.7*data1$ldlcm +30
  data1$ldlc12<-data1$ldlc12m+rnorm(n,0,sd = 6)
  data0$drug = 0
  data1$drug = 1
  data = rbind(data0,data1)
  data<-data%>%mutate(ldlc_change = ldlc12 - ldlcb)
  return(data)
}

#function23: reduced individual variation 14.2->5

simulation.function23.new<-function(sample.size = 20000,cutoff = 1000){
  n<-sample.size/2
  osp = 2
  n0<-0
  n1<-0
  data0<-NULL
  data1<-NULL
  #control arm
  for(i in 1:100){
    #print(i)
    if(n0<n){
      ldlcm=rnorm(n = n*osp,mean = 116,sd = 5)#mean distribution of LDLC
      d0<-data.frame(ldlcm,ldlcb=rnorm(n = n*osp,mean = ldlcm,sd = 16),
                     ldlc12=rnorm(n = n*osp,mean = ldlcm,sd = 16))#dataframe for control arm
      d0<-d0%>%filter(ldlcb<cutoff)
      data0<-rbind(data0,d0)
      n0<-nrow(data0)
    } else {
      break
    }
  }
  
  data0<-data0[1:n,]
  data0$ldlc12m<-data0$ldlcm
  
  #treatment arm
  for(i in 1:100){
    #print(i)
    if(n1<n){
      ldlcm=rnorm(n = n*osp,mean = 116,sd = 5)
      d1<-data.frame(ldlcm,ldlcb=rnorm(n = n*osp,mean = ldlcm,sd = 16))
      d1<-d1%>%filter(ldlcb<cutoff)
      data1<-rbind(data1,d1)
      n1<-nrow(data1)
    } else {
      break
    }
  }
  data1<-data1[1:n,]
  data1$ldlc12m<-data1$ldlcm - 0.7*data1$ldlcm +30
  data1$ldlc12<-data1$ldlc12m+rnorm(n,0,sd = 16)
  data0$drug = 0
  data1$drug = 1
  data = rbind(data0,data1)
  data<-data%>%mutate(ldlc_change = ldlc12 - ldlcb)
  return(data)
}

