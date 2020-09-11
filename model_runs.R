#LIBRARIES
library(ggplot2)
library(openxlsx)
library(FME)
library(DEoptim)
library(dplyr)
library(reshape)
library(deSolve)
library(minpack.lm)
library(ABCoptim)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
##GGPLOT THEME
theme_min<-theme(axis.text.x=element_text(vjust=0.2, size=18, colour="black"),
                 axis.text.y=element_text(hjust=0.2, size=18, colour="black"),
                 axis.title=element_text(size=18, colour="black"),
                 axis.line=element_line(size=0.5, colour="black"),
                 strip.text=element_text(size=18, face="bold"),
                 axis.ticks=element_line(size=1, colour="black"),
                 axis.ticks.length=unit(-0.05, "cm"),
                 panel.background=element_rect(colour="black", fill="white"),
                 panel.grid=element_line(linetype=0),
                 legend.text=element_text(size=14, colour="black"),
                 legend.title=element_text(size=14, colour="black"),
                 legend.position=c("right"),
                 legend.key.size=unit(1, "cm"),
                 strip.background=element_rect(fill="grey98", colour="black"),
                 legend.key=element_rect(fill="white", size=1.2),
                 legend.spacing=unit(0.5, "cm"),
                 plot.title=element_text(size=18, face="bold", hjust=-0.05))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
############################################################################################################################
############################################Marstorp and Witter, 1999#######################################################
############################################################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#DATA
mar<-read.csv("../Marstorp1999.csv", sep=',')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Sub-microbial model~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Model definition
SUBmodel<-function(time, state, pars){
    with(as.list(c(state, pars)),{
      #uptake
      uptake=v*G*Bs/(k+G)
      #Transfer function
      transfer=f*Br
      #death rate
      death=m*Bs
      #Respiration rate
      r=transfer*(1-Y)
      #Chloroform labile C and DNA
      CFC=fr*Br+(mar$Cmicinit[1]*fd/mar$DNAinit[1])*Bs
      CFC14=(fr*Br+(mar$Cmicinit[1]*fd/mar$DNAinit[1])*Bs)-mar$Cmicinit[1]
      DNA=fd*Bs
      
      
      #States
      dG<- - uptake
      dBr<- uptake - transfer
      dBs<- transfer*Y - death
      dCO2<- transfer*(1-Y)
      
      return(list(c(dG, dBr, dBs, dCO2), CFC=CFC, CFC14=CFC14, DNA=DNA, r=r))
    })
  }
  
#Goodness of fit
good<-function(x){
    p<-x
    names(p)<-c("v", "k", "f", "m", "Y", "fr", "fd")
    #Initial Br and Bs
    Bs_i<-mar$DNAinit[1]/p[["fd"]]
    #Br_i<-(mar$Cmicinit[1]-p[["fs"]]*Bs_i)/p[["fr"]]
    #Simulations
    yhat_all<-as.data.frame(ode(y=c(G=mar$Sinit[1], Br=0, Bs=Bs_i, CO2=0),
                                func = Twopool, parms=p,
                                times = as.numeric(mar$Time)))
    #Selecting measured variables
    yhat<-yhat_all[, c("time", "CO2", "DNA", "CFC", "G")]
    #Long format
    Yhat<-melt(yhat, id.vars=c("time"))
    #Observations
    Yhat$obs<-c(as.numeric(mar$CO212cumul), as.numeric(mar$DNA), #, as.numeric(mar$Cmic14)
                as.numeric(mar$Cmic12+mar$Cmic14), as.numeric(mar$S))
    Gfit<-Yhat %>% group_by(variable) %>% summarise(SSres=sum(((obs-value)^2), na.rm = T),
                                                    SStot=sum(((obs-mean(obs, na.rm = T))^2), na.rm = T),
                                                    ll=-sum(((obs-value)^2), na.rm = T)/2/(sd(obs, na.rm = T)^2))
    Gfit$R2<-with(Gfit, 1-SSres/SStot)
    Gfit$N<-length(p)
    Gfit$AIC<-with(Gfit, 2*N-2*ll)
    
    #Fine temporal scale fo r graphs
    yhat_all_fine<-as.data.frame(ode(y=c(G=mar$Sinit[1], Br=0, Bs=Bs_i, CO2=0),
                                     func = Twopool, parms=p,
                                     times = seq(0, 8.5, by=0.1)))
    Yhat_all_fine<-melt(yhat_all_fine, id.vars=c("time"))
    
    rsq_out<-list(Yhat=Yhat, Gfit=Gfit, Yhat_fine = Yhat_all_fine)
    
    return(rsq_out)
  }
  
#Read parameters estimated in python
opt_par<-as.numeric(read.csv("/home/capekp00/opt_parsa.csv", header = F))
    
Twopool_fit<-good(opt_par)
as.data.frame(Twopool_fit$Gfit)
    
#Figure
ggplot(Twopool_fit$Yhat, aes(time, obs))+geom_point(cex=6, pch=21, fill="grey")+
    geom_line(data=Twopool_fit$Yhat_fine, aes(time, value))+
    facet_wrap(~variable, scales="free")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Monod model~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Model definition
Monod<-function(time, state, pars){
  with(as.list(c(state, pars)),{
    #glucose uptake
    uptake=v*G*B/(k + G)
    #Decay rate
    decay=m*B
    #Respiration rate
    r=(1-Y)*uptake + decay
    #Chloroform labile C and DNA
    CFC=kec*B
    CFC14=kec*B-mar$Cmicinit[1]
    DNA=kd*B
    
    
    #States
    dB<- uptake*Y - decay
    dG<- - uptake
    
    return(list(c(dB, dG), CFC=CFC, CFC14=CFC14, DNA=DNA, r=r))
  })
}

#Goodness of fit
good<-function(x){
  p<-x
  names(p)<-c("v", "k", "m", "Y", "kec", "kd")
  #Initial Br and Bs
  B_i<-mar$Cmicinit[1]/p[["kec"]]
  
  #Simulations
  yhat_all<-as.data.frame(ode(y=c(B=B_i, G=mar$Sinit[1]),
                              func = Monod, parms=p,
                              times = as.numeric(mar$Time)))
  #Selecting measured variables
  yhat<-yhat_all[, c("time", "r", "DNA", "CFC", "CFC14", "G")]
  #Long format
  Yhat<-melt(yhat, id.vars=c("time"))
  #Observations
  Yhat$obs<-c(as.numeric(mar$CO212), as.numeric(mar$DNA), as.numeric(mar$Cmic12+mar$Cmic14),
               as.numeric(mar$Cmic14), as.numeric(mar$S))
  Gfit<-Yhat %>% group_by(variable) %>% summarise(SSres=sum(((obs-value)^2), na.rm = T),
                                                  SStot=sum(((obs-mean(obs, na.rm = T))^2), na.rm = T),
                                                  ll=-sum(((obs-value)^2), na.rm = T)/2/(sd(obs, na.rm = T)^2))
  Gfit$R2<-with(Gfit, 1-SSres/SStot)
  Gfit$N<-length(p)
  Gfit$AIC<-with(Gfit, 2*N-2*ll)
  
  #Fine temporal scale for graphs
  yhat_all_fine<-as.data.frame(ode(y=c(B=B_i, G=mar$Sinit[1]),
                                   func = Monod, parms=p,
                                   times = seq(0, 8.5, by=0.1)))
  Yhat_all_fine<-melt(yhat_all_fine, id.vars=c("time"))
  
  rsq_out<-list(Yhat=Yhat, Gfit=Gfit, Yhat_fine = Yhat_all_fine)
  
  return(rsq_out)
}

#Read parameters estimated in python
Monod_opt<-as.numeric(read.csv("../monod_opt.csv", header = F))

Monod_fit<-good(Monod_opt)
as.data.frame(Monod_fit$Gfit)

as.data.frame(Monod_fit$Gfit$R2)
as.data.frame(Twopool_fit$Gfit$R2)

#Figure
ggplot(Monod_fit$Yhat, aes(time, obs))+geom_point(cex=6, pch=21, fill="grey")+
  geom_line(data=Monod_fit$Yhat_fine[, c("time", "variable", "value")], 
            aes(time, value))+
  facet_wrap(~variable, scales="free")+
  geom_line(data=Twopool_fit$Yhat_fine[, c("time", "variable", "value")], 
            aes(time, value), color="blue")
        
######################################Statistics####################################
#F test - based on residual sum of squares, number of parameters 
#and number of measurements
##Monod model
###residual sum of squares
M1ss = sum((Monod_fit$Yhat$obs-Monod_fit$Yhat$value)^2, na.rm = T)
###number of parameters 
M1p = length(Monod_opt)

##Two pool model
###residual sum of squares
M2ss = sum((Twopool_fit$Yhat$obs-Twopool_fit$Yhat$value)^2, na.rm = T)
###number of parameters 
M2p = length(opt_par)

###total number of measurements
nt = nrow(Twopool_fit$Yhat[!is.na(Twopool_fit$Yhat), ])

####F value =  (M1ss - M2ss)*(nt - M2p)/M2ss/(M2p - M1p)
(M1ss - M2ss)*(nt - M2p)/M2ss/(M2p - M1p)

####associated p value
pf(q=(M1ss - M2ss)*(nt - M2p)/M2ss/(M2p - M1p), 
   df1=(M2p - M1p), 
   df2=(nt - M2p), 
   lower.tail=F)

