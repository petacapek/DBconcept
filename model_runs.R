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
mar<-read.csv("Data/Marstorp1999.csv", sep=',')

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
good_sub<-function(x){
    p<-x
    names(p)<-c("v", "k", "f", "m", "Y", "fr", "fd")
    #Initial Br and Bs
    Bs_i<-mar$DNAinit[1]/p[["fd"]]
    #Br_i<-(mar$Cmicinit[1]-p[["fs"]]*Bs_i)/p[["fr"]]
    #Simulations
    yhat_all<-as.data.frame(ode(y=c(G=mar$Sinit[1], Br=0, Bs=Bs_i, CO2=0),
                                func = SUBmodel, parms=p,
                                times = as.numeric(mar$Time)))
    #Selecting measured variables
    yhat<-yhat_all[, c("time", "CO2", "DNA", "CFC14", "CFC", "G")]
    #Long format
    Yhat<-melt(yhat, id.vars=c("time"))
    #Observations
    Yhat$obs<-c(as.numeric(mar$CO212cumul), as.numeric(mar$DNA), as.numeric(mar$Cmic14),
                as.numeric(mar$Cmic12+mar$Cmic14), as.numeric(mar$S))
    Gfit<-Yhat %>% group_by(variable) %>% summarise(SSres=sum(((obs-value)^2), na.rm = T),
                                                    SStot=sum(((obs-mean(obs, na.rm = T))^2), na.rm = T),
                                                    ll=-sum(((obs-value)^2), na.rm = T)/2/(sd(obs, na.rm = T)^2))
    Gfit$R2<-with(Gfit, 1-SSres/SStot)
    Gfit$N<-length(p)
    Gfit$AIC<-with(Gfit, 2*N-2*ll)
    
    #Fine temporal scale fo r graphs
    yhat_all_fine<-as.data.frame(ode(y=c(G=mar$Sinit[1], Br=0, Bs=Bs_i, CO2=0),
                                     func = SUBmodel, parms=p,
                                     times = seq(0, 8.5, by=0.1)))
    Yhat_all_fine<-melt(yhat_all_fine, id.vars=c("time"))
    
    rsq_out<-list(Yhat=Yhat, Gfit=Gfit, Yhat_fine = Yhat_all_fine)
    
    return(rsq_out)
  }
  
#Read parameters estimated in python
marstorp_optpar<-as.numeric(read.csv("parameters/marstorp_optpars.csv", header = F))
    
Marstorp_fit<-good_sub(marstorp_optpar)
as.data.frame(Marstorp_fit$Gfit)

#Figure
Marstorp_fit$Yhat$variable2<-Marstorp_fit$Yhat$variable
levels(Marstorp_fit$Yhat$variable2)<-c("CO[2]", "DNA", "MB^{14}~C", "MBC", "Glucose")

Marstorp_fit$Yhat_fine$variable2<-Marstorp_fit$Yhat_fine$variable
levels(Marstorp_fit$Yhat_fine$variable2)<-c("Glucose", "Br", "Bs", "CO[2]", "MBC", "MB^{14}~C", "DNA", "Respiration")

ggplot(subset(Marstorp_fit$Yhat, variable=="CO2" | variable=="DNA" | variable=="CFC14" | 
                variable=="CFC" | variable=="G"), aes(time, obs))+
  geom_point(cex=6, pch=21, fill="grey")+
  geom_line(data=subset(Marstorp_fit$Yhat_fine, variable=="CO2" | variable=="DNA" | variable=="CFC14" | 
                          variable=="CFC" | variable=="G"), aes(time, value), lwd=1.2, color="grey30")+theme_min+
  facet_wrap(~variable2, scales="free", labeller = label_parsed) + 
  ylab(expression(paste("Carbon pool (", mu, "mol ", g(DW)^{-1}, ")"))) +
  xlab("Time (days)")

    
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Monod model~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
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
    dCO2<- uptake*(1-Y)
    
    return(list(c(dB, dG, dCO2), CFC=CFC, CFC14=CFC14, DNA=DNA, r=r))
  })
}

#Goodness of fit
monodgood<-function(x){
  p<-x
  names(p)<-c("v", "k", "m", "Y", "kec", "kd")
  #Initial Br and Bs
  B_i<-mar$DNAinit[1]/p[["kd"]]
  
  #Simulations
  yhat_all<-as.data.frame(ode(y=c(B=B_i, G=mar$Sinit[1], CO2=0),
                              func = Monod, parms=p,
                              times = as.numeric(mar$Time)))
  #Selecting measured variables
  yhat<-yhat_all[, c("time", "CO2", "DNA", "CFC", "CFC14", "G")]
  #Long format
  Yhat<-melt(yhat, id.vars=c("time"))
  #Observations
  Yhat$obs<-c(as.numeric(mar$CO212cumul), as.numeric(mar$DNA), as.numeric(mar$Cmic12+mar$Cmic14),
               as.numeric(mar$Cmic14), as.numeric(mar$S))
  Gfit<-Yhat %>% group_by(variable) %>% summarise(SSres=sum(((obs-value)^2), na.rm = T),
                                                  SStot=sum(((obs-mean(obs, na.rm = T))^2), na.rm = T),
                                                  ll=-sum(((obs-value)^2), na.rm = T)/2/(sd(obs, na.rm = T)^2))
  Gfit$R2<-with(Gfit, 1-SSres/SStot)
  Gfit$N<-length(p)
  Gfit$AIC<-with(Gfit, 2*N-2*ll)
  
  #Fine temporal scale for graphs
  yhat_all_fine<-as.data.frame(ode(y=c(B=B_i, G=mar$Sinit[1], CO2=0),
                                   func = Monod, parms=p,
                                   times = seq(0, 8.5, by=0.1)))
  Yhat_all_fine<-melt(yhat_all_fine, id.vars=c("time"))
  
  rsq_out<-list(Yhat=Yhat, Gfit=Gfit, Yhat_fine = Yhat_all_fine)
  
  return(rsq_out)
}

#Read parameters estimated in python
marstorp_monodopt<-as.numeric(read.csv("parameters/marstorp_monodpars.csv", header = F))

Marstorp_monodfit<-monodgood(marstorp_monodopt)
as.data.frame(Marstorp_monodfit$Gfit)

#Figure
Marstorp_monodfit$Yhat_fine$variable2<-Marstorp_monodfit$Yhat_fine$variable
levels(Marstorp_monodfit$Yhat_fine$variable2)<-c("B", "Glucose",  "CO[2]", "MBC", "MB^{14}~C", "DNA", "Respiration")

Marstorp_fita<-subset(Marstorp_fit$Yhat_fine, variable=="CO2" | variable=="DNA" | variable=="CFC14" | 
                        variable=="CFC" | variable=="G")
Marstorp_fita$Model<-c("Sub-microbial")
Marstorp_fitb<-subset(Marstorp_monodfit$Yhat_fine, variable=="CO2" | variable=="DNA" | variable=="CFC14" | 
                        variable=="CFC" | variable=="G")
Marstorp_fitb$Model<-c("Monod")

Marstorp_fit_all<-rbind(Marstorp_fita, Marstorp_fitb)

ggplot(subset(Marstorp_fit$Yhat, variable=="CO2" | variable=="DNA" | variable=="CFC14" | 
                variable=="CFC" | variable=="G"), aes(time, obs))+
  geom_point(cex=6, pch=21, fill="grey")+
  geom_line(data=Marstorp_fit_all, aes(time, value, color=Model), lwd=1.2)+theme_min+
  facet_wrap(~variable2, scales="free", labeller = label_parsed) + 
  ylab(expression(paste("Carbon pool (", mu, "mol ", g(DW)^{-1}, ")"))) +
  xlab("Time (days)") + scale_color_manual(values = c("indianred", "grey30")) +
  theme(legend.position = c(0.85, 0.1))
        
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Statistics~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Log Likelihood ratio test
-2*(Marstorp_monodfit$Gfit$ll-Marstorp_fit$Gfit$ll)

round(pchisq(-2*(Marstorp_monodfit$Gfit$ll-Marstorp_fit$Gfit$ll), df=(length(marstorp_optpar)-length(marstorp_monodopt)),
       lower.tail = F), 3)

#F test - based on residual sum of squares, number of parameters and number of measurements
##Monod model
###residual sum of squares
M1ss = sum((Marstorp_monodfit$Yhat$obs-Marstorp_monodfit$Yhat$value)^2, na.rm = T)
###number of parameters 
M1p = length(marstorp_monodopt)

##Sub-microbial model
###residual sum of squares
M2ss = sum((Marstorp_fit$Yhat$obs-Marstorp_fit$Yhat$value)^2, na.rm = T)
###number of parameters 
M2p = length(marstorp_optpar)

###total number of measurements
nt = nrow(Marstorp_fit$Yhat[!is.na(Marstorp_fit$Yhat), ])

####F value =  (M1ss - M2ss)*(nt - M2p)/M2ss/(M2p - M1p)
(M1ss - M2ss)*(nt - M2p)/M2ss/(M2p - M1p)

####associated p value
pf(q=(M1ss - M2ss)*(nt - M2p)/M2ss/(M2p - M1p), 
   df1=(M2p - M1p), 
   df2=(nt - M2p), 
   lower.tail=F)

