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

#Model definition
DEBmodel<-function(time, state, pars){
  with(as.list(c(state, pars)),{
    #Define fluxes
    ##scaling function for substrate
    f=G/(500+G)
    ##growth rate
    growth=(v*e-m*g)/(m+g)
    ##CO2 yield
    Yco2=((v*f/Im)+(m*g/Im)+max(g*growth/Im,0))*ce/f
    
    #Chloroform labile C and DNA
    Cwcfc=mar$Cmicinit[1]*Cwdna/mar$DNAinit[1]
    CFC=(Cwcfc+Cecfc*e)*w
    CFC14=(Cwcfc+Cecfc*e)*w-mar$Cmicinit[1]
    DNA=Cwdna*w
    
    
    #States
    #Define derivatives
    dG=-f*w*Im
    de=v*f-v*e
    dw=growth*w
    dCO2=f*w*Im*Yco2
    
    return(list(c(dG, de, dw, dCO2), CFC=CFC, CFC14=CFC14, DNA=DNA))
  })
}

#Goodness of fit
good_DEB<-function(x){
  p<-x
  names(p)<-c("Im", "v", "m", "g", "ce", "Cwdna", "Cecfc")
  #Initial Br and Bs
  w_i<-mar$DNAinit[1]/p[["Cwdna"]]
  #Br_i<-(mar$Cmicinit[1]-p[["fs"]]*Bs_i)/p[["fr"]]
  #Simulations
  yhat_all<-as.data.frame(ode(y=c(G=mar$Sinit[1], e=0, w=w_i, CO2=0),
                              func = DEBmodel, parms=p,
                              times = as.numeric(mar$Time)*24))
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
  yhat_all_fine<-as.data.frame(ode(y=c(G=mar$Sinit[1], e=0, w=w_i, CO2=0),
                                   func = DEBmodel, parms=p,
                                   times = seq(0, 8.5, by=0.1)*24))
  Yhat_all_fine<-melt(yhat_all_fine, id.vars=c("time"))
  
  rsq_out<-list(Yhat=Yhat, Gfit=Gfit, Yhat_fine = Yhat_all_fine)
  
  return(rsq_out)
}

#Read parameters estimated in python
marstorp_optparDEB<-as.numeric(read.csv("parameters/marstorp_optparsDEB.csv", header = F))

Marstorp_fitDEB<-good_DEB(marstorp_optparDEB)
as.data.frame(Marstorp_fitDEB$Gfit)

#Figure
Marstorp_fitDEB$Yhat$variable2<-Marstorp_fitDEB$Yhat$variable
levels(Marstorp_fitDEB$Yhat$variable2)<-c("CO[2]", "DNA", "MB^{14}~C", "MBC", "Glucose")

Marstorp_fitDEB$Yhat_fine$variable2<-Marstorp_fitDEB$Yhat_fine$variable
levels(Marstorp_fitDEB$Yhat_fine$variable2)<-c("Glucose", "Br", "Bs", "CO[2]", "MBC", "MB^{14}~C", "DNA", "Respiration")

ggplot(subset(Marstorp_fitDEB$Yhat, variable=="CO2" | variable=="DNA" | variable=="CFC14" | 
                variable=="CFC" | variable=="G"), aes(time, obs))+
  geom_point(cex=6, pch=21, fill="grey")+
  geom_line(data=subset(Marstorp_fitDEB$Yhat_fine, variable=="CO2" | variable=="DNA" | variable=="CFC14" | 
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
-2*(Marstorp_monodfit$Gfit$ll-Marstorp_fitDEB$Gfit$ll)

round(pchisq(-2*(Marstorp_monodfit$Gfit$ll-Marstorp_fitDEB$Gfit$ll), df=(length(marstorp_optpar)-length(marstorp_monodopt)),
       lower.tail = F), 3)

#F test - based on residual sum of squares, number of parameters and number of measurements
##Monod model
###residual sum of squares
M1ss = sum((Marstorp_monodfit$Yhat$obs-Marstorp_monodfit$Yhat$value)^2, na.rm = T)
###number of parameters 
M1p = length(marstorp_monodopt)

##Sub-microbial model
###residual sum of squares
M2ss = sum((Marstorp_fitDEB$Yhat$obs-Marstorp_fitDEB$Yhat$value)^2, na.rm = T)
###number of parameters 
M2p = length(marstorp_optpar)

###total number of measurements
  nt = nrow(Marstorp_fitDEB$Yhat[!is.na(Marstorp_fit$Yhat), ])

####F value =  (M1ss - M2ss)*(nt - M2p)/M2ss/(M2p - M1p)
(M1ss - M2ss)*(nt - M2p)/M2ss/(M2p - M1p)

####associated p value
pf(q=(M1ss - M2ss)*(nt - M2p)/M2ss/(M2p - M1p), 
   df1=(M2p - M1p), 
   df2=(nt - M2p), 
   lower.tail=F)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
############################################################################################################################
#############################################Santruckova et al, 2004########################################################
############################################################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#DATA
has<-read.csv("Data/Santruckova2004.csv", sep=',')

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
    CFC=fr*Br+fs*Bs
    CFC14=(fr*Br+fs*Bs)-has$Cmicinit[1]
    
    #States
    dG<- - uptake
    dBr<- uptake - transfer
    dBs<- transfer*Y - death
    dCO2<- transfer*(1-Y)
    
    return(list(c(dG, dBr, dBs, dCO2), CFC=CFC, CFC14=CFC14, r=r))
  })
}

#Goodness of fit
good_sub<-function(x){
  p<-x
  names(p)<-c("v", "k", "f", "m", "Y", "fr", "fs")
  #Initial Br and Bs
  Bs_i<-has$Cmicinit[1]/p[["fs"]]
  #Br_i<-(mar$Cmicinit[1]-p[["fs"]]*Bs_i)/p[["fr"]]
  #Simulations
  yhat_all<-as.data.frame(ode(y=c(G=has$Sinit[1], Br=0, Bs=Bs_i, CO2=0),
                              func = SUBmodel, parms=p,
                              times = as.numeric(has$Time)))
  #Selecting measured variables
  yhat<-yhat_all[, c("time", "CO2", "CFC14", "G")]
  #Long format
  Yhat<-melt(yhat, id.vars=c("time"))
  #Observations
  Yhat$obs<-c(as.numeric(has$CO214cumul), as.numeric(has$Cmic14),
              as.numeric(has$S))
  Gfit<-Yhat %>% group_by(variable) %>% summarise(SSres=sum(((obs-value)^2), na.rm = T),
                                                  SStot=sum(((obs-mean(obs, na.rm = T))^2), na.rm = T),
                                                  ll=-sum(((obs-value)^2), na.rm = T)/2/(sd(obs, na.rm = T)^2))
  Gfit$R2<-with(Gfit, 1-SSres/SStot)
  Gfit$N<-length(p)
  Gfit$AIC<-with(Gfit, 2*N-2*ll)
  
  #Fine temporal scale fo r graphs
  yhat_all_fine<-as.data.frame(ode(y=c(G=has$Sinit[1], Br=0, Bs=Bs_i, CO2=0),
                                   func = SUBmodel, parms=p,
                                   times = seq(0, 3, by=0.1)))
  Yhat_all_fine<-melt(yhat_all_fine, id.vars=c("time"))
  
  rsq_out<-list(Yhat=Yhat, Gfit=Gfit, Yhat_fine = Yhat_all_fine)
  
  return(rsq_out)
}

#Read parameters estimated in python
hasan_optpar<-as.numeric(read.csv("parameters/hasan_optpars.csv", header = F))

Hasan_fit<-good_sub(hasan_optpar)
as.data.frame(Hasan_fit$Gfit)

#Figure
Hasan_fit$Yhat$variable2<-Hasan_fit$Yhat$variable
levels(Hasan_fit$Yhat$variable2)<-c("CO[2]", "MB^{14}~C", "Glucose")

Hasan_fit$Yhat_fine$variable2<-Hasan_fit$Yhat_fine$variable
levels(Hasan_fit$Yhat_fine$variable2)<-c("Glucose", "Br", "Bs", "CO[2]", "MBC", "MB^{14}~C", "Respiration")

ggplot(subset(Hasan_fit$Yhat, variable=="CO2" | variable=="CFC14" | 
                variable=="G"), aes(time, obs))+
  geom_point(cex=6, pch=21, fill="grey")+
  geom_line(data=subset(Hasan_fit$Yhat_fine, variable=="CO2" | variable=="CFC14" | 
                          variable=="G"), aes(time, value), lwd=1.2, color="grey30")+theme_min+
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
    CFC14=kec*B-has$Cmicinit[1]
    
    #States
    dB<- uptake*Y - decay
    dG<- - uptake
    dCO2<- uptake*(1-Y)
    
    return(list(c(dB, dG, dCO2), CFC=CFC, CFC14=CFC14, r=r))
  })
}

#Goodness of fit
monodgood<-function(x){
  p<-x
  names(p)<-c("v", "k", "m", "Y", "kec")
  #Initial Br and Bs
  B_i<-has$Cmicinit[1]/p[["kec"]]
  
  #Simulations
  yhat_all<-as.data.frame(ode(y=c(B=B_i, G=has$Sinit[1], CO2=0),
                              func = Monod, parms=p,
                              times = as.numeric(has$Time)))
  #Selecting measured variables
  yhat<-yhat_all[, c("time", "CO2", "CFC14", "G")]
  #Long format
  Yhat<-melt(yhat, id.vars=c("time"))
  #Observations
  Yhat$obs<-c(as.numeric(has$CO214cumul), 
              as.numeric(has$Cmic14), as.numeric(has$S))
  Gfit<-Yhat %>% group_by(variable) %>% summarise(SSres=sum(((obs-value)^2), na.rm = T),
                                                  SStot=sum(((obs-mean(obs, na.rm = T))^2), na.rm = T),
                                                  ll=-sum(((obs-value)^2), na.rm = T)/2/(sd(obs, na.rm = T)^2))
  Gfit$R2<-with(Gfit, 1-SSres/SStot)
  Gfit$N<-length(p)
  Gfit$AIC<-with(Gfit, 2*N-2*ll)
  
  #Fine temporal scale for graphs
  yhat_all_fine<-as.data.frame(ode(y=c(B=B_i, G=has$Sinit[1], CO2=0),
                                   func = Monod, parms=p,
                                   times = seq(0, 3, by=0.1)))
  Yhat_all_fine<-melt(yhat_all_fine, id.vars=c("time"))
  
  rsq_out<-list(Yhat=Yhat, Gfit=Gfit, Yhat_fine = Yhat_all_fine)
  
  return(rsq_out)
}

#Read parameters estimated in python
hasan_monodopt<-as.numeric(read.csv("parameters/hasan_monodpars.csv", header = F))

Hasan_monodfit<-monodgood(hasan_monodopt)
as.data.frame(Hasan_monodfit$Gfit)

#Figure
Hasan_monodfit$Yhat_fine$variable2<-Hasan_monodfit$Yhat_fine$variable
levels(Hasan_monodfit$Yhat_fine$variable2)<-c("B", "Glucose",  "CO[2]", "MBC", "MB^{14}~C", "Respiration")

Hasan_fita<-subset(Hasan_fit$Yhat_fine, variable=="CO2" | variable=="CFC14" | 
                        variable=="G")
Hasan_fita$Model<-c("Sub-microbial")
Hasan_fitb<-subset(Hasan_monodfit$Yhat_fine, variable=="CO2" | variable=="CFC14" | 
                        variable=="G")
Hasan_fitb$Model<-c("Monod")

Hasan_fit_all<-rbind(Hasan_fita, Hasan_fitb)

ggplot(subset(Hasan_fit$Yhat, variable=="CO2" | variable=="CFC14" | 
                variable=="G"), aes(time, obs))+
  geom_point(cex=6, pch=21, fill="grey")+
  geom_line(data=Hasan_fit_all, aes(time, value, color=Model), lwd=1.2)+theme_min+
  facet_wrap(~variable2, scales="free", labeller = label_parsed) + 
  ylab(expression(paste("Carbon pool (", mu, "mol ", g(DW)^{-1}, ")"))) +
  xlab("Time (days)") + scale_color_manual(values = c("indianred", "grey30")) +
  theme(legend.position = c(0.85, 0.8))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Statistics~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Log Likelihood ratio test
-2*(Hasan_monodfit$Gfit$ll-Hasan_fit$Gfit$ll)

round(pchisq(-2*(Hasan_monodfit$Gfit$ll-Hasan_fit$Gfit$ll), df=(length(hasan_optpar)-length(hasan_monodopt)),
             lower.tail = F), 3)

#F test - based on residual sum of squares, number of parameters and number of measurements
##Monod model
###residual sum of squares
M1ss = sum((Hasan_monodfit$Yhat$obs-Hasan_monodfit$Yhat$value)^2, na.rm = T)
###number of parameters 
M1p = length(hasan_monodopt)

##Sub-microbial model
###residual sum of squares
M2ss = sum((Hasan_fit$Yhat$obs-Hasan_fit$Yhat$value)^2, na.rm = T)
###number of parameters 
M2p = length(hasan_optpar)

###total number of measurements
nt = nrow(Hasan_fit$Yhat[!is.na(Hasan_fit$Yhat), ])

####F value =  (M1ss - M2ss)*(nt - M2p)/M2ss/(M2p - M1p)
(M1ss - M2ss)*(nt - M2p)/M2ss/(M2p - M1p)

####associated p value
pf(q=(M1ss - M2ss)*(nt - M2p)/M2ss/(M2p - M1p), 
   df1=(M2p - M1p), 
   df2=(nt - M2p), 
   lower.tail=F)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
############################################################################################################################
#################################################Ziegler et al, 2004########################################################
############################################################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#DATA
zieg<-read.csv("Data/Ziegler2005.csv", sep=',')

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
    CFC=fr*Br+fs*Bs
    
    #States
    dG<- - uptake
    dBr<- uptake - transfer
    dBs<- transfer*Y - death
    dCO2<- transfer*(1-Y)
    
    return(list(c(dG, dBr, dBs, dCO2), CFC=CFC, r=r))
  })
}

#Goodness of fit
good_sub<-function(x){
  p<-x
  names(p)<-c("v", "k", "f", "m", "Y", "fr", "fs")
  #Initial Br and Bs
  Bs_i<-zieg$Cmicinit[1]/p[["fs"]]
  #Br_i<-(mar$Cmicinit[1]-p[["fs"]]*Bs_i)/p[["fr"]]
  #Simulations
  yhat_all<-as.data.frame(ode(y=c(G=zieg$Sinit[1], Br=0, Bs=Bs_i, CO2=0),
                              func = SUBmodel, parms=p,
                              times = as.numeric(zieg$Time)))
  #Selecting measured variables
  yhat<-yhat_all[, c("time", "CO2", "CFC", "G")]
  #Long format
  Yhat<-melt(yhat, id.vars=c("time"))
  #Observations
  Yhat$obs<-c(as.numeric(zieg$CO2), as.numeric(zieg$Cmic),
              as.numeric(zieg$S))
  Gfit<-Yhat %>% group_by(variable) %>% summarise(SSres=sum(((obs-value)^2), na.rm = T),
                                                  SStot=sum(((obs-mean(obs, na.rm = T))^2), na.rm = T),
                                                  ll=-sum(((obs-value)^2), na.rm = T)/2/(sd(obs, na.rm = T)^2))
  Gfit$R2<-with(Gfit, 1-SSres/SStot)
  Gfit$N<-length(p)
  Gfit$AIC<-with(Gfit, 2*N-2*ll)
  
  #Fine temporal scale fo r graphs
  yhat_all_fine<-as.data.frame(ode(y=c(G=zieg$Sinit[1], Br=0, Bs=Bs_i, CO2=0),
                                   func = SUBmodel, parms=p,
                                   times = seq(0, 2, by=0.1)))
  Yhat_all_fine<-melt(yhat_all_fine, id.vars=c("time"))
  
  rsq_out<-list(Yhat=Yhat, Gfit=Gfit, Yhat_fine = Yhat_all_fine)
  
  return(rsq_out)
}

#Read parameters estimated in python
ziegler_optpar<-as.numeric(read.csv("parameters/ziegler_optpars.csv", header = F))

Ziegler_fit<-good_sub(ziegler_optpar)
as.data.frame(Ziegler_fit$Gfit)

#Figure
Ziegler_fit$Yhat$variable2<-Ziegler_fit$Yhat$variable
levels(Ziegler_fit$Yhat$variable2)<-c("CO[2]", "MBC", "Glucose")

Ziegler_fit$Yhat_fine$variable2<-Ziegler_fit$Yhat_fine$variable
levels(Ziegler_fit$Yhat_fine$variable2)<-c("Glucose", "Br", "Bs", "CO[2]", "MBC", "Respiration")

ggplot(subset(Ziegler_fit$Yhat, variable=="CO2" | variable=="CFC" | 
                variable=="G"), aes(time, obs))+
  geom_point(cex=6, pch=21, fill="grey")+
  geom_line(data=subset(Ziegler_fit$Yhat_fine, variable=="CO2" | variable=="CFC" | 
                          variable=="G"), aes(time, value), lwd=1.2, color="grey30")+theme_min+
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
    
    #States
    dB<- uptake*Y - decay
    dG<- - uptake
    dCO2<- uptake*(1-Y)
    
    return(list(c(dB, dG, dCO2), CFC=CFC, r=r))
  })
}

#Goodness of fit
monodgood<-function(x){
  p<-x
  names(p)<-c("v", "k", "m", "Y", "kec")
  #Initial Br and Bs
  B_i<-zieg$Cmicinit[1]/p[["kec"]]
  
  #Simulations
  yhat_all<-as.data.frame(ode(y=c(B=B_i, G=zieg$Sinit[1], CO2=0),
                              func = Monod, parms=p,
                              times = as.numeric(zieg$Time)))
  #Selecting measured variables
  yhat<-yhat_all[, c("time", "CO2", "CFC", "G")]
  #Long format
  Yhat<-melt(yhat, id.vars=c("time"))
  #Observations
  Yhat$obs<-c(as.numeric(zieg$CO2), 
              as.numeric(zieg$Cmic), as.numeric(zieg$S))
  Gfit<-Yhat %>% group_by(variable) %>% summarise(SSres=sum(((obs-value)^2), na.rm = T),
                                                  SStot=sum(((obs-mean(obs, na.rm = T))^2), na.rm = T),
                                                  ll=-sum(((obs-value)^2), na.rm = T)/2/(sd(obs, na.rm = T)^2))
  Gfit$R2<-with(Gfit, 1-SSres/SStot)
  Gfit$N<-length(p)
  Gfit$AIC<-with(Gfit, 2*N-2*ll)
  
  #Fine temporal scale for graphs
  yhat_all_fine<-as.data.frame(ode(y=c(B=B_i, G=zieg$Sinit[1], CO2=0),
                                   func = Monod, parms=p,
                                   times = seq(0, 2, by=0.1)))
  Yhat_all_fine<-melt(yhat_all_fine, id.vars=c("time"))
  
  rsq_out<-list(Yhat=Yhat, Gfit=Gfit, Yhat_fine = Yhat_all_fine)
  
  return(rsq_out)
}

#Read parameters estimated in python
ziegler_monodopt<-as.numeric(read.csv("parameters/ziegler_monodpars.csv", header = F))

Ziegler_monodfit<-monodgood(ziegler_monodopt)
as.data.frame(Ziegler_monodfit$Gfit)

#Figure
Ziegler_monodfit$Yhat_fine$variable2<-Ziegler_monodfit$Yhat_fine$variable
levels(Ziegler_monodfit$Yhat_fine$variable2)<-c("B", "Glucose",  "CO[2]", "MBC", "Respiration")

Ziegler_fita<-subset(Ziegler_fit$Yhat_fine, variable=="CO2" | variable=="CFC" | 
                     variable=="G")
Ziegler_fita$Model<-c("Sub-microbial")
Ziegler_fitb<-subset(Ziegler_monodfit$Yhat_fine, variable=="CO2" | variable=="CFC" | 
                     variable=="G")
Ziegler_fitb$Model<-c("Monod")

Ziegler_fit_all<-rbind(Ziegler_fita, Ziegler_fitb)

ggplot(subset(Ziegler_fit$Yhat, variable=="CO2" | variable=="CFC" | 
                variable=="G"), aes(time, obs))+
  geom_point(cex=6, pch=21, fill="grey")+
  geom_line(data=Ziegler_fit_all, aes(time, value, color=Model), lwd=1.2)+theme_min+
  facet_wrap(~variable2, scales="free", labeller = label_parsed) + 
  ylab(expression(paste("Carbon pool (", mu, "mol ", g(DW)^{-1}, ")"))) +
  xlab("Time (days)") + scale_color_manual(values = c("indianred", "grey30")) +
  theme(legend.position = c(0.85, 0.8))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Statistics~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Log Likelihood ratio test
-2*(Ziegler_monodfit$Gfit$ll-Ziegler_fit$Gfit$ll)

round(pchisq(-2*(Ziegler_monodfit$Gfit$ll-Ziegler_fit$Gfit$ll), df=(length(ziegler_optpar)-length(ziegler_monodopt)),
             lower.tail = F), 3)

#F test - based on residual sum of squares, number of parameters and number of measurements
##Monod model
###residual sum of squares
M1ss = sum((Ziegler_monodfit$Yhat$obs-Ziegler_monodfit$Yhat$value)^2, na.rm = T)
###number of parameters 
M1p = length(ziegler_monodopt)

##Sub-microbial model
###residual sum of squares
M2ss = sum((Ziegler_fit$Yhat$obs-Ziegler_fit$Yhat$value)^2, na.rm = T)
###number of parameters 
M2p = length(ziegler_optpar)

###total number of measurements
nt = nrow(Ziegler_fit$Yhat[!is.na(Ziegler_fit$Yhat), ])

####F value =  (M1ss - M2ss)*(nt - M2p)/M2ss/(M2p - M1p)
(M1ss - M2ss)*(nt - M2p)/M2ss/(M2p - M1p)

####associated p value
pf(q=(M1ss - M2ss)*(nt - M2p)/M2ss/(M2p - M1p), 
   df1=(M2p - M1p), 
   df2=(nt - M2p), 
   lower.tail=F)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
############################################################################################################################
#################################################Blagodatskaya et al, 2014########################################################
############################################################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#DATA
bla<-read.csv("Data/Blagodatskaya2014.csv", sep=',')

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
    CFC14=(fr*Br+fs*Bs)-bla$Cmicinit[1]
    
    #States
    dG<- - uptake
    dBr<- uptake - transfer
    dBs<- transfer*Y - death
    dCO2<- transfer*(1-Y)
    
    return(list(c(dG, dBr, dBs, dCO2), CFC14=CFC14, r=r))
  })
}

#Goodness of fit
good_sub<-function(x){
  p<-x
  names(p)<-c("v", "k", "f", "m", "Y", "fr", "fs")
  #Initial Br and Bs
  Bs_i<-bla$Cmicinit[1]/p[["fs"]]
  #Br_i<-(mar$Cmicinit[1]-p[["fs"]]*Bs_i)/p[["fr"]]
  #Simulations
  yhat_all<-as.data.frame(ode(y=c(G=bla$Sinit[1], Br=0, Bs=Bs_i, CO2=0),
                              func = SUBmodel, parms=p,
                              times = as.numeric(bla$Time)))
  #Selecting measured variables
  yhat<-yhat_all[, c("time", "CO2", "CFC14")]
  #Long format
  Yhat<-melt(yhat, id.vars=c("time"))
  #Observations
  Yhat$obs<-c(as.numeric(bla$CO214), as.numeric(bla$Cmic14))
  Gfit<-Yhat %>% group_by(variable) %>% summarise(SSres=sum(((obs-value)^2), na.rm = T),
                                                  SStot=sum(((obs-mean(obs, na.rm = T))^2), na.rm = T),
                                                  ll=-sum(((obs-value)^2), na.rm = T)/2/(sd(obs, na.rm = T)^2))
  Gfit$R2<-with(Gfit, 1-SSres/SStot)
  Gfit$N<-length(p)
  Gfit$AIC<-with(Gfit, 2*N-2*ll)
  
  #Fine temporal scale fo r graphs
  yhat_all_fine<-as.data.frame(ode(y=c(G=bla$Sinit[1], Br=0, Bs=Bs_i, CO2=0),
                                   func = SUBmodel, parms=p,
                                   times = seq(0, 103, by=1)))
  Yhat_all_fine<-melt(yhat_all_fine, id.vars=c("time"))
  
  rsq_out<-list(Yhat=Yhat, Gfit=Gfit, Yhat_fine = Yhat_all_fine)
  
  return(rsq_out)
}

#Read parameters estimated in python
blag_optpar<-as.numeric(read.csv("parameters/blag_optpars.csv", header = F))

Blag_fit<-good_sub(blag_optpar)
as.data.frame(Blag_fit$Gfit)

#Figure
Blag_fit$Yhat$variable2<-Blag_fit$Yhat$variable
levels(Blag_fit$Yhat$variable2)<-c("CO[2]", "MBC")

Blag_fit$Yhat_fine$variable2<-Blag_fit$Yhat_fine$variable
levels(Blag_fit$Yhat_fine$variable2)<-c("Glucose", "Br", "Bs", "CO[2]", "MBC", "Respiration")

ggplot(subset(Blag_fit$Yhat, variable=="CO2" | variable=="CFC14"), aes(time, obs))+
  geom_point(cex=6, pch=21, fill="grey")+
  geom_line(data=subset(Blag_fit$Yhat_fine, variable=="CO2" | variable=="CFC14"), aes(time, value), lwd=1.2, color="grey30")+theme_min+
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
    CFC14=kec*B-bla$Cmicinit[1]
    
    #States
    dB<- uptake*Y - decay
    dG<- - uptake
    dCO2<- uptake*(1-Y)
    
    return(list(c(dB, dG, dCO2), CFC14=CFC14, r=r))
  })
}

#Goodness of fit
monodgood<-function(x){
  p<-x
  names(p)<-c("v", "k", "m", "Y", "kec")
  #Initial Br and Bs
  B_i<-bla$Cmicinit[1]/p[["kec"]]
  
  #Simulations
  yhat_all<-as.data.frame(ode(y=c(B=B_i, G=bla$Sinit[1], CO2=0),
                              func = Monod, parms=p,
                              times = as.numeric(bla$Time)))
  #Selecting measured variables
  yhat<-yhat_all[, c("time", "CO2", "CFC14")]
  #Long format
  Yhat<-melt(yhat, id.vars=c("time"))
  #Observations
  Yhat$obs<-c(as.numeric(bla$CO214), 
              as.numeric(bla$Cmic14))
  Gfit<-Yhat %>% group_by(variable) %>% summarise(SSres=sum(((obs-value)^2), na.rm = T),
                                                  SStot=sum(((obs-mean(obs, na.rm = T))^2), na.rm = T),
                                                  ll=-sum(((obs-value)^2), na.rm = T)/2/(sd(obs, na.rm = T)^2))
  Gfit$R2<-with(Gfit, 1-SSres/SStot)
  Gfit$N<-length(p)
  Gfit$AIC<-with(Gfit, 2*N-2*ll)
  
  #Fine temporal scale for graphs
  yhat_all_fine<-as.data.frame(ode(y=c(B=B_i, G=bla$Sinit[1], CO2=0),
                                   func = Monod, parms=p,
                                   times = seq(0, 103, by=1)))
  Yhat_all_fine<-melt(yhat_all_fine, id.vars=c("time"))
  
  rsq_out<-list(Yhat=Yhat, Gfit=Gfit, Yhat_fine = Yhat_all_fine)
  
  return(rsq_out)
}

#Read parameters estimated in python
blag_monodopt<-as.numeric(read.csv("parameters/blag_monodpars.csv", header = F))

Blag_monodfit<-monodgood(blag_monodopt)
as.data.frame(Blag_monodfit$Gfit)

#Figure
Blag_monodfit$Yhat_fine$variable2<-Blag_monodfit$Yhat_fine$variable
levels(Blag_monodfit$Yhat_fine$variable2)<-c("B", "Glucose",  "CO[2]", "MBC", "Respiration")

Blag_fita<-subset(Blag_fit$Yhat_fine, variable=="CO2" | variable=="CFC14")
Blag_fita$Model<-c("Sub-microbial")
Blag_fitb<-subset(Blag_monodfit$Yhat_fine, variable=="CO2" | variable=="CFC14")
Blag_fitb$Model<-c("Monod")

Blag_fit_all<-rbind(Blag_fita, Blag_fitb)

ggplot(subset(Blag_fit$Yhat, variable=="CO2" | variable=="CFC14"), aes(time, obs))+
  geom_point(cex=6, pch=21, fill="grey")+
  geom_line(data=Blag_fit_all, aes(time, value, color=Model), lwd=1.2)+theme_min+
  facet_wrap(~variable2, scales="free", labeller = label_parsed) + 
  ylab(expression(paste("Carbon pool (", mu, "mol ", g(DW)^{-1}, ")"))) +
  xlab("Time (days)") + scale_color_manual(values = c("indianred", "grey30")) +
  theme(legend.position = c(0.85, 0.8))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Statistics~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Log Likelihood ratio test
-2*(Blag_monodfit$Gfit$ll-Blag_fit$Gfit$ll)

round(pchisq(-2*(Blag_monodfit$Gfit$ll-Blag_fit$Gfit$ll), df=(length(blag_optpar)-length(blag_monodopt)),
             lower.tail = F), 3)

#F test - based on residual sum of squares, number of parameters and number of measurements
##Monod model
###residual sum of squares
M1ss = sum((Blag_monodfit$Yhat$obs-Blag_monodfit$Yhat$value)^2, na.rm = T)
###number of parameters 
M1p = length(blag_monodopt)

##Sub-microbial model
###residual sum of squares
M2ss = sum((Blag_fit$Yhat$obs-Blag_fit$Yhat$value)^2, na.rm = T)
###number of parameters 
M2p = length(blag_optpar)

###total number of measurements
nt = nrow(Blag_fit$Yhat[!is.na(Blag_fit$Yhat), ])

####F value =  (M1ss - M2ss)*(nt - M2p)/M2ss/(M2p - M1p)
(M1ss - M2ss)*(nt - M2p)/M2ss/(M2p - M1p)

####associated p value
pf(q=(M1ss - M2ss)*(nt - M2p)/M2ss/(M2p - M1p), 
   df1=(M2p - M1p), 
   df2=(nt - M2p), 
   lower.tail=F)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
############################################################################################################################
#################################################Blagodatskaya et al, 2011########################################################
############################################################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#DATA
bla<-read.csv("Data/Blagodatskaya2011.csv", sep=',')

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
    CFC14=(fr*Br+fs*Bs)-bla$Cmicinit[1]
    
    #States
    dG<- - uptake
    dBr<- uptake - transfer
    dBs<- transfer*Y - death
    dCO2<- transfer*(1-Y)
    
    return(list(c(dG, dBr, dBs, dCO2), CFC14=CFC14, r=r))
  })
}

#Goodness of fit
good_sub<-function(x){
  p<-x
  names(p)<-c("v", "k", "f", "m", "Y", "fr", "fs")
  #Initial Br and Bs
  Bs_i<-bla$Cmicinit[1]/p[["fs"]]
  #Br_i<-(mar$Cmicinit[1]-p[["fs"]]*Bs_i)/p[["fr"]]
  #Simulations
  yhat_all<-as.data.frame(ode(y=c(G=bla$Sinit[1], Br=0, Bs=Bs_i, CO2=0),
                              func = SUBmodel, parms=p,
                              times = as.numeric(bla$Time)))
  #Selecting measured variables
  yhat<-yhat_all[, c("time", "CO2", "CFC14")]
  #Long format
  Yhat<-melt(yhat, id.vars=c("time"))
  #Observations
  Yhat$obs<-c(as.numeric(bla$CO214), as.numeric(bla$Cmic14))
  Gfit<-Yhat %>% group_by(variable) %>% summarise(SSres=sum(((obs-value)^2), na.rm = T),
                                                  SStot=sum(((obs-mean(obs, na.rm = T))^2), na.rm = T),
                                                  ll=-sum(((obs-value)^2), na.rm = T)/2/(sd(obs, na.rm = T)^2))
  Gfit$R2<-with(Gfit, 1-SSres/SStot)
  Gfit$N<-length(p)
  Gfit$AIC<-with(Gfit, 2*N-2*ll)
  
  #Fine temporal scale fo r graphs
  yhat_all_fine<-as.data.frame(ode(y=c(G=bla$Sinit[1], Br=0, Bs=Bs_i, CO2=0),
                                   func = SUBmodel, parms=p,
                                   times = seq(0, 55, by=1)))
  Yhat_all_fine<-melt(yhat_all_fine, id.vars=c("time"))
  
  rsq_out<-list(Yhat=Yhat, Gfit=Gfit, Yhat_fine = Yhat_all_fine)
  
  return(rsq_out)
}

#Read parameters estimated in python
blag11_optpar<-as.numeric(read.csv("parameters/blag11_optpars.csv", header = F))

Blag11_fit<-good_sub(blag11_optpar)
as.data.frame(Blag11_fit$Gfit)

#Figure
Blag11_fit$Yhat$variable2<-Blag11_fit$Yhat$variable
levels(Blag11_fit$Yhat$variable2)<-c("CO[2]", "MBC")

Blag11_fit$Yhat_fine$variable2<-Blag11_fit$Yhat_fine$variable
levels(Blag11_fit$Yhat_fine$variable2)<-c("Glucose", "Br", "Bs", "CO[2]", "MBC", "Respiration")

ggplot(subset(Blag11_fit$Yhat, variable=="CO2" | variable=="CFC14"), aes(time, obs))+
  geom_point(cex=6, pch=21, fill="grey")+
  geom_line(data=subset(Blag11_fit$Yhat_fine, variable=="CO2" | variable=="CFC14"), aes(time, value), lwd=1.2, color="grey30")+theme_min+
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
    CFC14=kec*B-bla$Cmicinit[1]
    
    #States
    dB<- uptake*Y - decay
    dG<- - uptake
    dCO2<- uptake*(1-Y)
    
    return(list(c(dB, dG, dCO2), CFC14=CFC14, r=r))
  })
}

#Goodness of fit
monodgood<-function(x){
  p<-x
  names(p)<-c("v", "k", "m", "Y", "kec")
  #Initial Br and Bs
  B_i<-bla$Cmicinit[1]/p[["kec"]]
  
  #Simulations
  yhat_all<-as.data.frame(ode(y=c(B=B_i, G=bla$Sinit[1], CO2=0),
                              func = Monod, parms=p,
                              times = as.numeric(bla$Time)))
  #Selecting measured variables
  yhat<-yhat_all[, c("time", "CO2", "CFC14")]
  #Long format
  Yhat<-melt(yhat, id.vars=c("time"))
  #Observations
  Yhat$obs<-c(as.numeric(bla$CO214), 
              as.numeric(bla$Cmic14))
  Gfit<-Yhat %>% group_by(variable) %>% summarise(SSres=sum(((obs-value)^2), na.rm = T),
                                                  SStot=sum(((obs-mean(obs, na.rm = T))^2), na.rm = T),
                                                  ll=-sum(((obs-value)^2), na.rm = T)/2/(sd(obs, na.rm = T)^2))
  Gfit$R2<-with(Gfit, 1-SSres/SStot)
  Gfit$N<-length(p)
  Gfit$AIC<-with(Gfit, 2*N-2*ll)
  
  #Fine temporal scale for graphs
  yhat_all_fine<-as.data.frame(ode(y=c(B=B_i, G=bla$Sinit[1], CO2=0),
                                   func = Monod, parms=p,
                                   times = seq(0, 55, by=1)))
  Yhat_all_fine<-melt(yhat_all_fine, id.vars=c("time"))
  
  rsq_out<-list(Yhat=Yhat, Gfit=Gfit, Yhat_fine = Yhat_all_fine)
  
  return(rsq_out)
}

#Read parameters estimated in python
blag11_monodopt<-as.numeric(read.csv("parameters/blag11_monodpars.csv", header = F))

Blag11_monodfit<-monodgood(blag11_monodopt)
as.data.frame(Blag11_monodfit$Gfit)

#Figure
Blag11_monodfit$Yhat_fine$variable2<-Blag11_monodfit$Yhat_fine$variable
levels(Blag11_monodfit$Yhat_fine$variable2)<-c("B", "Glucose",  "CO[2]", "MBC", "Respiration")

Blag11_fita<-subset(Blag11_fit$Yhat_fine, variable=="CO2" | variable=="CFC14")
Blag11_fita$Model<-c("Sub-microbial")
Blag11_fitb<-subset(Blag11_monodfit$Yhat_fine, variable=="CO2" | variable=="CFC14")
Blag11_fitb$Model<-c("Monod")

Blag11_fit_all<-rbind(Blag11_fita, Blag11_fitb)

ggplot(subset(Blag11_fit$Yhat, variable=="CO2" | variable=="CFC14"), aes(time, obs))+
  geom_point(cex=6, pch=21, fill="grey")+
  geom_line(data=Blag11_fit_all, aes(time, value, color=Model), lwd=1.2)+theme_min+
  facet_wrap(~variable2, scales="free", labeller = label_parsed) + 
  ylab(expression(paste("Carbon pool (", mu, "mol ", g(DW)^{-1}, ")"))) +
  xlab("Time (days)") + scale_color_manual(values = c("indianred", "grey30")) +
  theme(legend.position = c(0.85, 0.8))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Statistics~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Log Likelihood ratio test
-2*(Blag11_monodfit$Gfit$ll-Blag11_fit$Gfit$ll)

round(pchisq(-2*(Blag11_monodfit$Gfit$ll-Blag11_fit$Gfit$ll), df=(length(blag11_optpar)-length(blag11_monodopt)),
             lower.tail = F), 3)

#F test - based on residual sum of squares, number of parameters and number of measurements
##Monod model
###residual sum of squares
M1ss = sum((Blag11_monodfit$Yhat$obs-Blag11_monodfit$Yhat$value)^2, na.rm = T)
###number of parameters 
M1p = length(blag11_monodopt)

##Sub-microbial model
###residual sum of squares
M2ss = sum((Blag11_fit$Yhat$obs-Blag11_fit$Yhat$value)^2, na.rm = T)
###number of parameters 
M2p = length(blag11_optpar)

###total number of measurements
nt = nrow(Blag11_fit$Yhat[!is.na(Blag11_fit$Yhat), ])

####F value =  (M1ss - M2ss)*(nt - M2p)/M2ss/(M2p - M1p)
(M1ss - M2ss)*(nt - M2p)/M2ss/(M2p - M1p)

####associated p value
pf(q=(M1ss - M2ss)*(nt - M2p)/M2ss/(M2p - M1p), 
   df1=(M2p - M1p), 
   df2=(nt - M2p), 
   lower.tail=F)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
############################################################################################################################
#################################################Tsai et al, 1997###########################################################
############################################################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#DATA
tsai<-read.csv("Data/Tsai1997.csv", sep=',')
#High glucose treatment
tsai1<-subset(tsai, Treatment=="HighG")
tsai2<-subset(tsai, Treatment=="LowG")
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
    CFC=fr*Br+fs*Bs
    ATP=ar*Br+as*Bs
    
    #States
    dG<- - uptake
    dBr<- uptake - transfer
    dBs<- transfer*Y - death
    dCO2<- transfer*(1-Y)
    
    return(list(c(dG, dBr, dBs, dCO2), CFC=CFC, ATP=ATP, r=r))
  })
}

#Goodness of fit
good_sub<-function(x){
  p<-x
  names(p)<-c("v", "k", "f", "m", "Y", "fr", "fs", "ar", "as")
  #Initial Br and Bs
  Bs_i<-tsai1$Cmicinit[1]/p[["fs"]]
  
  #Simulations
  yhat_all1<-as.data.frame(ode(y=c(G=tsai1$Sinit[1], Br=0, Bs=Bs_i, CO2=0),
                              func = SUBmodel, parms=p,
                              times = as.numeric(tsai1$Time)))
  yhat_all2<-as.data.frame(ode(y=c(G=tsai2$Sinit[1], Br=0, Bs=Bs_i, CO2=0),
                               func = SUBmodel, parms=p,
                               times = as.numeric(tsai1$Time)))
  #Selecting measured variables
  yhat1<-yhat_all1[, c("time", "CO2", "CFC", "ATP")]
  #Long format
  Yhat1<-melt(yhat1, id.vars=c("time"))
  Yhat1$Treatment<-c("HighG")
  #Selecting measured variables
  yhat2<-yhat_all2[, c("time", "CO2", "CFC", "ATP")]
  #Long format
  Yhat2<-melt(yhat2, id.vars=c("time"))
  Yhat2$Treatment<-c("LowG")
  Yhat<-rbind(Yhat1, Yhat2)
  #Observations
  Yhat$obs<-c(as.numeric(tsai1$CO2cumul), as.numeric(tsai1$Cmic), as.numeric(tsai1$ATP),
              as.numeric(tsai2$CO2cumul), as.numeric(tsai2$Cmic), as.numeric(tsai2$ATP))
  Gfit<-Yhat %>% group_by(variable) %>% summarise(SSres=sum(((obs-value)^2), na.rm = T),
                                                  SStot=sum(((obs-mean(obs, na.rm = T))^2), na.rm = T),
                                                  ll=-sum(((obs-value)^2), na.rm = T)/2/(sd(obs, na.rm = T)^2))
  Gfit$R2<-with(Gfit, 1-SSres/SStot)
  Gfit$N<-length(p)
  Gfit$AIC<-with(Gfit, 2*N-2*ll)
  
  #Fine temporal scale fo r graphs
  yhat_all_fine1<-as.data.frame(ode(y=c(G=tsai1$Sinit[1], Br=0, Bs=Bs_i, CO2=0),
                                   func = SUBmodel, parms=p,
                                   times = seq(0, 9, by=0.1)))
  yhat_all_fine2<-as.data.frame(ode(y=c(G=tsai2$Sinit[1], Br=0, Bs=Bs_i, CO2=0),
                                    func = SUBmodel, parms=p,
                                    times = seq(0, 9, by=0.1)))
  Yhat_all_fine1<-melt(yhat_all_fine1, id.vars=c("time"))
  Yhat_all_fine1$Treatment<-c("HighG")
  Yhat_all_fine2<-melt(yhat_all_fine2, id.vars=c("time"))
  Yhat_all_fine2$Treatment<-c("LowG")
  
  Yhat_all_fine<-rbind(Yhat_all_fine1, Yhat_all_fine2)
  
  rsq_out<-list(Yhat=Yhat, Gfit=Gfit, Yhat_fine = Yhat_all_fine)
  
  return(rsq_out)
}

#Read parameters estimated in python
tsai_optpar<-as.numeric(read.csv("parameters/tsai_optpars.csv", header = F))

Tsai_fit<-good_sub(tsai_optpar)
as.data.frame(Tsai_fit$Gfit)

#Figure
Tsai_fit$Yhat$variable2<-Tsai_fit$Yhat$variable
levels(Tsai_fit$Yhat$variable2)<-c("CO[2]", "MBC", "ATP")

Tsai_fit$Yhat_fine$variable2<-Tsai_fit$Yhat_fine$variable
levels(Tsai_fit$Yhat_fine$variable2)<-c("Glucose", "Br", "Bs", "CO[2]", "MBC", "ATP", "Respiration")

ggplot(subset(Tsai_fit$Yhat, variable=="CO2" | variable=="CFC" | variable=="ATP"), aes(time, obs))+
  geom_point(cex=6, aes(pch=Treatment), fill="grey")+
  scale_shape_manual(values = c(21, 22)) + 
  geom_line(data=subset(Tsai_fit$Yhat_fine, variable=="CO2" | variable=="CFC" | variable=="ATP"), 
            aes(time, value, linetype=Treatment), lwd=1.2, color="grey30")+theme_min+
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
    ATP=ka*B
    
    #States
    dB<- uptake*Y - decay
    dG<- - uptake
    dCO2<- uptake*(1-Y)
    
    return(list(c(dB, dG, dCO2), CFC=CFC, ATP=ATP, r=r))
  })
}

#Goodness of fit
monodgood<-function(x){
  p<-x
  names(p)<-c("v", "k", "m", "Y", "kec", "ka")
  #Initial Br and Bs
  B_i<-tsai$Cmicinit[1]/p[["kec"]]
  
  #Simulations
  yhat_all1<-as.data.frame(ode(y=c(B=B_i, G=tsai1$Sinit[1], CO2=0),
                              func = Monod, parms=p,
                              times = as.numeric(tsai1$Time)))
  yhat_all2<-as.data.frame(ode(y=c(B=B_i, G=tsai2$Sinit[1], CO2=0),
                               func = Monod, parms=p,
                               times = as.numeric(tsai1$Time)))
  #Selecting measured variables
  yhat1<-yhat_all1[, c("time", "CO2", "CFC", "ATP")]
  #Long format
  Yhat1<-melt(yhat1, id.vars=c("time"))
  Yhat1$Treatment<-c("HighG")
  #Selecting measured variables
  yhat2<-yhat_all2[, c("time", "CO2", "CFC", "ATP")]
  #Long format
  Yhat2<-melt(yhat2, id.vars=c("time"))
  Yhat2$Treatment<-c("LowG")
  Yhat<-rbind(Yhat1, Yhat2)
  #Observations
  Yhat$obs<-c(as.numeric(tsai1$CO2cumul), as.numeric(tsai1$Cmic), as.numeric(tsai1$ATP),
              as.numeric(tsai2$CO2cumul), as.numeric(tsai2$Cmic), as.numeric(tsai2$ATP))
  Gfit<-Yhat %>% group_by(variable) %>% summarise(SSres=sum(((obs-value)^2), na.rm = T),
                                                  SStot=sum(((obs-mean(obs, na.rm = T))^2), na.rm = T),
                                                  ll=-sum(((obs-value)^2), na.rm = T)/2/(sd(obs, na.rm = T)^2))
  Gfit$R2<-with(Gfit, 1-SSres/SStot)
  Gfit$N<-length(p)
  Gfit$AIC<-with(Gfit, 2*N-2*ll)
  
  #Fine temporal scale for graphs
  yhat_all_fine1<-as.data.frame(ode(y=c(B=B_i, G=tsai1$Sinit[1], CO2=0),
                                   func = Monod, parms=p,
                                   times = seq(0, 9, by=0.1)))
  yhat_all_fine2<-as.data.frame(ode(y=c(B=B_i, G=tsai2$Sinit[1], CO2=0),
                                    func = Monod, parms=p,
                                    times = seq(0, 9, by=0.1)))
  Yhat_all_fine1<-melt(yhat_all_fine1, id.vars=c("time"))
  Yhat_all_fine1$Treatment<-c("HighG")
  Yhat_all_fine2<-melt(yhat_all_fine2, id.vars=c("time"))
  Yhat_all_fine2$Treatment<-c("LowG")
  
  Yhat_all_fine<-rbind(Yhat_all_fine1, Yhat_all_fine2)
  
  rsq_out<-list(Yhat=Yhat, Gfit=Gfit, Yhat_fine = Yhat_all_fine)
  
  return(rsq_out)
}

#Read parameters estimated in python
tsai_monodopt<-as.numeric(read.csv("parameters/tsai_monodpars.csv", header = F))

Tsai_monodfit<-monodgood(tsai_monodopt)
as.data.frame(Tsai_monodfit$Gfit)

#Figure
Tsai_monodfit$Yhat_fine$variable2<-Tsai_monodfit$Yhat_fine$variable
levels(Tsai_monodfit$Yhat_fine$variable2)<-c("B", "Glucose",  "CO[2]", "MBC", "ATP","Respiration")

Tsai_fita<-subset(Tsai_fit$Yhat_fine, variable=="CO2" | variable=="CFC" | variable=="ATP")
Tsai_fita$Model<-c("Sub-microbial")
Tsai_fitb<-subset(Tsai_monodfit$Yhat_fine, variable=="CO2" | variable=="CFC" | variable=="ATP")
Tsai_fitb$Model<-c("Monod")

Tsai_fit_all<-rbind(Tsai_fita, Tsai_fitb)

ggplot(subset(Tsai_fit$Yhat, variable=="CO2" | variable=="CFC" | variable=="ATP"), aes(time, obs))+
  geom_point(cex=6, aes(pch=Treatment), fill="grey")+
  scale_shape_manual(values=c(21, 22))+
  geom_line(data=Tsai_fit_all, aes(time, value, color=Model, linetype=Treatment), lwd=1.2)+theme_min+
  facet_wrap(~variable2, scales="free", labeller = label_parsed) + 
  ylab(expression(paste("Carbon pool (", mu, "mol ", g(DW)^{-1}, ")"))) +
  xlab("Time (days)") + scale_color_manual(values = c("indianred", "grey30")) +
  theme(legend.position = c(0.85, 0.8))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Statistics~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Log Likelihood ratio test
-2*(Tsai_monodfit$Gfit$ll-Tsai_fit$Gfit$ll)

round(pchisq(-2*(Tsai_monodfit$Gfit$ll-Tsai_fit$Gfit$ll), df=(length(tsai_optpar)-length(tsai_monodopt)),
             lower.tail = F), 3)

#F test - based on residual sum of squares, number of parameters and number of measurements
##Monod model
###residual sum of squares
M1ss = sum((Tsai_monodfit$Yhat$obs-Tsai_monodfit$Yhat$value)^2, na.rm = T)
###number of parameters 
M1p = length(tsai_monodopt)

##Sub-microbial model
###residual sum of squares
M2ss = sum((Tsai_fit$Yhat$obs-Tsai_fit$Yhat$value)^2, na.rm = T)
###number of parameters 
M2p = length(tsai_optpar)

###total number of measurements
nt = nrow(Tsai_fit$Yhat[!is.na(Tsai_fit$Yhat), ])

####F value =  (M1ss - M2ss)*(nt - M2p)/M2ss/(M2p - M1p)
(M1ss - M2ss)*(nt - M2p)/M2ss/(M2p - M1p)

####associated p value
pf(q=(M1ss - M2ss)*(nt - M2p)/M2ss/(M2p - M1p), 
   df1=(M2p - M1p), 
   df2=(nt - M2p), 
   lower.tail=F)
  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
############################################################################################################################
###################################################Wu et al, 1993###########################################################
############################################################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#DATA
wu<-read.csv("Data/Wu1993.csv", sep=',')
#Glucose
wuG<-subset(wu, Substrate=='Glucose')
wuG1<-subset(wuG, Treatment=="Large")
wuG2<-subset(wuG, Treatment=="Small")
#Ryegrass
wuR<-subset(wu, Substrate=='Ryegrass')
wuR1<-subset(wuR, Treatment=="Large")
wuR2<-subset(wuR, Treatment=="Small")
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
    CFC14=(fr*Br+fs*Bs)-wu$Cmicinit[1]
    
    #States
    dG<- - uptake
    dBr<- uptake - transfer
    dBs<- transfer*Y - death
    dCO2<- transfer*(1-Y)
    
    return(list(c(dG, dBr, dBs, dCO2), CFC14=CFC14))
  })
}

#Goodness of fit
good_sub<-function(xG, xR){
  
  names(xG)<-c("v", "k", "f", "m", "Y", "fr", "fs")
  names(xR)<-c("v", "k", "f", "m", "Y", "fr", "fs")
  #Initial Br and Bs
  Bs_iG<-wu$Cmicinit[1]/xG[["fs"]]
  Bs_iR<-wu$Cmicinit[1]/xR[["fs"]]
  
  #Simulations
  yhat_allG1<-as.data.frame(ode(y=c(G=wuG1$Sinit[1], Br=0, Bs=Bs_iG, CO2=0),
                               func = SUBmodel, parms=xG, 
                               times = as.numeric(wuG1$Time)))
  yhat_allG2<-as.data.frame(ode(y=c(G=wuG2$Sinit[1], Br=0, Bs=Bs_iG, CO2=0),
                               func = SUBmodel, parms=xG,
                               times = as.numeric(wuG2$Time)))
  yhat_allR1<-as.data.frame(ode(y=c(G=wuR1$Sinit[1], Br=0, Bs=Bs_iR, CO2=0),
                                func = SUBmodel, parms=xR,
                                times = as.numeric(wuR1$Time)))
  yhat_allR2<-as.data.frame(ode(y=c(G=wuR2$Sinit[1], Br=0, Bs=Bs_iR, CO2=0),
                                func = SUBmodel, parms=xR,
                                times = as.numeric(wuR2$Time)))
  #Selecting measured variables
  yhatG1<-yhat_allG1[, c("time", "CO2", "CFC14")]
  yhatG2<-yhat_allG2[, c("time", "CO2", "CFC14")]
  yhatR1<-yhat_allR1[, c("time", "CO2", "CFC14")]
  yhatR2<-yhat_allR2[, c("time", "CO2", "CFC14")]
  #Long format
  YhatG1<-melt(yhatG1, id.vars=c("time"))
  YhatG1$Substrate<-c("Glucose")
  YhatG1$Treatment<-c("Large")
  YhatG2<-melt(yhatG2, id.vars=c("time"))
  YhatG2$Substrate<-c("Glucose")
  YhatG2$Treatment<-c("Small")
  YhatR1<-melt(yhatR1, id.vars=c("time"))
  YhatR1$Substrate<-c("Ryegrass")
  YhatR1$Treatment<-c("Large")
  YhatR2<-melt(yhatR2, id.vars=c("time"))
  YhatR2$Substrate<-c("Ryegrass")
  YhatR2$Treatment<-c("Small")
  Yhat<-rbind(YhatG1, YhatG2, YhatR1, YhatR2)
  #Observations
  Yhat$obs<-c(as.numeric(wuG1$CO2cumul), as.numeric(wuG1$Cmic14),
              as.numeric(wuG2$CO2cumul), as.numeric(wuG2$Cmic14),
              as.numeric(wuR1$CO2cumul), as.numeric(wuR1$Cmic14),
              as.numeric(wuR2$CO2cumul), as.numeric(wuR2$Cmic14))
  Gfit<-Yhat %>% group_by(variable, Substrate) %>% summarise(SSres=sum(((obs-value)^2), na.rm = T),
                                                  SStot=sum(((obs-mean(obs, na.rm = T))^2), na.rm = T),
                                                  ll=-sum(((obs-value)^2), na.rm = T)/2/(sd(obs, na.rm = T)^2))
  Gfit$R2<-with(Gfit, 1-SSres/SStot)
  Gfit$N<-length(xG)
  Gfit$AIC<-with(Gfit, 2*N-2*ll)
  
  #Fine temporal scale fo r graphs
  yhat_all_fineG1<-as.data.frame(ode(y=c(G=wuG1$Sinit[1], Br=0, Bs=Bs_iG, CO2=0),
                                    func = SUBmodel, parms=xG,
                                    times = seq(0, 105)))
  yhat_all_fineG2<-as.data.frame(ode(y=c(G=wuG2$Sinit[1], Br=0, Bs=Bs_iG, CO2=0),
                                    func = SUBmodel, parms=xG,
                                    times = seq(0, 105)))
  yhat_all_fineR1<-as.data.frame(ode(y=c(G=wuR1$Sinit[1], Br=0, Bs=Bs_iR, CO2=0),
                                    func = SUBmodel, parms=xR,
                                    times = seq(0, 150)))
  yhat_all_fineR2<-as.data.frame(ode(y=c(G=wuR2$Sinit[1], Br=0, Bs=Bs_iR, CO2=0),
                                    func = SUBmodel, parms=xR,
                                    times = seq(0, 150)))
  Yhat_all_fineG1<-melt(yhat_all_fineG1, id.vars=c("time"))
  Yhat_all_fineG1$Substrate<-c("Glucose")
  Yhat_all_fineG1$Treatment<-c("Large")
  Yhat_all_fineG2<-melt(yhat_all_fineG2, id.vars=c("time"))
  Yhat_all_fineG2$Substrate<-c("Glucose")
  Yhat_all_fineG2$Treatment<-c("Small")
  Yhat_all_fineR1<-melt(yhat_all_fineR1, id.vars=c("time"))
  Yhat_all_fineR1$Substrate<-c("Ryegrass")
  Yhat_all_fineR1$Treatment<-c("Large")
  Yhat_all_fineR2<-melt(yhat_all_fineR2, id.vars=c("time"))
  Yhat_all_fineR2$Substrate<-c("Ryegrass")
  Yhat_all_fineR2$Treatment<-c("Small")
  
  Yhat_all_fine<-rbind(Yhat_all_fineG1, Yhat_all_fineG2,
                       Yhat_all_fineR1, Yhat_all_fineR2)
  
  rsq_out<-list(Yhat=Yhat, Gfit=Gfit, Yhat_fine = Yhat_all_fine)
  
  return(rsq_out)
}

#Read parameters estimated in python
wu_optparG<-as.numeric(read.csv("parameters/wu_optparsG.csv", header = F))
wu_optparR<-as.numeric(read.csv("parameters/wu_optparsR.csv", header = F))

Wu_fit<-good_sub(wu_optparG, wu_optparR)
as.data.frame(Wu_fit$Gfit)

#Figure
Wu_fit$Yhat$variable2<-Wu_fit$Yhat$variable
levels(Wu_fit$Yhat$variable2)<-c("CO[2]", "MBC")

Wu_fit$Yhat_fine$variable2<-Wu_fit$Yhat_fine$variable
levels(Wu_fit$Yhat_fine$variable2)<-c("Glucose", "Br", "Bs", "CO[2]", "MBC")

ggplot(subset(Wu_fit$Yhat, variable=="CO2" | variable=="CFC14"), aes(time, obs))+
  geom_point(cex=6, aes(pch=Treatment), fill="grey")+
  scale_shape_manual(values = c(21, 22)) + 
  geom_line(data=subset(Wu_fit$Yhat_fine, variable=="CO2" | variable=="CFC14"), 
            aes(time, value, linetype=Treatment), lwd=1.2, color="grey30")+theme_min+
  facet_wrap(Substrate~variable2, scales="free", labeller = label_parsed) + 
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
    CFC14=(kec*B)-wu$Cmicinit[1]
    
    #States
    dB<- uptake*Y - decay
    dG<- - uptake
    dCO2<- uptake*(1-Y)
    
    return(list(c(dB, dG, dCO2), CFC14=CFC14))
  })
}

#Goodness of fit
monodgood<-function(xG, xR){
  names(xG)<-c("v", "k", "m", "Y", "kec")
  names(xR)<-c("v", "k", "m", "Y", "kec")
  
  #Initial Br and Bs
  B_iG<-wu$Cmicinit[1]/xG[["kec"]]
  B_iR<-wu$Cmicinit[1]/xR[["kec"]]
  
  #Simulations
  yhat_allG1<-as.data.frame(ode(y=c(B=B_iG, G=wuG1$Sinit[1], CO2=0),
                               func = Monod, parms=xG,
                               times = as.numeric(wuG1$Time)))
  yhat_allG2<-as.data.frame(ode(y=c(B=B_iG, G=wuG2$Sinit[1], CO2=0),
                               func = Monod, parms=xG,
                               times = as.numeric(wuG2$Time)))
  yhat_allR1<-as.data.frame(ode(y=c(B=B_iR, G=wuR1$Sinit[1], CO2=0),
                                func = Monod, parms=xR,
                                times = as.numeric(wuR1$Time)))
  yhat_allR2<-as.data.frame(ode(y=c(B=B_iR, G=wuR2$Sinit[1], CO2=0),
                                func = Monod, parms=xR,
                                times = as.numeric(wuR2$Time)))
  #Selecting measured variables
  yhatG1<-yhat_allG1[, c("time", "CO2", "CFC14")]
  yhatG2<-yhat_allG2[, c("time", "CO2", "CFC14")]
  yhatR1<-yhat_allR1[, c("time", "CO2", "CFC14")]
  yhatR2<-yhat_allR2[, c("time", "CO2", "CFC14")]
  #Long format
  YhatG1<-melt(yhatG1, id.vars=c("time"))
  YhatG1$Substrate<-c("Glucose")
  YhatG1$Treatment<-c("Large")
  YhatG2<-melt(yhatG2, id.vars=c("time"))
  YhatG2$Substrate<-c("Glucose")
  YhatG2$Treatment<-c("Small")
  YhatR1<-melt(yhatR1, id.vars=c("time"))
  YhatR1$Substrate<-c("Ryegrass")
  YhatR1$Treatment<-c("Large")
  YhatR2<-melt(yhatR2, id.vars=c("time"))
  YhatR2$Substrate<-c("Ryegrass")
  YhatR2$Treatment<-c("Small")
  Yhat<-rbind(YhatG1, YhatG2, YhatR1, YhatR2)
  #Observations
  Yhat$obs<-c(as.numeric(wuG1$CO2cumul), as.numeric(wuG1$Cmic14),
              as.numeric(wuG2$CO2cumul), as.numeric(wuG2$Cmic14),
              as.numeric(wuR1$CO2cumul), as.numeric(wuR1$Cmic14),
              as.numeric(wuR2$CO2cumul), as.numeric(wuR2$Cmic14))
  Gfit<-Yhat %>% group_by(variable, Substrate) %>% summarise(SSres=sum(((obs-value)^2), na.rm = T),
                                                                        SStot=sum(((obs-mean(obs, na.rm = T))^2), na.rm = T),
                                                                        ll=-sum(((obs-value)^2), na.rm = T)/2/(sd(obs, na.rm = T)^2))
  Gfit$R2<-with(Gfit, 1-SSres/SStot)
  Gfit$N<-length(xG)
  Gfit$AIC<-with(Gfit, 2*N-2*ll)
  
  #Fine temporal scale for graphs
  yhat_all_fineG1<-as.data.frame(ode(y=c(B=B_iG, G=wuG1$Sinit[1], CO2=0),
                                    func = Monod, parms=xG,
                                    times = seq(0, 105)))
  yhat_all_fineG2<-as.data.frame(ode(y=c(B=B_iG, G=wuG2$Sinit[1], CO2=0),
                                    func = Monod, parms=xG,
                                    times = seq(0, 105)))
  yhat_all_fineR1<-as.data.frame(ode(y=c(B=B_iR, G=wuR1$Sinit[1], CO2=0),
                                     func = Monod, parms=xR,
                                     times = seq(0, 150)))
  yhat_all_fineR2<-as.data.frame(ode(y=c(B=B_iR, G=wuR2$Sinit[1], CO2=0),
                                     func = Monod, parms=xR,
                                     times = seq(0, 150)))
  Yhat_all_fineG1<-melt(yhat_all_fineG1, id.vars=c("time"))
  Yhat_all_fineG1$Substrate<-c("Glucose")
  Yhat_all_fineG1$Treatment<-c("Large")
  Yhat_all_fineG2<-melt(yhat_all_fineG2, id.vars=c("time"))
  Yhat_all_fineG2$Substrate<-c("Glucose")
  Yhat_all_fineG2$Treatment<-c("Small")
  Yhat_all_fineR1<-melt(yhat_all_fineR1, id.vars=c("time"))
  Yhat_all_fineR1$Substrate<-c("Ryegrass")
  Yhat_all_fineR1$Treatment<-c("Large")
  Yhat_all_fineR2<-melt(yhat_all_fineR2, id.vars=c("time"))
  Yhat_all_fineR2$Substrate<-c("Ryegrass")
  Yhat_all_fineR2$Treatment<-c("Small")
  
  Yhat_all_fine<-rbind(Yhat_all_fineG1, Yhat_all_fineG2,
                       Yhat_all_fineR1, Yhat_all_fineR2)
  
  rsq_out<-list(Yhat=Yhat, Gfit=Gfit, Yhat_fine = Yhat_all_fine)
  
  return(rsq_out)
}

#Read parameters estimated in python
wu_monodoptG<-as.numeric(read.csv("parameters/wu_monodparsG.csv", header = F))
wu_monodoptR<-as.numeric(read.csv("parameters/wu_monodparsR.csv", header = F))

Wu_monodfit<-monodgood(wu_monodoptG, wu_monodoptR)
as.data.frame(Wu_monodfit$Gfit)

#Figure
Wu_monodfit$Yhat_fine$variable2<-Wu_monodfit$Yhat_fine$variable
levels(Wu_monodfit$Yhat_fine$variable2)<-c("B", "Glucose",  "CO[2]", "MBC")

Wu_fita<-subset(Wu_fit$Yhat_fine, variable=="CO2" | variable=="CFC14")
Wu_fita$Model<-c("Sub-microbial")
Wu_fitb<-subset(Wu_monodfit$Yhat_fine, variable=="CO2" | variable=="CFC14")
Wu_fitb$Model<-c("Monod")

Wu_fit_all<-rbind(Wu_fita, Wu_fitb)

ggplot(subset(Wu_fit$Yhat, variable=="CO2" | variable=="CFC14"), aes(time, obs))+
  geom_point(cex=6, aes(pch=Treatment), fill="grey")+
  scale_shape_manual(values=c(21, 22))+
  geom_line(data=Wu_fit_all, aes(time, value, color=Model, linetype=Treatment), lwd=1.2)+theme_min+
  facet_wrap(Substrate~variable2, scales="free", labeller = label_parsed) + 
  ylab(expression(paste("Carbon pool (", mu, "mol ", g(DW)^{-1}, ")"))) +
  xlab("Time (days)") + scale_color_manual(values = c("indianred", "grey30")) +
  theme(legend.position = c(0.85, 0.8))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Statistics~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Log Likelihood ratio test
-2*(Wu_monodfit$Gfit$ll-Wu_fit$Gfit$ll)

round(pchisq(-2*(Wu_monodfit$Gfit$ll-Wu_fit$Gfit$ll), df=4,
             lower.tail = F), 3)

#F test - based on residual sum of squares, number of parameters and number of measurements
##Monod model
###residual sum of squares
M1ss = sum((Wu_monodfit$Yhat$obs-Wu_monodfit$Yhat$value)^2, na.rm = T)
###number of parameters 
M1p = 10

##Sub-microbial model
###residual sum of squares
M2ss = sum((Wu_fit$Yhat$obs-Wu_fit$Yhat$value)^2, na.rm = T)
###number of parameters 
M2p = 14

###total number of measurements
nt = nrow(Wu_fit$Yhat[!is.na(Wu_fit$Yhat), ])

####F value =  (M1ss - M2ss)*(nt - M2p)/M2ss/(M2p - M1p)
(M1ss - M2ss)*(nt - M2p)/M2ss/(M2p - M1p)

####associated p value
pf(q=(M1ss - M2ss)*(nt - M2p)/M2ss/(M2p - M1p), 
   df1=(M2p - M1p), 
   df2=(nt - M2p), 
   lower.tail=F)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
############################################################################################################################
###################################################Wu et al, 2011###########################################################
############################################################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#DATA
wu<-read.csv("Data/Wu2011.csv", sep=',')
#Glucose
wuG<-subset(wu, Substrate=='Glucose')
wuG1<-subset(wuG, Soil=="Paddy")
wuG2<-subset(wuG, Soil=="Upland")
#Ryegrass
wuR<-subset(wu, Substrate=='Rice')
wuR1<-subset(wuR, Soil=="Paddy")
wuR2<-subset(wuR, Soil=="Upland")
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
    #CFC14=(fr*Br+fs*Bs)-wu$Cmicinit[1]
    
    #States
    dG<- - uptake
    dBr<- uptake - transfer
    dBs<- transfer*Y - death
    dCO2<- transfer*(1-Y)
    
    return(list(c(dG, dBr, dBs, dCO2)))
  })
}

#Goodness of fit
good_sub<-function(xG1, xG2, xR1, xR2){
  
  names(xG1)<-c("v", "k", "f", "m", "Y", "fr", "fs")
  names(xG2)<-c("v", "k", "f", "m", "Y", "fr", "fs")
  names(xR1)<-c("v", "k", "f", "m", "Y", "fr", "fs")
  names(xR2)<-c("v", "k", "f", "m", "Y", "fr", "fs")
  #Initial Br and Bs
  Bs_iG1<-wuG1$Cmicinit[1]/xG1[["fs"]]
  Bs_iG2<-wuG2$Cmicinit[1]/xG2[["fs"]]
  Bs_iR1<-wuR1$Cmicinit[1]/xR1[["fs"]]
  Bs_iR2<-wuR2$Cmicinit[1]/xR2[["fs"]]
  
  #Simulations
  yhat_allG1<-as.data.frame(ode(y=c(G=wuG1$Sinit[1], Br=0, Bs=Bs_iG1, CO2=0),
                                func = SUBmodel, parms=xG1, 
                                times = as.numeric(wuG1$Time)))
  yhat_allG1$CFC14<-with(yhat_allG1, (xG1[["fr"]]*Br+xG1[["fs"]]*Bs)-wuG1$Cmicinit[1])
  yhat_allG2<-as.data.frame(ode(y=c(G=wuG2$Sinit[1], Br=0, Bs=Bs_iG2, CO2=0),
                                func = SUBmodel, parms=xG2,
                                times = as.numeric(wuG2$Time)))
  yhat_allG2$CFC14<-with(yhat_allG2, (xG2[["fr"]]*Br+xG2[["fs"]]*Bs)-wuG2$Cmicinit[1])
  yhat_allR1<-as.data.frame(ode(y=c(G=wuR1$Sinit[1], Br=0, Bs=Bs_iR1, CO2=0),
                                func = SUBmodel, parms=xR1,
                                times = as.numeric(wuR1$Time)))
  yhat_allR1$CFC14<-with(yhat_allR1, (xR1[["fr"]]*Br+xR1[["fs"]]*Bs)-wuR1$Cmicinit[1])
  yhat_allR2<-as.data.frame(ode(y=c(G=wuR2$Sinit[1], Br=0, Bs=Bs_iR2, CO2=0),
                                func = SUBmodel, parms=xR2,
                                times = as.numeric(wuR2$Time)))
  yhat_allR2$CFC14<-with(yhat_allR2, (xR2[["fr"]]*Br+xR2[["fs"]]*Bs)-wuR2$Cmicinit[1])
  #Selecting measured variables
  yhatG1<-yhat_allG1[, c("time", "CO2", "CFC14")]
  yhatG2<-yhat_allG2[, c("time", "CO2", "CFC14")]
  yhatR1<-yhat_allR1[, c("time", "CO2", "CFC14")]
  yhatR2<-yhat_allR2[, c("time", "CO2", "CFC14")]
  #Long format
  YhatG1<-melt(yhatG1, id.vars=c("time"))
  YhatG1$Substrate<-c("Glucose")
  YhatG1$Soil<-c("Paddy")
  YhatG2<-melt(yhatG2, id.vars=c("time"))
  YhatG2$Substrate<-c("Glucose")
  YhatG2$Soil<-c("Upland")
  YhatR1<-melt(yhatR1, id.vars=c("time"))
  YhatR1$Substrate<-c("Rice")
  YhatR1$Soil<-c("Paddy")
  YhatR2<-melt(yhatR2, id.vars=c("time"))
  YhatR2$Substrate<-c("Rice")
  YhatR2$Soil<-c("Upland")
  Yhat<-rbind(YhatG1, YhatG2, YhatR1, YhatR2)
  #Observations
  Yhat$obs<-c(as.numeric(wuG1$CO2cumul), as.numeric(wuG1$Cmic14),
              as.numeric(wuG2$CO2cumul), as.numeric(wuG2$Cmic14),
              as.numeric(wuR1$CO2cumul), as.numeric(wuR1$Cmic14),
              as.numeric(wuR2$CO2cumul), as.numeric(wuR2$Cmic14))
  Gfit<-Yhat %>% group_by(variable, Substrate, Soil) %>% summarise(SSres=sum(((obs-value)^2), na.rm = T),
                                                                        SStot=sum(((obs-mean(obs, na.rm = T))^2), na.rm = T),
                                                                        ll=-sum(((obs-value)^2), na.rm = T)/2/(sd(obs, na.rm = T)^2))
  Gfit$R2<-with(Gfit, 1-SSres/SStot)
  Gfit$N<-length(xG1)
  Gfit$AIC<-with(Gfit, 2*N-2*ll)
  
  #Fine temporal scale fo r graphs
  yhat_all_fineG1<-as.data.frame(ode(y=c(G=wuG1$Sinit[1], Br=0, Bs=Bs_iG1, CO2=0),
                                     func = SUBmodel, parms=xG1,
                                     times = seq(0, 100)))
  yhat_all_fineG1$CFC14<-with(yhat_all_fineG1, (xG1[["fr"]]*Br+xG1[["fs"]]*Bs)-wuG1$Cmicinit[1])
  yhat_all_fineG2<-as.data.frame(ode(y=c(G=wuG2$Sinit[1], Br=0, Bs=Bs_iG2, CO2=0),
                                     func = SUBmodel, parms=xG2,
                                     times = seq(0, 100)))
  yhat_all_fineG2$CFC14<-with(yhat_all_fineG2, (xG2[["fr"]]*Br+xG2[["fs"]]*Bs)-wuG2$Cmicinit[1])
  yhat_all_fineR1<-as.data.frame(ode(y=c(G=wuR1$Sinit[1], Br=0, Bs=Bs_iR1, CO2=0),
                                     func = SUBmodel, parms=xR1,
                                     times = seq(0, 100)))
  yhat_all_fineR1$CFC14<-with(yhat_all_fineR1, (xR1[["fr"]]*Br+xR1[["fs"]]*Bs)-wuR1$Cmicinit[1])
  yhat_all_fineR2<-as.data.frame(ode(y=c(G=wuR2$Sinit[1], Br=0, Bs=Bs_iR2, CO2=0),
                                     func = SUBmodel, parms=xR2,
                                     times = seq(0, 100)))
  yhat_all_fineR2$CFC14<-with(yhat_all_fineR2, (xR2[["fr"]]*Br+xR2[["fs"]]*Bs)-wuR2$Cmicinit[1])
  
  Yhat_all_fineG1<-melt(yhat_all_fineG1, id.vars=c("time"))
  Yhat_all_fineG1$Substrate<-c("Glucose")
  Yhat_all_fineG1$Soil<-c("Paddy")
  Yhat_all_fineG2<-melt(yhat_all_fineG2, id.vars=c("time"))
  Yhat_all_fineG2$Substrate<-c("Glucose")
  Yhat_all_fineG2$Soil<-c("Upland")
  Yhat_all_fineR1<-melt(yhat_all_fineR1, id.vars=c("time"))
  Yhat_all_fineR1$Substrate<-c("Rice")
  Yhat_all_fineR1$Soil<-c("Paddy")
  Yhat_all_fineR2<-melt(yhat_all_fineR2, id.vars=c("time"))
  Yhat_all_fineR2$Substrate<-c("Rice")
  Yhat_all_fineR2$Soil<-c("Upland")
  
  Yhat_all_fine<-rbind(Yhat_all_fineG1, Yhat_all_fineG2,
                       Yhat_all_fineR1, Yhat_all_fineR2)
  
  rsq_out<-list(Yhat=Yhat, Gfit=Gfit, Yhat_fine = Yhat_all_fine)
  
  return(rsq_out)
}

#Read parameters estimated in python
wu11_optparG1<-as.numeric(read.csv("parameters/wu11_optparsG1.csv", header = F))
wu11_optparG2<-as.numeric(read.csv("parameters/wu11_optparsG2.csv", header = F))
wu11_optparR1<-as.numeric(read.csv("parameters/wu11_optparsR1.csv", header = F))
wu11_optparR2<-as.numeric(read.csv("parameters/wu11_optparsR2.csv", header = F))

Wu11_fit<-good_sub(wu11_optparG1, wu11_optparG2, wu11_optparR1, wu11_optparR2)
as.data.frame(Wu11_fit$Gfit)

#Figure
Wu11_fit$Yhat$variable2<-Wu11_fit$Yhat$variable
levels(Wu11_fit$Yhat$variable2)<-c("CO[2]", "MBC")

Wu11_fit$Yhat_fine$variable2<-Wu11_fit$Yhat_fine$variable
levels(Wu11_fit$Yhat_fine$variable2)<-c("Glucose", "Br", "Bs", "CO[2]", "MBC")

ggplot(subset(Wu11_fit$Yhat, variable=="CO2" | variable=="CFC14"), aes(time, obs))+
  geom_point(cex=6, aes(pch=Soil), fill="grey")+
  scale_shape_manual(values = c(21, 22)) + 
  geom_line(data=subset(Wu11_fit$Yhat_fine, variable=="CO2" | variable=="CFC14"), 
            aes(time, value, linetype=Soil), lwd=1.2, color="grey30")+theme_min+
  facet_wrap(Substrate~variable2, scales="free", labeller = label_parsed) + 
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
    #CFC14=(kec*B)-wu$Cmicinit[1]
    
    #States
    dB<- uptake*Y - decay
    dG<- - uptake
    dCO2<- uptake*(1-Y)
    
    return(list(c(dB, dG, dCO2)))
  })
}

#Goodness of fit
monodgood<-function(xG1, xG2, xR1, xR2){
  names(xG1)<-c("v", "k", "m", "Y", "kec")
  names(xG2)<-c("v", "k", "m", "Y", "kec")
  names(xR1)<-c("v", "k", "m", "Y", "kec")
  names(xR2)<-c("v", "k", "m", "Y", "kec")
  
  #Initial Br and Bs
  B_iG1<-wuG1$Cmicinit[1]/xG1[["kec"]]
  B_iG2<-wuG2$Cmicinit[1]/xG2[["kec"]]
  B_iR1<-wuR1$Cmicinit[1]/xR1[["kec"]]
  B_iR2<-wuR2$Cmicinit[1]/xR2[["kec"]]
  
  #Simulations
  yhat_allG1<-as.data.frame(ode(y=c(B=B_iG1, G=wuG1$Sinit[1], CO2=0),
                                func = Monod, parms=xG1,
                                times = as.numeric(wuG1$Time)))
  yhat_allG1$CFC14<-with(yhat_allG1, (xG1[["kec"]]*B)-wuG1$Cmicinit[1])
  yhat_allG2<-as.data.frame(ode(y=c(B=B_iG2, G=wuG2$Sinit[1], CO2=0),
                                func = Monod, parms=xG2,
                                times = as.numeric(wuG2$Time)))
  yhat_allG2$CFC14<-with(yhat_allG2, (xG2[["kec"]]*B)-wuG2$Cmicinit[1])
  yhat_allR1<-as.data.frame(ode(y=c(B=B_iR1, G=wuR1$Sinit[1], CO2=0),
                                func = Monod, parms=xR1,
                                times = as.numeric(wuR1$Time)))
  yhat_allR1$CFC14<-with(yhat_allR1, (xR1[["kec"]]*B)-wuR1$Cmicinit[1])
  yhat_allR2<-as.data.frame(ode(y=c(B=B_iR2, G=wuR2$Sinit[1], CO2=0),
                                func = Monod, parms=xR2,
                                times = as.numeric(wuR2$Time)))
  yhat_allR2$CFC14<-with(yhat_allR2, (xR2[["kec"]]*B)-wuR2$Cmicinit[1])
  #Selecting measured variables
  yhatG1<-yhat_allG1[, c("time", "CO2", "CFC14")]
  yhatG2<-yhat_allG2[, c("time", "CO2", "CFC14")]
  yhatR1<-yhat_allR1[, c("time", "CO2", "CFC14")]
  yhatR2<-yhat_allR2[, c("time", "CO2", "CFC14")]
  #Long format
  YhatG1<-melt(yhatG1, id.vars=c("time"))
  YhatG1$Substrate<-c("Glucose")
  YhatG1$Soil<-c("Paddy")
  YhatG2<-melt(yhatG2, id.vars=c("time"))
  YhatG2$Substrate<-c("Glucose")
  YhatG2$Soil<-c("Upland")
  YhatR1<-melt(yhatR1, id.vars=c("time"))
  YhatR1$Substrate<-c("Rice")
  YhatR1$Soil<-c("Paddy")
  YhatR2<-melt(yhatR2, id.vars=c("time"))
  YhatR2$Substrate<-c("Rice")
  YhatR2$Soil<-c("Upland")
  Yhat<-rbind(YhatG1, YhatG2, YhatR1, YhatR2)
  #Observations
  Yhat$obs<-c(as.numeric(wuG1$CO2cumul), as.numeric(wuG1$Cmic14),
              as.numeric(wuG2$CO2cumul), as.numeric(wuG2$Cmic14),
              as.numeric(wuR1$CO2cumul), as.numeric(wuR1$Cmic14),
              as.numeric(wuR2$CO2cumul), as.numeric(wuR2$Cmic14))
  Gfit<-Yhat %>% group_by(variable, Substrate, Soil) %>% summarise(SSres=sum(((obs-value)^2), na.rm = T),
                                                                        SStot=sum(((obs-mean(obs, na.rm = T))^2), na.rm = T),
                                                                        ll=-sum(((obs-value)^2), na.rm = T)/2/(sd(obs, na.rm = T)^2))
  Gfit$R2<-with(Gfit, 1-SSres/SStot)
  Gfit$N<-length(xG1)
  Gfit$AIC<-with(Gfit, 2*N-2*ll)
  
  #Fine temporal scale for graphs
  yhat_all_fineG1<-as.data.frame(ode(y=c(B=B_iG1, G=wuG1$Sinit[1], CO2=0),
                                     func = Monod, parms=xG1,
                                     times = seq(0, 100)))
  yhat_all_fineG1$CFC14<-with(yhat_all_fineG1, (xG1[["kec"]]*B)-wuG1$Cmicinit[1])
  yhat_all_fineG2<-as.data.frame(ode(y=c(B=B_iG2, G=wuG2$Sinit[1], CO2=0),
                                     func = Monod, parms=xG2,
                                     times = seq(0, 100)))
  yhat_all_fineG2$CFC14<-with(yhat_all_fineG2, (xG2[["kec"]]*B)-wuG2$Cmicinit[1])
  yhat_all_fineR1<-as.data.frame(ode(y=c(B=B_iR1, G=wuR1$Sinit[1], CO2=0),
                                     func = Monod, parms=xR1,
                                     times = seq(0, 100)))
  yhat_all_fineR1$CFC14<-with(yhat_all_fineR1, (xR1[["kec"]]*B)-wuR1$Cmicinit[1])
  yhat_all_fineR2<-as.data.frame(ode(y=c(B=B_iR2, G=wuR2$Sinit[1], CO2=0),
                                     func = Monod, parms=xR2,
                                     times = seq(0, 100)))
  yhat_all_fineR2$CFC14<-with(yhat_all_fineR2, (xR2[["kec"]]*B)-wuR2$Cmicinit[1])
  Yhat_all_fineG1<-melt(yhat_all_fineG1, id.vars=c("time"))
  Yhat_all_fineG1$Substrate<-c("Glucose")
  Yhat_all_fineG1$Soil<-c("Paddy")
  Yhat_all_fineG2<-melt(yhat_all_fineG2, id.vars=c("time"))
  Yhat_all_fineG2$Substrate<-c("Glucose")
  Yhat_all_fineG2$Soil<-c("Upland")
  Yhat_all_fineR1<-melt(yhat_all_fineR1, id.vars=c("time"))
  Yhat_all_fineR1$Substrate<-c("Rice")
  Yhat_all_fineR1$Soil<-c("Paddy")
  Yhat_all_fineR2<-melt(yhat_all_fineR2, id.vars=c("time"))
  Yhat_all_fineR2$Substrate<-c("Rice")
  Yhat_all_fineR2$Soil<-c("Upland")
  
  Yhat_all_fine<-rbind(Yhat_all_fineG1, Yhat_all_fineG2,
                       Yhat_all_fineR1, Yhat_all_fineR2)
  
  rsq_out<-list(Yhat=Yhat, Gfit=Gfit, Yhat_fine = Yhat_all_fine)
  
  return(rsq_out)
}

#Read parameters estimated in python
wu11_monodoptG1<-as.numeric(read.csv("parameters/wu11_monodparsG1.csv", header = F))
wu11_monodoptG2<-as.numeric(read.csv("parameters/wu11_monodparsG2.csv", header = F))
wu11_monodoptR1<-as.numeric(read.csv("parameters/wu11_monodparsR1.csv", header = F))
wu11_monodoptR2<-as.numeric(read.csv("parameters/wu11_monodparsR2.csv", header = F))

Wu11_monodfit<-monodgood(wu11_monodoptG1, wu11_monodoptG2, wu11_monodoptR1, wu11_monodoptR2)
as.data.frame(Wu11_monodfit$Gfit)

#Figure
Wu11_monodfit$Yhat_fine$variable2<-Wu11_monodfit$Yhat_fine$variable
levels(Wu11_monodfit$Yhat_fine$variable2)<-c("B", "Glucose",  "CO[2]", "MBC")

Wu11_fita<-subset(Wu11_fit$Yhat_fine, variable=="CO2" | variable=="CFC14")
Wu11_fita$Model<-c("Sub-microbial")
Wu11_fitb<-subset(Wu11_monodfit$Yhat_fine, variable=="CO2" | variable=="CFC14")
Wu11_fitb$Model<-c("Monod")

Wu11_fit_all<-rbind(Wu11_fita, Wu11_fitb)

ggplot(subset(Wu11_fit$Yhat, variable=="CO2" | variable=="CFC14"), aes(time, obs))+
  geom_point(cex=6, aes(pch=Soil), fill="grey")+
  scale_shape_manual(values=c(21, 22))+
  geom_line(data=Wu11_fit_all, aes(time, value, color=Model, linetype=Soil), lwd=1.2)+theme_min+
  facet_wrap(Substrate~variable2, scales="free", labeller = label_parsed) + 
  ylab(expression(paste("Carbon pool (", mu, "mol ", g(DW)^{-1}, ")"))) +
  xlab("Time (days)") + scale_color_manual(values = c("indianred", "grey30")) +
  theme(legend.position = c(0.85, 0.8))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Statistics~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Log Likelihood ratio test
-2*(Wu11_monodfit$Gfit$ll-Wu11_fit$Gfit$ll)

round(pchisq(-2*(Wu11_monodfit$Gfit$ll-Wu11_fit$Gfit$ll), df=8,
             lower.tail = F), 3)

#F test - based on residual sum of squares, number of parameters and number of measurements
##Monod model
###residual sum of squares
M1ss = sum((Wu11_monodfit$Yhat$obs-Wu11_monodfit$Yhat$value)^2, na.rm = T)
###number of parameters 
M1p = 20

##Sub-microbial model
###residual sum of squares
M2ss = sum((Wu11_fit$Yhat$obs-Wu11_fit$Yhat$value)^2, na.rm = T)
###number of parameters 
M2p = 28

###total number of measurements
nt = nrow(Wu11_fit$Yhat[!is.na(Wu11_fit$Yhat), ])

####F value =  (M1ss - M2ss)*(nt - M2p)/M2ss/(M2p - M1p)
(M1ss - M2ss)*(nt - M2p)/M2ss/(M2p - M1p)

####associated p value
pf(q=(M1ss - M2ss)*(nt - M2p)/M2ss/(M2p - M1p), 
   df1=(M2p - M1p), 
   df2=(nt - M2p), 
   lower.tail=F)
  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
############################################################################################################################
###################################################Chander and Brookes, 1991################################################
############################################################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#DATA
ChB<-read.csv("Data/ChB1991.csv", sep=',')
#Glucose
ChBG<-subset(ChB, Substrate=='Glucose')
#Maize
ChBR<-subset(ChB, Substrate=='Maize')
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
    #CFC14=(fr*Br+fs*Bs)-wu$Cmicinit[1]
    
    #States
    dG<- - uptake
    dBr<- uptake - transfer
    dBs<- transfer*Y - death
    dCO2<- transfer*(1-Y)
    
    return(list(c(dG, dBr, dBs, dCO2)))
  })
}

#Goodness of fit
good_sub<-function(xG, xR){
  
  names(xG)<-c("v", "k", "f", "m", "Y", "fr", "fs")
  names(xR)<-c("v", "k", "f", "m", "Y", "fr", "fs")
  #Initial Br and Bs
  Bs_iG<-ChBG$Cmicinit[1]/xG[["fs"]]
  Bs_iR<-ChBR$Cmicinit[1]/xR[["fs"]]
  
  #Simulations
  yhat_allG<-as.data.frame(ode(y=c(G=ChBG$Sinit[1], Br=0, Bs=Bs_iG, CO2=0),
                                func = SUBmodel, parms=xG, 
                                times = as.numeric(ChBG$Time)))
  yhat_allG$CFC14<-with(yhat_allG, (xG[["fr"]]*Br+xG[["fs"]]*Bs)-ChBG$Cmicinit[1])
  yhat_allR<-as.data.frame(ode(y=c(G=ChBR$Sinit[1], Br=0, Bs=Bs_iR, CO2=0),
                                func = SUBmodel, parms=xR, 
                                times = as.numeric(ChBR$Time)))
  yhat_allR$CFC14<-with(yhat_allR, (xR[["fr"]]*Br+xR[["fs"]]*Bs)-ChBR$Cmicinit[1])
  #Selecting measured variables
  yhatG<-yhat_allG[, c("time", "CO2", "CFC14")]
  yhatR<-yhat_allR[, c("time", "CO2", "CFC14")]
  #Long format
  YhatG<-melt(yhatG, id.vars=c("time"))
  YhatG$Substrate<-c("Glucose")
  YhatR<-melt(yhatR, id.vars=c("time"))
  YhatR$Substrate<-c("Maize")
  Yhat<-rbind(YhatG, YhatR)
  #Observations
  Yhat$obs<-c(as.numeric(ChBG$CO2cumul), as.numeric(ChBG$Cmic14),
              as.numeric(ChBR$CO2cumul), as.numeric(ChBR$Cmic14))
  Gfit<-Yhat %>% group_by(variable, Substrate) %>% summarise(SSres=sum(((obs-value)^2), na.rm = T),
                                                                   SStot=sum(((obs-mean(obs, na.rm = T))^2), na.rm = T),
                                                                   ll=-sum(((obs-value)^2), na.rm = T)/2/(sd(obs, na.rm = T)^2))
  Gfit$R2<-with(Gfit, 1-SSres/SStot)
  Gfit$N<-length(xG)
  Gfit$AIC<-with(Gfit, 2*N-2*ll)
  
  #Fine temporal scale fo r graphs
  yhat_all_fineG<-as.data.frame(ode(y=c(G=ChBG$Sinit[1], Br=0, Bs=Bs_iG, CO2=0),
                                     func = SUBmodel, parms=xG, 
                                     times = seq(0, 50)))
  yhat_all_fineG$CFC14<-with(yhat_all_fineG, (xG[["fr"]]*Br+xG[["fs"]]*Bs)-ChBG$Cmicinit[1])
  yhat_all_fineR<-as.data.frame(ode(y=c(G=ChBR$Sinit[1], Br=0, Bs=Bs_iR, CO2=0),
                                     func = SUBmodel, parms=xR, 
                                     times = seq(0, 100)))
  yhat_all_fineR$CFC14<-with(yhat_all_fineR, (xR[["fr"]]*Br+xR[["fs"]]*Bs)-ChBR$Cmicinit[1])
  
  Yhat_all_fineG<-melt(yhat_all_fineG, id.vars=c("time"))
  Yhat_all_fineG$Substrate<-c("Glucose")
  Yhat_all_fineR<-melt(yhat_all_fineR, id.vars=c("time"))
  Yhat_all_fineR$Substrate<-c("Maize")
  
  Yhat_all_fine<-rbind(Yhat_all_fineG, Yhat_all_fineR)
  
  rsq_out<-list(Yhat=Yhat, Gfit=Gfit, Yhat_fine = Yhat_all_fine)
  
  return(rsq_out)
}

#Read parameters estimated in python
ChB_optparG<-as.numeric(read.csv("parameters/ChB1991_optparsG.csv", header = F))
ChB_optparR<-as.numeric(read.csv("parameters/ChB1991_optparsR.csv", header = F))

ChB_fit<-good_sub(ChB_optparG, ChB_optparR)
as.data.frame(ChB_fit$Gfit)

#Figure
ChB_fit$Yhat$variable2<-ChB_fit$Yhat$variable
levels(ChB_fit$Yhat$variable2)<-c("CO[2]", "MBC")

ChB_fit$Yhat_fine$variable2<-ChB_fit$Yhat_fine$variable
levels(ChB_fit$Yhat_fine$variable2)<-c("Glucose", "Br", "Bs", "CO[2]", "MBC")

ggplot(subset(ChB_fit$Yhat, variable=="CO2" | variable=="CFC14"), aes(time, obs))+
  geom_point(cex=6, fill="grey")+
  geom_line(data=subset(ChB_fit$Yhat_fine, variable=="CO2" | variable=="CFC14"), 
            aes(time, value), lwd=1.2, color="grey30")+theme_min+
  facet_wrap(Substrate~variable2, scales="free", labeller = label_parsed) + 
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
    #CFC14=(kec*B)-wu$Cmicinit[1]
    
    #States
    dB<- uptake*Y - decay
    dG<- - uptake
    dCO2<- uptake*(1-Y)
    
    return(list(c(dB, dG, dCO2)))
  })
}

#Goodness of fit
monodgood<-function(xG, xR){
  names(xG)<-c("v", "k", "m", "Y", "kec")
  names(xR)<-c("v", "k", "m", "Y", "kec")
  
  #Initial Br and Bs
  B_iG<-ChBG$Cmicinit[1]/xG[["kec"]]
  B_iR<-ChBR$Cmicinit[1]/xR[["kec"]]
  
  #Simulations
  yhat_allG<-as.data.frame(ode(y=c(B=B_iG, G=ChBG$Sinit[1], CO2=0),
                                func = Monod, parms=xG,
                                times = as.numeric(ChBG$Time)))
  yhat_allG$CFC14<-with(yhat_allG, (xG[["kec"]]*B)-ChBG$Cmicinit[1])
  yhat_allR<-as.data.frame(ode(y=c(B=B_iR, G=ChBR$Sinit[1], CO2=0),
                                func = Monod, parms=xR,
                                times = as.numeric(ChBR$Time)))
  yhat_allR$CFC14<-with(yhat_allR, (xR[["kec"]]*B)-ChBR$Cmicinit[1])
  #Selecting measured variables
  yhatG<-yhat_allG[, c("time", "CO2", "CFC14")]
  yhatR<-yhat_allR[, c("time", "CO2", "CFC14")]
  #Long format
  YhatG<-melt(yhatG, id.vars=c("time"))
  YhatG$Substrate<-c("Glucose")
  YhatR<-melt(yhatR, id.vars=c("time"))
  YhatR$Substrate<-c("Maize")
  Yhat<-rbind(YhatG, YhatR)
  #Observations
  Yhat$obs<-c(as.numeric(ChBG$CO2cumul), as.numeric(ChBG$Cmic14),
              as.numeric(ChBR$CO2cumul), as.numeric(ChBR$Cmic14))
  Gfit<-Yhat %>% group_by(variable, Substrate) %>% summarise(SSres=sum(((obs-value)^2), na.rm = T),
                                                                   SStot=sum(((obs-mean(obs, na.rm = T))^2), na.rm = T),
                                                                   ll=-sum(((obs-value)^2), na.rm = T)/2/(sd(obs, na.rm = T)^2))
  Gfit$R2<-with(Gfit, 1-SSres/SStot)
  Gfit$N<-length(xG)
  Gfit$AIC<-with(Gfit, 2*N-2*ll)
  
  #Fine temporal scale for graphs
  yhat_all_fineG<-as.data.frame(ode(y=c(B=B_iG, G=ChBG$Sinit[1], CO2=0),
                                     func = Monod, parms=xG,
                                     times = seq(0, 50)))
  yhat_all_fineG$CFC14<-with(yhat_all_fineG, (xG[["kec"]]*B)-ChBG$Cmicinit[1])
  yhat_all_fineR<-as.data.frame(ode(y=c(B=B_iR, G=ChBR$Sinit[1], CO2=0),
                                     func = Monod, parms=xR,
                                     times = seq(0, 100)))
  yhat_all_fineR$CFC14<-with(yhat_all_fineR, (xR[["kec"]]*B)-ChBR$Cmicinit[1])
  Yhat_all_fineG<-melt(yhat_all_fineG, id.vars=c("time"))
  Yhat_all_fineG$Substrate<-c("Glucose")
  Yhat_all_fineR<-melt(yhat_all_fineR, id.vars=c("time"))
  Yhat_all_fineR$Substrate<-c("Maize")
  
  Yhat_all_fine<-rbind(Yhat_all_fineG, Yhat_all_fineR)
  
  rsq_out<-list(Yhat=Yhat, Gfit=Gfit, Yhat_fine = Yhat_all_fine)
  
  return(rsq_out)
}

#Read parameters estimated in python
ChB_monodoptG<-as.numeric(read.csv("parameters/ChB1991_monodparsG.csv", header = F))
ChB_monodoptR<-as.numeric(read.csv("parameters/ChB1991_monodparsR.csv", header = F))

ChB_monodfit<-monodgood(ChB_monodoptG, ChB_monodoptR)
as.data.frame(ChB_monodfit$Gfit)

#Figure
ChB_monodfit$Yhat_fine$variable2<-ChB_monodfit$Yhat_fine$variable
levels(ChB_monodfit$Yhat_fine$variable2)<-c("B", "Glucose",  "CO[2]", "MBC")

ChB_fita<-subset(ChB_fit$Yhat_fine, variable=="CO2" | variable=="CFC14")
ChB_fita$Model<-c("Sub-microbial")
ChB_fitb<-subset(ChB_monodfit$Yhat_fine, variable=="CO2" | variable=="CFC14")
ChB_fitb$Model<-c("Monod")

ChB_fit_all<-rbind(ChB_fita, ChB_fitb)

ggplot(subset(ChB_fit$Yhat, variable=="CO2" | variable=="CFC14"), aes(time, obs))+
  geom_point(cex=6, fill="grey")+
  geom_line(data=ChB_fit_all, aes(time, value, color=Model), lwd=1.2)+theme_min+
  facet_wrap(Substrate~variable2, scales="free", labeller = label_parsed) + 
  ylab(expression(paste("Carbon pool (", mu, "mol ", g(DW)^{-1}, ")"))) +
  xlab("Time (days)") + scale_color_manual(values = c("indianred", "grey30")) +
  theme(legend.position = c(0.85, 0.8))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Statistics~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Log Likelihood ratio test
-2*(ChB_monodfit$Gfit$ll-ChB_fit$Gfit$ll)

round(pchisq(-2*(ChB_monodfit$Gfit$ll-ChB_fit$Gfit$ll), df=8,
             lower.tail = F), 3)

#F test - based on residual sum of squares, number of parameters and number of measurements
##Monod model
###residual sum of squares
M1ss = sum((ChB_monodfit$Yhat$obs-ChB_monodfit$Yhat$value)^2, na.rm = T)
###number of parameters 
M1p = 10

##Sub-microbial model
###residual sum of squares
M2ss = sum((ChB_fit$Yhat$obs-ChB_fit$Yhat$value)^2, na.rm = T)
###number of parameters 
M2p = 14

###total number of measurements
nt = nrow(ChB_fit$Yhat[!is.na(ChB_fit$Yhat), ])

####F value =  (M1ss - M2ss)*(nt - M2p)/M2ss/(M2p - M1p)
(M1ss - M2ss)*(nt - M2p)/M2ss/(M2p - M1p)

####associated p value
pf(q=(M1ss - M2ss)*(nt - M2p)/M2ss/(M2p - M1p), 
   df1=(M2p - M1p), 
   df2=(nt - M2p), 
   lower.tail=F)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
############################################################################################################################
###################################################Bremer and Kessel, 1990##################################################
############################################################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#DATA
B1990<-read.csv("Data/Bremer1990.csv", sep=',')
B1990HC1<-subset(B1990, Treatment=="HC1")
B1990HC2<-subset(B1990, Treatment=="HC2")
B1990LC<-subset(B1990, Treatment=="LC")
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
    CFC14=(fr*Br+fs*Bs)-B1990$Cmicinit[1]
    
    #States
    dG<- - uptake
    dBr<- uptake - transfer
    dBs<- transfer*Y - death
    dCO2<- transfer*(1-Y)
    
    return(list(c(dG, dBr, dBs, dCO2), CFC14=CFC14))
  })
}
#Goodness of fit
good_sub<-function(x){
  
  names(x)<-c("v", "k", "f", "m", "Y", "fr", "fs")
  #Initial Br and Bs
  Bs_i<-B1990HC1$Cmicinit[1]/x[["fs"]]
  
  
  #Simulations
  ##HC1
  yhat_allHC1<-as.data.frame(ode(y=c(G=B1990HC1$Sinit[1], Br=0, Bs=Bs_i, CO2=0),
                               func = SUBmodel, parms=x, 
                               times = as.numeric(B1990HC1$Time)))
  ##HC2
  yhat_allHC2<-as.data.frame(ode(y=c(G=B1990HC1$Sinit[1], Br=0, Bs=Bs_i, CO2=0),
                                func = SUBmodel, parms=x, 
                                times = as.numeric(B1990HC1$Time)))
  ##LC
  yhat_allLC<-as.data.frame(ode(y=c(G=B1990LC$Sinit[1], Br=0, Bs=Bs_i, CO2=0),
                                func = SUBmodel, parms=x, 
                                times = as.numeric(B1990HC1$Time)))
  #Selecting measured variables
  yhatHC1<-yhat_allHC1[, c("time", "CO2", "CFC14")]
  yhatHC1$Treatment<-c("HC")
  yhatHC2<-yhat_allHC2[, c("time", "CO2", "CFC14")]
  yhatHC2$Treatment<-c("HC")
  yhatLC<-yhat_allLC[, c("time", "CO2", "CFC14")]
  yhatLC$Treatment<-c("LC")
  #Long format
  Yhat<-rbind(melt(yhatHC1, id.vars=c("time", "Treatment")), 
              melt(yhatHC2, id.vars=c("time", "Treatment")),
              melt(yhatLC, id.vars=c("time", "Treatment")))
  
  #Observations
  Yhat$obs<-c(as.numeric(B1990HC1$CO2cumul), as.numeric(B1990HC1$Cmic14),
              as.numeric(B1990HC2$CO2cumul), as.numeric(B1990HC2$Cmic14),
              as.numeric(B1990LC$CO2cumul), as.numeric(B1990LC$Cmic14))
  Gfit<-Yhat %>% group_by(variable) %>% summarise(SSres=sum(((obs-value)^2), na.rm = T),
                                                  SStot=sum(((obs-mean(obs, na.rm = T))^2), na.rm = T),
                                                  ll=-sum(((obs-value)^2), na.rm = T)/2/(sd(obs, na.rm = T)^2))
  Gfit$R2<-with(Gfit, 1-SSres/SStot)
  Gfit$N<-length(x)
  Gfit$AIC<-with(Gfit, 2*N-2*ll)
  
  #Fine temporal scale fo r graphs
  yhat_all_fineHC1<-as.data.frame(ode(y=c(G=B1990HC1$Sinit[1], Br=0, Bs=Bs_i, CO2=0),
                                    func = SUBmodel, parms=x, 
                                    times = seq(0, 7, by=0.1)))
  yhat_all_fineHC2<-as.data.frame(ode(y=c(G=B1990HC1$Sinit[1], Br=0, Bs=Bs_i, CO2=0),
                                      func = SUBmodel, parms=x, 
                                      times = seq(0, 7, by=0.1)))
  yhat_all_fineLC<-as.data.frame(ode(y=c(G=B1990LC$Sinit[1], Br=0, Bs=Bs_i, CO2=0),
                                      func = SUBmodel, parms=x, 
                                      times = seq(0, 7, by=0.1)))
  yhat_all_fineHC1$Treatment<-c("HC")
  yhat_all_fineHC2$Treatment<-c("HC")
  yhat_all_fineLC$Treatment<-c("LC")
  Yhat_all_fine<-rbind(melt(yhat_all_fineHC1, id.vars=c("time", "Treatment")),
                       melt(yhat_all_fineHC2, id.vars=c("time", "Treatment")),
                       melt(yhat_all_fineLC, id.vars=c("time", "Treatment")))
  
  rsq_out<-list(Yhat=Yhat, Gfit=Gfit, Yhat_fine = Yhat_all_fine)
  
  return(rsq_out)
}

#Read parameters estimated in python
B1990_optpar<-as.numeric(read.csv("parameters/Bremer1990_optpars.csv", header = F))

B1990_fit<-good_sub(B1990_optpar)
as.data.frame(B1990_fit$Gfit)

#Figure
B1990_fit$Yhat$variable2<-B1990_fit$Yhat$variable
levels(B1990_fit$Yhat$variable2)<-c("CO[2]", "MBC")

B1990_fit$Yhat_fine$variable2<-B1990_fit$Yhat_fine$variable
levels(B1990_fit$Yhat_fine$variable2)<-c("Glucose", "Br", "Bs", "CO[2]", "MBC")

ggplot(subset(B1990_fit$Yhat, variable=="CO2" | variable=="CFC14"), aes(time, obs))+
  geom_point(cex=6, aes(color=Treatment))+
  scale_color_manual(values=c("black", "grey40"))+
  geom_line(data=subset(B1990_fit$Yhat_fine, variable=="CO2" | variable=="CFC14"), 
            aes(time, value, color=Treatment), lwd=1.2)+theme_min+
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
    CFC14=(kec*B)-B1990$Cmicinit[1]
    
    #States
    dB<- uptake*Y - decay
    dG<- - uptake
    dCO2<- uptake*(1-Y)
    
    return(list(c(dB, dG, dCO2), CFC14=CFC14))
  })
}

#Goodness of fit
monodgood<-function(x){
  names(x)<-c("v", "k", "m", "Y", "kec")
  
  #Initial Br and Bs
  B_i<-B1990$Cmicinit[1]/x[["kec"]]
  
  ##HC1
  yhat_allHC1<-as.data.frame(ode(y=c(B=B_i, G=B1990HC1$Sinit[1], CO2=0),
                                 func = Monod, parms=x, 
                                 times = as.numeric(B1990HC1$Time)))
  ##HC2
  yhat_allHC2<-as.data.frame(ode(y=c(B=B_i, G=B1990HC1$Sinit[1], CO2=0),
                                 func = Monod, parms=x, 
                                 times = as.numeric(B1990HC1$Time)))
  ##LC
  yhat_allLC<-as.data.frame(ode(y=c(B=B_i, G=B1990LC$Sinit[1], CO2=0),
                                func = Monod, parms=x, 
                                times = as.numeric(B1990HC1$Time)))
  #Selecting measured variables
  yhatHC1<-yhat_allHC1[, c("time", "CO2", "CFC14")]
  yhatHC1$Treatment<-c("HC")
  yhatHC2<-yhat_allHC2[, c("time", "CO2", "CFC14")]
  yhatHC2$Treatment<-c("HC")
  yhatLC<-yhat_allLC[, c("time", "CO2", "CFC14")]
  yhatLC$Treatment<-c("LC")
  #Long format
  Yhat<-rbind(melt(yhatHC1, id.vars=c("time", "Treatment")), 
              melt(yhatHC2, id.vars=c("time", "Treatment")),
              melt(yhatLC, id.vars=c("time", "Treatment")))
  
  #Observations
  Yhat$obs<-c(as.numeric(B1990HC1$CO2cumul), as.numeric(B1990HC1$Cmic14),
              as.numeric(B1990HC2$CO2cumul), as.numeric(B1990HC2$Cmic14),
              as.numeric(B1990LC$CO2cumul), as.numeric(B1990LC$Cmic14))
  Gfit<-Yhat %>% group_by(variable) %>% summarise(SSres=sum(((obs-value)^2), na.rm = T),
                                                             SStot=sum(((obs-mean(obs, na.rm = T))^2), na.rm = T),
                                                             ll=-sum(((obs-value)^2), na.rm = T)/2/(sd(obs, na.rm = T)^2))
  Gfit$R2<-with(Gfit, 1-SSres/SStot)
  Gfit$N<-length(x)
  Gfit$AIC<-with(Gfit, 2*N-2*ll)
  
  #Fine temporal scale fo r graphs
  yhat_all_fineHC1<-as.data.frame(ode(y=c(B=B_i, G=B1990HC1$Sinit[1], CO2=0),
                                      func = Monod, parms=x, 
                                      times = seq(0, 7, by=0.1)))
  yhat_all_fineHC2<-as.data.frame(ode(y=c(B=B_i, G=B1990HC1$Sinit[1], CO2=0),
                                      func = Monod, parms=x, 
                                      times = seq(0, 7, by=0.1)))
  yhat_all_fineLC<-as.data.frame(ode(y=c(B=B_i, G=B1990LC$Sinit[1], CO2=0),
                                     func = Monod, parms=x, 
                                     times = seq(0, 7, by=0.1)))
  yhat_all_fineHC1$Treatment<-c("HC")
  yhat_all_fineHC2$Treatment<-c("HC")
  yhat_all_fineLC$Treatment<-c("LC")
  Yhat_all_fine<-rbind(melt(yhat_all_fineHC1, id.vars=c("time", "Treatment")),
                       melt(yhat_all_fineHC2, id.vars=c("time", "Treatment")),
                       melt(yhat_all_fineLC, id.vars=c("time", "Treatment")))
  
  rsq_out<-list(Yhat=Yhat, Gfit=Gfit, Yhat_fine = Yhat_all_fine)
  
  return(rsq_out)
}

#Read parameters estimated in python
B1990_monodopt<-as.numeric(read.csv("parameters/Bremer1990_monodpars.csv", header = F))

B1990_monodfit<-monodgood(B1990_monodopt)
as.data.frame(ChB_monodfit$Gfit)

#Figure
B1990_monodfit$Yhat_fine$variable2<-B1990_monodfit$Yhat_fine$variable
levels(B1990_monodfit$Yhat_fine$variable2)<-c("B", "Glucose",  "CO[2]", "MBC")

B1990_fita<-subset(B1990_fit$Yhat_fine, variable=="CO2" | variable=="CFC14")
B1990_fita$Model<-c("Sub-microbial")
B1990_fitb<-subset(B1990_monodfit$Yhat_fine, variable=="CO2" | variable=="CFC14")
B1990_fitb$Model<-c("Monod")

B1990_fit_all<-rbind(B1990_fita, B1990_fitb)

ggplot(subset(B1990_fit$Yhat, variable=="CO2" | variable=="CFC14"), aes(time, obs))+
  geom_point(cex=6, aes(pch=Treatment))+
  geom_line(data=B1990_fit_all, aes(time, value, color=Model, linetype=Treatment), lwd=1.2)+theme_min+
  facet_wrap(~variable2, scales="free", labeller = label_parsed) + 
  ylab(expression(paste("Carbon pool (", mu, "mol ", g(DW)^{-1}, ")"))) +
  xlab("Time (days)") + scale_color_manual(values = c("indianred", "grey30")) +
  theme(legend.position = c(0.85, 0.8))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Statistics~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Log Likelihood ratio test
-2*(B1990_monodfit$Gfit$ll-B1990_fit$Gfit$ll)

round(pchisq(-2*(B1990_monodfit$Gfit$ll-B1990_fit$Gfit$ll), df=2,
             lower.tail = F), 3)

#F test - based on residual sum of squares, number of parameters and number of measurements
##Monod model
###residual sum of squares
M1ss = sum((B1990_monodfit$Yhat$obs-B1990_monodfit$Yhat$value)^2, na.rm = T)
###number of parameters 
M1p = 5

##Sub-microbial model
###residual sum of squares
M2ss = sum((B1990_fit$Yhat$obs-B1990_fit$Yhat$value)^2, na.rm = T)
###number of parameters 
M2p = 7

###total number of measurements
nt = nrow(B1990_fit$Yhat[!is.na(B1990_fit$Yhat), ])

####F value =  (M1ss - M2ss)*(nt - M2p)/M2ss/(M2p - M1p)
(M1ss - M2ss)*(nt - M2p)/M2ss/(M2p - M1p)

####associated p value
pf(q=(M1ss - M2ss)*(nt - M2p)/M2ss/(M2p - M1p), 
   df1=(M2p - M1p), 
   df2=(nt - M2p), 
   lower.tail=F)
  