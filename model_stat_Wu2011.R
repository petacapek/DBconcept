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
############################################Santruckova et al. 2004#######################################################
############################################################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#DATA
d<-read.csv("Data/Wu2011.csv", sep=',')
d<-subset(d, Substrate=="Glucose")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Monod model~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
MonodA<-function(time, state, pars){
  with(as.list(c(state, pars)),{
    #glucose uptake
    uptake=v*G*B/(k + G)
    #Decay rate
    decay=m*B
    #Chloroform labile C
    CFC14=kec*B-d[d$Soil=="Upland", "Cmicinit"][1]
    
    #States
    dB<- uptake*Y - decay
    dG<- - uptake
    dCO2<- uptake*(1-Y)
    
    return(list(c(dB, dG, dCO2), CFC14=CFC14))
  })
}
MonodB<-function(time, state, pars){
  with(as.list(c(state, pars)),{
    #glucose uptake
    uptake=v*G*B/(k + G)
    #Decay rate
    decay=m*B
    #Chloroform labile C
    CFC14=kec*B-d[d$Soil=="Paddy", "Cmicinit"][1]
    
    #States
    dB<- uptake*Y - decay
    dG<- - uptake
    dCO2<- uptake*(1-Y)
    
    return(list(c(dB, dG, dCO2), CFC14=CFC14))
  })
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Mend model~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
MendA<-function(time, state, pars){
  with(as.list(c(state, pars)),{
    #Define fluxes
    ##glucose uptake
    uptake=1/Y*(v+h)*B*G/(k+G)
    ##growth respiration
    gresp=(1/Y-1)*v*B*G/(k+G)
    ##maintenance respiration
    maintenance=(1/Y-1)*h*B*G/(k+G)
    ##decay
    decay=B*h
    #Chloroform labile C
    CFC14=kec*B-d[d$Soil=="Upland", "Cmicinit"][1]
    
    #Define derivatives
    dB=uptake-gresp-maintenance-decay
    dG=-uptake
    dCO2=gresp+maintenance
    
    return(list(c(dB, dG, dCO2), CFC14=CFC14))
  })
}
MendB<-function(time, state, pars){
  with(as.list(c(state, pars)),{
    #Define fluxes
    ##glucose uptake
    uptake=1/Y*(v+h)*B*G/(k+G)
    ##growth respiration
    gresp=(1/Y-1)*v*B*G/(k+G)
    ##maintenance respiration
    maintenance=(1/Y-1)*h*B*G/(k+G)
    ##decay
    decay=B*h
    #Chloroform labile C
    CFC14=kec*B-d[d$Soil=="Paddy", "Cmicinit"][1]
    
    #Define derivatives
    dB=uptake-gresp-maintenance-decay
    dG=-uptake
    dCO2=gresp+maintenance
    
    return(list(c(dB, dG, dCO2), CFC14=CFC14))
  })
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Pirt model~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
PirtA<-function(time, state, pars){
  with(as.list(c(state, pars)),{
    #glucose uptake
    uptake=v*G*B/(k + G)
    #Decay rate
    decay=m*B
    #Maintenance
    maintenance=B*h
    #Chloroform labile C
    CFC14=kec*B-d[d$Soil=="Upland", "Cmicinit"][1]
    
    #States
    dB<- uptake*Y - decay - maintenance
    dG<- - uptake
    dCO2<- uptake*(1-Y) + maintenance
    
    return(list(c(dB, dG, dCO2), CFC14=CFC14))
  })
}
PirtB<-function(time, state, pars){
  with(as.list(c(state, pars)),{
    #glucose uptake
    uptake=v*G*B/(k + G)
    #Decay rate
    decay=m*B
    #Maintenance
    maintenance=B*h
    #Chloroform labile C
    CFC14=kec*B-d[d$Soil=="Paddy", "Cmicinit"][1]
    
    #States
    dB<- uptake*Y - decay - maintenance
    dG<- - uptake
    dCO2<- uptake*(1-Y) + maintenance
    
    return(list(c(dB, dG, dCO2), CFC14=CFC14))
  })
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~DEB model~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
DEBmodelA<-function(time, state, pars){
  with(as.list(c(state, pars)),{
    #Define fluxes
    ##scaling function for substrate
    f=G/(500+G)
    ##growth rate
    growth=(v*e-m*g)/(m+g)
    ##CO2 yield
    Yco2=((v*f/Im)+(m*g/Im)+max(g*growth/Im,0))*ce/f
    
    #Chloroform labile C
    CFC14=(Cwcfc+Cecfc*e)*w-d[d$Soil=="Upland", "Cmicinit"][1]
    
    #States
    #Define derivatives
    dG=-f*w*Im
    de=v*f-v*e
    dw=growth*w
    dCO2=f*w*Im*Yco2
    
    return(list(c(dG, de, dw, dCO2), CFC14=CFC14))
  })
}
DEBmodelB<-function(time, state, pars){
  with(as.list(c(state, pars)),{
    #Define fluxes
    ##scaling function for substrate
    f=G/(500+G)
    ##growth rate
    growth=(v*e-m*g)/(m+g)
    ##CO2 yield
    Yco2=((v*f/Im)+(m*g/Im)+max(g*growth/Im,0))*ce/f
    
    #Chloroform labile C
    CFC14=(Cwcfc+Cecfc*e)*w-d[d$Soil=="Paddy", "Cmicinit"][1]
    
    #States
    #Define derivatives
    dG=-f*w*Im
    de=v*f-v*e
    dw=growth*w
    dCO2=f*w*Im*Yco2
    
    return(list(c(dG, de, dw, dCO2), CFC14=CFC14))
  })
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Goodness of fit
good_all<-function(xmonodA, xmendA, xpirtA, xdebA,
                   xmonodB, xmendB, xpirtB, xdebB){
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Monod~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  pmonodA<-xmonodA
  pmonodB<-xmonodB
  names(pmonodA)<-c("v", "k", "m", "Y", "kec")
  names(pmonodB)<-c("v", "k", "m", "Y", "kec")
  #Initial B
  B_iA<-d[d$Soil=="Upland", "Cmicinit"][1]/pmonodA[["kec"]]
  B_iB<-d[d$Soil=="Paddy", "Cmicinit"][1]/pmonodB[["kec"]]
  #Simulations
  yhat_monodA<-as.data.frame(ode(y=c(B=B_iA, G=d$Sinit[1], CO2=0),
                                func = MonodA, parms=pmonodA, 
                                times = as.numeric(d[d$Soil=="Upland", "Time"])*24))
  
  #Selecting measured variables
  yhat_monodA<-yhat_monodA[, c("time", "CO2", "CFC14")]
  #Long format
  Yhat_monodA<-melt(yhat_monodA, id.vars=c("time"))
  #Observations
  Yhat_monodA$obs<-c(as.numeric(d[d$Soil=="Upland", "CO2cumul"]), as.numeric(d[d$Soil=="Upland", "Cmic14"]))
  Yhat_monodA$Soil<-"Upland"
  
  yhat_monodB<-as.data.frame(ode(y=c(B=B_iB, G=d$Sinit[1], CO2=0),
                                 func = MonodB, parms=pmonodB, 
                                 times = as.numeric(d[d$Soil=="Paddy", "Time"])*24))
  
  #Selecting measured variables
  yhat_monodB<-yhat_monodB[, c("time", "CO2", "CFC14")]
  #Long format
  Yhat_monodB<-melt(yhat_monodB, id.vars=c("time"))
  #Observations
  Yhat_monodB$obs<-c(as.numeric(d[d$Soil=="Paddy", "CO2cumul"]), as.numeric(d[d$Soil=="Paddy", "Cmic14"]))
  Yhat_monodB$Soil<-"Paddy"
  
  
  Yhat_monod<-rbind(Yhat_monodA, Yhat_monodB)
  
  Gfit_monod<-Yhat_monod %>% group_by(Soil, variable) %>% summarise(SSres=sum(((obs-value)^2), na.rm = T),
                                                  SStot=sum(((obs-mean(obs, na.rm = T))^2), na.rm = T),
                                                  ll=-sum(((obs-value)^2), na.rm = T)/2/(sd(obs, na.rm = T)^2))
  Gfit_monod$R2<-with(Gfit_monod, 1-SSres/SStot)
  Gfit_monod$N<-length(pmonod)
  Gfit_monod$AIC<-with(Gfit_monod, 2*N-2*ll)
  Gfit_monod$Model<-"Monod"
  
  #Fine temporal scale fo r graphs
  yhat_monod_finea<-as.data.frame(ode(y=c(B=B_iA, G=d$Sinit[1], CO2=0),
                                   func = MonodA, parms=pmonodA,
                                   times = seq(0, 100, by=0.1)*24))
  Yhat_monod_finea<-melt(yhat_monod_finea, id.vars=c("time"))
  Yhat_monod_finea$Model<-"Monod"
  Yhat_monod_finea$Soil<-"Upland"
  yhat_monod_fineb<-as.data.frame(ode(y=c(B=B_iB, G=d$Sinit[1], CO2=0),
                                      func = MonodB, parms=pmonodB,
                                      times = seq(0, 100, by=0.1)*24))
  Yhat_monod_fineb<-melt(yhat_monod_fineb, id.vars=c("time"))
  Yhat_monod_fineb$Model<-"Monod"
  Yhat_monod_fineb$Soil<-"Paddy"
  
  Yhat_monod_fine<-rbind(Yhat_monod_finea, Yhat_monod_fineb)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Mend~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  pmendA<-xmendA
  pmendB<-xmendB
  names(pmendA)<-c("v", "k", "h", "Y", "kec")
  names(pmendB)<-c("v", "k", "h", "Y", "kec")
  #Initial B
  B_iA<-d[d$Soil=="Upland", "Cmicinit"][1]/pmendA[["kec"]]
  B_iB<-d[d$Soil=="Paddy", "Cmicinit"][1]/pmendB[["kec"]]
  #Simulations
  yhat_mendA<-as.data.frame(ode(y=c(B=B_iA, G=d$Sinit[1], CO2=0),
                               func = MendA, parms=pmendA,
                               times = as.numeric(d[d$Soil=="Upland", "Time"])*24))
  
  #Selecting measured variables
  yhat_mendA<-yhat_mendA[, c("time", "CO2", "CFC14")]
  #Long format
  Yhat_mendA<-melt(yhat_mendA, id.vars=c("time"))
  #Observations
  Yhat_mendA$obs<-c(as.numeric(d[d$Soil=="Upland", "CO2cumul"]), as.numeric(d[d$Soil=="Upland", "Cmic14"]))
  Yhat_mendA$Soil<-"Upland"
  
  yhat_mendB<-as.data.frame(ode(y=c(B=B_iB, G=d$Sinit[1], CO2=0),
                                func = MendB, parms=pmendB,
                                times = as.numeric(d[d$Soil=="Paddy", "Time"])*24))
  
  #Selecting measured variables
  yhat_mendB<-yhat_mendB[, c("time", "CO2", "CFC14")]
  #Long format
  Yhat_mendB<-melt(yhat_mendB, id.vars=c("time"))
  #Observations
  Yhat_mendB$obs<-c(as.numeric(d[d$Soil=="Paddy", "CO2cumul"]), as.numeric(d[d$Soil=="Paddy", "Cmic14"]))
  Yhat_mendB$Soil<-"Paddy"
  
  
  Yhat_mend<-rbind(Yhat_mendA, Yhat_mendB)
  
  Gfit_mend<-Yhat_mend %>% group_by(Soil, variable) %>% summarise(SSres=sum(((obs-value)^2), na.rm = T),
                                                              SStot=sum(((obs-mean(obs, na.rm = T))^2), na.rm = T),
                                                              ll=-sum(((obs-value)^2), na.rm = T)/2/(sd(obs, na.rm = T)^2))
  Gfit_mend$R2<-with(Gfit_mend, 1-SSres/SStot)
  Gfit_mend$N<-length(pmend)
  Gfit_mend$AIC<-with(Gfit_mend, 2*N-2*ll)
  Gfit_mend$Model<-"Mend"
  
  #Fine temporal scale fo r graphs
  yhat_mend_finea<-as.data.frame(ode(y=c(B=B_iA, G=d$Sinit[1], CO2=0),
                                     func = MendA, parms=pmendA,
                                     times = seq(0, 100, by=0.1)*24))
  Yhat_mend_finea<-melt(yhat_mend_finea, id.vars=c("time"))
  Yhat_mend_finea$Model<-"Mend"
  Yhat_mend_finea$Soil<-"Upland"
  yhat_mend_fineb<-as.data.frame(ode(y=c(B=B_iB, G=d$Sinit[1], CO2=0),
                                     func = MendB, parms=pmendB,
                                     times = seq(0, 100, by=0.1)*24))
  Yhat_mend_fineb<-melt(yhat_mend_fineb, id.vars=c("time"))
  Yhat_mend_fineb$Model<-"Mend"
  Yhat_mend_fineb$Soil<-"Paddy"
  
  Yhat_mend_fine<-rbind(Yhat_mend_finea, Yhat_mend_fineb)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Pirt~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  ppirtA<-xpirtA
  ppirtB<-xpirtB
  names(ppirtA)<-c("v", "k", "m", "h", "Y", "kec")
  names(ppirtB)<-c("v", "k", "m", "h", "Y", "kec")
  #Initial B
  B_iA<-d[d$Soil=="Upland", "Cmicinit"][1]/ppirtA[["kec"]]
  B_iB<-d[d$Soil=="Paddy", "Cmicinit"][1]/ppirtB[["kec"]]
  #Simulations
  yhat_pirtA<-as.data.frame(ode(y=c(B=B_iA, G=d$Sinit[1], CO2=0),
                               func = PirtA, parms=ppirtA,
                               times = as.numeric(d[d$Soil=="Upland", "Time"])*24))
  
  #Selecting measured variables
  yhat_pirtA<-yhat_pirtA[, c("time", "CO2", "CFC14")]
  #Long format
  Yhat_pirtA<-melt(yhat_pirtA, id.vars=c("time"))
  #Observations
  Yhat_pirtA$obs<-c(as.numeric(d[d$Soil=="Upland", "CO2cumul"]), as.numeric(d[d$Soil=="Upland", "Cmic14"]))
  Yhat_pirtA$Soil<-"Upland"
  
  yhat_pirtB<-as.data.frame(ode(y=c(B=B_iB, G=d$Sinit[1], CO2=0),
                                func = PirtB, parms=ppirtB,
                                times = as.numeric(d[d$Soil=="Paddy", "Time"])*24))
  
  #Selecting measured variables
  yhat_pirtB<-yhat_pirtB[, c("time", "CO2", "CFC14")]
  #Long format
  Yhat_pirtB<-melt(yhat_pirtB, id.vars=c("time"))
  #Observations
  Yhat_pirtB$obs<-c(as.numeric(d[d$Soil=="Paddy", "CO2cumul"]), as.numeric(d[d$Soil=="Paddy", "Cmic14"]))
  Yhat_pirtB$Soil<-"Paddy"
  
  Yhat_pirt<-rbind(Yhat_pirtA, Yhat_pirtB)
  
  Gfit_pirt<-Yhat_pirt %>% group_by(Soil, variable) %>% summarise(SSres=sum(((obs-value)^2), na.rm = T),
                                                            SStot=sum(((obs-mean(obs, na.rm = T))^2), na.rm = T),
                                                            ll=-sum(((obs-value)^2), na.rm = T)/2/(sd(obs, na.rm = T)^2))
  Gfit_pirt$R2<-with(Gfit_pirt, 1-SSres/SStot)
  Gfit_pirt$N<-length(ppirt)
  Gfit_pirt$AIC<-with(Gfit_pirt, 2*N-2*ll)
  Gfit_pirt$Model<-"Pirt"
  
  #Fine temporal scale fo r graphs
  yhat_pirt_finea<-as.data.frame(ode(y=c(B=B_iA, G=d$Sinit[1], CO2=0),
                                    func = PirtA, parms=ppirtA,
                                    times = seq(0, 100, by=0.1)*24))
  Yhat_pirt_finea<-melt(yhat_pirt_finea, id.vars=c("time"))
  Yhat_pirt_finea$Model<-"Pirt"
  Yhat_pirt_finea$Soil<-"Upland"
  
  yhat_pirt_fineb<-as.data.frame(ode(y=c(B=B_iB, G=d$Sinit[1], CO2=0),
                                     func = PirtB, parms=ppirtB,
                                     times = seq(0, 100, by=0.1)*24))
  Yhat_pirt_fineb<-melt(yhat_pirt_fineb, id.vars=c("time"))
  Yhat_pirt_fineb$Model<-"Pirt"
  Yhat_pirt_fineb$Soil<-"Paddy"
  
  Yhat_pirt_fine<-rbind(Yhat_pirt_finea, Yhat_pirt_fineb)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~DEB~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  pdebA<-xdebA
  pdebB<-xdebB
  names(pdebA)<-c("Im", "v", "m", "g", "ce", "Cwcfc", "Cecfc")
  names(pdebB)<-c("Im", "v", "m", "g", "ce", "Cwcfc", "Cecfc")
  #Initial w 
  w_iA = d[d$Soil=="Upland", "Cmicinit"][1]/pdebA[["Cwcfc"]]
  w_iB = d[d$Soil=="Paddy", "Cmicinit"][1]/pdebB[["Cwcfc"]]
  #Br_i<-(mar$Cmicinit[1]-p[["fs"]]*Bs_i)/p[["fr"]]
  #Simulations
  yhat_debA<-as.data.frame(ode(y=c(G=d$Sinit[1], e=0, w=w_iA, CO2=0),
                              func = DEBmodelA, parms=pdebA,
                              times = as.numeric(d[d$Soil=="Upland", "Time"])*24))
  #Selecting measured variables
  yhat_debA<-yhat_debA[, c("time", "CO2", "CFC14")]
  #Long format
  Yhat_debA<-melt(yhat_debA, id.vars=c("time"))
  #Observations
  Yhat_debA$obs<-c(as.numeric(d[d$Soil=="Upland", "CO2cumul"]), as.numeric(d[d$Soil=="Upland", "Cmic14"]))
  Yhat_debA$Soil<-"Upland"
  
  yhat_debB<-as.data.frame(ode(y=c(G=d$Sinit[1], e=0, w=w_iB, CO2=0),
                               func = DEBmodelB, parms=pdebB,
                               times = as.numeric(d[d$Soil=="Paddy", "Time"])*24))
  #Selecting measured variables
  yhat_debB<-yhat_debB[, c("time", "CO2", "CFC14")]
  #Long format
  Yhat_debB<-melt(yhat_debB, id.vars=c("time"))
  #Observations
  Yhat_debB$obs<-c(as.numeric(d[d$Soil=="Paddy", "CO2cumul"]), as.numeric(d[d$Soil=="Paddy", "Cmic14"]))
  Yhat_debB$Soil<-"Paddy"
  
  Yhat_deb<-rbind(Yhat_debA, Yhat_debB)
  
  Gfit_deb<-Yhat_deb %>% group_by(Soil, variable) %>% summarise(SSres=sum(((obs-value)^2), na.rm = T),
                                                  SStot=sum(((obs-mean(obs, na.rm = T))^2), na.rm = T),
                                                  ll=-sum(((obs-value)^2), na.rm = T)/2/(sd(obs, na.rm = T)^2))
  Gfit_deb$R2<-with(Gfit_deb, 1-SSres/SStot)
  Gfit_deb$N<-length(pdeb)
  Gfit_deb$AIC<-with(Gfit_deb, 2*N-2*ll)
  Gfit_deb$Model<-"DEB"
  
  #Fine temporal scale fo r graphs
  yhat_deb_finea<-as.data.frame(ode(y=c(G=d$Sinit[1], e=0, w=w_iA, CO2=0),
                                   func = DEBmodelA, parms=pdebA, method = "daspk",
                                   times = seq(0, 100, by=0.1)*24))
  Yhat_deb_finea<-melt(yhat_deb_finea, id.vars=c("time"))
  Yhat_deb_finea$Model<-"DEB"
  Yhat_deb_finea$Soil<-"Upland"
  
  yhat_deb_fineb<-as.data.frame(ode(y=c(G=d$Sinit[1], e=0, w=w_iB, CO2=0),
                                    func = DEBmodelB, parms=pdebB, method = "daspk",
                                    times = seq(0, 100, by=0.1)*24))
  Yhat_deb_fineb<-melt(yhat_deb_fineb, id.vars=c("time"))
  Yhat_deb_fineb$Model<-"DEB"
  Yhat_deb_fineb$Soil<-"Paddy"
  
  Yhat_deb_fine<-rbind(Yhat_deb_finea, Yhat_deb_fineb)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Statistics~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #For each variable
  #Ftest
  ##Monod vs DEB
  nt<-as.numeric(as.data.frame(subset(Yhat_deb, !is.na(obs)) %>% group_by(Soil, variable) %>% summarize(length(obs)))[,3])
  
  stat1ll<-(Gfit_monod$SSres - Gfit_deb$SSres)*(nt - length(pdebA))/Gfit_deb$SSres/(length(pdebA) - length(pmonodA))
  stat1p<-pf(q=stat1ll, df1=(length(pdebA) - length(pmonodA)), df2=(nt - length(pdeb)), lower.tail=F)
  ##Mend vs DEB
  stat2ll<-(Gfit_mend$SSres - Gfit_deb$SSres)*(nt - length(pdebA))/Gfit_deb$SSres/(length(pdebA) - length(pmendA))
  stat2p<-pf(q=stat2ll, df1=(length(pdebA) - length(pmendA)), df2=(nt - length(pdebA)), lower.tail=F)
  ##Pirt vs DEB
  stat3ll<-(Gfit_pirt$SSres - Gfit_deb$SSres)*(nt - length(pdebA))/Gfit_deb$SSres/(length(pdebA) - length(ppirtA))
  stat3p<-pf(q=stat3ll, df1=(length(pdebA) - length(ppirtA)), df2=(nt - length(pdebA)), lower.tail=F)
  
  stat_eachF<-data.frame(Variable=c("CO2", "CFC14"),
                         MvsDEB_F=stat1ll, MvsDEB_p=stat1p,
                         MdvsDEB_F=stat2ll, MdvsDEB_p=stat2p,
                         PvsDEB_F=stat3ll, PvsDEB_p=stat3p)
  
  #LRtest
  ##Monod vs DEB
  stat1ll<--2*(Gfit_monod$ll-Gfit_deb$ll)
  stat1p<-round(pchisq(-2*(Gfit_monod$ll-Gfit_deb$ll), df=(length(pdebA)-length(pmonodA)),
                       lower.tail = F), 3)
  ##Mend vs DEB
  stat2ll<--2*(Gfit_mend$ll-Gfit_deb$ll)
  stat2p<-round(pchisq(-2*(Gfit_mend$ll-Gfit_deb$ll), df=(length(pdebA)-length(pmendA)),
                       lower.tail = F), 3)
  
  ##Pirt vs DEB
  stat3ll<--2*(Gfit_pirt$ll-Gfit_deb$ll)
  stat3p<-round(pchisq(-2*(Gfit_pirt$ll-Gfit_deb$ll), df=(length(pdebA)-length(ppirtA)),
                       lower.tail = F), 3)
  
  stat_eachLR<-data.frame(Variable=c("CO2", "CFC14"),
                          MvsDEB_ll=stat1ll, MvsDEB_p=stat1p,
                          MdvsDEB_ll=stat2ll, MdvsDEB_p=stat2p,
                          PvsDEB_ll=stat3ll, PvsDEB_p=stat3p)
  
  #Across variables
  #F test
  ###residual sum of squares Monod, Mend, Pirt
  M1ssmonod = sum((Yhat_monod$obs-Yhat_monod$value)^2, na.rm = T)
  M1ssmend = sum((Yhat_mend$obs-Yhat_mend$value)^2, na.rm = T)
  M1sspirt = sum((Yhat_pirt$obs-Yhat_pirt$value)^2, na.rm = T)
  ###number of parameters Monod, Mend, Pirt
  M1pmonod = length(pmonodA)*2
  M1pmend = length(pmendA)*2
  M1ppirt = length(ppirtA)*2
  
  ##DEB model
  ###residual sum of squares
  M2ss = sum((Yhat_deb$obs-Yhat_deb$value)^2, na.rm = T)
  ###number of parameters 
  M2p = length(pdebA)*2
  
  ###total number of measurements
  nt = nrow(Yhat_deb[!is.na(Yhat_deb), ])
  
  ####F value =  (M1ss - M2ss)*(nt - M2p)/M2ss/(M2p - M1p)
  fstat1=(M1ssmonod - M2ss)*(nt - M2p)/M2ss/(M2p - M1pmonod)
  fstat2=(M1ssmend - M2ss)*(nt - M2p)/M2ss/(M2p - M1pmend)
  fstat3=(M1sspirt - M2ss)*(nt - M2p)/M2ss/(M2p - M1ppirt)
  
  ####associated p value
  fstatp1<-pf(q=fstat1, df1=(M2p - M1pmonod), df2=(nt - M2p), lower.tail=F)
  fstatp2<-pf(q=fstat2, df1=(M2p - M1pmend), df2=(nt - M2p), lower.tail=F)
  fstatp3<-pf(q=fstat3, df1=(M2p - M1ppirt), df2=(nt - M2p), lower.tail=F)
  
  stat_allF<-data.frame(Model=c("Monod", "Mend", "Pirt"),
                        F=c(fstat1, fstat2, fstat3),
                        p=c(fstatp1, fstatp2, fstatp3))
  ##LR test
  ##Monod vs DEB
  stat1llall<--2*(with(Yhat_monod, -sum(((obs-value)^2), na.rm = T)/2/(sd(obs, na.rm = T)^2))-
                    with(Yhat_deb, -sum(((obs-value)^2), na.rm = T)/2/(sd(obs, na.rm = T)^2)))
  stat1pall<-round(pchisq(-2*(with(Yhat_monod, -sum(((obs-value)^2), na.rm = T)/2/(sd(obs, na.rm = T)^2))-
                                with(Yhat_deb, -sum(((obs-value)^2), na.rm = T)/2/(sd(obs, na.rm = T)^2))), 
                          df=(length(pdeb)*2-length(pmonod)*2),
                          lower.tail = F), 3)
  ##Mend vs DEB
  stat2llall<--2*(with(Yhat_mend, -sum(((obs-value)^2), na.rm = T)/2/(sd(obs, na.rm = T)^2))-
                    with(Yhat_deb, -sum(((obs-value)^2), na.rm = T)/2/(sd(obs, na.rm = T)^2)))
  stat2pall<-round(pchisq(-2*(with(Yhat_mend, -sum(((obs-value)^2), na.rm = T)/2/(sd(obs, na.rm = T)^2))-
                                with(Yhat_deb, -sum(((obs-value)^2), na.rm = T)/2/(sd(obs, na.rm = T)^2))), 
                          df=(length(pdeb)*2-length(pmonod)*2),
                          lower.tail = F), 3)
  
  ##Pirt vs DEB
  stat3llall<--2*(with(Yhat_pirt, -sum(((obs-value)^2), na.rm = T)/2/(sd(obs, na.rm = T)^2))-
                    with(Yhat_deb, -sum(((obs-value)^2), na.rm = T)/2/(sd(obs, na.rm = T)^2)))
  stat3pall<-round(pchisq(-2*(with(Yhat_pirt, -sum(((obs-value)^2), na.rm = T)/2/(sd(obs, na.rm = T)^2))-
                                with(Yhat_deb, -sum(((obs-value)^2), na.rm = T)/2/(sd(obs, na.rm = T)^2))), 
                          df=(length(pdeb)*2-length(pmonod)*2),
                          lower.tail = F), 3)
  
  stat_allLR<-data.frame(Model=c("Monod", "Mend", "Pirt"),
                         LR=c(stat1llall, stat2llall, stat3llall),
                         p=c(stat1pall, stat2pall, stat3pall))
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~mering data~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  Yhat_all_fine<-rbind(Yhat_monod_fine, Yhat_mend_fine, Yhat_pirt_fine, Yhat_deb_fine)
  Gfit<-rbind(Gfit_monod, Gfit_mend, Gfit_pirt, Gfit_deb)
  
  good_all_out<-list(Yhat=Yhat_monod, Gfit=Gfit, Yhat_fine = Yhat_all_fine, 
                     stat_eachF=stat_eachF, stat_eachLR=stat_eachLR, 
                     stat_allF=stat_allF, stat_allLR=stat_allLR)
  
  return(good_all_out)
}
  
#Read parameters estimated in python
##Monod
monod_parA<-as.numeric(read.csv("parameters/wu11_monodparsA.csv", header = F))
monod_parB<-as.numeric(read.csv("parameters/wu11_monodparsB.csv", header = F))
##Mend
mend_parA<-as.numeric(read.csv("parameters/wu11_mendparsA.csv", header = F))
mend_parB<-as.numeric(read.csv("parameters/wu11_mendparsB.csv", header = F))
##Pirt
pirt_parA<-as.numeric(read.csv("parameters/wu11_pirtparsA.csv", header = F))
pirt_parB<-as.numeric(read.csv("parameters/wu11_pirtparsB.csv", header = F))
##DEB
deb_parA<-as.numeric(read.csv("parameters/wu11_debparsA.csv", header = F))
deb_parB<-as.numeric(read.csv("parameters/wu11_debparsB.csv", header = F))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Models evaluation~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
Wu11_fit<-good_all(monod_parA, mend_parA, pirt_parA, deb_parA,
                   monod_parB, mend_parB, pirt_parB, deb_parB)
Wu11_fit$Gfit
Wu11_fit$stat_eachF
Wu11_fit$stat_eachLR
Wu11_fit$stat_allF
Wu11_fit$stat_allLR
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Figure~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
Rawdata<-subset(Wu11_fit$Yhat, variable=="CO2" | variable=="CFC14")

Fits<-subset(Wu11_fit$Yhat_fine, variable=="CO2" |  variable=="CFC14")
ggplot(subset(Rawdata), aes(time, obs))+
  geom_point(cex=6, pch=21, aes(fill=Soil))+
  scale_fill_manual(values = c("grey60", "white"))+
  geom_line(data=subset(Fits), aes(time, value, color=Model, lty=Soil), lwd=1.2)+theme_min+
  facet_wrap(~variable, scales="free", labeller = label_parsed) + 
  ylab(expression(paste("Carbon pool (", mu, "mol ", g(DW)^{-1}, ")"))) +
  xlab("Time (days)")
    


