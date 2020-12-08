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
d<-read.csv("Data/Blagodatskaya2011.csv", sep=',')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Monod model~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
Monod<-function(time, state, pars){
  with(as.list(c(state, pars)),{
    #glucose uptake
    uptake=v*G*B/(k + G)
    #Decay rate
    decay=m*B
    #Chloroform labile C
    CFC14=kec*B-d$Cmicinit[1]
    
    #States
    dB<- uptake*Y - decay
    dG<- - uptake
    dCO2<- uptake*(1-Y)
    
    return(list(c(dB, dG, dCO2), CFC14=CFC14))
  })
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Mend model~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
Mend<-function(time, state, pars){
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
    CFC14=kec*B-d$Cmicinit[1]
    
    #Define derivatives
    dB=uptake-gresp-maintenance-decay
    dG=-uptake
    dCO2=gresp+maintenance
    
    return(list(c(dB, dG, dCO2), CFC14=CFC14))
  })
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Pirt model~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
Pirt<-function(time, state, pars){
  with(as.list(c(state, pars)),{
    #glucose uptake
    uptake=v*G*B/(k + G)
    #Decay rate
    decay=m*B
    #Maintenance
    maintenance=B*h
    #Chloroform labile C
    CFC14=kec*B-d$Cmicinit[1]
    
    #States
    dB<- uptake*Y - decay - maintenance
    dG<- - uptake
    dCO2<- uptake*(1-Y) + maintenance
    
    return(list(c(dB, dG, dCO2), CFC14=CFC14))
  })
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~DEB model~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
DEBmodel<-function(time, state, pars){
  with(as.list(c(state, pars)),{
    #Define fluxes
    ##scaling function for substrate
    f=G/(500+G)
    ##growth rate
    growth=(v*e-m*g)/(m+g)
    ##CO2 yield
    Yco2=((v*f/Im)+(m*g/Im)+max(g*growth/Im,0))*ce/f
    
    #Chloroform labile C
    CFC14=(Cwcfc+Cecfc*e)*w-d$Cmicinit[1]
    
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
good_all<-function(xmonod, xmend, xpirt, xdeb){
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Monod~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  pmonod<-xmonod
  names(pmonod)<-c("v", "k", "m", "Y", "kec")
  #Initial B
  B_i<-d$Cmicinit[1]/pmonod[["kec"]]
  #Simulations
  yhat_monod<-as.data.frame(ode(y=c(B=B_i, G=d$Sinit[1], CO2=0),
                              func = Monod, parms=pmonod,
                              times = as.numeric(d$Time)*24))
    
  #Selecting measured variables
  yhat_monod<-yhat_monod[, c("time", "CO2", "CFC14")]
  #Long format
  Yhat_monod<-melt(yhat_monod, id.vars=c("time"))
  #Observations
  Yhat_monod$obs<-c(as.numeric(d$CO214), as.numeric(d$Cmic14))
  Gfit_monod<-Yhat_monod %>% group_by(variable) %>% summarise(SSres=sum(((obs-value)^2), na.rm = T),
                                                  SStot=sum(((obs-mean(obs, na.rm = T))^2), na.rm = T),
                                                  ll=-sum(((obs-value)^2), na.rm = T)/2/(sd(obs, na.rm = T)^2))
  Gfit_monod$R2<-with(Gfit_monod, 1-SSres/SStot)
  Gfit_monod$N<-length(pmonod)
  Gfit_monod$AIC<-with(Gfit_monod, 2*N-2*ll)
  Gfit_monod$Model<-"Monod"
  
  #Fine temporal scale fo r graphs
  yhat_monod_fine<-as.data.frame(ode(y=c(B=B_i, G=d$Sinit[1], CO2=0),
                                   func = Monod, parms=pmonod,
                                   times = seq(0, 55, by=0.1)*24))
  Yhat_monod_fine<-melt(yhat_monod_fine, id.vars=c("time"))
  Yhat_monod_fine$Model<-"Monod"
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Mend~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  pmend<-xmend
  names(pmend)<-c("v", "k", "h", "Y", "kec")
  #Initial B
  B_i<-d$Cmicinit[1]/pmend[["kec"]]
  #Simulations
  yhat_mend<-as.data.frame(ode(y=c(B=B_i, G=d$Sinit[1], CO2=0),
                                func = Mend, parms=pmend,
                                times = as.numeric(d$Time)*24))
  
  #Selecting measured variables
  yhat_mend<-yhat_mend[, c("time", "CO2", "CFC14")]
  #Long format
  Yhat_mend<-melt(yhat_mend, id.vars=c("time"))
  #Observations
  Yhat_mend$obs<-c(as.numeric(d$CO214), as.numeric(d$Cmic14))
  Gfit_mend<-Yhat_mend %>% group_by(variable) %>% summarise(SSres=sum(((obs-value)^2), na.rm = T),
                                                              SStot=sum(((obs-mean(obs, na.rm = T))^2), na.rm = T),
                                                              ll=-sum(((obs-value)^2), na.rm = T)/2/(sd(obs, na.rm = T)^2))
  Gfit_mend$R2<-with(Gfit_mend, 1-SSres/SStot)
  Gfit_mend$N<-length(pmend)
  Gfit_mend$AIC<-with(Gfit_mend, 2*N-2*ll)
  Gfit_mend$Model<-"Mend"
  
  #Fine temporal scale fo r graphs
  yhat_mend_fine<-as.data.frame(ode(y=c(B=B_i, G=d$Sinit[1], CO2=0),
                                     func = Mend, parms=pmend,
                                     times = seq(0, 55, by=0.1)*24))
  Yhat_mend_fine<-melt(yhat_mend_fine, id.vars=c("time"))
  Yhat_mend_fine$Model<-"Mend"
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Pirt~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  ppirt<-xpirt
  names(ppirt)<-c("v", "k", "m", "h", "Y", "kec")
  #Initial B
  B_i<-d$Cmicinit[1]/ppirt[["kec"]]
  #Simulations
  yhat_pirt<-as.data.frame(ode(y=c(B=B_i, G=d$Sinit[1], CO2=0),
                               func = Pirt, parms=ppirt,
                               times = as.numeric(d$Time)*24))
  
  #Selecting measured variables
  yhat_pirt<-yhat_pirt[, c("time", "CO2", "CFC14")]
  #Long format
  Yhat_pirt<-melt(yhat_pirt, id.vars=c("time"))
  #Observations
  Yhat_pirt$obs<-c(as.numeric(d$CO214), as.numeric(d$Cmic14))
  Gfit_pirt<-Yhat_pirt %>% group_by(variable) %>% summarise(SSres=sum(((obs-value)^2), na.rm = T),
                                                            SStot=sum(((obs-mean(obs, na.rm = T))^2), na.rm = T),
                                                            ll=-sum(((obs-value)^2), na.rm = T)/2/(sd(obs, na.rm = T)^2))
  Gfit_pirt$R2<-with(Gfit_pirt, 1-SSres/SStot)
  Gfit_pirt$N<-length(ppirt)
  Gfit_pirt$AIC<-with(Gfit_pirt, 2*N-2*ll)
  Gfit_pirt$Model<-"Pirt"
  
  #Fine temporal scale fo r graphs
  yhat_pirt_fine<-as.data.frame(ode(y=c(B=B_i, G=d$Sinit[1], CO2=0),
                                    func = Pirt, parms=ppirt,
                                    times = seq(0, 2, by=0.1)*24))
  Yhat_pirt_fine<-melt(yhat_pirt_fine, id.vars=c("time"))
  Yhat_pirt_fine$Model<-"Pirt"
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~DEB~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  pdeb<-xdeb
  names(pdeb)<-c("Im", "v", "m", "g", "ce", "Cwcfc", "Cecfc")
  #Initial w 
  w_i<-d$Cmicinit[1]/pdeb[["Cwcfc"]]
  #Br_i<-(mar$Cmicinit[1]-p[["fs"]]*Bs_i)/p[["fr"]]
  #Simulations
  yhat_deb<-as.data.frame(ode(y=c(G=d$Sinit[1], e=0, w=w_i, CO2=0),
                              func = DEBmodel, parms=pdeb,
                              times = as.numeric(d$Time)*24))
  #Selecting measured variables
  yhat_deb<-yhat_deb[, c("time", "CO2", "CFC14")]
  #Long format
  Yhat_deb<-melt(yhat_deb, id.vars=c("time"))
  #Observations
  Yhat_deb$obs<-c(as.numeric(d$CO214), as.numeric(d$Cmic14))
  Gfit_deb<-Yhat_deb %>% group_by(variable) %>% summarise(SSres=sum(((obs-value)^2), na.rm = T),
                                                  SStot=sum(((obs-mean(obs, na.rm = T))^2), na.rm = T),
                                                  ll=-sum(((obs-value)^2), na.rm = T)/2/(sd(obs, na.rm = T)^2))
  Gfit_deb$R2<-with(Gfit_deb, 1-SSres/SStot)
  Gfit_deb$N<-length(pdeb)
  Gfit_deb$AIC<-with(Gfit_deb, 2*N-2*ll)
  Gfit_deb$Model<-"DEB"
  
  #Fine temporal scale fo r graphs
  yhat_deb_fine<-as.data.frame(ode(y=c(G=d$Sinit[1], e=0, w=w_i, CO2=0),
                                   func = DEBmodel, parms=pdeb,
                                   times = seq(0, 55, by=0.1)*24))
  Yhat_deb_fine<-melt(yhat_deb_fine, id.vars=c("time"))
  Yhat_deb_fine$Model<-"DEB"
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Statistics~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #For each variable
  #Ftest
  ##Monod vs DEB
  nt<-as.numeric(as.data.frame(subset(Yhat_deb, !is.na(obs)) %>% group_by(variable) %>% summarize(length(obs)))[,2])
  
  stat1ll<-(Gfit_monod$SSres - Gfit_deb$SSres)*(nt - length(pdeb))/Gfit_deb$SSres/(length(pdeb) - length(pmonod))
  stat1p<-pf(q=stat1ll, df1=(length(pdeb) - length(pmonod)), df2=(nt - length(pdeb)), lower.tail=F)
  ##Mend vs DEB
  stat2ll<-(Gfit_mend$SSres - Gfit_deb$SSres)*(nt - length(pdeb))/Gfit_deb$SSres/(length(pdeb) - length(pmend))
  stat2p<-pf(q=stat2ll, df1=(length(pdeb) - length(pmend)), df2=(nt - length(pdeb)), lower.tail=F)
  ##Pirt vs DEB
  stat3ll<-(Gfit_pirt$SSres - Gfit_deb$SSres)*(nt - length(pdeb))/Gfit_deb$SSres/(length(pdeb) - length(ppirt))
  stat3p<-pf(q=stat3ll, df1=(length(pdeb) - length(ppirt)), df2=(nt - length(pdeb)), lower.tail=F)
  
  stat_eachF<-data.frame(Variable=c("CO2", "CFC"),
                         MvsDEB_F=stat1ll, MvsDEB_p=stat1p,
                         MdvsDEB_F=stat2ll, MdvsDEB_p=stat2p,
                         PvsDEB_F=stat3ll, PvsDEB_p=stat3p)
  
  #LRtest
  ##Monod vs DEB
  stat1ll<--2*(Gfit_monod$ll-Gfit_deb$ll)
  stat1p<-round(pchisq(-2*(Gfit_monod$ll-Gfit_deb$ll), df=(length(pdeb)-length(pmonod)),
                       lower.tail = F), 3)
  ##Mend vs DEB
  stat2ll<--2*(Gfit_mend$ll-Gfit_deb$ll)
  stat2p<-round(pchisq(-2*(Gfit_mend$ll-Gfit_deb$ll), df=(length(pdeb)-length(pmend)),
                       lower.tail = F), 3)
  
  ##Pirt vs DEB
  stat3ll<--2*(Gfit_pirt$ll-Gfit_deb$ll)
  stat3p<-round(pchisq(-2*(Gfit_pirt$ll-Gfit_deb$ll), df=(length(pdeb)-length(ppirt)),
                       lower.tail = F), 3)
  
  stat_eachLR<-data.frame(Variable=c("CO2", "CFC"),
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
  M1pmonod = length(pmonod)
  M1pmend = length(pmend)
  M1ppirt = length(ppirt)
  
  ##DEB model
  ###residual sum of squares
  M2ss = sum((Yhat_deb$obs-Yhat_deb$value)^2, na.rm = T)
  ###number of parameters 
  M2p = length(pdeb)
  
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
                          df=(length(pdeb)-length(pmonod)),
                          lower.tail = F), 3)
  ##Mend vs DEB
  stat2llall<--2*(with(Yhat_mend, -sum(((obs-value)^2), na.rm = T)/2/(sd(obs, na.rm = T)^2))-
                    with(Yhat_deb, -sum(((obs-value)^2), na.rm = T)/2/(sd(obs, na.rm = T)^2)))
  stat2pall<-round(pchisq(-2*(with(Yhat_mend, -sum(((obs-value)^2), na.rm = T)/2/(sd(obs, na.rm = T)^2))-
                                with(Yhat_deb, -sum(((obs-value)^2), na.rm = T)/2/(sd(obs, na.rm = T)^2))), 
                          df=(length(pdeb)-length(pmend)),
                          lower.tail = F), 3)
  
  ##Pirt vs DEB
  stat3llall<--2*(with(Yhat_pirt, -sum(((obs-value)^2), na.rm = T)/2/(sd(obs, na.rm = T)^2))-
                    with(Yhat_deb, -sum(((obs-value)^2), na.rm = T)/2/(sd(obs, na.rm = T)^2)))
  stat3pall<-round(pchisq(-2*(with(Yhat_pirt, -sum(((obs-value)^2), na.rm = T)/2/(sd(obs, na.rm = T)^2))-
                                with(Yhat_deb, -sum(((obs-value)^2), na.rm = T)/2/(sd(obs, na.rm = T)^2))), 
                          df=(length(pdeb)-length(ppirt)),
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
monod_par<-as.numeric(read.csv("parameters/b11_monodpars.csv", header = F))
##Mend
mend_par<-as.numeric(read.csv("parameters/b11_mendpars.csv", header = F))
##Pirt
pirt_par<-as.numeric(read.csv("parameters/b11_pirtpars.csv", header = F))
##DEB
deb_par<-as.numeric(read.csv("parameters/b11_debpars.csv", header = F))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Models evaluation~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
B11_fit<-good_all(monod_par, mend_par, pirt_par, deb_par)
B11_fit$Gfit
B11_fit$stat_eachF
B11_fit$stat_eachLR
B11_fit$stat_allF
B11_fit$stat_allLR
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Figure~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
Rawdata<-subset(B11_fit$Yhat, variable=="CO2" | variable=="CFC14")
Fits<-subset(B11_fit$Yhat_fine, variable=="CO2" |  variable=="CFC14")
ggplot(subset(Rawdata), aes(time, obs))+
  geom_point(cex=6, pch=21, fill="grey")+
  geom_line(data=subset(Fits), aes(time, value, color=Model), lwd=1.2)+theme_min+
  facet_wrap(~variable, scales="free", labeller = label_parsed) + 
  ylab(expression(paste("Carbon pool (", mu, "mol ", g(DW)^{-1}, ")"))) +
  xlab("Time (days)")
  
