##################################################################
### COMPARE IC's tot & sampling on 50, 100, 150, 200 years samples
##################################################################
rm(list=ls())
dev.off()
source("C://Users/mathieu.lucas/Desktop/GitMat/PropagMaxAn/Codes/dirs.R")
source(paste0(dir.codes,"module_BaRatin.r"))
source(paste0(dir.codes,"Fun_SPD.r"))
### dir for plots
dir.plots.gev = paste0(dir.plots,"/GeV_ShapeFlatPrior")
dir.create(dir.plots.gev,showWarnings = F)
plotonly = F

######################################
############ DATA LOADING ############
######################################
  SpagsHydro = read.table(paste0(dir.res,"spags_uTot_Amax.txt"))
  QuantHydro = read.table(paste0(dir.res,"Quantiles_Amax.txt"),header=T)
######################################
######### MODEL PREPARATION ##########
######################################
  #### 4 cases of sample size
  nyears = c(50,100,150, length(QuantHydro$an))
  nspag = dim(SpagsHydro)[2] ; nsim = 5000
  #### SToods parameters, common GEV priors to all cases
  Pos <- parameter(name='Pos',init = 6000) 
  Ech =  parameter(name='Ech',init = 1000) 
  Form = parameter(name='Form',init = 0)#,priorDist='Gaussian',priorPar=c(0,0.2))
  Param = list(Pos,Ech,Form)
  #### Quantiles extracted up to Q1000 & return period associated
  prob = seq(0.01,0.999,0.001) ; Pr = 1/(1-prob)
  #### Initializing results list
  QuantAll = list()
  FreqAll = list()
  FormHyd = list()
  FormTot = list()
  #### Initializing plots list
  GgGeV = list()
  GgU = list()
  #### Creating results folder
  dir.res.stoods = paste0(dir.res,"/GeV_Amax")
  dir.create(dir.res.stoods,showWarnings = F)
######################################
######### MODEL PREPARATION ##########
######################################
  if(plotonly==F){
  for(case in 1: length(nyears)){
    message("***************************************************************")
    message(paste0("Case = ",case))
    message("***************************************************************")
    #### Creating case folder
    dir.case = paste0(dir.res.stoods,"/",nyears[case],"y")
    dir.create(dir.case,showWarnings = F)
    #### Maxpost for case years
    Mp.case = tail(QuantHydro$mp,nyears[case])
    Q.case = tail(QuantHydro,nyears[case])
    #### Spags for case years
    Spag.case = tail(SpagsHydro,nyears[case])
    #### Initializing DF for MCMC results
    MegaSpagTot = data.frame()
    MegaSpagHyd = data.frame()
    ######### ESTIMATING GEV PARAMETERS FOR EACH HYDRO SPAG ##########
    for(spag in 1:nspag){
      print(paste0("Spag = ",spag))
      #### GEV model definition
      dat <- dataset(Y=data.frame(Q=Spag.case[,spag]))
      mod <- model(dataset=dat, parentDist ='GEV', par=Param)
      #### Create 1 sub-folder per hydro spag 
      dir.spag = file.path(dir.case,"Spag")#paste0("Spag",i)) 
      dir.create(dir.spag,showWarnings = F)
      STooDs(mod,workspace = dir.spag, mcmcOptions = mcmc(Nsim = nsim, Nslim = 2))
      #### Read MCMC results for each spag
      MCMCres = readMCMC(file = file.path(dir.spag,"mcmc.txt"), burnFactor = 0.5, slimFactor = 2)
      #### Streamflow uncertainty : Maxpost GeV for the i'th hydro spag
      MegaSpagHyd = rbind(MegaSpagHyd, MCMCres[which.max(MCMCres$post),(1:3)])
      #### Total uncertainty : Nspag GeV for the i'th hydro spag
      MegaSpagTot = rbind(MegaSpagTot, MCMCres[,(1:3)])
    }
    #### Create matrix to stock N spags quantiles
    MegaGevHyd = matrix(nrow = nrow(MegaSpagHyd),ncol = length(prob))
    MegaGevTot = matrix(nrow = nrow(MegaSpagTot),ncol = length(prob))
    #### Compute the quantiles for all the GeV realisations spags
    #### WARNING, HERE SHAPE HAS TO BE THE OPPOSITE THE SHAPE GIVEN BY STOODS (BY CONVENTION)
    for ( i in 1 : nrow(MegaSpagHyd)) {MegaGevHyd[i,] = qgev(p = prob,loc = MegaSpagHyd$Pos[i],
        scale = MegaSpagHyd$Ech[i],shape = -1*MegaSpagHyd$Form[i]) }
    for ( i in 1 : nrow(MegaSpagTot)) {MegaGevTot[i,] = qgev(p = prob,loc = MegaSpagTot$Pos[i],
        scale = MegaSpagTot$Ech[i], shape = -1*MegaSpagTot$Form[i]) }
    #### Stock the form parameters
    FormHyd[[case]] = MegaSpagHyd$Form
    FormTot[[case]] = MegaSpagTot$Form
    #### Compute the true maxpost quantiles : maxpost of hydro sample x maxpost of GeV estim
    dat <- dataset(Y = data.frame(Q = Mp.case))
    mod <- model(dataset=dat, parentDist ='GEV', par=Param)
    #### Create MP sub-folder & run Stoods
    dir.mp = file.path(dir.case,"/Mp") ; dir.create(dir.mp,showWarnings = F)
    STooDs(mod,workspace = dir.mp, mcmcOptions = mcmc(Nsim = nsim, Nslim = 2))
    MCMCres = readMCMC(file = file.path(dir.spag,"mcmc.txt"), burnFactor = 0.5, slimFactor = 1)
    Mp.GeV = MCMCres[which.max(MCMCres$post),(1:3)]
    Mp.Quant = qgev(p = prob,loc = Mp.GeV$Pos, scale = Mp.GeV$Ech,shape = -1*Mp.GeV$Form)
    quanthyd = apply(MegaGevHyd,MARGIN = 2,FUN = quantile, probs = c(0.025,0.975))
    quanttot = apply(MegaGevTot,MARGIN = 2,FUN = quantile, probs = c(0.025,0.975))
    ## floods PR
    Freq = Q.case[order(Q.case$mp),]
    Freq$Fr = (seq(1:length(Q.case$an))-0.5)/length(Q.case$an)
    Freq$Pr = 1/(1-Freq$Fr)
    
    Quants = data.frame(Pr=Pr, Mp=Mp.Quant, Qhyd_2=quanthyd[1,] ,Qhyd_9=quanthyd[2,],
                        Qtot_2=quanttot[1,],Qtot_9=quanttot[2,])
    
    QuantBoth = ggplot()+
      geom_ribbon(data=Quants,aes(x=Pr, ymin = Qtot_2, ymax=Qtot_9,
                                  fill="2-Total uncertainty : \n Streamflow + sampling"),alpha=0.8)+
      geom_ribbon(data=Quants,aes(x=Pr, ymin = Qhyd_2, ymax=Qhyd_9, 
                                  fill="1-Streamflow uncertainty"),alpha=1)+
      geom_line(data=Quants,aes(x=Pr,y=Mp,col="Maxpost"))+
      geom_point(data = Freq, aes(x=Pr, y = mp))+
      geom_errorbar(data=Freq, aes(x=Pr,ymin = tot2.5, ymax = tot97.5))+
      scale_x_continuous(trans="log10")+
      xlab("Return period [years]")+
      ylab(expression(paste("Discharge [",m^3,".",s^-1,"]",sep="")))+
      scale_fill_manual(name = element_blank(),
                        values = c("#fec44f","#67a9cf"))+
      scale_color_manual(values = "royalblue")+
      theme_bw(base_size=15)+
      labs(title = paste0(head(Q.case$an,1)," - ",tail(Q.case$an,1)))+
      coord_cartesian(xlim=c(1,1000))+
      theme(legend.title=element_blank(),
            plot.title = element_text(hjust = 0.01, vjust = -7),
            legend.position = c(0.8,0.2),
            axis.text=element_text(size=20),
            axis.title=element_text(size=20),
            legend.text=element_text(size=15),
            legend.key.size=unit(1, "cm"))
    
    # QuantBoth
    # ggsave(paste0(dir.plots.gev,"/GeV_",nyears[case],"years.pdf"),
    #        device = "pdf", width = 10,height = 7,units = "in")
    
    QuantICs2 = ggplot()+
      geom_ribbon(data=Quants,aes(x=Pr, ymin = 0, ymax = 100, fill = "2-Sampling U") )+
      geom_ribbon(data=Quants,aes(x=Pr, ymin = 0, 
                                  ymax = ((Qhyd_9-Qhyd_2)/(Qtot_9-Qtot_2))*100,
                                  fill = "1-Streamflow U" ))+
      scale_x_continuous(trans="log10")+
      xlab("Return period [years]")+
      ylab("Percentage of total uncertainty [%]")+
      scale_fill_manual(name = element_blank(),
                        values = c("#fec44f","#67a9cf"))+
      theme_bw(base_size=15)+
      coord_cartesian(xlim=c(1,1000))
    # QuantICs2
    
    # ggarrange(QuantBoth + theme(axis.title.x=element_blank()),
    #           QuantICs2, nrow = 2, common.legend = T, align = "v", legend = "right")
    # 
    # ggsave(paste0(dir.plots.gev,"/GeV_IC_",nyears[case],"years.pdf"),
    #        device = "pdf", width = 10,height = 8,units = "in")
    
    QuantAll[[case]] = as.matrix(Quants)
    FreqAll[[case]] = Freq
    
    # GgU[[case]]  = QuantICs2
    # if(case==1 | case==2){GgU[[case]]  = QuantICs2 + theme(axis.title.x=element_blank())} 
    # if(case==2 | case==4) {GgU[[case]]  = QuantICs2 + theme(axis.title.y=element_blank())}
    # if(case==1){ GgU[[case]]  = QuantICs2 + theme(legend.position = c(0.8, 0.8),
    #                                           legend.background = element_rect(fill="lightgrey", 
    #                                           size=0.5, linetype="solid",colour = "black"))}
    # 
    write.table(Quants,file = paste0(dir.case,"/Quants",case,".txt"),col.names = F,row.names = F)
    write.table(Freq,file = paste0(dir.case,"/Freqs",case,".txt"),col.names = F,row.names = F)
    write.table(FormTot[[case]],file = paste0(dir.case,"/FormTot",case,".txt"),
                row.names = F)
    write.table(FormHyd[[case]],file = paste0(dir.case,"/FormHyd",case,".txt"),
                row.names = F)
  
  }
  }
    
    
  # ggarrange(GgU[[1]] ,GgU[[2]] ,GgU[[4]] ,GgU[[4]] ,ncol=2,nrow=2,align = "hv")
  # ggsave("Ukplot4cases.pdf",path = dir.plots.gev,device = "pdf", width = 14,height = 8,units = "in")
    
  # write.table(QuantAll,file = paste0(dir.res,"/QuantsGEVAll.txt"),col.names = F,row.names = F)
  # write.table(FreqAll,file = paste0(dir.res,"/FreqsGEVAll.txt"),col.names = F,row.names = F)
  
  
######################################
######### PLOTS QUANTILES ############
######################################
  
  
  QuantAll = list()
  QuantGG = data.frame()
  Guk = list()
  FreqAll = list()
  for(case in 1:4){
    dir.case = paste0(dir.res.stoods,"/",nyears[case],"y")
    QuantAll[[case]] = read.table(paste0(dir.case,"/Quants",case,".txt"),header = F)
    names(QuantAll[[case]]) = c("Pr","Mp","Qhyd_2","Qhyd_9","Qtot_2","Qtot_9")
    QuantGG = rbind(QuantGG, data.frame(QuantAll[[case]],case = nyears[case]))
    FreqAll[[case]] = read.table(paste0(dir.case,"/Freqs",case,".txt"),header = F)
    names(FreqAll[[case]]) = c(names(QuantHydro),"Fr","Pr")
    
    FormTot[[case]] = read.table(paste0(dir.case,"/FormTot",case,".txt"), header = T)
    FormHyd[[case]] = read.table(paste0(dir.case,"/FormHyd",case,".txt"), header = T)
    
    
    QuantBoth = ggplot()+
      geom_ribbon(data=QuantAll[[case]],aes(x=Pr, ymin = Qtot_2, ymax=Qtot_9,
                                  fill="2-Total uncertainty : \n Streamflow + sampling"),alpha=0.8)+
      geom_ribbon(data=QuantAll[[case]],aes(x=Pr, ymin = Qhyd_2, ymax=Qhyd_9, 
                                  fill="1-Streamflow uncertainty"),alpha=1)+
      geom_line(data=QuantAll[[case]],aes(x=Pr,y=Mp,col="Maxpost"))+
      geom_point(data = FreqAll[[case]], aes(x=Pr, y = mp))+
      geom_errorbar(data=FreqAll[[case]], aes(x=Pr,ymin = tot2.5, ymax = tot97.5))+
      scale_x_continuous(trans="log10")+
      xlab("Return period [years]")+
      ylab(expression(paste("Discharge [",m^3,".",s^-1,"]",sep="")))+
      scale_fill_manual(name = element_blank(),
                        values = c("#fec44f","#67a9cf"))+
      scale_color_manual(values = "royalblue")+
      theme_bw(base_size=15)+
      labs(title = paste0(min(FreqAll[[case]]$an)," - ",max(FreqAll[[case]]$an)))+
      coord_cartesian(xlim=c(1,1000))+
      theme(legend.title=element_blank(),
            plot.title = element_text(hjust = 0.01, vjust = -7),
            legend.position = c(0.8,0.2),
            axis.text=element_text(size=20),
            axis.title=element_text(size=20),
            legend.text=element_text(size=15),
            legend.key.size=unit(1, "cm"))
    
    ggsave(paste0(dir.plots.gev,"/GeV_",nyears[case],"years.pdf"),
           device = "pdf", width = 18,height = 8,units = "in")
    
    QuantICs2 = ggplot()+
      geom_ribbon(data=data.frame(QuantAll[[case]]),aes(x=Pr, ymin = 0, ymax = 100,
                                                        fill = "2-Sampling U") )+
      geom_ribbon(data=data.frame(QuantAll[[case]]),aes(x=Pr, ymin = 0, 
                                  ymax = ((Qhyd_9-Qhyd_2)/(Qtot_9-Qtot_2))*100,
                                  fill = "1-Streamflow U" ))+
      scale_x_continuous(trans="log10")+
      xlab("Return period [years]")+
      ylab("Perc. of tot. uncertainty [%]")+
      scale_fill_manual(name = element_blank(),
                        values = c("#fec44f","#67a9cf"))+
      theme_bw(base_size=10)+
      coord_cartesian(xlim=c(1,1000))+
      theme(legend.position = "null",
            axis.text=element_text(size=25),
            axis.title=element_text(size=25),
            legend.text=element_text(size=18),
            plot.title = element_text(hjust = 0.08, vjust = -8,size = 25))+
      labs(title=paste0(nyears[case], " years"))
    
    Guk[[case]] = QuantICs2
    
    if(case==1 | case==2){Guk[[case]]  = QuantICs2 + theme(axis.title.x=element_blank())} 
    if(case==2 | case==4) {Guk[[case]]  = QuantICs2 + theme(axis.title.y=element_blank())}
    if(case==2){ Guk[[case]]  = QuantICs2 + theme(legend.position = c(0.72, 0.75),
                                                  legend.background = element_rect(fill="lightgrey", 
                                                  size=0.5, linetype="solid",colour = "black"),
                                                  axis.title.y=element_blank(),
                                                  axis.text.y = element_blank())}
  }
  
  ggarrange(Guk[[1]]+theme(axis.title.x=element_blank(), axis.text.x = element_blank()),
            Guk[[2]]+theme(axis.title.x=element_blank(), axis.text.x = element_blank()),
            Guk[[3]],
            Guk[[4]]+theme(axis.text.y = element_blank()) ,
            ncol=2,nrow=2,align = "hv")
  ggsave("Ukplot4cases.pdf",path = dir.plots.gev,device = "pdf", width = 16,height = 9,units = "in")

  QuantGG$Pr = round(QuantGG$Pr,2)
  QuantGG$Pr=factor(QuantGG$Pr)
  QuantGG$case = factor(QuantGG$case)
  axis.text.size = 20
  titles.size = 25
  ymax = round(max(QuantGG$Qtot_9))+100
  
  palgrey = RColorBrewer::brewer.pal(5,"Greys")[2:5]
  
  BpQuant10 = ggplot(QuantGG[which(QuantGG$Pr==10),]) +
    geom_bar( aes(x= case, y=Mp, fill=case), stat="identity", alpha=0.8, width = 0.8) +
    geom_errorbar( aes(x=case, ymin=Qtot_2, ymax=Qtot_9), width=0.2,
                   colour="black"#"orange"
                   , alpha=0.9, size=1)+
    ylim(c(0,ymax))+
    theme_bw(base_size=15)+
    theme(legend.position = "none",
          axis.title.x=element_blank(),
          axis.text=element_text(size=axis.text.size),
          axis.title=element_text(size=titles.size))+
    scale_fill_manual(name = element_blank(),
                      values = palgrey)+
    ylab(label=expression(paste("Discharge [",m^3,".",s^-1,"]",sep="")))
  
  BpQuant100 = ggplot(QuantGG[which(QuantGG$Pr==100),]) +
    geom_bar( aes(x= case, y=Mp,fill=case), stat="identity",  alpha=0.7, width = 0.8) +
    geom_errorbar( aes(x=case, ymin=Qtot_2, ymax=Qtot_9), width=0.2,
                   colour="black", alpha=0.9, size=1)+
    ylim(c(0,ymax))+
    theme_bw(base_size=15)+
    theme(legend.position = "none",
          axis.text.y=element_blank(), #remove x axis labels
          axis.ticks.y=element_blank(),
          axis.title.y=element_blank(),
          #axis.title.x=element_blank(),
          axis.title=element_text(size=titles.size),
          axis.text=element_text(size=axis.text.size))+
    xlab("Sample size [years]")+
    scale_fill_manual(name = element_blank(),
                      values = palgrey)#c("#bdd7e7","#6baed6","#3182bd","#08519c"))
  
  BpQuant1000 = ggplot(QuantGG[which(QuantGG$Pr==1000),]) +
    geom_bar( aes(x= case, y=Mp,fill=case), stat="identity",  alpha=0.7,  width = 0.8) +
    geom_errorbar( aes(x=case, ymin=Qtot_2, ymax=Qtot_9), width=0.2, colour="black",
                   alpha=0.9, size=1)+
    ylim(c(0,ymax))+
    theme_bw(base_size=15)+
    theme(legend.position = "none",
          axis.text.y=element_blank(), #remove y axis labels
          axis.ticks.y=element_blank(),
          axis.title.y=element_blank(),
          axis.title.x=element_blank(),
          axis.text=element_text(size=axis.text.size))+
    scale_fill_manual(name = element_blank(),
                      values = palgrey)#c("#bdd7e7","#6baed6","#3182bd","#08519c"))
  
  #write plot
  Ggall = align_plots(BpQuant10, BpQuant100, BpQuant1000,align="hv")
  ggarrange(Ggall[[1]],Ggall[[2]],Ggall[[3]],
            labels = c("Q10","Q100","Q1000"),label.y = 0.98,
            label.x = 0.62,font.label = list(size=25),
            ncol = 3, nrow = 1)
  ggsave(paste0(dir.plots.gev,"/BarplotQuantiles4cases.pdf"),device = "pdf",
         width = 14,height = 8,units = "in")
  
  
  #### PLOT FORM PARAMS
  require(ggridges)
 
  # Ggforms=list()
  # 
  #  for(i in 1:4){
  # 
  #   Ggforms[[i]] = ggplot(data = data.frame( tot = sample(FormTot[[i]][,1],length(FormHyd[[i]][,1])),
  #                                            hyd = FormHyd[[i]][,1], case = i ) )+
  #     geom_density(aes(x=tot, fill = "2-Total uncertainty"),alpha = 0.8)+
  #     geom_density(aes(x=hyd, fill = "1-Streamflow uncertainty"),alpha = 0.8)+
  #     scale_fill_manual(values = c("#fec44f","#67a9cf"))+
  #     xlab("")+
  #     theme_light()+
  #     theme(legend.position = c(0.25,0.8), #legend.title = element_blank(),
  #           axis.text=element_text(size=20),
  #           axis.title=element_text(size=25),
  #           legend.text=element_text(size=20),
  #           legend.title=element_text(size=15))+
  #     ggtitle(paste0(nyears[[i]]," years"))+
  #     coord_cartesian(xlim = c(-0.3,0.3),ylim = c(0,40))
  #   
  # }
  # ggarrange(Ggforms[[1]]+theme(legend.position="none"),
  #           Ggforms[[2]]+theme(axis.title.y = element_blank()),
  #           Ggforms[[3]]+theme(legend.position="none"),
  #           Ggforms[[4]]+theme(axis.title.y = element_blank(),legend.position="none"),
  #             align = "hv", ncol = 2, nrow=2)
  # 
  # ggsave(path = dir.plots.gev, filename = "FormParam_4cases.pdf",
  #        width = 10, height = 8, units = "in")
  # 
  #### ridges try
  
  Hyd = data.frame(t(matrix(unlist(FormHyd), nrow=length(FormHyd), byrow=TRUE)))
  names(Hyd) = nyears
  Hyd = melt(Hyd)
  Hyd$type = "1-Streamflow uncertainty"
  
  Tot = data.frame(t(matrix(unlist(FormTot), nrow=length(FormTot), byrow=TRUE)))
  names(Tot) = nyears
  # Tot = Tot[sample(1:nrow(Tot),500,F),]
  Tot = melt(Tot)
  Tot$type = "2-Total uncertainty"
  
  Forms = rbind(Hyd,Tot)
  Forms$variable = as.factor(Forms$variable)
  
  ggplot(Forms, aes(y = variable,
                    x = value,
                    fill = type)) +
    scale_fill_manual(values = c("#fec44f","#67a9cf"))+
    geom_density_ridges(alpha=0.7)+
    scale_x_continuous(expand = c(0, 0))+
    scale_y_discrete(expand = c(0, 0)) +
    coord_cartesian(xlim = c(-0.25,0.41))+
    ylab("Sample size [years]")+
    xlab("Shape parameter []")+
    theme_light()+
    theme(legend.position = c(0.22,0.93), legend.title = element_blank(),
          #axis.title.x = element_blank(),
          axis.text=element_text(size=20),
          axis.title=element_text(size=25),
          legend.text=element_text(size=20))

  ggsave(path = dir.plots.gev, filename = "Form_4cases_ridgeline.pdf",
         width = 10, height = 7, units = "in")
    
beepr::beep()
  
  