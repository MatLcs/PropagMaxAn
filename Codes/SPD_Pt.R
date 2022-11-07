rm(list=ls())
dev.off()
source("C://Users/mathieu.lucas/Desktop/GitMat/PropagMaxAn/Codes/dirs.R")
source(paste0(dir.codes,"module_BaRatin.r"))
source(paste0(dir.codes,"Fun_SPD.r"))

start.time = Sys.time()
case = "Pt"
dir.spd = paste0(dir.bam,"SPD_",case)
dir.bar = paste0(dir.bam,"BaR_",case)

######################################
############ DATA LOADING ############
######################################
  Amax = read.table(paste0(dir.res,"Limni/amaxH_Pt.txt"),header = T)
  Noisy = as.matrix(read.table(paste0(dir.res,"Limni/NoisyHmaxAn_Pt.txt")))
  Gau = read.table(paste0(dir.data,"GauPtBcrWperiods.txt"),header=T)
  #Shifts from segmentation : col 2 = greatest flood within post. dates. 
  #As indices from 1st day of limni
  TshiftsPt = as.Date("1816-05-15")+read.table(paste0(dir.data,"shift_timesPt.txt"),header=T)[,2]
  
######################################
########## DATA FORMATTING ###########
######################################
  # add 1840 flood as a timeshift 
  Tshifts = c(as.Date("1840-11-02"),TshiftsPt)
  nperiod = length(Tshifts) + 1
  # identifiying the period of maxAn stages
  Amax$Date = ymd(Amax$Date) ; nyears = nrow(Amax) ; Amax$period = 1
  for(shift in 1:length(Tshifts)){ Amax$period[which(Amax$Date >= Tshifts[shift])] = shift + 1 }
  #no gaugings in period 1 : before 1840
  Gau$Period = as.numeric(Gau$Period)
  Gau$Period = Gau$Period + 1
######################################
######### PRIORS DEFINITION ##########
######################################
  # 1. Define number of controls, periods and MonteCarlo samples
  ncontrol = 2; nsim = 1000 ; RunSPD = F ; nspag = ncol(Noisy)
  # 2. Define which parameters are varying in parameterization (b1,a1,c1,b2,a2,c2)
  isVar=c(T,F,F,T,F,F)
  # 3. Define parameter names and prior distributions for b,a,c (of any control)
  names.bac=c('b','a','c')
  margins.bac=c('Gaussian','LogNormal','Gaussian')
  #---------------------------------------
  # Define priors on "physical" parameters
  #---------------------------------------
  # Priors are defined through Monte Carlo samples
  # 1. Priors for low-flow control (control 1)
  b1<-rnorm(nsim,mean=-4,sd=0.5) # weir activation stage
  c1<-rnorm(nsim,mean=5/3,sd=0.025) # exponent
  Bc1<-rlnorm(nsim,meanlog=log(300),sdlog=as.numeric(Transf_Gauss_lognorm(300,50)[2])) # channel width
  KS1<-rlnorm(nsim,meanlog=log(35),sdlog=as.numeric(Transf_Gauss_lognorm(35,5)[2])) # Strickler coefficient
  S01<-rlnorm(nsim,meanlog=log(1.5e-4),sdlog=as.numeric(Transf_Gauss_lognorm(1.5e-4,9e-5)[2])) # slope
  # 2. Priors for main channel (control 2)
  b2<-rnorm(nsim,mean=2,sd=0.5) # channel bed stage
  c2<-rnorm(nsim,mean=5/3,sd=0.025) # exponent
  Bc2<-rlnorm(nsim,meanlog=log(500),sdlog=as.numeric(Transf_Gauss_lognorm(500,50)[2])) # channel width
  KS2<-rlnorm(nsim,meanlog=log(30),sdlog=as.numeric(Transf_Gauss_lognorm(30,5)[2])) # Strickler coefficient
  S02<-rlnorm(nsim,meanlog=log(2.6e-4),sdlog=as.numeric(Transf_Gauss_lognorm(2.6e-4,9e-5)[2])) # slope
  # 4. Priors for incremental global changes
  d.g<-matrix(rnorm(nsim*(nperiod-1),mean=0,sd=0.3),nrow=nsim,ncol=nperiod-1)
  # 5. Priors for incremental local changes
  d.l<-matrix(rnorm(nsim*(nperiod-1),mean=0,sd=0.3),nrow=nsim,ncol=nperiod-1)
  # 6. Starting point for all parameters
  start=list(
    b1=-4,c1=5/3,Bc1=300,KS1=35,S01=1.5e-4,  # control 1
    b2= 2,c2=5/3,Bc2=500,KS2=30,S02=2.6e-4, # control 2
    d.l=rep(0,nperiod-1),     # incremental local changes d.l
    d.g=rep(0,nperiod-1)      # incremental global changes d.g
  )
  M                 = matrix(0, ncontrol, ncontrol)  #matrix of controls:
  M[1,]             = c(1,0)   # control section (rectangular weir in critical condition)
  M[2,]             = c(1,1)   # main control channel (rectangular wide channel in uniform condition)
  Hmax = 10 # by default, 10 but be careful
  #---------------------------------------
  # Perform Monte-Carlo propagation
  #---------------------------------------
  MC=propagate_Bcr_b1b2reversed(b1,c1,Bc1,KS1,S01,  # control 1: rectangular weir
                                b2,c2,Bc2,KS2,S02, # control 2: lrectangular weir
                                d.g,d.l,           # global and local incremental changes
                                start              # starting vector
  )
  #plot priors verif
  # boxplot(MC$sim,ylim=c(-10,10))
  #---------------------------------------
  # Fit prior distribution on MC samples
  #---------------------------------------
  # define marginal prior distributions and parameter names
  margins=c();names=c();k=0
  for(i in 1:ncontrol){
    for(j in 1:3){
      k=k+1
      if(isVar[k]){
        margins=c(margins,rep(margins.bac[j],nperiod))
        names=c(names,paste(names.bac[j],i,'_',1:nperiod,sep=''))      
      }else{
        margins=c(margins,margins.bac[j])
        names=c(names,paste(names.bac[j],i,sep=''))
      }  
    }
  }
  # fit multivariate prior distribution
  prior=fit(MC$sim,margins,names)
  #Remnant error model:
  #--------------------
  #LINEAR: err = g1+g2*Q   or    CONSTANT: err = g1):
  remnant.err.model = "Linear"              #"Linear" or "Constant"
  g1.prior          = c(0, 1000, 0.1)       #c(min,max, starting point)
  g2.prior          = c(0, 100, 0.1)        #c(min,max, starting point) or c(mean, stdev, starting point) 
  #if model = "Constant" then put FALSE !.
  g1.distr.type     = "Uniform"             #  "Uniform.
  g2.distr.type     = "Uniform"             # "Lognormal" or "Uniform

######################################
############ BARATIN SPD #############
######################################
  message("--------------------------")
  message("--- RUNNING BARATIN SPD---")
  message("--------------------------")
  #dir for SPD 
  dir.create(dir.spd,showWarnings = F) ; setwd(dir.spd)
  #write gau
  write.table(Gau,"Gaugings_data.txt",row.names = F)
  # write config files
  colperiod=4 # columns where period will be stored in BaM data file
  writeConfigFilesSPD(prior = prior,start = MC$start,ncontrol = ncontrol,
                      nperiod = nperiod,colperiod = colperiod,isVar = isVar)
  ##write config BAM
  BaRatin_SPD_config(dir.BaM = dir.bam,
                     dir.SPD.config = dir.spd,
                     dir.spd.short = paste0("SPD_",case,"/"),
                     pred = F,
                     nobs = nrow(Gau),
                     M = M,
                     remnant = remnant.err.model,
                     g1.prior = g1.prior, g2.prior = g2.prior,
                     g1.distr.type = g1.distr.type, g2.distr.type = g2.distr.type, 
                     Ncycles = 500,
                     Hmax = Hmax)
  setwd(dir.bam) 
  if(RunSPD == T){system2(paste0(dir.bam,"BaM.exe"))}
  #dir for results
  dir.res.case = paste0(dir.res,case,"/")
  dir.create(dir.res.case,showWarnings = F)
  dir.create(paste0(dir.res.case,"/SPD/"),showWarnings = F)
  #copy SPD results in dir res
  file.copy(from = file.path(dir.spd,"/Results_MCMC_Cooked.txt"),
            to = file.path(dir.res.case,"/SPD/Results_MCMC_Cooked.txt"),overwrite = T)
  file.copy(from = file.path(dir.spd,"Results_Summary.txt"),
            to = file.path(dir.res.case,"/SPD/Results_Summary.txt"),overwrite = T)
  
  model = ReadModel(ModelFile = paste0(dir.spd,"/Config_Model.txt"))
  #### MCMC plots
  vertical.length = length(names)*4
  mcmc = MCMCplot(doLogPost = T,
                  doPar     = T,
                  doDPar    = T, 
                  MCMCfile  = paste0(dir.res.case,"/SPD/Results_MCMC_Cooked.txt") , 
                  type      = "trace",  #="trace", # "histogram","density","scatterplot"
                  xlab      = '',
                  ylab      = '',
                  ncol      = 1, 
                  prior     = NULL,
                  burn      = 0, 
                  slim      = 1,
                  theme_size= 15)
  # ggsave(mcmc, filename =paste(workspace,"/mcmc_it",iter,".png", sep=""), 
  #        bg = "transparent", width = 12, height =vertical.length, dpi = 100)
  mcmc2 = MCMCplot(doLogPost = T,
                   doPar     = T,
                   doDPar    = T, 
                   MCMCfile  =  paste0(dir.res.case,"/SPD/Results_MCMC_Cooked.txt")  , 
                   type      = 'density', 
                   prior     = model$par,
                   xlab      = '',
                   ylab      = '',
                   ncol      = 1, 
                   burn      = 0, 
                   slim      = 1,
                   theme_size= 15)  
  # ggsave(mcmc2, filename =paste(workspace,"/mcmc2_it",iter,".png", sep=""), 
  # bg = "transparent", width = 12, height =vertical.length, dpi = 100)
  
  dir.create(paste0(dir.plots,"/",case,"/"),showWarnings = F)
  dir.plot.case = paste0(dir.plots,"/",case,"/")
  dir.create(paste0(dir.plot.case,"SPD/"),showWarnings = F)
  
  pdf(paste0(dir.plot.case,"/SPD/mcmc.pdf"), 17, vertical.length, useDingbats=F)
    print(plot_grid(mcmc, mcmc2, nrow=1, ncol = 2, 
                    labels = c("Trace plots", "Density plots"),label_size = 20) )
  dev.off()
  rm(mcmc);  rm(mcmc2); rm(model)

######################################
############ PLOT RC's ###############
######################################
  message("--------------------------")
  message("PLOTTING SPD RC's RESULTS")
  message("--------------------------")
  # read SPD results
  data.MCMC.spd = as.matrix(read.table(paste0(dir.spd,"/Results_MCMC_Cooked.txt"),
                                   header=TRUE,dec=".", sep=""))
  data.MP.spd = tail(data.frame(read.table(paste0(dir.spd,"/Results_Summary.txt"),
                                                 row.names = 1, dec=".",sep="")),1)
  nsample = length(data.MCMC.spd[,1])
  min.grid = min(data.MCMC.spd[,1]) 
  Hinf = -3;  Hsup = 13
  Hgrid  = seq(Hinf, Hsup, 0.01)
  ngrid = length(Hgrid)
  #colors
  RCcol = c("#a50026","#d73027","#f46d43","#fdae61","#fee090","#e0f3f8","#abd9e9","#74add1",
            "#4575b4","#313695")
  #init
  RCs = list()
  DatasRC = list()
  for(periode in 1:nperiod){
    GauPer = Gau[which(Gau$Period==periode),]
    if (periode != nperiod){
      cook=data.MCMC.spd[,c(paste0("b1_",periode),"a1","c1", #bac1
                        paste0("b2_",periode),"a2","c2", #bac2
                        #paste0("b3_",periode),"a3","c3", #bac3
                        "Y1_gamma1","Y1_gamma2","LogPost",#y & logpost
                        paste0("k1_",periode),paste0("k2_",periode),#paste0("k3_",periode),#ks
                        paste0("ktype1_",periode),paste0("ktype2_",periode)#,paste0("ktype3_",periode)
      )]
      cookmp  = data.MP.spd[,c(paste0("b1_",periode),"a1","c1", #bac1
                                     paste0("b2_",periode),"a2","c2", #bac2
                                     #paste0("b3_",periode),"a3","c3", #bac3
                                     "Y1_gamma1","Y1_gamma2",#y & logpost
                                     paste0("k1_",periode),paste0("k2_",periode),
                                     #paste0("k3_",periode),#ks
                                     paste0("ktype1_",periode),paste0("ktype2_",periode)
                                     #paste0("ktype3_",periode) #ktypes
      )]     
    } else {
      cook=data.MCMC.spd[,c(paste0("b1_",periode),"a1","c1", #bac1
                        paste0("b2_",periode),"a2","c2", #bac2
                        #paste0("b3_",periode),"a3","c3", #bac3
                        "Y1_gamma1","Y1_gamma2","LogPost",#y & logpost
                        paste0("k1_",periode-1),paste0("k2_",periode-1),#paste0("k3_",periode),#ks
                        paste0("ktype1_",periode-1),paste0("ktype2_",periode-1)
                        #,paste0("ktype3_",periode) #ktypes
      )]
      cookmp  = data.MP.spd[,c(paste0("b1_",periode),"a1","c1", #bac1
                                     paste0("b2_",periode),"a2","c2", #bac2
                                     #paste0("b3_",periode),"a3","c3", #bac3
                                     "Y1_gamma1","Y1_gamma2",#y & logpost
                                     paste0("k1_",periode-1),paste0("k2_",periode-1),
                                     #paste0("k3_",periode),#ks
                                     paste0("ktype1_",periode-1),paste0("ktype2_",periode-1)
                                     #,paste0("ktype3_",periode) #ktypes
      )]     
    }
    
    MCMC.save    =  matrix(NA, nrow = nsample, ncol = ncol(cook))  
    MaxPost.save =  matrix(NA, nrow = 1,  ncol = ncol(cook)-1)  
    MCMC.save = cbind(cook, rep(1, nsample))
    MaxPost.save = as.numeric(c(cookmp,1))
    RC.Post    =  apply(MCMC.save,MARGIN = 1,  RC_controls,h = Hgrid,  M = M, ncontrols = ncontrol)
    RC.MaxPost =  RC_controls_mp(theta = as.numeric(MaxPost.save),h=Hgrid,M=M,ncontrols = ncontrol)
    data.tmp            = apply(RC.Post, MARGIN=1, quantile, probs=c(0.025,0.975), na.rm=TRUE)
    data.tmp            = apply(data.tmp, MARGIN=c(1,2), function(x){ifelse(x<0,0,x)})
    List.RC.quants      = data.frame(cbind(Hgrid, t(data.tmp),  RC.MaxPost))
    rm(data.tmp)
    colnames(List.RC.quants) = c("h", "inf", "sup", "maxpost")
    # data.RC = List.RC.quants
    # DatasRC[[periode]] = data.RC
    DatasRC[[periode]] = List.RC.quants
  }
  ## writing data for plot
  DatPlot = data.frame(List.RC.quants[0,])
  for(per in 1:nperiod){
    DatPlot = rbind(DatPlot,cbind(round(DatasRC[[per]],2),rep(per,length(DatasRC[[per]]$h))))  }
  names(DatPlot) = c('h','inf','sup','maxpost','period')
  write.table(DatPlot,paste0(dir.res.case,"/SPD/DataRCs.txt"))
  ### PLOTS RC
  DataRC = DatPlot
  rm(DatPlot)
  DataRC = DataRC[which(is.na(DataRC$maxpost)==F  & DataRC$maxpost > 0),]
  DataRC$period=as.factor(DataRC$period)
  GauPlot = Gau
  GauPlot$Period=as.factor(GauPlot$Period)
  
  plot.RC=ggplot(DataRC)+
    #RC's
    geom_smooth(aes(x=h,y=maxpost,ymax=sup,ymin=inf,colour=period),size=1,stat='identity',alpha=0.2)+
    geom_path(aes(x=h,y=maxpost,colour=period),size=1.5)+
    ### Gaugings
    geom_linerange(aes(x=h,ymax=Q+2*uQ,ymin=Q-2*uQ,colour=Period),data=GauPlot,size=1)+
    geom_point(aes(x=h,y=Q,colour=Period),data=GauPlot,na.rm=T,shape=16,size=3)+
    ### Labels
    xlab(expression(paste("Stage [m]",sep="")))+
    ylab(expression(paste("Discharge [",m^3,".",s^-1,"]",sep="")))+
    labs(colour = "Period")+
    scale_fill_brewer(type = "div",palette = 5,aesthetics = "color")+
    coord_cartesian(ylim=c(500,15000),xlim=c(-2,10))+
    ### Theme
    theme_bw(base_size=20)+
    theme(axis.text=element_text(size=15),axis.title=element_text(size=20,face="bold")
          ,panel.grid.major=element_line(size=1.2),panel.grid.minor=element_line(size=0.8)
          ,legend.text=element_text(size=15),legend.title=element_text(size=20)
          ,legend.key.size=unit(1.5, "cm"),legend.position="right")
  ## Normal scale
  pdf(paste(dir.plot.case,"SPD/RC.pdf",sep=""),18,10,useDingbats=F)
    print(plot.RC)
  dev.off()
  ## Logarithmic scale
  ylim.wind=c(200,14000);xlim.wind=c(-2,10)
  
  pdf(paste(dir.plot.case,"SPD/RClog.pdf",sep=""),18,10,useDingbats=F)
    print(plot.RC+
          scale_y_log10()+
          coord_cartesian(ylim=ylim.wind,xlim=xlim.wind))
  dev.off()
  
  ## with IC
  IC=ggplot(DataRC)+
    geom_ribbon(aes(x=h,ymax=((sup-maxpost)/maxpost)*100,ymin= ((inf-maxpost)/maxpost)*100,
                    colour=period,fill = period),size=0.6,stat='identity',alpha=0.2)+
    scale_fill_brewer(type = "div",palette = 5,aesthetics = c("colour","fill"))+
    coord_cartesian(ylim=c(-100,100), xlim = xlim.wind)+
    ### Theme
    theme_bw(base_size=20)+
    theme(axis.text=element_text(size=15),axis.title=element_text(size=20,face="plain")
          ,panel.grid.major=element_line(size=1.2),panel.grid.minor=element_line(size=0.8)
          ,legend.text=element_text(size=15),legend.title=element_text(size=20)
          ,legend.key.size=unit(1.5, "cm"),legend.position="right")+
    ylab("Relative uncertainty [%]")+
    xlab("Stage [m]")
  
  ggarrange(plot.RC+
              scale_y_log10()+
              coord_cartesian(ylim=ylim.wind,xlim=xlim.wind)+
              theme(axis.title.x = element_blank()),
            IC,
            ncol = 1,align = "v",common.legend = T, legend = "right",
            heights = c(2,1))
  ggsave(filename = "RClog_ICdown.pdf",path = paste0(dir.plot.case,"SPD/"), width = 12, height = 7)
  
  #### all param boxplot
  prior.bs = MC$sim[,c(1:10,13:22)]
  post.bs = data.frame(tail(data.MCMC.spd[,c(1:10,13:22)],1000))
  names(prior.bs) = names(post.bs) #= names(MC$sim)[c(1:10,13:22)]
  prior= melt(prior.bs)
  post = melt(post.bs)#[,(2:3)]
  #ind for periods
  periodsim= rep(NA,(nsim*nperiod)*2)
  periodsim[1:1000] = 1
  for(per in 2:(nperiod*2)) {  
    if(per <= nperiod) {periodsim[(((per-1)*nsim)+1):(per*nsim)] = per
    } else {periodsim[(((per-1)*nsim)+1):(per*nsim)] = per-nperiod}  }
  #factors
  prior$per = as.factor(periodsim)
  post$per = as.factor(periodsim)
  prior$type = as.factor("Prior")
  post$type = as.factor("Posterior")
  # names(post) = names(prior)
  ParamsBoxplot = rbind(prior,post)
  ### B's
  boxbs=ggplot(data=ParamsBoxplot)+
    geom_boxplot(aes(x=variable,y=value,fill=type),colour="black",width=1)+
    facet_grid(~per,scales="free",labeller = label_parsed)+
    ylab("Offset [m]")+
    theme_bw(base_size=20)+
    labs(fill = "")+
    theme(axis.title.x = element_blank(),
          axis.text=element_text(size=20),axis.title=element_text(size=25,face="bold")
          ,legend.text=element_text(size=5),legend.title=element_text(size=20)
          ,legend.key.size=unit(1.5, "cm"),legend.position="none" #"right"
          ,strip.text.x=element_text(size = 25,face="bold"))+
    scale_x_discrete(breaks=c("b1_1","b2_1",
                              "b1_2","b2_2",
                              "b1_3","b2_3",
                              "b1_4","b2_4",
                              "b1_5","b2_5",
                              "b1_6","b2_6",
                              "b1_7","b2_7",
                              "b1_8","b2_8",
                              "b1_9","b2_9",
                              "b1_10","b2_10"),
                     labels=c(bquote(b[1]^(1)),bquote(b[2]^(1)),
                              bquote(b[1]^(2)),bquote(b[2]^(2)),
                              bquote(b[1]^(3)),bquote(b[2]^(3)),
                              bquote(b[1]^(4)),bquote(b[2]^(4)),
                              bquote(b[1]^(5)),bquote(b[2]^(5)),
                              bquote(b[1]^(6)),bquote(b[2]^(6)),
                              bquote(b[1]^(7)),bquote(b[2]^(7)),
                              bquote(b[1]^(8)),bquote(b[2]^(8)),
                              bquote(b[1]^(9)),bquote(b[2]^(9)),
                              bquote(b[1]^(10)),bquote(b[2]^(10))
                     ))
  
  pdf(paste(dir.plot.case,"SPD/bs.pdf",sep=""),20,9,useDingbats=F)
    print(boxbs)
  dev.off()
  
  ###### A's
  priors.as = melt(MC$sim[,c("a1","a2")])
  priors.as$type = as.factor("Prior")
  post.as = melt(data.frame(tail(data.MCMC.spd[,c("a1","a2")],1000)))#[,(2:3)]
  post.as$type =as.factor("Posterior")
  # names(post.as) = names(priors.as)
  as = rbind(priors.as,post.as)
  boxas=ggplot(data=as)+
    geom_boxplot(aes(x=variable,y=value,fill=type),colour="black")+
    theme_bw(base_size=15)+
    ylab("Value [-]")+
    theme(axis.text=element_text(size=20)#,axis.title=element_blank()
          ,panel.grid.major=element_line(size=1),panel.grid.minor=element_line(size=0.8)
          ,legend.text=element_text(size=20),legend.title=element_blank()
          ,legend.key.size=unit(1.5, "cm"),legend.position="bottom"
          ,axis.title.x=element_blank(), axis.title = element_text(size = 20))+
    scale_x_discrete(breaks=c("a1","a2"),
                     labels=c(bquote(a[1]),bquote(a[2])))
  ######## C's
  priors.cs = melt(MC$sim[,c("c1","c2")])
  priors.cs$type = as.factor("Prior")
  post.cs = melt(data.frame(tail(data.MCMC.spd[,c("c1","c2")],1000)))#[,(2:3)]
  post.cs$type = as.factor("Posterior")
  # names(post.cs) = names(priors.cs)
  cs = rbind(priors.cs,post.cs)
  boxcs=ggplot(data=cs)+
    geom_boxplot(aes(x=variable,y=value,fill=type),colour="black")+
    theme_bw(base_size=15)+
    # ylab("Exponent value [-]")+
    theme(axis.text=element_text(size=25),axis.title=element_blank()
          ,panel.grid.major=element_line(size=1),panel.grid.minor=element_line(size=0.8)
          ,legend.text=element_text(size=15),legend.title=element_blank()
          ,legend.key.size=unit(1.5, "cm"),legend.position="right"
          ,axis.title.x=element_blank())+
    scale_x_discrete(breaks=c("c1","c2"),
                     labels=c(bquote(c[1]),bquote(c[2])))

  pdf(paste(dir.plot.case,"SPD/a&cs.pdf",sep=""),12,7,useDingbats=F,onefile = F)
    ggarrange(boxas,boxcs,legend = 'right',ncol = 2,align = "hv",common.legend = T)
  dev.off()
  
  
######################################
### PROPAG U LIMNI VIA BAM PREDICT ###
######################################
  message("----------------------------------------")
  message("PROPAGATING U LIMNI VIA BAM PREDICTIONS")
  message("---------------------------------------")
  
  dir.create(dir.bar,showWarnings = F)
  write.table(x = Hgrid,
              file = paste0(dir.bar,"/Hgrid.txt"), col.names = F ,row.names = F, sep = "    ",)

  
  #init matrix
  Tot = matrix(nrow = length(Amax$Y),ncol = nspag)
  Param = Tot
  Mp = Tot
  TrueMp = rep(NA,length(Amax$Y))
  
  for(period in 1:nperiod){
    # period = 1
    message(paste0("Processing period ",period," ....."))
    
    bufferAn = which(Amax$period == period)
    ######### SPD results cleaning
    #Folder period n
    dir.create(file.path(dir.res.case,paste0("P",period)),showWarnings = F)
    dir.per = file.path(dir.res.case,paste0("P",period))
    # 
    # #Copy and rename SPD results file to period dir
    # file.copy(from = file.path(dir.spd,list.files(dir.spd,pattern = "Cooked")),
    #           to = file.path(dir.per,list.files(dir.spd,pattern = "Cooked")),overwrite = T)
    # file.rename(file.path(dir.per,"Results_MCMC_Cooked.txt"),
    #             file.path(dir.per,"Results_MCMC_Cooked_OLD.txt"))
    
    cook = read.table(file.path(dir.res.case,"/SPD/Results_MCMC_Cooked.txt"),header=T)
    #Cleaning SPD results to keep only the b of period n
    ### WARNING, to adapt to the different possibility of par var in SPD
    
    ## why no kx_10??
    if (period != nperiod){
      cook=cook[c(paste0("b1_",period),"a1","c1", #bac1
                  paste0("b2_",period),"a2","c2", #bac2
                  #paste0("b3_",period),"a3","c3", #bac3
                  "Y1_gamma1","Y1_gamma2","LogPost",#y & logpost
                  paste0("k1_",period),paste0("k2_",period),#paste0("k3_",period),#ks
                  paste0("ktype1_",period),paste0("ktype2_",period)#,paste0("ktype3_",period) #ktypes
      )]
    } else {
      cook=cook[c(paste0("b1_",period),"a1","c1", #bac1
                  paste0("b2_",period),"a2","c2", #bac2
                  #paste0("b3_",period),"a3","c3", #bac3
                  "Y1_gamma1","Y1_gamma2","LogPost",#y & logpost
                  paste0("k1_",period-1),paste0("k2_",period-1),#paste0("k3_",period),#ks
                  paste0("ktype1_",period-1),paste0("ktype2_",period-1)#,paste0("ktype3_",period) #ktypes
      )]
    }
    
    #Write the new cleaned result file
    write.table(cook,file.path(dir.bar,"Results_MCMC_Cooked.txt"),
                row.names = F,sep= "     ")
    write.table(cook,file.path(dir.per,"Results_MCMC_Cooked.txt"),
                row.names = F,sep= "     ")
    
    ############ Apply BaM prediction on the cleaned results file
    ######## writing data in BaM folder for the current period
    ## Gaugings
    GauPer = Gau[which(Gau$Period == period),]
    GauPer$Period = as.numeric(GauPer$Period)
    write.table(x = GauPer,
                file = paste0(dir.bar,"/Gaugings_data.txt"), row.names = F, sep = "    ",)
    nobs_gaug = nrow(GauPer)
    
    buffer=c(paste0("b1.",period),"a1","c1", #bac1
             paste0("b2.",period),"a2","c2")#, #bac2
    
    priorP =  fit(MC$sim[c(paste0("b1.",period),"a1","c1", #bac1 #simulations
                           paste0("b2.",period),"a2","c2")],#, #bac2
                  #paste0("b3.",period),"a3","c3")], #bac3
                  rep(margins.bac,ncontrol), #margins
                  c(paste0("b1_",period),"a1","c1", #bac1 #names
                    paste0("b2_",period),"a2","c2")#, #bac2
                  )
    
    #### Limni period n 
    AmaxP = Amax[which(Amax$period == period),3]
    MatSpag = Noisy[which(Amax$period == period),]
    
    Limnisize=length(AmaxP)
    write.table(x = list(MatSpag),
                file = paste0(dir.bar,"/limni.txt"), row.names = F,sep = "    ",col.names = F)
    
      #################################################
      # setwd(dir.bar)
      #Baratin Config corrected to allow spagetti predict files & do not print BaM state
      BaRatin_configNspagPeriod( dir               = dir.bar,
                                 period            = period,
                                 nperiod           = nperiod,
                                 nsim              = nspag,#1000
                                 prior             = priorP,
                                 ncontrol          = ncontrol,
                                 M                 = M,
                                 nobs              = nobs_gaug,
                                 Ncycles           = 1000,#100
                                 ngrid             = ngrid,
                                 nlimni            = length(AmaxP),
                                 predictionRC      = T,#T
                                 predictionQt      = T,#T
                                 predictionPrior   = F,#T
                                 simMCMC           = F,#T
                                 mcmc.prior        = 100,#1000
                                 remnant.err.model = remnant.err.model,
                                 g1.prior          = g1.prior,
                                 g2.prior          = g2.prior,
                                 g1.distr.type     = g1.distr.type,
                                 g2.distr.type     = g2.distr.type,
                                 nspag = nspag, Hmax = Hmax, dir_code = dir.bam,
                                 dir.bam.short = paste0("BaR_",case,"/"))
      #Run BaM
      message("***************************************************************"); flush.console()
      message(c("Applying BaRatin to period!!!  Wait ... ")); flush.console()
      message("***************************************************************"); flush.console()
      setwd(dir.bam)
      system2(paste0(dir.bam,"/BaM.exe")) 
      
      # save results in period folder
      # #Mcmc files
      # file.copy(from = file.path(dir.bar,list.files(dir.BaM,pattern = "Cooked" )),
      #           to = file.path(dirSub,paste0(list.files(dir.BaM,pattern = "Cooked"),subp)), overwrite = T)
      #pred config files
      # file.copy(from = file.path(dir.BaM,list.files(dir.BaM,pattern = "Pred" )),
      #           to = file.path(dirSub,paste0(list.files(dir.BaM,pattern = "Pred"),subp)), overwrite = T)
      #Qt predictions
      file.copy(from = file.path(dir.bar,list.files(dir.bar,pattern = "Qt" )),
                to = file.path(dir.per,list.files(dir.bar,pattern = "Qt")), overwrite = T)
      #RC predictions
      file.copy(from = file.path(dir.bar,list.files(dir.bar,pattern = "Qrc" )),
                to = file.path(dir.per,list.files(dir.bar,pattern = "Qrc")), overwrite = T)
      #model
      # file.copy(from = file.path(dir.bar,list.files(dir.bar,pattern = "Model" )),
      #           to = file.path(dirSub,paste0(list.files(dir.BaM,pattern = "Model"),subp)), overwrite = T)
    
    ### period done
    message(paste0("period : ",period,"    "))
    
    ### TOT U
    Tot[bufferAn,] = as.matrix(round(read.table(paste0(dir.per,"/Qt_TotalU.spag"),
                                                header = F,),2)[,c(1:nspag)])   
    ### PARAM U
    Param[bufferAn,] = as.matrix(round(read.table(paste0(dir.per,"/Qt_ParamU.spag"),
                                                  header = F,),2)[,c(1:nspag)])  
    ### MAXPOST U (NOT THE REAL MAXPOST, 500 spag limni x 1 spag MAXPOST CT)
    Mp[bufferAn,] = as.matrix(round(read.table(paste0(dir.per,"/Qt_Maxpost.spag"),
                                               header = F,),2)[,c(1:nspag)])  
    ### REAL MAXPOST ???? 1 spag MP limni x 1 spag MP CT & fun
    ### maxpost limni = median h
    if(length(AmaxP)>1){TrueMp[bufferAn] = approx(x = Hgrid, y = DatasRC[[period]]$maxpost, 
                              xout = apply(MatSpag,MARGIN = 1, median) )$y
    } else {TrueMp[bufferAn] = approx(x = Hgrid, y = DatasRC[[period]]$maxpost,
                                      xout = median(MatSpag)    )$y}
    
  }

  message(paste0("Computation duration :", (Sys.time()-start.time), " min"))
  
  
  AmaxQuants = data.frame(
    an = Amax$Y,
    mp = TrueMp,
    mp2.5 = apply(Mp,MARGIN = 1,quantile,probs=0.025,na.rm=T),
    mp97.5 = apply(Mp,MARGIN = 1,quantile,probs=0.975,na.rm=T),
    param2.5 = apply(Param,MARGIN = 1,quantile,probs=0.025,na.rm=T),
    param97.5 = apply(Param,MARGIN = 1,quantile,probs=0.975,na.rm=T),
    tot2.5 = apply(Tot,MARGIN = 1,quantile,probs=0.025,na.rm=T),
    tot97.5 = apply(Tot,MARGIN = 1,quantile,probs=0.975,na.rm=T)
  )
  
  ## Number of gaugings by each year
  Gaudat = data.frame(Date = ymd(read.csv2(paste0(dir.data,"JauBeaucaireOLD.csv"))[,1]),
                      Y = NA)
  Gaudat$Y = year(Gaudat$Date)
  Njau = 1816 : 1967
  for(i in 1816:1967){Njau[i-1815] = length(which(Gaudat$Y==i))}
  Njau[which(Njau==0)] = NA
  
  ICreldif = ggplot()+
    scale_x_continuous(name = expression("Year"), expand=c(0,2))+
    scale_y_continuous(name = expression(paste("Annual maximum discharge [",m^3,".",s^-1,"]",sep="")),
                       expand=c(0.1,0))+
    labs(x = "Time [day]", y = expression(paste("Annual maximum discharge [",m^3,".",s^-1,"]",sep=""))) +
    theme_bw(base_size=15)+
    geom_ribbon(data = AmaxQuants,
                aes(x = an, ymin = ((tot2.5-mp)/mp)*100, ymax= ((tot97.5-mp)/mp)*100,#),
                    fill= "Remnant uncertainty"))+#, alpha = 0.5,show.legend = T)+
    geom_ribbon(data = AmaxQuants,
                aes(x = an, ymin = ((param2.5-mp)/mp)*100, ymax= ((param97.5-mp)/mp)*100,#),
                    fill= "Parametric uncertainty"), alpha = 0.9)+#, show.legend = T)+
    geom_ribbon(data = AmaxQuants,
                aes(x = an, ymin = ((mp2.5-mp)/mp)*100, ymax= ((mp97.5-mp)/mp)*100,#),
                    fill="Limnimetric uncertainty"))+#, alpha = 0.5,show.legend = T)+
    geom_line(data = AmaxQuants, aes(x = an, y = 0))+
    geom_vline(xintercept = year(Tshifts),show.legend = T)+
    theme(legend.position = "right")+#,legend.text = element_text(size=8))+
    scale_fill_manual(name = element_blank(),
                      values = c("#fec44f","#fa9fb5","#f03b20"))+
    geom_point(aes(x =1816:1967,y = 0,size = Njau),color='blue',alpha = 0.5)+
    theme( axis.text=element_text(size=12)
           ,axis.title=element_text(size=14)
           ,legend.text=element_text(size=15)
           ,legend.title=element_text(size=15)
           #,plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm")
           ,legend.key.size=unit(1, "cm"))+
    scale_size(range = c(1, 10), name="Gaugings / year")+
    coord_cartesian(ylim = c(-30,70))
  
  
  pdf(paste0(dir.plots,"/",case,"/ICrel_AMAX.pdf"),14)
    print(ICreldif)
  dev.off()
    
  
  IC_Amax = ggplot()+
    scale_x_continuous(name = expression("Year"), expand=c(0,2))+
    scale_y_continuous(name = expression("Relative uncertainty for annual maximum discharges [%]"),
                       expand=c(0.1,0))+
    labs(x = "Time [day]", y = "Relative uncertainty for annual maximum discharges [%]") +
    theme_bw(base_size=15)+
    geom_ribbon(data = AmaxQuants,
                aes(x = an, ymin = tot2.5, ymax = tot97.5,#),
                    fill= "Remnant uncertainty"))+  #, alpha = 0.5,show.legend = T)+
    geom_ribbon(data = AmaxQuants,
                aes(x = an, ymin = param2.5, ymax = param97.5,#),
                    fill= "Parametric uncertainty"), alpha = 0.9)+#, show.legend = T)+
    geom_ribbon(data = AmaxQuants,
                aes(x = an, ymin = mp2.5, ymax = mp97.5,#),
                    fill="Limnimetric uncertainty"))+#, alpha = 0.5,show.legend = T)+
    geom_line(data = AmaxQuants, aes(x = an, y = mp))+
    geom_vline(xintercept = year(Tshifts),show.legend = T)+
    theme(legend.position = "right")+#,legend.text = element_text(size=8))+
    scale_fill_manual(name = element_blank(),
                values = c("#fec44f","#fa9fb5","#f03b20"))+
    theme( axis.text=element_text(size=12)
           ,axis.title=element_text(size=14)
           ,legend.text=element_text(size=15)
           ,legend.title=element_text(size=15)
           ,legend.key.size=unit(1, "cm"))
    
  
  pdf(paste0(dir.plots,"/",case,"/IC_AMAX.pdf"),10)
    print(IC_Amax)
  dev.off()
  
######################################
#####  SAVE PROPAGATION RESULTS ######
######################################
    
    write.table(round(AmaxQuants,2), paste0(dir.res.case,"Quantiles_AmaxPt.txt"),row.names = F)
    write.table(round(Mp,2), paste0(dir.res.case,"MpSpagsAMAX_Pt.txt"),row.names = F)
    write.table(round(Param,2), paste0(dir.res.case,"ParamSpagsAMAX_Pt.txt"),row.names = F)
    write.table(round(Tot,2), paste0(dir.res.case,"TotSpagsAMAX_Pt.txt"),row.names = F)
    
    
    
    
    