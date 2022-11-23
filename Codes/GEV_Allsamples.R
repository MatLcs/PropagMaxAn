##################################################################
### COMPARE IC's tot & sampling on 50, 100, 150, 200 years samples
##################################################################
rm(list=ls())
dev.off()
source("C://Users/mathieu.lucas/Desktop/GitMat/PropagMaxAn/Codes/dirs.R")
source(paste0(dir.codes,"module_BaRatin.r"))
source(paste0(dir.codes,"Fun_SPD.r"))
### dir for plots
dir.plots.gev = paste0(dir.plots,"/GeV_all")
dir.create(dir.plots.gev,showWarnings = F)

######################################
############ DATA LOADING ############
######################################
SpagsHydro = read.table(paste0(dir.res,"spags_uTot_Amax.txt"))
QuantHydro = read.table(paste0(dir.res,"Quantiles_Amax.txt"),header=T)
######################################
######### MODEL PREPARATION ##########
######################################
#### 205 cases of sample size
nyears = seq(20, length(QuantHydro$an), 2)#c(20 : length(QuantHydro$an))
nspag = 200#dim(SpagsHydro)[2] 
nsim = 4000
#### SToods parameters, common GEV priors to all cases
Pos <- parameter(name='Pos',init = 6000) 
Ech =  parameter(name='Ech',init = 1000) 
Form = parameter(name='Form',init = 0,priorDist='Gaussian',priorPar=c(0,0.2))
Param = list(Pos,Ech,Form)
#### Quantiles extracted up to Q1000 & return period associated
prob = c(0.99,0.999) ; Pr = c(100,1000)#1/(1-prob)
#### Initializing results list
QuantAll = data.frame()
#### Initializing plots list
GgGeV = list()
GgU = list()
#### Creating results folder
dir.res.stoods = paste0(dir.res,"/GeV_Amax_Allsamples")
dir.create(dir.res.stoods,showWarnings = F)
######################################
######### MODEL PREPARATION ##########
######################################
for(case in 1: length(nyears)){
# case = 1
  message("***************************************************************")
  message(paste0("Case = ",case))
  message("***************************************************************")
  #### Creating case folder
  dir.case = paste0(dir.res.stoods,"/dircase")
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
  for(spag in sample(1:500,nspag)){
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
                                                           scale = MegaSpagTot$Ech[i],shape = -1*MegaSpagTot$Form[i]) }
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
  # Freq = Q.case[order(Q.case$mp),]
  # Freq$Fr = (seq(1:length(Q.case$an))-0.5)/length(Q.case$an)
  # Freq$Pr = 1/(1-Freq$Fr)
  
  Quants = data.frame(Pr=Pr, Mp=Mp.Quant, Qhyd_2=quanthyd[1,] ,Qhyd_9=quanthyd[2,],
                      Qtot_2=quanttot[1,],Qtot_9=quanttot[2,],nyears = nyears[case])
  
  QuantAll = rbind(QuantAll,Quants)

}



# QuantAll$Pr = rep(c(100,1000),53)

# Quant1000 = QuantAll[which(QuantAll$Pr == 1000),]

GGQ1000 = ggplot(QuantAll[which(QuantAll$Pr == 1000),])+
  geom_ribbon(aes(x=nyears, ymin=0, ymax = Qtot_9 - Qtot_2, fill = "1-Total unc."))+
  geom_ribbon(aes(x=nyears, ymin=0, ymax = Qhyd_9 - Qhyd_2, fill = "2-Streamflow unc."))+
  geom_line(aes(x=nyears, y = Mp, color = "1000-year flood maxpost value"),lwd=1)+
  xlab("Sample size [years]")+
  ylab("95% confidence interval width [m3/s]")+
  scale_fill_manual(name = element_blank(),
                    values = c("#67a9cf","#fec44f"))+
  theme_bw(base_size=15)+
  labs(title="1000-year flood")+
  theme(legend.title = element_blank())

GGQ1000
ggsave(filename = "WidthQ1000.pdf",path = dir.plots.gev,width = 10, height = 7)

GGQ100 = ggplot(QuantAll[which(QuantAll$Pr == 100),])+
  geom_ribbon(aes(x=nyears, ymin=0, ymax = Qtot_9 - Qtot_2, fill = "1-Total unc."))+
  geom_ribbon(aes(x=nyears, ymin=0, ymax = Qhyd_9 - Qhyd_2, fill = "2-Streamflow unc."))+
  geom_line(aes(x=nyears, y = Mp, color = "1000-year flood maxpost value"),lwd=1)+
  xlab("Sample size [years]")+
  ylab("95% confidence interval width [m3/s]")+
  scale_fill_manual(name = element_blank(),
                    values = c("#67a9cf","#fec44f"))+
  theme_bw(base_size=15)+
  labs(title="100-year flood")+
  theme(legend.title = element_blank())


GGQ100
ggsave(filename = "WidthQ100.pdf",path = dir.plots.gev,width = 10, height = 7)


write.table(round(QuantAll,3),file = paste0(dir.res,"QuantWidth.txt"),row.names = F)



####### PLOTS ONLY



QuantAll = read.table(paste0(dir.res,"QuantWidth.txt"), header = T)


GGQ1000 = ggplot(QuantAll[which(QuantAll$Pr == 1000),])+
  # geom_line(data = rev(QuantHydro[1:(205-20),]), aes(x = QuantAll$nyears, y = mp), col ="lightgrey")+
  geom_ribbon(aes(x=nyears, ymin=Qtot_9, ymax =  Qtot_2, fill = "1-Total unc."))+
  geom_ribbon(aes(x=nyears, ymin=Qhyd_9, ymax =  Qhyd_2, fill = "2-Streamflow unc."))+
  geom_line(aes(x=nyears, y = Mp, color = "maxpost"),lwd=1)+
  geom_vline(xintercept = 2020-1970, col ="lightgrey", lwd=1, lty = 2)+
  geom_vline(xintercept = 2020-1856, col ="lightgrey", lwd=1, lty = 2)+
  geom_vline(xintercept = 2020-1840, col ="lightgrey", lwd=1, lty = 2)+
  xlab("Sample size [years]")+
  ylab("Discharge [m3/s]")+
  scale_fill_manual(name = element_blank(),
                    values = c("#67a9cf","#fec44f"))+
  theme_bw(base_size=15)+
  labs(title="1000-year flood")+
  theme(legend.title = element_blank())+
  coord_cartesian(ylim=c(10000,50000))+
  scale_y_log10()

GGQ1000

ggsave(filename = "Q1000SSize.pdf",path = dir.plots.gev,width = 10, height = 7)


GGQ100 = ggplot(QuantAll[which(QuantAll$Pr == 100),])+
  geom_ribbon(aes(x=nyears, ymin=Qtot_9, ymax =  Qtot_2, fill = "1-Total unc."))+
  geom_ribbon(aes(x=nyears, ymin=Qhyd_9, ymax =  Qhyd_2, fill = "2-Streamflow unc."))+
  geom_line(aes(x=nyears, y = Mp, color = "1000-year flood maxpost value"),lwd=1)+
  xlab("Sample size [years]")+
  ylab("Discharge [m3/s]")+
  scale_fill_manual(name = element_blank(),
                    values = c("#67a9cf","#fec44f"))+
  theme_bw(base_size=15)+
  labs(title="1000-year flood")+
  theme(legend.title = element_blank())

GGQ100







