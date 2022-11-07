rm(list=ls())
dev.off()
source("C://Users/mathieu.lucas/Desktop/GitMat/PropagMaxAn/Codes/dirs.R")
source(paste0(dir.codes,"module_BaRatin.r"))
source(paste0(dir.codes,"Fun_SPD.r"))

case1 = "Pt"
case2 = "Res"

######################################
############ DATA LOADING ############
######################################
  Spagtot1 = read.table(paste0(dir.res,case1,"/TotSpagsAMAX_",case1,".txt"),header = T)
  Quants1 = read.table(paste0(dir.res,case1,"/Quantiles_Amax",case1,".txt"), header = T)
  Spagtot2 = read.table(paste0(dir.res,case2,"/TotSpagsAMAX_",case2,".txt"), header = T)
  Quants2 = read.table(paste0(dir.res,case2,"/Quantiles_Amax",case2,".txt"), header = T)
  Gau1 = data.frame(Date = ymd(read.csv2(paste0(dir.data,"JauBeaucaireOLD.csv"))[,1]), Y = NA)
  Gau2 = data.frame(Date = ymd(read.csv2(paste0(dir.data,"JauBeaucaireCNR.csv"),header=T)[,1]), Y = NA)
  Tshifts1 = c(as.Date("1840-11-02"),
               as.Date("1816-05-15")+read.table(paste0(dir.data,"shift_timesPt.txt"),header=T)[,2])
  Tshifts2 = as.Date(c("1976-11-11","1994-01-09","2002-11-28","2003-12-05","2004-07-15","2004-12-03",
                      "2005-12-24","2008-06-28","2009-10-20","2016-11-24"))

######################################
######### DATA FORMATTING ############
######################################
  ### 3 missing years = annuaire hydro + ic +/- 10%
  Miss = data.frame(an = 1968:1970, mp = c(4760,4995,5510), mp2.5=NA,mp97.5=NA,param2.5=NA,param97.5=NA,
                    tot2.5 = NA , tot97.5 =NA)
  Miss$tot2.5 = qnorm(0.025,Miss$mp,(0.05*Miss$mp)) ; Miss$tot97.5 = qnorm(0.975,Miss$mp,(0.05*Miss$mp))
  MissNoisy = rbind(  rnorm(n = 500, mean = Miss$mp[1], sd = 0.05*Miss$mp[1]),
                      rnorm(n = 500, mean = Miss$mp[2], sd = 0.05*Miss$mp[2]),
                      rnorm(n = 500, mean = Miss$mp[3], sd = 0.05*Miss$mp[3]) )
  # combining both stations data
  Tshifts = c(Tshifts1,Tshifts2)
  Quants = rbind(Quants1,Miss,Quants2)
  Spags = rbind(Spagtot1, MissNoisy, Spagtot2)
  Gau = rbind(Gau1, Gau2)
  Gau$Y = year(Gau$Date)
  
  # number of gaugings by year for plots
  Njau = 1816 : 2020
  for(i in 1816 : 2020){Njau[i-1815] = length(which(Gau$Y==i))}
  Njau[which(Njau==0)] = NA

######################################
############## PLOTS #################
######################################
ICreldif = ggplot()+
  scale_x_continuous(name = expression("Year"), expand=c(0,2))+
  scale_y_continuous(name = expression("Relative uncertainty for annual maximum discharges [%]"),
                     expand=c(0.1,0))+
  labs(x = "Time [day]", y = "Relative uncertainty for annual maximum discharges [%]") +
  theme_bw(base_size=15)+
  geom_ribbon(data = Quants,
              aes(x = an, ymin = ((tot2.5-mp)/mp)*100, ymax= ((tot97.5-mp)/mp)*100,#),
                  fill= "3 - Remnant uncertainty"))+#, alpha = 0.5,show.legend = T)+
  geom_ribbon(data = Quants,
              aes(x = an, ymin = ((param2.5-mp)/mp)*100, ymax= ((param97.5-mp)/mp)*100,#),
                  fill= "2 - Parametric uncertainty"), alpha = 0.9)+#, show.legend = T)+
  geom_ribbon(data = Quants,
              aes(x = an, ymin = ((mp2.5-mp)/mp)*100, ymax= ((mp97.5-mp)/mp)*100,#),
                  fill="1 - Limnimetric uncertainty"))+#, alpha = 0.5,show.legend = T)+
  geom_ribbon(data = Quants[which(Quants$an >= 1967 & Quants$an <= 1971),], aes(x = an, 
              ymin = ((tot2.5-mp)/mp)*100, ymax= ((tot97.5-mp)/mp)*100, 
              fill = "4 - Reconstructed")   )+
  geom_line(data = Quants, aes(x = an, y = 0))+
  geom_vline(xintercept = year(Tshifts),show.legend = T, linetype=2)+
  theme(legend.position = "right")+#,legend.text = element_text(size=8))+
  scale_fill_manual(name = element_blank(),
                    values = c("#fec44f","#fa9fb5","#f03b20","grey"))+
  geom_point(aes(x = 1816 : 2020,y = 0,size = Njau),color='blue',alpha = 0.5)+
  theme( axis.text=element_text(size=12)
         ,axis.title=element_text(size=14)
         ,legend.text=element_text(size=15)
         ,legend.title=element_text(size=15)
         #,plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm")
         ,legend.key.size=unit(1, "cm"))+
  scale_size(range = c(1, 10), name="Gaugings / year")

pdf(paste0(dir.plots,"/ICrel_AMAX_Both.pdf"),14)
  print(ICreldif)
dev.off()


IC_Amax = ggplot()+
  scale_x_continuous(name = expression("Year"), expand=c(0,2))+
  scale_y_continuous(name = expression(paste("Annual maximum discharge [",m^3,".",s^-1,"]",sep="")),
                     expand=c(0.1,0))+
  labs(x = "Time [day]", y = expression(paste("Annual maximum discharge [",m^3,".",s^-1,"]",sep=""))) +
  theme_bw(base_size=15)+
  geom_ribbon(data = Quants,
              aes(x = an, ymin = tot2.5, ymax = tot97.5,#),
                  fill= "3-Remnant uncertainty"))+  #, alpha = 0.5,show.legend = T)+
  geom_ribbon(data = Quants,
              aes(x = an, ymin = param2.5, ymax = param97.5,#),
                  fill= "2-Parametric uncertainty"), alpha = 0.9)+#, show.legend = T)+
  geom_ribbon(data = Quants,
              aes(x = an, ymin = mp2.5, ymax = mp97.5,#),
                  fill="1-Limnimetric uncertainty"))+#, alpha = 0.5,show.legend = T)+
  geom_ribbon(data = Quants[which(Quants$an >= 1967 & Quants$an <= 1971),], aes(x = an, 
                         ymin = tot2.5, ymax= tot97.5, fill = "4 - Reconstructed") )+
  geom_line(data = Quants, aes(x = an, y = mp))+
  geom_point(aes(x = 1816 : 2020,y = 2000,size = Njau),color='blue',alpha = 0.5)+
  geom_vline(xintercept = year(Tshifts),show.legend = T, linetype = 2, size = 0.8, color = "lightgrey")+
  theme(legend.position = "right")+#,legend.text = element_text(size=8))+
  scale_fill_manual(name = element_blank(),
                    values = c("#fec44f","#fa9fb5","#f03b20","grey"))+
  theme( axis.text=element_text(size=12)
         ,axis.title=element_text(size=14)
         ,legend.text=element_text(size=15)
         ,legend.title=element_text(size=15)
         ,legend.key.size=unit(1, "cm"))+
  scale_size_continuous(range = c(3, 15),name="Gaugings / year")+
  coord_cartesian(ylim=c(2000,16000))

pdf(paste0(dir.plots,"IC_AMAX_Both.pdf"),14)
  print(IC_Amax)
dev.off()

######################################
##### WRITING FINAL HYDRO DATA #######
######################################

  ## SPAGS
  write.table(round(Spags,2), paste0(dir.res,"Spags_uTot_Amax.txt"),row.names = F, col.names = F)
  ## QUANTILES
  write.table(round(Quants,2), paste0(dir.res,"Quantiles_Amax.txt"),row.names = F)
  
  
######################################
######## HOMOGENEITY TESTS ###########
######################################
  nspag = ncol(Spags)
  pet = rep(NA,nspag)
  rupt = rep(NA,nspag)
  mann = rep(NA,nspag)

  for(spag in 1:nspag){
      pet[spag] = as.numeric(pettitt.test(Spags[,spag])[4])
      rupt[spag] = as.numeric(pettitt.test(Spags[,spag])$estimate[1])
      mann[spag] = mk.test(Spags[,spag])$p.value
  }

  hist(pet,breaks = 40)
  hist(mann, breaks = 40)
  hist(rupt[which(pet < 0.05)])
  
  length(which(pet <0.05))/500
  length(which(mann <0.05))/500

  plot(x = Quants$an, y = Quants$mp,type='l')
  lines(x = 1816:1886, y = rep( mean(Quants$mp[which(Quants$an <= 1886)]),
                                length(Quants$mp[which(Quants$an <= 1886)] ) ) )   
  lines(x = 1887:2020, y = rep( mean(Quants$mp[which(Quants$an > 1886)]),
                                length(Quants$mp[which(Quants$an > 1886)] ) ) )   
  
  # Quants$mp[1:70] = Quants$mp[1:70]-600
  
  pettitt.test(Quants$mp[41:70])
  mean(Quants$mp[1:41])
  mean(Quants$mp[41:70])
  
  bard = read.table("C://Users/mathieu.lucas/Desktop/MaxAnBard.txt",header = T)[,c(1,3)]
  bard$Date = dmy(bard$Date)
  bard$y = year(bard$Date)
  plot(x=1816:2016, y = bard$Valeur - Quants$mp[1:201], ylim = c(-2000,2000), ylab = "Q_Bard - Q_Mathieu [m3/s]")
  abline(v = year(Tshifts), col = "grey", lty = 2)
  
  mk.test(bard$Valeur)
  mk.test(Quants$mp)        
  pettitt.test(bard$Valeur)
  pettitt.test(Quants$mp)  
  
  plot(x=1816:2016, 
       y = ((bard$Valeur - Quants$mp[1:201])/Quants$mp[1:201])*100,# ylim = c(-2000,2000), 
       ylab = "(Q_Bard - Q_Mathieu)/Q_Mathieu * 100 [%]")
  abline(v = year(Tshifts), col = "grey", lty = 2)
  
  

    