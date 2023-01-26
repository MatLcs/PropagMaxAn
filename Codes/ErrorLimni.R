rm(list=ls())
dev.off()
source("C://Users/mathieu.lucas/Desktop/GitMat/PropagMaxAn/Codes/dirs.R")
# set.seed(42)

####### 
####### DATA LOADING AND FORMATTING
# Data Bcr Restit
H_hor=read.csv2(paste0(dir.data,"/BcrResHorCor.csv"),header=T)
H_hor$Date = ymd_hms(H_hor$Date)
H_hor$Day = as_date(H_hor$Date)
H_hor$Time = hour(H_hor$Date)
H_hor$y = year(H_hor$Date)
# without 01-01-2021
H_hor = H_hor[-(nrow(H_hor)),]
# Data Pt Bcr
LimniPt = read.csv2(paste0(dir.data,"PtBcrJournalier.csv"))
LimniPt$Date = ymd(LimniPt$Date)
LimniPt$H = as.numeric(LimniPt$H)

#GradPt = LimniPt$H[-1] - LimniPt$H[-nrow(LimniPt)]

####### 
### ERROR MODEL FITTING BASED ON BCR RESTITUTION STATION
# True daily max of bcr restit
MaxJ = aggregate(x = H_hor$H, by = list(H_hor$Day), FUN = max)
names(MaxJ) = c("day","h")
# 12h stage maxima
Midi = H_hor[which(H_hor$Time==12),c(3,2)]
# 3 stages/day maxima (7,12,18)
Three = H_hor[which(H_hor$Time==7 | H_hor$Time==12 | H_hor$Time==18),c(3,2)]
MaxThree = aggregate(x=Three$H, by = list(Three$Day), FUN=max)
names(MaxThree) = c("day","h")

# Residuals between measured and true maxima
DifMid = MaxJ$h - Midi$H
DifThree = MaxJ$h - MaxThree$h 

####### 
#### Annual maxima selection
#init DF
MaxARestit = data.frame(Y = unique(year(H_hor$Date)), h = NA, 
                        Date = as.Date("1970-01-01"),Emid = NA, Ethree = NA)
#loop on all years
for(years in unique(H_hor$y)){
  #indice of the max
  indmax = which(H_hor$y == years)[which.max(H_hor$H[which(H_hor$y == years)])]
  #stage max
  MaxARestit$h[which(MaxARestit$Y == years)] = H_hor$H[indmax]
  #stage day
  MaxARestit$Date[which(MaxARestit$Y == years)] = H_hor$Day[indmax]
  #1/d error of the max (true daily max - 12h stage)
  MaxARestit$Emid[which(MaxARestit$Y == years)] = max(H_hor$H[which(H_hor$Day == H_hor$Day[indmax])]) -
                                                     H_hor$H[which(H_hor$Day == H_hor$Day[indmax])][13]
  #3/d error of the max (true daily max - max(7,12,18) stage)
  MaxARestit$Ethree[which(MaxARestit$Y == years)] = max(H_hor$H[which(H_hor$Day == H_hor$Day[indmax])]) -
                                          max(H_hor$H[which(H_hor$Day == H_hor$Day[indmax])][c(8,13,19)]) 
}

####### 
## LOGNORM or EXP fit on annual max stage errors + plot
# pdf(file = paste0(dir.plots,"FitsMaxAstage_BcrRestit.pdf"),width = 6,height = 10)
  # layout(matrix(1:8, nrow = 4, ncol = 2))
  xtheo = seq(0,4,0.01)
  br = 20
  #EXP FIT 1/D
  E1 = mean(MaxARestit$Emid,na.rm=T)
  lamb1 = 1/E1
  hist(MaxARestit$Emid,prob = T, main = "Exponential 1/d",ylim = c(0,3),
       breaks = br, xlab = "Stage error [m]")
  ytheo = dexp(xtheo,lamb1);lines(xtheo, ytheo,col="deepskyblue",lwd=3)
  q_theo=qexp((1:length(MaxARestit$Emid))/length(MaxARestit$Emid),lamb1)
  qqplot(q_theo,MaxARestit$Emid,xlim = c(0,3),ylim=c(0,4),col="deepskyblue",ylab="Stage error [m]")
  abline(a=0,b=1, col="red")
  #EXP FIT 3/D
  E3 = mean(MaxARestit$Ethree,na.rm=T)
  lamb3 = 1/E3
  hist(MaxARestit$Ethree,prob = T, main = "Exponential 3/d",ylim = c(0,15),
       breaks = br,xlab = "Stage error [m]")
  ytheo = dexp(xtheo,lamb3);lines(xtheo, ytheo,col="deepskyblue",lwd=3)
  q_theo=qexp((1:length(MaxARestit$Ethree))/length(MaxARestit$Ethree),lamb3)
  qqplot(q_theo,MaxARestit$Ethree,xlim = c(0,1),ylim=c(0,1),col="deepskyblue",ylab="Stage error [m]")
  abline(a=0,b=1, col="red")
  #LN FIT 1/D
  mu1 = mean(log(MaxARestit$Emid)[which(log(MaxARestit$Emid)!=-Inf)])
  sigma1 = sd(log(MaxARestit$Emid)[which(log(MaxARestit$Emid)!=-Inf)])
  hist(MaxARestit$Emid,prob = T,ylim=c(0,3),main = "LogNorm 1/d",
       breaks = br,xlab = "Stage error [m]")
  ytheo = dlnorm(xtheo,mu1,sigma1)
  lines(xtheo, ytheo,col="blue",lwd=3)
  q_theo=qlnorm((1:length(MaxARestit$Emid))/length(MaxARestit$Emid),mu1,sigma1)
  qqplot(q_theo,MaxARestit$Emid,xlim = c(0,3),ylim=c(0,4),col="blue",ylab="Stage error [m]")
  abline(a=0,b=1, col="red")
  #LN FIT 3/D
  mu3 = mean(log(MaxARestit$Ethree)[which(log(MaxARestit$Ethree)!=-Inf)])
  sigma3 = sd(log(MaxARestit$Ethree)[which(log(MaxARestit$Ethree)!=-Inf)])
  hist(MaxARestit$Ethree,prob = T,ylim=c(0,15),main = "LogNorm 3/d",
       breaks = br,xlab = "Stage error [m]")
  ytheo = dlnorm(xtheo,mu3,sigma3)
  lines(xtheo, ytheo,col="blue",lwd=3)
  q_theo=qlnorm((1:length(MaxARestit$Ethree))/length(MaxARestit$Ethree),mu3,sigma3)
  qqplot(q_theo,MaxARestit$Ethree,xlim = c(0,1),ylim=c(0,1),col="blue",ylab="Stage error [m]")
  abline(a=0,b=1, col="red")

# dev.off()


####### 
####### APPLY THE FITTED ERROR TO PT BCR ANNUAL MAX STAGES
LimniPt$y = year(LimniPt$Date)
#init matrix
nspag = 500
years = 1816:1967
NoisyLimni = matrix(ncol = nspag, nrow = length(years))
HmaxA = data.frame(y = years,Date = as.Date("1900-01-01"), H = NA)
countrec = 1

## datum reference error init delta4
d4_pt = rnorm(1,0,1)*0.06

for(i in 1 : nrow(NoisyLimni) ){
  
  if(years[i] != 1873){
    indmax = which(LimniPt$y == years[i])[which.max(LimniPt$H[which(LimniPt$y == years[i])])]
    Hmax = LimniPt$H[indmax]
    Dmax = LimniPt$Date[indmax]
  ## annuaire hydrologique 1873
  } else {Hmax = 5.7; Dmax = as.Date("1873-03-19")}

  ## draw a new datum reference err each 25 years
  if(countrec == 25){
    countrec = 1
    d4_pt = rnorm(1,0,1)*0.06  }
  
  HmaxA$H[i] = Hmax
  HmaxA$Date[i] = Dmax
  
  ##### before 1841
  if(years[i] < 1841){
    NoisyLimni[i,] = Hmax  +
      ## measurement frequency err delta5 1/d
      rexp(n = nspag, rate = lamb1)  +
      ## gauge reading error delta 1
      rnorm(nspag,0, 0.10/1.96) +
      ## datum reference error
      d4_pt
    
  ##### after 1841
  } else {
    #### < 5m
    if(Hmax < 5){
      NoisyLimni[i,] = Hmax  +
        ## measurement frequency err delta5 3/d
        rexp(n = nspag, rate = lamb3)  +
        ## gauge reading error delta 1
        rnorm(nspag,0, 0.10/1.96) +
        ## datum reference error
        d4_pt
      
      #### > 5m
    } else {
      NoisyLimni[i,] = Hmax +
        ## gauge reading error delta 1
        rnorm(nspag,0,1)*(0.10/1.96) +
        ## datum reference error
        d4_pt
    }
  }
  countrec = countrec+1
}
  

## Bard & Lang dike breaks corrections
# 1840 +/- 0.94
NoisyLimni[which(HmaxA$y == 1840),] = HmaxA$H[which(HmaxA$y == 1840)] + rnorm(nspag,0,0.94/1.96)
# 1856 +/- 0.4
NoisyLimni[which(HmaxA$y == 1856),] = HmaxA$H[which(HmaxA$y == 1856)] + rnorm(nspag,0,0.4/1.96)
# 1841 +/- 0.4
NoisyLimni[which(HmaxA$y == 1841),] = HmaxA$H[which(HmaxA$y == 1841)] + rnorm(nspag,0,0.94/1.96)

## Condition C4 flood Pichard before 1841
# for(i in 1 : (1841-1816) ){
#     if ( HmaxA$H[i] < 7 & any(NoisyLimni[i,] > 7) ){
#       NoisyLimni[i,which(NoisyLimni[i,] > 7)] = 7
#     }
# }

### Quantiles
HmaxA$min = round(apply(NoisyLimni,MARGIN = 1,quantile,probs = 0.025),2)
HmaxA$max = round(apply(NoisyLimni,MARGIN = 1,quantile,probs = 0.975),2)
HmaxA$med = round(apply(NoisyLimni,MARGIN = 1,median),2)
NoisyLimni = round(NoisyLimni,2)


### Pettitt 
# pettitt.test(HmaxA$med)
# pettitt.test(HmaxA$max)
# pettitt.test(HmaxA$H)
# mk.test(HmaxA$med)
# mk.test(HmaxA$max)

### IC abs plot
IC_an = ggplot()+
  labs(x = "Year", y = "AMAX stage & uncertainty interval [m]") +
  theme_bw(base_size=15)+
  geom_ribbon(data = HmaxA,
              aes(x = y, ymin = min, ymax=max,fill= "Uncertainty interval"))+
  geom_line(data = HmaxA, aes(x = y, y = H),lwd=1)+
  # geom_line(data = HmaxA, aes(x = y, y = med),lwd=1,col="red")+
  theme(legend.position = "NULL")+
  theme( axis.text=element_text(size=12)
         ,axis.title=element_text(size=14)
         ,legend.text=element_text(size=15)
         ,legend.title=element_text(size=15)
         ,legend.key.size=unit(1, "cm"))+
  geom_line(data=HmaxA, aes(x=y, y = 5 ),col="lightgrey",lwd=1,lty=2)  +
  # geom_vline(aes(xintercept=1841),col="grey" ,lty=2,lwd = 1)+
  coord_cartesian(ylim = c(2.5,10.5))  
  
# IC_an


### IC rel plot
IC_perc = ggplot()+
  labs(x = "Year", y = "Uncertainty interval diff. to measured stage [m]") +
  theme_bw(base_size=15)+
  geom_ribbon(data = HmaxA,
              aes(x = y, ymin = min-H, ymax=max-H, fill = "Uncertainty interval"))+
  geom_line(data = HmaxA, aes(x = y, y = 0), lwd=1)+
  geom_line(data = HmaxA, aes(x = y, y = med-H), col="red",lwd=1)+
  theme(legend.position = "NULL")+
  theme( axis.text=element_text(size=12)
         ,axis.title=element_text(size=14)
         ,legend.text=element_text(size=15)
         ,legend.title=element_text(size=15)
         ,legend.key.size=unit(1, "cm"))  
  # geom_vline(aes(xintercept=1841),lwd=1, lty =2 , col="grey")
  # coord_cartesian(ylim = c(2,10))  

# IC_perc

ggarrange(IC_an,IC_perc,ncol = 1,align = "hv")
# ggsave(path = dir.plots, filename = "StageErrorAMAX_bcr.pdf",width = 10, height = 8)

###### NOISY RESTIT

NoisyRestit = matrix(ncol = nspag, nrow = nrow(MaxARestit))

for(i in 1:nrow(NoisyRestit)){
  NoisyRestit[i,] = MaxARestit$h[i] +
    ## sensor precision
    rnorm(nspag,0, 0.01/sqrt(3))+
    #calibration every 6 month, new draw every year
    rnorm(nspag,0, 0.12/1.96) }

# NoisyRestit = round(NoisyRestit,2)

HmaxA_Res = MaxARestit[,c(1,3,2)]
HmaxA_Res$min = round(apply(NoisyRestit,MARGIN = 1,quantile,probs = 0.025),3)
HmaxA_Res$max = round(apply(NoisyRestit,MARGIN = 1,quantile,probs = 0.975),3)
HmaxA_Res$med = round(apply(NoisyRestit,MARGIN = 1,median),2)
names(HmaxA_Res) = c("y","Date","H","min","max","med")



### IC abs plot
IC_an_res = ggplot()+
  labs(x = "Year", y = "AMAX stage & uncertainty interval [m]") +
  theme_bw(base_size=15)+
  geom_ribbon(data = HmaxA_Res,
              aes(x = y, ymin = min, ymax=max,fill= "Uncertainty interval"))+
  geom_line(data = HmaxA_Res, aes(x = y, y = H),lwd=1)+
  # geom_line(data = HmaxA, aes(x = y, y = med),lwd=1,col="red")+
  theme(legend.position = "NULL")+
  theme( axis.text=element_text(size=12)
         ,axis.title=element_text(size=14)
         ,legend.text=element_text(size=15)
         ,legend.title=element_text(size=15)
         ,legend.key.size=unit(1, "cm"))
  #coord_cartesian(ylim = c(2.5,10.5))  

# IC_an_res

HmaxA$st = as.factor("Pt")
HmaxA_Res$st = as.factor("Res")
BothMax = rbind(HmaxA,
                c(1968,NA,NA,NA,NA,NA,NA),
                c(1969,NA,NA,NA,NA,NA,NA),
                c(1970,NA,NA,NA,NA,NA,NA),
                HmaxA_Res)


### IC abs plot
IC_an_both = ggplot()+
  labs(x = "Year", y = "Stage [m]") +
  theme_bw(base_size=15)+
  geom_ribbon(data = BothMax,
              aes(x = y, ymin = min, ymax=max,fill= "Stage uncertainty"))+
  geom_line(data = BothMax, aes(x = y, y = H, color = "Measured Stage"),lwd=0.7)+
  # geom_vline(aes(xintercept=1841), lwd = 1, lty = 2)+
  geom_vline(aes(xintercept=c(1968,1969,1970), color = "Vallabrègues Works"), lwd = 4)+
  geom_line(data = BothMax, aes(x = y, y = 0#med
                                , col = "Median stage"),lwd = 0.5)+
  scale_fill_manual(values = c("#fec44f"))+
  scale_color_manual(values = c("black","red","lightgrey"))+
  theme( axis.text=element_text(size=20)
         ,axis.title=element_text(size=20)
         ,legend.text=element_text(size=20)
         ,legend.title=element_blank()#element_text(size=18)
         ,legend.key.size=unit(1, "cm"))+
  coord_cartesian(ylim = c(2.5,11.5))
  
# IC_an_both

### IC rel both
IC_rel_both = ggplot()+
  labs(x = "Year", y = "Difference to measured stage [m]") +
  theme_bw(base_size=15)+
  geom_ribbon(data = BothMax,
              aes(x = y, ymin = min-H, ymax=max-H, fill = "Stage uncertainty"))+
  geom_line(data = BothMax, aes(x = y, y = 0, col = "Measured stage"), lwd=1)+
  geom_line(data = BothMax, aes(x = y, y = med-H, color = "Median stage"),lwd=1)+
  geom_vline(aes(xintercept=c(1968,1969,1970), color = "Vallabrègues Works"), lwd = 4)+
  # geom_vline(aes(xintercept=1841), lwd = 1, lty = 2)+
  scale_fill_manual(values = c("#fec44f"))+
  scale_color_manual(values = c("black","red","lightgrey"))+
  theme( axis.text=element_text(size=20)
         ,axis.title=element_text(size=20)
         ,legend.text=element_text(size=20)
         ,legend.title=element_blank()#element_text(size=18)
         ,legend.key.size=unit(1, "cm"))  

# IC_rel_both

ggarrange(IC_an_both+
            theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x.bottom = element_blank()),
          IC_rel_both+theme(axis.text.x = element_text(size=20), axis.title.x = element_text(size=25))
          ,ncol = 1,align = "v", common.legend = T, legend = "right")
ggsave(path = dir.plots, filename = "StageErrorAMAX_BOTH.pdf",width = 15, height = 10)

### Save
dir.create(paste0(dir.res,"Limni/"))
dir.limni = paste0(dir.res,"Limni/")
names(HmaxA) = c("Y","Date","Hmes","Hinf","Hsup","Hmed")
names(HmaxA_Res) = c("Y","Date","Hmes","Hinf","Hsup","Hmed")
write.table(HmaxA[,(1:6)], paste0(dir.limni,"amaxH_Pt.txt"),row.names = F)
write.table(HmaxA_Res[,(1:6)], paste0(dir.limni,"amaxH_Res.txt"),row.names = F)

write.table(NoisyLimni,paste0(dir.limni,"NoisyHmaxAn_Pt.txt"), row.names = F, col.names = F)
write.table(NoisyRestit,paste0(dir.limni,"NoisyHmaxAn_Res.txt"), row.names = F, col.names = F)


### width IC

bef40 = which(BothMax$y<1841)
aft40 = which(BothMax$y>=1841 & BothMax$y <= 1967)
res = which(BothMax$y > 1970)

mean(BothMax$max[bef40] - BothMax$min[bef40])

mean(BothMax$max[aft40] - BothMax$min[aft40])

mean(BothMax$max[res] - BothMax$min[res])

