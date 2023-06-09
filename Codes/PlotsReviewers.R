rm(list=ls())
source("C://Users/mathieu.lucas/Desktop/GitMat/PropagMaxAn/Codes/dirs.R")

#### stage old
PtStage = read.csv2(paste0(dir.data,"PtBcrJournalier.csv"), header = T)[,(1:2)]
PtStage$Date = ymd(PtStage$Date)
PtStage$H = as.numeric(PtStage$H)
# PtStage$H = PtStage$H + 3.37

#### stage restit
RestitStage = read.csv2(paste0(dir.data,"BcrResHorCor.csv"), header = T)[,(1:2)]
RestitStage$Date = ymd_hms(RestitStage$Date)
RestitStage$H = as.numeric(RestitStage$H)
# RestitStage$H = RestitStage$H - 0.06

### gau restit
RestitGau = read.csv2(paste0(dir.data,"JauBeaucaireCNR.csv"),header = T)

#### flood dates
date1840 = ymd("1840-11-03")
date1856 = ymd("1856-05-31")
date2003 = ymd_hms("2003-12-03 12:00:00")

delta = 5
ylimmax = 12
Sys.setlocale("LC_ALL","English")

gg1840 = ggplot(data = PtStage[(which(PtStage$Date==(date1840-delta)) : 
                                which(PtStage$Date==(date1840+delta))),], aes(x=Date, y = H))+
  geom_segment(aes(x=Date, xend=Date, y=0, yend=H))+
  geom_point(color="orange", size=3)+
  ylab(label = "Stage [m]")+
  coord_cartesian(ylim = c(0,ylimmax))+
  theme_bw(base_size = 20)+
  geom_text(aes(x = date1840-4, y = 12, label = "1840"),size = 10)

gg1856 = ggplot(data = PtStage[(which(PtStage$Date==(date1856-delta)) : 
                                which(PtStage$Date==(date1856+delta))),], aes(x=Date, y = H))+
  geom_segment(aes(x=Date, xend=Date, y=0, yend=H))+
  geom_point(color="orange", size=3)+
  ylab(label = "Stage [m]")+
  coord_cartesian(ylim = c(0,ylimmax))+
  geom_text(aes(x = date1856-4, y = 12, label = "1856"),size = 10)+
  theme_bw(base_size = 20)

gg2003 = ggplot(data = RestitStage[(which(RestitStage$Date==(date2003-(delta*86400))) : 
                                    which(RestitStage$Date==(date2003+(delta*86400)))),], aes(x=Date, y = H))+
  geom_segment(aes(x=Date, xend=Date, y=0, yend=H), size = .1, color = "black", lty = 2)+
  geom_point(color="royalblue", size=2)+
  ylab(label = "Stage [m]")+
  coord_cartesian(ylim = c(0,ylimmax))+
  theme_bw(base_size = 20)+
  geom_text(aes(x = date2003-(4*86400), y = 12, label = "2003"),size = 10)

ggarrange(gg1840,
          gg1856+theme(axis.title.y = element_blank()),
          gg2003+theme(axis.title.y = element_blank()),
          ncol = 3, nrow = 1, align = "hv")

ggsave(path = dir.plots,filename = "Stages.pdf", width = 16, height = 10)





