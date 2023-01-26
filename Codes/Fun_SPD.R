require(ggplot2);require(cowplot);require(lubridate);require(reshape2)
##### FONCTIONS POUR BaRatin SPD

### Fonctions Valou
propagate_Meyras<-function(b1,c1,Bw1,Cr1,g1,  # control 1: rectangular weir
                           b2,c2,Bc2,KS2,S02, # control 2: large rectangular channel
                           b3,c3,Bc3,KS3,S03, # control 3: large rectangular channel
                           d.g,d.l,           # global and local incremental changes
                           start              # starting vector
){
  #^******************************************************************************
  #^* PURPOSE: Monte-Carlo propagation to move from "physical" parameterization P1
  #^*          to "inference" parameterization P2:
  #^*          P1=(b1,c1,Bw1,Cr1,g1,  # 1st control: offset, exponent, weir width, discharge coeff., gravity
  #^*              b2,c2,Bc2,KS2,S02, # 2nd control: offset, exponent, channel width, Strickler, slope
  #^*              b3,c3,Bc3,KS3,S03, # 3rd control: offset, exponent, channel width, Strickler, slope
  #^*              d.g,               # incremental overall changes affecting both the weir and the channel (size nperiod-1)
  #^*              d.l)               # incremental local changes affecting only the weir (size nperiod-1)
  #^*          P2=(b1,       # PERIOD-SPECIFIC: parameters b1 for all periods (size nperiod)
  #^*              a1,c1,    # STATIC: coefficient and exponent of 1st control
  #^*              b2,       # PERIOD-SPECIFIC: parameters b2 for all periods (size nperiod)
  #^*              a2,c2,    # STATIC: coefficient and exponent of 2nd control
  #^*              b3,a3,c3) # STATIC: all parameters of 3rd control
  #^******************************************************************************
  #^* PROGRAMMERS: Benjamin Renard, Irstea Lyon
  #^******************************************************************************
  #^* CREATED/MODIFIED: 09/01/2018
  #^******************************************************************************
  #^* IN
  #^*    1. [numeric vector] b1, weir activation stage (1st control)
  #^*    2. [numeric vector] c1, exponent for a rectangualar weir (1st control)
  #^*    3. [numeric vector] Bw1, weir width (1st control)
  #^*    4. [numeric vector] Cr1, discharge coefficienr for the weir (1st control)
  #^*    5. [numeric vector] g1, gravity (1st control)
  #^*    6. [numeric vector] b2, main channel bed stage (2nd control)
  #^*    7. [numeric vector] c2, exponent for a large rectangular channel (2nd control)
  #^*    8. [numeric vector] Bc2, channel width (2nd control)
  #^*    9. [numeric vector] KS2, Strickler coefficient (2nd control)
  #^*    10. [numeric vector] S02, Channel slope (2nd control)
  #^*    11. [numeric vector] b3, floodway channel bed stage (3rd control)
  #^*    12. [numeric vector] c3, exponent for a large rectangular channel (3rd control)
  #^*    13. [numeric vector] Bc3, channel width (3rd control)
  #^*    14. [numeric vector] KS3, Strickler coefficient (3rd control)
  #^*    15. [numeric vector] S03, Channel slope (3rd control)
  #^*    16. [numeric matrix] d.g, incremental global changes; size nsim*(nperiod-1)
  #^*    17. [numeric matrix] d.l, incremental local changes; size nsim*(nperiod-1)
  #^*    18. [list] start, starting vector for all parameters above
  #^* OUT
  #^*    A list containing:
  #^*    1. [data.frame] sim: Monte-Carlo samples in "inference" parameterization P2   
  #^*    2. [numeric vector] start: starting point in "inference" parameterization P2   
  #^******************************************************************************
  #^* REF.: Mansanarez et al. (2019) Shift happens! Adjusting stage-discharge 
  #^*       rating curves to morphological changes at known times. Water
  #^*       Resources Research.
  #^******************************************************************************
  #^* TO DO: /
  #^******************************************************************************
  #^* COMMENTS: /
  #^******************************************************************************
  
  nsim=length(b1);nperiod=NCOL(d.g)+1
  # initialise parameters that will be computed by propagation
  a1=vector(mode='double',length=nsim)
  a2=vector(mode='double',length=nsim)
  a3=vector(mode='double',length=nsim)
  allb1=matrix(NA,nsim,nperiod)
  allb2=matrix(NA,nsim,nperiod)
  #------------------------------
  # Monte-Carlo propagation
  #------------------------------
  for(i in 1:nsim){
    # Compute a1, a2 and a3
    a1[i]=Cr1[i]*Bw1[i]*sqrt(2*g1[i])
    a2[i]=KS2[i]*Bc2[i]*sqrt(S02[i])
    a3[i]=KS3[i]*Bc3[i]*sqrt(S03[i])
    # compute b1 and b2 for all periods
    allb1[i,1]=b1[i]
    allb1[i,2:nperiod]=b1[i]-cumsum(d.g[i,]+d.l[i,])
    allb2[i,1]=b2[i]
    allb2[i,2:nperiod]=b2[i]-cumsum(d.g[i,])
  }
  #-------------------------------------------------
  # transform starting vector in parameterization P2
  #-------------------------------------------------
  allb1.s=vector(mode='double',length=nperiod)
  allb2.s=vector(mode='double',length=nperiod)
  # All b1's
  allb1.s[1]=start$b1
  allb1.s[2:nperiod]=start$b1-cumsum(start$d.g+start$d.l)
  # All b2's
  allb2.s[1]=start$b2
  allb2.s[2:nperiod]=start$b2-cumsum(start$d.g)
  # a's
  a1.s=start$Cr1*start$Bw1*sqrt(2*start$g1)
  a2.s=start$KS2*start$Bc2*sqrt(start$S02)
  a3.s=start$KS3*start$Bc3*sqrt(start$S03)
  # return results
  return(list(
    sim=data.frame(b1=allb1,a1=a1,c1=c1,b2=allb2,a2=a2,c2=c2,b3=b3,a3=a3,c3=c3),
    start=c(allb1.s,a1.s,start$c1,allb2.s,a2.s,start$c2,start$b3,a3.s,start$c3)
  )
  )
}

propagate_Wairau<-function(b1,c1,Bc1,KS1,S01, # control 1: low-flow channel
                           b2,c2,Bc2,KS2,S02, # control 2: high-flow channel
                           d.b1,d.Bc1,        # incremental changes on low-flow channel bed and width
                           start              # starting vector
){
  #^******************************************************************************
  #^* PURPOSE: Monte-Carlo propagation to move from "physical" parameterization P1
  #^*          to "inference" parameterization P2:
  #^*          P1=(b1,c1,Bc1,KS1,S01, # 1st control: offset, exponent, channel width, Strickler, slope
  #^*              b2,c2,Bc2,KS2,S02, # 2nd control: offset, exponent, channel width, Strickler, slope
  #^*              d.b1,              # incremental changes on low-flow channel bed (size nperiod-1)
  #^*              d.Bc1)             # incremental changes on low-flow channel width (size nperiod-1)
  #^*          P2=(b1,       # PERIOD-SPECIFIC: parameters b1 for all periods (size nperiod)
  #^*              a1,       # PERIOD-SPECIFIC: parameters a1 for all periods (size nperiod)
  #^*              c1,       # STATIC: exponent of 1st control
  #^*              b2,a2,c2) # STATIC: all parameters of 2nd control
  #^******************************************************************************
  #^* PROGRAMMERS: Benjamin Renard, Irstea Lyon
  #^******************************************************************************
  #^* CREATED/MODIFIED: 16/01/2018
  #^******************************************************************************
  #^* IN
  #^*    1. [numeric vector] b1, low-flow channel bed stage (1st control)
  #^*    2. [numeric vector] c1, exponent for a large rectangular channel (1st control)
  #^*    3. [numeric vector] Bc1, channel width (1st control)
  #^*    4. [numeric vector] KS1, Strickler coefficient (1st control)
  #^*    5. [numeric vector] S01, Channel slope (1st control)
  #^*    6. [numeric vector] b2, high-flow channel bed stage (2nd control)
  #^*    7. [numeric vector] c2, exponent for a large rectangular channel (2nd control)
  #^*    8. [numeric vector] Bc2, channel width (2nd control)
  #^*    9. [numeric vector] KS2, Strickler coefficient (2nd control)
  #^*    10. [numeric vector] S02, Channel slope (2nd control)
  #^*    11. [numeric matrix] d.b1, incremental changes on low-flow channel bed; size nsim*(nperiod-1)
  #^*    12. [numeric matrix] d.Bc1, incremental changes on low-flow channel width; size nsim*(nperiod-1)
  #^*    13. [list] start, starting vector for all parameters above
  #^* OUT
  #^*    A list containing:
  #^*    1. [data.frame] sim: Monte-Carlo samples in "inference" parameterization P2   
  #^*    2. [numeric vector] start: starting point in "inference" parameterization P2   
  #^******************************************************************************
  #^* REF.: Mansanarez et al. (2019) Shift happens! Adjusting stage-discharge 
  #^*       rating curves to morphological changes at known times. Water
  #^*       Resources Research.
  #^******************************************************************************
  #^* TO DO: /
  #^******************************************************************************
  #^* COMMENTS: /
  #^******************************************************************************
  
  nsim=length(b1);nperiod=NCOL(d.b1)+1
  #------------------------------
  # Monte-Carlo propagation
  #------------------------------
  a2=vector(mode='double',length=nsim)
  allb1=matrix(NA,nsim,nperiod)
  alla1=matrix(NA,nsim,nperiod)
  # Start Monte-Carlo propagation
  for(i in 1:nsim){
    # Compute a2
    a2[i]=KS2[i]*Bc2[i]*sqrt(S02[i])
    # compute b1 and a1 for all periods
    allb1[i,1]=b1[i]
    allb1[i,2:nperiod]=b1[i]-cumsum(d.b1[i,])
    alla1[i,1]=KS1[i]*Bc1[i]*sqrt(S01[i])
    alla1[i,2:nperiod]=alla1[i,1]/cumprod(d.Bc1[i,])
  }
  #-------------------------------------------------
  # transform starting vector in parameterization P2
  #-------------------------------------------------
  allb1.s=vector(mode='double',length=nperiod)
  alla1.s=vector(mode='double',length=nperiod)
  # All b1's
  allb1.s[1]=start$b1
  allb1.s[2:nperiod]=start$b1-cumsum(start$d.b1)
  # All a1's
  alla1.s[1]=start$KS1*start$Bc1*sqrt(start$S01)
  alla1.s[2:nperiod]=alla1.s[1]/cumprod(start$d.Bc1)
  # a2
  a2.s=start$KS2*start$Bc2*sqrt(start$S02)
  
  # return results
  return(list(
    sim=data.frame(b1=allb1,a1=alla1,c1=c1,b2=b2,a2=a2,c2=c2),
    start=c(allb1.s,alla1.s,start$c1,start$b2,a2.s,start$c2)
  )
  )
}

propagate_Wairau_V2<-function(b1,c1,Bc1,KS1,S01, # control 1: low-flow channel
                              b2,c2,Bc2,KS2,S02, # control 2: high-flow channel
                              d.b1,d.Bc1,d.b2,   # incremental changes on low-flow channel bed and width and high-flow channel bed
                              start              # starting vector
){
  #^******************************************************************************
  #^* PURPOSE: Monte-Carlo propagation to move from "physical" parameterization P1
  #^*          to "inference" parameterization P2:
  #^*          P1=(b1,c1,Bc1,KS1,S01, # 1st control: offset, exponent, channel width, Strickler, slope
  #^*              b2,c2,Bc2,KS2,S02, # 2nd control: offset, exponent, channel width, Strickler, slope
  #^*              d.b1,              # incremental changes on low-flow channel bed (size nperiod-1)
  #^*              d.Bc1,             # incremental changes on low-flow channel width (size nperiod-1)
  #^*              d.b2)              # incremental changes on high-flow channel bed (size nperiod-1)
  #^*          P2=(b1,       # PERIOD-SPECIFIC: parameters b1 for all periods (size nperiod)
  #^*              a1,       # PERIOD-SPECIFIC: parameters a1 for all periods (size nperiod)
  #^*              c1,       # STATIC: exponent of 1st control
  #^*              b2,       # PERIOD-SPECIFIC: parameters b2 for all periods (size nperiod)
  #^*              a2,c2)    # STATIC: a/c parameters of 2nd control
  #^******************************************************************************
  #^* PROGRAMMERS: Benjamin Renard, Irstea Lyon
  #^******************************************************************************
  #^* CREATED/MODIFIED: 16/01/2018
  #^******************************************************************************
  #^* IN
  #^*    1. [numeric vector] b1, low-flow channel bed stage (1st control)
  #^*    2. [numeric vector] c1, exponent for a large rectangular channel (1st control)
  #^*    3. [numeric vector] Bc1, channel width (1st control)
  #^*    4. [numeric vector] KS1, Strickler coefficient (1st control)
  #^*    5. [numeric vector] S01, Channel slope (1st control)
  #^*    6. [numeric vector] b2, high-flow channel bed stage (2nd control)
  #^*    7. [numeric vector] c2, exponent for a large rectangular channel (2nd control)
  #^*    8. [numeric vector] Bc2, channel width (2nd control)
  #^*    9. [numeric vector] KS2, Strickler coefficient (2nd control)
  #^*    10. [numeric vector] S02, Channel slope (2nd control)
  #^*    11. [numeric matrix] d.b1, incremental changes on low-flow channel bed; size nsim*(nperiod-1)
  #^*    12. [numeric matrix] d.Bc1, incremental changes on low-flow channel width; size nsim*(nperiod-1)
  #^*    13. [numeric matrix] d.b2, incremental changes on high-flow channel bed; size nsim*(nperiod-1)
  #^*    14. [list] start, starting vector for all parameters above
  #^* OUT
  #^*    A list containing:
  #^*    1. [data.frame] sim: Monte-Carlo samples in "inference" parameterization P2   
  #^*    2. [numeric vector] start: starting point in "inference" parameterization P2   
  #^******************************************************************************
  #^* REF.: Mansanarez et al. (2019) Shift happens! Adjusting stage-discharge 
  #^*       rating curves to morphological changes at known times. Water
  #^*       Resources Research.
  #^******************************************************************************
  #^* TO DO: /
  #^******************************************************************************
  #^* COMMENTS: /
  #^******************************************************************************
  
  nsim=length(b1);nperiod=NCOL(d.b1)+1
  #------------------------------
  # Monte-Carlo propagation
  #------------------------------
  a2=vector(mode='double',length=nsim)
  allb1=matrix(NA,nsim,nperiod)
  alla1=matrix(NA,nsim,nperiod)
  allb2=matrix(NA,nsim,nperiod)
  # Start Monte-Carlo propagation
  for(i in 1:nsim){
    # Compute a2
    a2[i]=KS2[i]*Bc2[i]*sqrt(S02[i])
    # compute b1, a1 and b2 for all periods
    allb1[i,1]=b1[i]
    allb1[i,2:nperiod]=b1[i]-cumsum(d.b1[i,])
    alla1[i,1]=KS1[i]*Bc1[i]*sqrt(S01[i])
    alla1[i,2:nperiod]=alla1[i,1]/cumprod(d.Bc1[i,])
    allb2[i,1]=b2[i]
    allb2[i,2:nperiod]=b2[i]-cumsum(d.b2[i,])
  }
  #-------------------------------------------------
  # transform starting vector in parameterization P2
  #-------------------------------------------------
  allb1.s=vector(mode='double',length=nperiod)
  alla1.s=vector(mode='double',length=nperiod)
  allb2.s=vector(mode='double',length=nperiod)
  # All b1's
  allb1.s[1]=start$b1
  allb1.s[2:nperiod]=start$b1-cumsum(start$d.b1)
  # All a1's
  alla1.s[1]=start$KS1*start$Bc1*sqrt(start$S01)
  alla1.s[2:nperiod]=alla1.s[1]/cumprod(start$d.Bc1)
  # All b2's
  allb2.s[1]=start$b2
  allb2.s[2:nperiod]=start$b2-cumsum(start$d.b2)
  # a2
  a2.s=start$KS2*start$Bc2*sqrt(start$S02)
  
  # return results
  return(list(
    sim=data.frame(b1=allb1,a1=alla1,c1=c1,b2=allb2,a2=a2,c2=c2),
    start=c(allb1.s,alla1.s,start$c1,allb2.s,a2.s,start$c2)
  )
  )
}

propagate_Rakaia<-function(b1,c1,Bc1,KS1,S01, # control 1: low-flow channel
                           b2,c2,Bc2,KS2,S02, # control 2: high-flow channel
                           d.b1,d.Bc1,        # incremental changes on low-flow channel bed and width
                           start              # starting vector
){
  #^******************************************************************************
  #^* PURPOSE: Monte-Carlo propagation to move from "physical" parameterization P1
  #^*          to "inference" parameterization P2:
  #^*          P1=(b1,c1,Bc1,KS1,S01, # 1st control: offset, exponent, channel width, Strickler, slope
  #^*              b2,c2,Bc2,KS2,S02, # 2nd control: offset, exponent, channel width, Strickler, slope
  #^*              d.b1,              # incremental changes on low-flow channel bed (size nperiod-1)
  #^*              d.Bc1)             # incremental changes on low-flow channel width (size nperiod-1)
  #^*          P2=(b1,       # PERIOD-SPECIFIC: parameters b1 for all periods (size nperiod)
  #^*              a1,       # PERIOD-SPECIFIC: parameters a1 for all periods (size nperiod)
  #^*              c1,       # STATIC: exponent of 1st control
  #^*              b2,a2,c2) # STATIC: all parameters of 2nd control
  #^******************************************************************************
  #^* PROGRAMMERS: Benjamin Renard, Irstea Lyon
  #^******************************************************************************
  #^* CREATED/MODIFIED: 16/01/2018
  #^******************************************************************************
  #^* IN
  #^*    1. [numeric vector] b1, low-flow channel bed stage (1st control)
  #^*    2. [numeric vector] c1, exponent for a large rectangular channel (1st control)
  #^*    3. [numeric vector] Bc1, channel width (1st control)
  #^*    4. [numeric vector] KS1, Strickler coefficient (1st control)
  #^*    5. [numeric vector] S01, Channel slope (1st control)
  #^*    6. [numeric vector] b2, high-flow channel bed stage (2nd control)
  #^*    7. [numeric vector] c2, exponent for a large rectangular channel (2nd control)
  #^*    8. [numeric vector] Bc2, channel width (2nd control)
  #^*    9. [numeric vector] KS2, Strickler coefficient (2nd control)
  #^*    10. [numeric vector] S02, Channel slope (2nd control)
  #^*    11. [numeric matrix] d.b1, incremental changes on low-flow channel bed; size nsim*(nperiod-1)
  #^*    12. [numeric matrix] d.Bc1, incremental changes on low-flow channel width; size nsim*(nperiod-1)
  #^*    13. [list] start, starting vector for all parameters above
  #^* OUT
  #^*    A list containing:
  #^*    1. [data.frame] sim: Monte-Carlo samples in "inference" parameterization P2   
  #^*    2. [numeric vector] start: starting point in "inference" parameterization P2   
  #^******************************************************************************
  #^* REF.: Mansanarez et al. (2019) Shift happens! Adjusting stage-discharge 
  #^*       rating curves to morphological changes at known times. Water
  #^*       Resources Research.
  #^******************************************************************************
  #^* TO DO: /
  #^******************************************************************************
  #^* COMMENTS: /
  #^******************************************************************************
  
  nsim=length(b1);nperiod=NCOL(d.b1)+1
  #------------------------------
  # Monte-Carlo propagation
  #------------------------------
  a2=vector(mode='double',length=nsim)
  allb1=matrix(NA,nsim,nperiod)
  alla1=matrix(NA,nsim,nperiod)
  # Start Monte-Carlo propagation
  for(i in 1:nsim){
    # Compute a2
    a2[i]=KS2[i]*Bc2[i]*sqrt(S02[i])
    # compute b1 and a1 for all periods
    allb1[i,1]=b1[i]
    allb1[i,2:nperiod]=b1[i]-cumsum(d.b1[i,])
    alla1[i,1]=KS1[i]*Bc1[i]*sqrt(S01[i])
    alla1[i,2:nperiod]=alla1[i,1]/cumprod(d.Bc1[i,])
  }
  #-------------------------------------------------
  # transform starting vector in parameterization P2
  #-------------------------------------------------
  allb1.s=vector(mode='double',length=nperiod)
  alla1.s=vector(mode='double',length=nperiod)
  # All b1's
  allb1.s[1]=start$b1
  allb1.s[2:nperiod]=start$b1-cumsum(start$d.b1)
  # All a1's
  alla1.s[1]=start$KS1*start$Bc1*sqrt(start$S01)
  alla1.s[2:nperiod]=alla1.s[1]/cumprod(start$d.Bc1)
  # a2
  a2.s=start$KS2*start$Bc2*sqrt(start$S02)
  
  # return results
  return(list(
    sim=data.frame(b1=allb1,a1=alla1,c1=c1,b2=b2,a2=a2,c2=c2),
    start=c(allb1.s,alla1.s,start$c1,start$b2,a2.s,start$c2)
  )
  )
}


fit<-function(sim, # Monte Carlo simulations produced by function propagate_XXX$sim
              margins=rep('Gaussian',NCOL(sim)), # marginal distributions
              names=paste('par',1:NCOL(sim),sep=''), # parameter names
              ncol=7,rowmax=4,plots = F # graphical parameters (here 7*4 parameters shown per figure)
){
  #^******************************************************************************
  #^* PURPOSE: fit a multivariate prior distribution to Monte-Carlo samples
  #^******************************************************************************
  #^* PROGRAMMERS: Benjamin Renard, Irstea Lyon
  #^******************************************************************************
  #^* CREATED/MODIFIED: 16/01/2018
  #^******************************************************************************
  #^* IN
  #^*    1. [numeric matrix] sim, Monte Carlo simulations (nsim*npar)
  #^*    2. [character vector] margins, marginal distributions (npar)
  #^*    3. [character vector] names, parameter names (npar)
  #^*    4. [integer] ncol, rowmax: graphical parameters controlling the number of subplots per figure
  #^* OUT
  #^*    A list containing:
  #^*    1. [character vector] margins, marginal distributions (same as the one used as input)
  #^*    2. [list of vectors] priorpar, prior parameters for each marginal prior    
  #^*    3. [numeric vector] start, starting points (e.g. prior mean)    
  #^*    4. [numeric matrix] cor, prior correlation (in Gaussian space)    
  #^******************************************************************************
  #^* REF.: Mansanarez et al. (2019) Shift happens! Adjusting stage-discharge 
  #^*       rating curves to morphological changes at known times. Water
  #^*       Resources Research.
  #^******************************************************************************
  #^* COMMENTS: only normal and lognormal margins available for the moment
  #^******************************************************************************
  #^* TO DO: can probably be generalized to any marginal distribution
  #^******************************************************************************
  
  # check sizes match and if so, initialize
  n=NCOL(sim)
  if(length(margins)!=n){
    message('Fatal:size mismatch [sim,margins]');return(NA)
  }
  priorpar=vector("list",n)
  transform=matrix(NA,nsim,n)
  # start computations for each margin
  for(i in 1:n){
    # marginal prior parameters
    priorpar[[i]]=switch(margins[i],
                         'Gaussian'={c(mean(sim[,i]),sd(sim[,i]))},
                         'LogNormal'={c(mean(log(sim[,i])),sd(log(sim[,i])))},
                         NA)
    # Transform simulations into Gaussian space to estimate the correlation of the Gaussian copula
    transform[,i]=switch(margins[i],
                         'Gaussian'={qnorm(pnorm(sim[,i],mean=priorpar[[i]][1],sd=priorpar[[i]][2]))},
                         'LogNormal'={qnorm(plnorm(sim[,i],meanlog=priorpar[[i]][1],sdlog=priorpar[[i]][2]))},
                         NA)
  }
  # prior correlation
  corel=cor(transform)
  # Verify that marginal fit is acceptable graphically
  if(plots == T){
    nrow=min(ceiling(n/ncol),rowmax)
    for(i in 1:n){
      if(i%%(nrow*ncol)==1) {X11();par(mfrow=c(nrow,ncol))}
      x=seq(min(sim[,i]),max(sim[,i]),,100)
      y=switch(margins[i],
               'Gaussian'={dnorm(x,mean=priorpar[[i]][1],sd=priorpar[[i]][2])},
               'LogNormal'={dlnorm(x,meanlog=priorpar[[i]][1],sdlog=priorpar[[i]][2])},
               NA)
      hist(sim[,i],breaks=50,freq=F,main=NULL,xlab=names[i],ylab='density')
      lines(x,y,col='red')
    }
  
    # show correlation matrix
    X11();image(1:n,1:n,corel)
  }
  # Return
  return(list(margins=margins,priorpar=priorpar,cor=corel))
}

writeConfigFilesSPD<-function(prior, # prior object produced by function 'fit'
                           start, # starting vector produced by function propagate_XXX$start
                           ncontrol, # number of hydraulic controls
                           nperiod,colperiod, # number of periods and column where it is stored in data file
                           isVar=rep(F,ncontrol*nperiod), # flag varying (period-specific) parameters
                           model='Config_Model.txt', # Model configuration file
                           correl='PriorCorrelation.txt', #prior correlation file
                           names.bac=c('b','a','c') # base name of the 3 parameters for one control
){
  #^******************************************************************************
  #^* PURPOSE: write configuration files used by BaM
  #^******************************************************************************
  #^* PROGRAMMERS: Benjamin Renard, Irstea Lyon
  #^******************************************************************************
  #^* CREATED/MODIFIED: 17/01/2018
  #^******************************************************************************
  #^* IN
  #^*    1. [list] prior, object produced by function 'fit'
  #^*    2. [numeric vector] start, starting vector produced by function propagate_XXX$start
  #^*    3. [integer] ncontrol, number of hydraulic controls
  #^*    4. [integer] nperiod, number of periods 
  #^*    5. [integer] colperiod, column where period is stored in data file
  #^*    6. [logical vector] isVar, flag varying (period-specific) parameters
  #^*    7. [character] model, Model configuration file
  #^*    8. [character] correl, prior correlation file
  #^*    9. [character vector]  base name of the 3 parameters for one control
  #^* OUT
  #^*    Nothing - just write files
  #^******************************************************************************
  #^* REF.: Mansanarez et al. (2019) Shift happens! Adjusting stage-discharge 
  #^*       rating curves to morphological changes at known times. Water
  #^*       Resources Research.
  #^******************************************************************************
  #^* COMMENTS: /
  #^******************************************************************************
  #^* TO DO: /
  #^******************************************************************************
  
  #---------------------------------------
  # Config_Model file
  #---------------------------------------
  # Start assembling values to be written in file and associated comments
  comment=c(
    'Model ID',
    'nX: number of input variables',
    'nY: number of output variables',
    'nPar: number of parameters theta'
  )
  val=list('BaRatinBAC',1,1,3*ncontrol)
  # specification for each parameter
  k=0;m=1
  for(i in 1:ncontrol){ # loop on each hydraulic control
    for(j in 1:length(names.bac)){ # loop on each b-a-c parameter of the control
      k=k+1
      # comments
      comment=c(comment,
                'Parameter name',
                'Initial guess',
                'Prior distribution',
                'Prior parameters')
      # parname
      pname=paste(names.bac[j],i,sep='')
      val=c(val,pname)
      if(!isVar[k]){ # Static parameter
        # starting point
        val=c(val,start[m])
        # prior distribution
        val=c(val,prior$margins[m])
        # prior parameters
        val=c(val,list(prior$priorpar[[m]]))
        m=m+1                
      } else { # period-specific parameter
        # starting point
        val=c(val,-666.666) # unused, can use a dummy number
        # prior distribution: use 'VAR'
        val=c(val,'VAR')
        # Name of the configuration file for this varying parameter
        fname=paste('Config_',pname,'_VAR.txt',sep='')
        val=c(val,fname)
        # Write this additional configuration file
        comment0=c(
          'Number of periods',
          'Column where period is written in data file'
        )
        val0=list(nperiod,colperiod)
        for(ii in 1:nperiod){ # one parameter block per period
          # comments
          comment0=c(comment0,
                     'Parameter name',
                     'Initial guess',
                     'Prior distribution',
                     'Prior parameters')
          # par name
          val0=c(val0,paste(pname,'_',ii,sep=''))
          # starting point
          val0=c(val0,start[m])
          # prior distribution
          val0=c(val0,prior$margins[m])
          # prior parameters
          val0=c(val0,list(prior$priorpar[[m]]))
          m=m+1
        }
        writeConfig(val0,comment0,getwd(),fname)
      }       
    }
  }  
  writeConfig(val,comment,getwd(),model)
  
  #---------------------------------------
  # Prior correlations
  #---------------------------------------
  # add 2 lines/columns, uncorrelated with all others, for structural uncertainty parameters gamma1 and gamma2
  n=NCOL(prior$cor)
  foo=matrix(0,n+2,n+2)
  foo[1:n,1:n]=prior$cor
  foo[n+1,n+1]=1
  foo[n+2,n+2]=1
  # write resulting matrix
  write(foo,file=correl,ncolumns=n+2)
}

writeConfig<-function(val,comment,dir,fname,addQuote=T){
  #^******************************************************************************
  #^* PURPOSE: Generic function to write a configuration file
  #^******************************************************************************
  #^* PROGRAMMERS: Benjamin Renard, Irstea Lyon
  #^******************************************************************************
  #^* CREATED/MODIFIED: 12/10/2018
  #^******************************************************************************
  #^* IN
  #^*    1. [list] val, list of values to be written in the configuration file 
  #^*    2. [string vector] comment, comments (same length as val) 
  #^*    3. [string] dir, destination directory 
  #^*    4. [string] fname, file name 
  #^*    5. [logical] addQuote, add double quotes to any string encountered in val? 
  #^* OUT
  #^*    1. nothing  
  #^******************************************************************************
  #^* REF.: 
  #^******************************************************************************
  #^* TO DO: 
  #^******************************************************************************
  #^* COMMENTS: Overwrites any pre-existing file. Directory created if needed
  #^******************************************************************************
  
  n=length(val)
  txt=vector("character",n)
  for (i in 1:n){
    foo=val[[i]]
    # add double quotes to strings
    if(is.character(foo) & addQuote){foo=addQuotes(foo)}
    # transform R logicals to Fortran logicals
    if(is.logical(foo)){foo=paste(ifelse(foo,'.true.','.false.'))}
    # stitch values and comment
    if(is.null(comment[i])){
      txt[i]=paste(foo,collapse=',')
    } else {
      txt[i]=paste(paste(foo,collapse=','),'!',comment[i])
    }
  }
  
  if(!dir.exists(dir)){dir.create(dir,recursive=T)}
  file=file.path(dir,fname)
  
  write.table(matrix(txt, ncol = 1, nrow = length(txt)), file = file,
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}

addQuotes<-function(txt){
  #^******************************************************************************
  #^* PURPOSE: add double quotes to a string or to each element of a string vector
  #^******************************************************************************
  #^* PROGRAMMERS: Benjamin Renard, Irstea Lyon
  #^******************************************************************************
  #^* CREATED/MODIFIED: 12/10/2018
  #^******************************************************************************
  #^* IN
  #^*    1. [string vector] txt
  #^* OUT
  #^*    1. [string vector] double-quoted string vector 
  #^******************************************************************************
  #^* REF.: 
  #^******************************************************************************
  #^* TO DO: 
  #^******************************************************************************
  #^* COMMENTS: 
  #^******************************************************************************
  
  n=length(txt)
  out=txt
  for(i in 1:n){
    out[i]=paste('"',txt[i],'"',sep='')
  }
  return(out)
}


propagate_Bcr_b1b2<-function(b1,c1,Bc1,KS1,S01, # control 1: low-flow channel
                              b2,c2,Bc2,KS2,S02, # control 2: high-flow channel
                              d.g,d.l,   # incremental changes on low-flow channel bed and high-flow channel bed
                              start              # starting vector
){

  #^* OUT
  #^*    A list containing:
  #^*    1. [data.frame] sim: Monte-Carlo samples in "inference" parameterization P2   
  #^*    2. [numeric vector] start: starting point in "inference" parameterization P2   
  #^******************************************************************************
  #^* REF.: Mansanarez et al. (2019) Shift happens! Adjusting stage-discharge 
  #^*       rating curves to morphological changes at known times. Water
  #^*       Resources Research.
  #^******************************************************************************
  #^* TO DO: /
  #^******************************************************************************
  #^* COMMENTS: /
  #^******************************************************************************
  
  nsim=length(b1);nperiod=NCOL(d.g)+1
  #------------------------------
  # Monte-Carlo propagation
  #------------------------------
  a1=vector(mode='double',length=nsim)
  a2=vector(mode='double',length=nsim)
  allb1=matrix(NA,nsim,nperiod)
  allb2=matrix(NA,nsim,nperiod)
  # Start Monte-Carlo propagation
  for(i in 1:nsim){
    # Compute a's
    a1[i]=KS1[i]*Bc1[i]*sqrt(S01[i])
    a2[i]=KS2[i]*Bc2[i]*sqrt(S02[i])
    # compute b1 and b2 for all periods
    allb1[i,1]=b1[i]
    allb1[i,2:nperiod]=b1[i]-cumsum(d.g[i,]+d.l[i,])
    allb2[i,1]=b2[i]
    allb2[i,2:nperiod]=b2[i]-cumsum(d.g[i,])
  }
  #-------------------------------------------------
  # transform starting vector in parameterization P2
  #-------------------------------------------------
  allb1.s=vector(mode='double',length=nperiod)
  allb2.s=vector(mode='double',length=nperiod)
  # All b1's
  allb1.s[1]=start$b1
  allb1.s[2:nperiod]=start$b1-cumsum(start$d.g+start$d.l)
  # All b2's
  allb2.s[1]=start$b2
  allb2.s[2:nperiod]=start$b2-cumsum(start$d.g)
  # a's
  a1.s=start$KS1*start$Bc1*sqrt(start$S01)
  a2.s=start$KS2*start$Bc2*sqrt(start$S02)
  
  # return results
  return(list(
    sim=data.frame(b1=allb1,a1=a1,c1=c1,b2=allb2,a2=a2,c2=c2),
    start=c(allb1.s,a1.s,start$c1,allb2.s,a2.s,start$c2)
  )
  )
}

propagate_Bcr_b1b2reversed<-function(b1,c1,Bc1,KS1,S01, # control 1: low-flow channel
                             b2,c2,Bc2,KS2,S02, # control 2: high-flow channel
                             d.g,d.l,   # incremental changes on low-flow channel bed and high-flow channel bed
                             start              # starting vector
){
  
  #^* OUT
  #^*    A list containing:
  #^*    1. [data.frame] sim: Monte-Carlo samples in "inference" parameterization P2   
  #^*    2. [numeric vector] start: starting point in "inference" parameterization P2   
  #^******************************************************************************
  #^* REF.: Mansanarez et al. (2019) Shift happens! Adjusting stage-discharge 
  #^*       rating curves to morphological changes at known times. Water
  #^*       Resources Research.
  #^******************************************************************************
  #^* TO DO: /
  #^******************************************************************************
  #^* COMMENTS: /
  #^******************************************************************************
  
  nsim=length(b1);nperiod=NCOL(d.g)+1
  #------------------------------
  # Monte-Carlo propagation
  #------------------------------
  a1=vector(mode='double',length=nsim)
  a2=vector(mode='double',length=nsim)
  allb1=matrix(NA,nsim,nperiod)
  allb2=matrix(NA,nsim,nperiod)
  # Start Monte-Carlo propagation
  for(i in 1:nsim){
    # Compute a's
    a1[i]=KS1[i]*Bc1[i]*sqrt(S01[i])
    a2[i]=KS2[i]*Bc2[i]*sqrt(S02[i])
    # compute b1 and b2 for all periods
    allb1[i,1]=b1[i]
    allb1[i,2:nperiod]=b1[i]-cumsum(d.g[i,]+d.l[i,])
    allb2[i,1]=b2[i]
    allb2[i,2:nperiod]=b2[i]-cumsum(d.g[i,])
  }
  
  ## REVERSE THE ORDER
  allb1 = allb1[,(nperiod:1)]
  allb2 = allb2[,(nperiod:1)]
  #-------------------------------------------------
  # transform starting vector in parameterization P2
  #-------------------------------------------------
  allb1.s=vector(mode='double',length=nperiod)
  allb2.s=vector(mode='double',length=nperiod)
  # All b1's
  allb1.s[1]=start$b1
  allb1.s[2:nperiod]=start$b1-cumsum(start$d.g+start$d.l)
  # All b2's
  allb2.s[1]=start$b2
  allb2.s[2:nperiod]=start$b2-cumsum(start$d.g)
  # a's
  a1.s=start$KS1*start$Bc1*sqrt(start$S01)
  a2.s=start$KS2*start$Bc2*sqrt(start$S02)
  
  # return results
  return(list(
    sim=data.frame(b1=allb1,a1=a1,c1=c1,b2=allb2,a2=a2,c2=c2),
    start=c(allb1.s,a1.s,start$c1,allb2.s,a2.s,start$c2)
  )
  )
}


propagate_BcrRes_b1b2b3<-function(b1,c1,Bc1,KS1,S01, # control 1: low-flow channel
                                  b2,c2,Bc2,KS2,S02, # control 2: mid-flow channel
                                  b3,c3,Bc3,KS3,S03, # 3 : high flow
                                  d.g,d.l,   # incremental changes on low-flow channel bed and high-flow channel bed
                                  start              # starting vector
){
  
  #^* OUT
  #^*    A list containing:
  #^*    1. [data.frame] sim: Monte-Carlo samples in "inference" parameterization P2   
  #^*    2. [numeric vector] start: starting point in "inference" parameterization P2   
  #^******************************************************************************
  #^* REF.: Mansanarez et al. (2019) Shift happens! Adjusting stage-discharge 
  #^*       rating curves to morphological changes at known times. Water
  #^*       Resources Research.

  
  nsim=length(b1);nperiod=NCOL(d.g)+1
  #------------------------------
  # Monte-Carlo propagation
  #------------------------------
  a1=vector(mode='double',length=nsim)
  a2=vector(mode='double',length=nsim)
  a3=vector(mode='double',length=nsim)
  allb1=matrix(NA,nsim,nperiod)
  allb2=matrix(NA,nsim,nperiod)
  allb3=matrix(NA,nsim,nperiod)
  # Start Monte-Carlo propagation
  for(i in 1:nsim){
    # Compute a's
    a1[i]=KS1[i]*Bc1[i]*sqrt(S01[i])
    a2[i]=KS2[i]*Bc2[i]*sqrt(S02[i])
    a3[i]=KS3[i]*Bc3[i]*sqrt(S03[i])
    # compute b1 and b2 for all periods
    allb1[i,1]=b1[i]
    allb1[i,2:nperiod]=b1[i]-cumsum(d.g[i,]+d.l[i,])
    allb2[i,1]=b2[i]
    allb2[i,2:nperiod]=b2[i]-cumsum(d.g[i,])
    allb3[i,1]=b3[i]
    allb3[i,2:nperiod]=b3[i]-cumsum(d.g[i,])
    
  }
  #-------------------------------------------------
  # transform starting vector in parameterization P2
  #-------------------------------------------------
  allb1.s=vector(mode='double',length=nperiod)
  allb2.s=vector(mode='double',length=nperiod)
  allb3.s=vector(mode='double',length=nperiod)
  # All b1's
  allb1.s[1]=start$b1
  allb1.s[2:nperiod]=start$b1-cumsum(start$d.g+start$d.l)
  # All b2's
  allb2.s[1]=start$b2
  allb2.s[2:nperiod]=start$b2-cumsum(start$d.g)
  # All b3's
  allb3.s[1]=start$b3
  allb3.s[2:nperiod]=start$b3-cumsum(start$d.g)
  # a's
  a1.s=start$KS1*start$Bc1*sqrt(start$S01)
  a2.s=start$KS2*start$Bc2*sqrt(start$S02)
  a3.s=start$KS3*start$Bc3*sqrt(start$S03)
  
  # return results
  return(list(
    sim=data.frame(b1=allb1,a1=a1,c1=c1,b2=allb2,a2=a2,c2=c2,b3=allb3,a3=a3,c3=c3),
    start=c(allb1.s,a1.s,start$c1,allb2.s,a2.s,start$c2,allb3.s,a3.s,start$c3)
  )
  )
}



propagate_BcrRes_varb1b2_b3<-function(b1,c1,Bc1,KS1,S01, # control 1: low-flow channel
                                      b2,c2,Bc2,KS2,S02, # control 2: mid-flow channel
                                      b3,c3,Bc3,KS3,S03, # 3 : high flow
                                      d.g,d.l,   # incremental changes on low-flow channel bed and high-flow channel bed
                                      start){              # starting vector

  
    nsim=length(b1);nperiod=NCOL(d.g)+1
    # initialise parameters that will be computed by propagation
    a1=vector(mode='double',length=nsim)
    a2=vector(mode='double',length=nsim)
    a3=vector(mode='double',length=nsim)
    allb1=matrix(NA,nsim,nperiod)
    allb2=matrix(NA,nsim,nperiod)
    #------------------------------
    # Monte-Carlo propagation
    #------------------------------
    for(i in 1:nsim){
      # Compute a1, a2 and a3
      a1[i]=KS1[i]*Bc1[i]*sqrt(S01[i])
      a2[i]=KS2[i]*Bc2[i]*sqrt(S02[i])
      a3[i]=KS3[i]*Bc3[i]*sqrt(S03[i])
      # compute b1 and b2 for all periods
      allb1[i,1]=b1[i]
      allb1[i,2:nperiod]=b1[i]-cumsum(d.g[i,]+d.l[i,])
      allb2[i,1]=b2[i]
      allb2[i,2:nperiod]=b2[i]-cumsum(d.g[i,])
    }
    #-------------------------------------------------
    # transform starting vector in parameterization P2
    #-------------------------------------------------
    allb1.s=vector(mode='double',length=nperiod)
    allb2.s=vector(mode='double',length=nperiod)
    # All b1's
    allb1.s[1]=start$b1
    allb1.s[2:nperiod]=start$b1-cumsum(start$d.g+start$d.l)
    # All b2's
    allb2.s[1]=start$b2
    allb2.s[2:nperiod]=start$b2-cumsum(start$d.g)
    # a's
    a1.s=start$KS1*start$Bc1*sqrt(start$S01)
    a2.s=start$KS2*start$Bc2*sqrt(start$S02)
    a3.s=start$KS3*start$Bc3*sqrt(start$S03)
    # return results
    return(list(
      sim=data.frame(b1=allb1,a1=a1,c1=c1,b2=allb2,a2=a2,c2=c2,b3=b3,a3=a3,c3=c3),
      start=c(allb1.s,a1.s,start$c1,allb2.s,a2.s,start$c2,start$b3,a3.s,start$c3)
    )
    )
  }

propagate_BcrRes_varb1b2_b3reverse<-function(b1,c1,Bc1,KS1,S01, # control 1: low-flow channel
                                      b2,c2,Bc2,KS2,S02, # control 2: mid-flow channel
                                      b3,c3,Bc3,KS3,S03, # 3 : high flow
                                      d.g,d.l,   # incremental changes on low-flow channel bed and high-flow channel bed
                                      start){              # starting vector
  
  
  nsim=length(b1);nperiod=NCOL(d.g)+1
  # initialise parameters that will be computed by propagation
  a1=vector(mode='double',length=nsim)
  a2=vector(mode='double',length=nsim)
  a3=vector(mode='double',length=nsim)
  allb1=matrix(NA,nsim,nperiod)
  allb2=matrix(NA,nsim,nperiod)
  #------------------------------
  # Monte-Carlo propagation
  #------------------------------
  for(i in 1:nsim){
    # Compute a1, a2 and a3
    a1[i]=KS1[i]*Bc1[i]*sqrt(S01[i])
    a2[i]=KS2[i]*Bc2[i]*sqrt(S02[i])
    a3[i]=KS3[i]*Bc3[i]*sqrt(S03[i])
    # compute b1 and b2 for all periods
    allb1[i,1]=b1[i]
    allb1[i,2:nperiod]=b1[i]-cumsum(d.g[i,]+d.l[i,])
    allb2[i,1]=b2[i]
    allb2[i,2:nperiod]=b2[i]-cumsum(d.g[i,])
  }
  
  allb1 = allb1[,nperiod:1]
  allb2 = allb2[,nperiod:1]
  
  #-------------------------------------------------
  # transform starting vector in parameterization P2
  #-------------------------------------------------
  allb1.s=vector(mode='double',length=nperiod)
  allb2.s=vector(mode='double',length=nperiod)
  # All b1's
  allb1.s[1]=start$b1
  allb1.s[2:nperiod]=start$b1-cumsum(start$d.g+start$d.l)
  # All b2's
  allb2.s[1]=start$b2
  allb2.s[2:nperiod]=start$b2-cumsum(start$d.g)
  # a's
  a1.s=start$KS1*start$Bc1*sqrt(start$S01)
  a2.s=start$KS2*start$Bc2*sqrt(start$S02)
  a3.s=start$KS3*start$Bc3*sqrt(start$S03)
  # return results
  return(list(
    sim=data.frame(b1=allb1,a1=a1,c1=c1,b2=allb2,a2=a2,c2=c2,b3=b3,a3=a3,c3=c3),
    start=c(allb1.s,a1.s,start$c1,allb2.s,a2.s,start$c2,start$b3,a3.s,start$c3)
  )
  )
}

propagate_BcrRes_localb1b2_reverse<-function(b1,c1,Bc1,KS1,S01, # control 1: low-flow channel
                                             b2,c2,Bc2,KS2,S02, # control 2: mid-flow channel
                                             b3,c3,Bc3,KS3,S03, # 3 : high flow
                                             d.g,d.l,   # incremental changes on low-flow channel bed and high-flow channel bed
                                             start){              # starting vector
  
  
  nsim=length(b1);nperiod=NCOL(d.l)+1
  # initialise parameters that will be computed by propagation
  a1=vector(mode='double',length=nsim)
  a2=vector(mode='double',length=nsim)
  a3=vector(mode='double',length=nsim)
  allb1=matrix(NA,nsim,nperiod)
  allb2=matrix(NA,nsim,nperiod)
  #------------------------------
  # Monte-Carlo propagation
  #------------------------------
  for(i in 1:nsim){
    # Compute a1, a2 and a3
    a1[i]=KS1[i]*Bc1[i]*sqrt(S01[i])
    a2[i]=KS2[i]*Bc2[i]*sqrt(S02[i])
    a3[i]=KS3[i]*Bc3[i]*sqrt(S03[i])
    # compute b1 and b2 for all periods
    allb1[i,1]=b1[i]
    allb1[i,2:nperiod]=b1[i]-cumsum(d.l[i,])
    allb2[i,1]=b2[i]
    allb2[i,2:nperiod]=b2[i]-cumsum(d.l[i,])
  }
  
  allb1 = allb1[,nperiod:1]
  allb2 = allb2[,nperiod:1]
  
  #-------------------------------------------------
  # transform starting vector in parameterization P2
  #-------------------------------------------------
  allb1.s=vector(mode='double',length=nperiod)
  allb2.s=vector(mode='double',length=nperiod)
  # All b1's
  allb1.s[1]=start$b1
  allb1.s[2:nperiod]=start$b1-cumsum(start$d.l)
  # All b2's
  allb2.s[1]=start$b2
  allb2.s[2:nperiod]=start$b2-cumsum(start$d.l)
  # a's
  a1.s=start$KS1*start$Bc1*sqrt(start$S01)
  a2.s=start$KS2*start$Bc2*sqrt(start$S02)
  a3.s=start$KS3*start$Bc3*sqrt(start$S03)
  # return results
  return(list(
    sim=data.frame(b1=allb1,a1=a1,c1=c1,b2=allb2,a2=a2,c2=c2,b3=b3,a3=a3,c3=c3),
    start=c(allb1.s,a1.s,start$c1,allb2.s,a2.s,start$c2,start$b3,a3.s,start$c3)
  )
  )
}

###############################################################################
Transf_Gauss_lognorm = function(E, stdev){
  ###############################################################################  
  # IN:
  # E = mean of the Gaussian distribution
  # stdev = standanrd deviartion of the gaussian distribution
  # OUT:
  # mu = mean of the Lognormal distribution
  # sd = standard deviation of the lognormal distribution
  mu = log((E^2)/(stdev^2+E^2)^0.5)
  sd = (log(1+ (stdev^2)/(E^2)))^0.5
  return(list(mu = mu,
              sd = sd))
}


###########################################################################################
BaRatin_SPD_config <- function(dir.BaM, dir.SPD.config, dir.spd.short, pred, nobs, M,     #
                               remnant, g1.prior, g2.prior,                               #
                               g1.distr.type, g2.distr.type,                              # 
                               Ncycles, Hmax=10) {                                        #
  ###########################################################################################
  # creation of Config_BaM.txt
  file.BaM <- paste(dir.BaM,"/Config_BaM.txt",sep="")
  # cat('"BaRatin_SPD/"', file = file.BaM, sep="\n", append = FALSE)
  cat(paste0("'",dir.spd.short,"'"), file = file.BaM, sep="\n", append = FALSE)
  cat('"Config_RunOptions.txt"', file = file.BaM , sep="\n", append = TRUE)
  cat('"Config_Model.txt"', file = file.BaM , sep="\n", append = TRUE)
  cat('"Config_ControlMatrix.txt"', file = file.BaM , sep="\n", append = TRUE)
  cat('"Config_Data.txt"', file = file.BaM , sep="\n", append = TRUE)
  cat('"Config_RemnantSigma.txt"', file = file.BaM , sep="\n", append = TRUE)
  cat('"Config_MCMC.txt"', file = file.BaM , sep="\n", append = TRUE)
  cat('"Config_Cooking.txt"', file = file.BaM , sep="\n", append = TRUE)
  cat('"Config_Summary.txt"', file = file.BaM , sep="\n", append = TRUE)
  cat('"Config_Residuals.txt"', file = file.BaM , sep="\n", append = TRUE)
  if (pred == TRUE) {
    cat('"Config_Pred_Master.txt"', file = file.BaM , sep="\n", append = TRUE)
  } else {
    cat('""', file = file.BaM , sep="\n", append = TRUE)
  }
  # creation of Config_Data.txt
  file.name2 = paste0(dir.SPD.config,"/Config_Data.txt")
  # cat("'BaRatin_SPD\\Gaugings_data.txt'", file =file.name2, sep="\n")
  cat(paste0("'",dir.spd.short,"Gaugings_data.txt'"), file =file.name2, sep="\n")
  cat(1, file = file.name2, append = TRUE,sep="\n")
  cat(nobs, file = file.name2, append = TRUE,sep="\n")
  cat(4, file =file.name2, append = TRUE,sep="\n")
  cat(1, file =file.name2, append = TRUE,sep="\n")
  cat(0, file =file.name2, append = TRUE,sep="\n")
  cat(0, file =file.name2, append = TRUE,sep="\n")
  cat(0, file =file.name2, append = TRUE,sep="\n")
  cat(2, file =file.name2, append = TRUE,sep="\n")
  cat(3, file =file.name2, append = TRUE,sep="\n")
  cat(0, file =file.name2, append = TRUE,sep="\n")
  cat(0, file =file.name2, append = TRUE,sep="\n")
  # creation of Config_MCMC.txt
  file.mcmc = paste(dir.SPD.config,"/Config_MCMC.txt",sep="")
  cat('"Results_MCMC.txt"', file =file.mcmc,sep="\n")
  cat(100, file = file.mcmc, append = TRUE,sep="\n")   #Nadapt
  cat(Ncycles, file = file.mcmc, append = TRUE,sep="\n")  #Ncycles
  cat(0.1, file =file.mcmc, append = TRUE,sep="\n")    #minMoveRate
  cat(0.5, file =file.mcmc, append = TRUE,sep="\n")    #maxMoveRate
  cat(0.9, file =file.mcmc, append = TRUE,sep="\n")    #DownMult
  cat(1.1, file =file.mcmc, append = TRUE,sep="\n")    #UpMult
  cat(0, file =file.mcmc, append = TRUE,sep="\n")      #mode for init jump distr
  cat("****", file =file.mcmc, append = TRUE,sep="\n") 
  cat(0.1, file =file.mcmc, append = TRUE,sep="\n")    #MultFact
  cat(0.1,file =file.mcmc, append = TRUE, sep=",")     #RC MultiFact
  cat(0.1,file =file.mcmc, append = TRUE, sep=",")     
  cat(0.1, file =file.mcmc, append = TRUE,sep="\n")
  cat(0.1,file =file.mcmc, append = TRUE, sep=",")      #Remnant MultiFact
  cat(0.1, file =file.mcmc, append = TRUE,sep="\n")
  # creation of Config_ControlMatrix.txt
  file.matrix = paste(dir.SPD.config,"/Config_ControlMatrix.txt",sep="")
  write.table(M, file =file.matrix, row.names = FALSE, col.names = FALSE)   #M control matrix
  cat(paste0(Hmax,"."), file = file.matrix,  append = TRUE, sep="\n")  #hmax
  # creation of Config_RemnantSigma.txt
  file.remnant = paste(dir.SPD.config,"/Config_RemnantSigma.txt",sep="")
  if (remnant == "Linear") {
    cat("'Linear'", file = file.remnant, sep="\n")                    #! Function f used in sdev=f(Qrc) 
    cat(2, file = file.remnant, append = TRUE, sep="\n")              #! Number of parameters gamma for f
    cat("gamma1", file = file.remnant, append = TRUE, sep="\n")       #! Parameter Name
    cat(1, file = file.remnant, append = TRUE, sep="\n")              #! Initial Guess
    cat('Uniform', file = file.remnant, append = TRUE, sep="\n")      #! Prior distribution
    cat(0,file =file.remnant, append = TRUE, sep=",")
    cat(",",file =file.remnant, append = TRUE, sep=",")
    cat(1000,file =file.remnant, append = TRUE, sep="\n")
    #
    cat("gamma2", file = file.remnant, append = TRUE, sep="\n")       #! Parameter Name
    cat(0.1, file = file.remnant, append = TRUE, sep="\n")            #! Initial Guess
    cat('Uniform', file = file.remnant, append = TRUE, sep="\n")      #! Initial Guess
    cat(0,file =file.remnant, append = TRUE, sep=",")
    cat(",",file =file.remnant, append = TRUE, sep=",")
    cat(100,file =file.remnant, append = TRUE, sep="\n")
  } else {
    cat("'Constant'", file = file.remnant, sep="\n")                     #! Function f used in sdev=f(Qrc) 
    cat(1, file = file.remnant, append = TRUE, sep="\n")                 #! Number of parameters gamma for f
    cat("gamma1", file = file.remnant, append = TRUE, sep="\n")          #! Parameter Name
    cat(1, file = file.remnant, append = TRUE, sep="\n")                 #! Initial Guess
    cat('Uniform', file = file.remnant, append = TRUE, sep="\n")         #! Prior distribution
    cat(0,file =file.remnant, append = TRUE, sep=",")
    cat(",",file =file.remnant, append = TRUE, sep=",")
    cat(1000,file =file.remnant, append = TRUE, sep="\n")
  }
  # creation of Config_RunOptions.txt
  file.run = paste(dir.SPD.config,"/Config_RunOptions.txt",sep="")
  cat(".true.", file = file.run, sep="\n")                             # Do MCMC?
  cat(".true.", file = file.run, append = TRUE, sep="\n")              # Do MCMC summary?
  cat(".true.", file = file.run, append = TRUE, sep="\n")              # Do Residual diagnostics?
  if (pred == TRUE) {
    cat(".true.", file = file.run, append = TRUE, sep="\n")         # Do Predictions?
  } else {
    cat(".false.", file = file.run, append = TRUE, sep="\n")        # Do Predictions?
  }
  # creation of Config_Residuals.txt
  file.residuals = paste(dir.SPD.config,"/Config_Residuals.txt",sep="")
  cat('"Results_Residuals.txt"' , file =file.residuals ,sep="\n")     # Result file
  # creation of Config_Summary.txt
  file.summary = paste(dir.SPD.config,"/Config_Summary.txt", sep="")
  cat('"Results_Summary.txt"' , file =file.summary, sep="\n")         # Result file
  # creation of Config_Cooking.txt
  file.cooking = paste(dir.SPD.config,"/Config_Cooking.txt", sep="")
  cat('"Results_MCMC_Cooked.txt"' , file =file.cooking, sep="\n")     # Result file
  cat(0.5, file =file.cooking, append = TRUE, sep="\n")               # Burn factor
  cat(10, file =file.cooking, append = TRUE, sep="\n")                # Nslim
}



MCMCplot<-function(MCMCfile,
                   doPar, #=T or F
                   doLogPost, # =T or F
                   doDPar,  # =T or F
                   type,  #="trace", # "histogram","density","scatterplot"
                   xlab, # =''
                   ylab, #=''
                   ncol, # number of columns in the plot
                   prior, # = NULL if traceplot type
                   burn,  # = 0 , fraction of mcmc to burn
                   slim,   # = 1, mcmc to slim
                   theme_size){
  #################################################################################################
  # Plot MCMC samples
  #------------------------------------------------------------
  # read MCMC simulation, burn and slim if required
  X = read.table(MCMCfile, header=T)
  X = X[seq(burn+1, nrow(X), slim), ]
  # determine logpot column
  lpcol=which(names(X)=='LogPost')
  # determine which columns should be kept
  keep=c()
  if(doPar){keep=c(keep,1:(lpcol-1))}
  if(doLogPost){keep=c(keep,lpcol)}
  if(doDPar){if(ncol(X)>lpcol){keep=c(keep,(lpcol+1):ncol(X))}}
  X = X[,keep]
  if(type=="scatterplot"){
    m=ggpairs(X,diag=list(continuous='density'),
              lower=list(continuous='points'),upper=list(continuous='cor'),axisLabels='show')
  } else {
    X=cbind(X, indx = 1:nrow(X))
    Xm=melt(X, id.vars = 'indx')
    m=ggplot(Xm)
    m=m+switch(type,
               trace     = geom_line(aes(x = indx, y = value, colour = variable)),
               histogram = geom_histogram(aes(value, fill = variable, ..density..)),
               density   = geom_density(aes(value, fill = variable), colour=NA, alpha = 0.8))
    m=m+facet_wrap(~variable, scales='free', ncol=ncol)+
      xlab(xlab)+
      ylab(ylab)+
      theme(legend.position="none")
    
    
    # Plot priors if required
    if(doPar & (type=='density' | type=='histogram') & !is.null(prior)){
      #priorDF = data.frame(matrix(NA, length(prior), 100))
      y= matrix(NA, 100, length(prior))
      grid = matrix(NA, 100, length(prior)) 
      for(i in 1:length(prior)){
        grid[,i]=seq(min(X[,i]), max(X[,i]),length.out=100)
        if(prior[[i]]$prior$dist=='Gaussian'){
          y[,i]=dnorm(grid[,i],mean=prior[[i]]$prior$par[1],sd=prior[[i]]$prior$par[2])
        }
        if(prior[[i]]$prior$dist=='LogNormal'){
          y[,i]=dlnorm(grid[,i], meanlog = prior[[i]]$prior$par[1], sdlog =prior[[i]]$prior$par[2])
        }
        if(prior[[i]]$prior$dist=='Uniform'){
          y[,i]=dunif(grid[,i],min=prior[[i]]$prior$par[1],max=prior[[i]]$prior$par[2])
        }
        if(prior[[i]]$prior$dist=='FlatPrior'){
          y[,i]=dunif(grid[,i],min=-99999999,max=99999999)
        }
        if(prior[[i]]$prior$dist=='FlatPrior+'){
          y[,i]=dunif(grid[,i],min=0,max=99999999)
        }
        if(prior[[i]]$prior$dist=='FlatPrior-'){
          y[,i]= dunif(grid[,i],min=-99999999,max=0)
        }
        if(prior[[i]]$prior$dist=='Gaussian'  | prior[[i]]$prior$dist=='Uniform'   |
           prior[[i]]$prior$dist=='FlatPrior' | prior[[i]]$prior$dist=='FlatPrior+' | prior[[i]]$prior$dist=='FlatPrior-' |
           prior[[i]]$prior$dist=='LogNormal'){
          # priorDF = rbind(priorDF, data.frame(variable = prior[[i]]$name,
          #                                     x        = grid,
          #                                     y        = y,
          #                                     ind      = i))
        }
      }
      
      priorDF = data.frame(y)
      priorDF.grid = data.frame(grid)
      names(priorDF) = unlist(lapply(prior, '[[', 1))
      names(priorDF.grid) = unlist(lapply(prior, '[[', 1))
      priorDF.grid.bis    = cbind(priorDF.grid, indx = 1:nrow(priorDF.grid))
      priorDF.grid.bis.m  = melt(priorDF.grid.bis, id.vars = 'indx')
      priorDF.bis    = cbind(priorDF, indx = 1:nrow(priorDF))
      priorDF.bis.m  = cbind(melt(priorDF.bis, id.vars = 'indx'), grid = priorDF.grid.bis.m$value)
      
      m=ggplot(Xm)
      m=m+switch(type,
                 trace     = geom_line(aes(x = indx, y = value, colour = variable)),
                 histogram = geom_histogram( aes(value, fill = variable, ..density..)),
                 density   = geom_density(aes(value, fill = variable), colour=NA, alpha = 0.8)) +
        geom_line(data = priorDF.bis.m, aes(x=grid, y=value), colour="black")
      m=m+
        facet_wrap(~variable, scales='free', ncol= ncol)+
        xlab(xlab)+
        ylab(ylab)+
        theme(legend.position="none")
      
    }
  }
  m= m +
    theme_bw(base_size      = theme_size)+
    theme(legend.position   = "none"
          ,axis.text        = element_text(size=10)
          ,plot.margin      = unit(c(1, 0.5, 0.5, 0.5),"cm")
          ,panel.grid.minor = element_blank())
  return(m)
}




#################################################################################################
qscan<-function(f,k,sep=' '){
  #################################################################################################
  # quick scan - just a convenience wrapper
  return(scan(f,what='character',skip=k,nlines=1,quiet=TRUE,sep=sep))
}




ReadModel<-function(ModelFile){
  #------------------------------------------------------------
  # Read information on the fitted model in Config_Model.txt
  #------------------------------------------------------------
  # read model info
  k=0;ID=qscan(ModelFile,k)[1]
  k=k+1;nX=as.numeric(qscan(ModelFile,k)[1])
  k=k+1;nY=as.numeric(qscan(ModelFile,k)[1])
  k=k+1;nPar=as.numeric(qscan(ModelFile,k)[1])
  # parameter blocs
  par=list();m=0
  for(i in 1:nPar){
    k=k+1;name=qscan(ModelFile,k)[1]
    k=k+1;init=as.numeric(qscan(ModelFile,k)[1])
    k=k+1;dist=qscan(ModelFile,k)[1]
    k=k+1;foo=qscan(ModelFile,k,sep='?')
    if(dist=='Gaussian' | dist=='Uniform' | dist=='FlatPrior' | dist=='FlatPrior+' | dist=='LogNormal' | dist=='FlatPrior-'){
      ppar=as.numeric(strsplit(gsub(pattern=",",replacement=" ",foo),split=" ")[[1]][1:2])
    } else {
      ppar=c(NA,NA)
    }        
    if(dist!='FIX'){
      m=m+1;par[[m]]=list(name=name, init=init, prior=list(dist=dist, par=ppar))
    }
  }
  return(list(ID=ID,nX=nX,nY=nY,nPar=nPar,par=par))
}

#########################################################################################################
BaRatin_configNspag <- function(dir, dir_code,
                                nsim, 
                                propagat,
                                b.distr, a.distr, c.distr,  
                                a.prior, st_a.prior, 
                                c.prior, st_c.prior, 
                                b.prior, st_b.prior, #priors
                                Bw.prior, Cr.prior, g.prior, 
                                Bc.prior, KS.prior, S0.prior,
                                st_Bw.prior, st_Cr.prior, st_g.prior,
                                st_Bc.prior, st_KS.prior, st_S0.prior,
                                ncontrol, 
                                M,Hmax=10,
                                nobs, 
                                Ncycles,
                                ngrid,
                                nlimni, 
                                predictionRC, 
                                predictionQt, 
                                predictionPrior,
                                simMCMC, 
                                mcmc.prior,
                                remnant.err.model, g1.prior, g2.prior, g1.distr.type, g2.distr.type,
                                nspag) {
  #########################################################################################################
  #Hmax = 10  ##### CAREFUL !!!!!!!!!!!!! PARAM
  
  
  #prior progation of RC parameters:
  prior = BaRatin.propagation(dir, 
                              nsim, 
                              propagat, 
                              b.distr,     a.distr,     c.distr,  
                              a.prior,     st_a.prior,  c.prior, 
                              st_c.prior,  b.prior,     st_b.prior,
                              Bw.prior,    Cr.prior,    g.prior,
                              Bc.prior,    KS.prior,    S0.prior,
                              st_Bw.prior, st_Cr.prior, st_g.prior, 
                              st_Bc.prior, st_KS.prior, st_S0.prior,
                              ncontrol, 
                              M)
  
  writeConfigFiles(prior     = prior[[1]],        # Prior object produced by function 'fit'
                   start     = prior[[2]],        # starting vector produced by function propagate_XXX$start
                   ncontrol,                      # number of hydraulic controls
                   model     ='Config_Model.txt', # Model configuration file
                   names.bac = c('b','a','c'),    # base name of the 3 parameters for one control
                   dir)
  
  npar = 3*ncontrol
  
  #**********************************************************************************************   DATA
  file.name2 = paste0(dir,"/Config_Data.txt")
  cat("'BaM_BaRatin_2\\Gaugings_data.txt'", file =file.name2,sep="\n")
  cat(1,      file = file.name2, append = TRUE,sep="\n")
  cat(nobs,   file = file.name2, append = TRUE,sep="\n")      #nobs in the file  
  cat(4,      file = file.name2, append = TRUE,sep="\n")      #number of columns
  cat(1,      file = file.name2, append = TRUE,sep="\n")      #column with the input variable obs
  cat(0,      file = file.name2, append = TRUE,sep="\n")      
  cat(0,      file = file.name2, append = TRUE,sep="\n")
  cat(0,      file = file.name2, append = TRUE,sep="\n")
  cat(2,      file = file.name2, append = TRUE,sep="\n")      #column with output variable obs
  cat(3,      file = file.name2, append = TRUE,sep="\n")      #column with uQ (=uQrel/100*Q*0.5)
  cat(0,      file = file.name2, append = TRUE,sep="\n")
  cat(0,      file = file.name2, append = TRUE,sep="\n")
  ###################################################################   MCMC
  file.mcmc = paste(dir,"/Config_MCMC.txt",sep="")
  cat('"Results_MCMC.txt"', file =file.mcmc,sep="\n")
  cat(100, file = file.mcmc, append = TRUE,sep="\n")   #Nadapt
  cat(Ncycles, file = file.mcmc, append = TRUE,sep="\n")  #Ncycles
  cat(0.1, file =file.mcmc, append = TRUE,sep="\n")    #minMoveRate
  cat(0.5, file =file.mcmc, append = TRUE,sep="\n")    #maxMoveRate
  cat(0.9, file =file.mcmc, append = TRUE,sep="\n")    #DownMult
  cat(1.1, file =file.mcmc, append = TRUE,sep="\n")    #UpMult
  cat(0, file =file.mcmc, append = TRUE,sep="\n")      #mode for init jump distr
  cat("****", file =file.mcmc, append = TRUE,sep="\n") 
  cat(0.1, file =file.mcmc, append = TRUE,sep="\n")    #MultFact
  cat(0.1,file =file.mcmc, append = TRUE, sep=",")     #RC MultiFact
  cat(0.1,file =file.mcmc, append = TRUE, sep=",")     
  cat(0.1, file =file.mcmc, append = TRUE,sep="\n")
  cat(0.1,file =file.mcmc, append = TRUE, sep=",")      #Remnant MultiFact
  cat(0.1, file =file.mcmc, append = TRUE,sep="\n")
  ###################################################################   MATRIX
  file.matrix = paste(dir,"/Config_ControlMatrix.txt",sep="")
  write.table(M, file =file.matrix, row.names = FALSE, col.names = FALSE)
  cat(Hmax, file =file.matrix, append = TRUE,sep="\n")
  #---------------------------------------------------------------------- 
  file.remnant = paste(dir,"/Config_RemnantSigma.txt",sep="")
  if (remnant.err.model == "Linear") {
    cat("'Linear'", file = file.remnant, sep="\n")                         #! Function f used in sdev=f(Qrc) 
    cat(2, file = file.remnant, append = TRUE, sep="\n")                    #! Number of parameters gamma for f
    cat("gamma1", file = file.remnant, append = TRUE, sep="\n")             #! Parameter Name
    cat(g1.prior[3], file = file.remnant, append = TRUE, sep="\n")                   #! Initial Guess
    cat(g1.distr.type, file = file.remnant, append = TRUE, sep="\n")            #! Prior distribution
    cat(g1.prior[1],file =file.remnant, append = TRUE, sep=",")
    cat(",",file =file.remnant, append = TRUE, sep=",")
    cat(g1.prior[2],file =file.remnant, append = TRUE, sep="\n")
    
    cat("gamma2", file = file.remnant, append = TRUE, sep="\n")             #! Parameter Name
    cat(g2.prior[3], file = file.remnant, append = TRUE, sep="\n")                   #! Initial Guess
    cat(g2.distr.type, file = file.remnant, append = TRUE, sep="\n")             #! Initial Guess
    cat(g2.prior[1],file =file.remnant, append = TRUE, sep=",")
    cat(",",file =file.remnant, append = TRUE, sep=",")
    cat(g2.prior[2],file =file.remnant, append = TRUE, sep="\n")
    
  } else if (remnant.err.model == "Constant") {
    cat("'Constant'", file = file.remnant, sep="\n")                       #! Function f used in sdev=f(Qrc) 
    cat(1, file = file.remnant, append = TRUE, sep="\n")                    #! Number of parameters gamma for f
    cat("gamma1", file = file.remnant, append = TRUE, sep="\n")             #! Parameter Name
    cat(g1.prior[3], file = file.remnant, append = TRUE, sep="\n")                   #! Initial Guess
    cat(g1.distr.type, file = file.remnant, append = TRUE, sep="\n")            #! Prior distribution
    cat(g1.prior[1],file =file.remnant, append = TRUE, sep=",")
    cat(",",file =file.remnant, append = TRUE, sep=",")
    cat(g1.prior[2],file =file.remnant, append = TRUE, sep="\n")
  }
  ###################################################################  Config_BaM
  file_BaM <- paste(dir_code,"/Config_BaM.txt",sep="")
  #creation of Config_BaM.txt
  cat('"BaM_BaRatin_2/"', file =file_BaM , sep="\n", append = FALSE)
  cat('"Config_RunOptions.txt"', file = file_BaM , sep="\n", append = TRUE)    
  cat('"Config_Model.txt"', file = file_BaM , sep="\n", append = TRUE)
  cat('"Config_ControlMatrix.txt"', file = file_BaM , sep="\n", append = TRUE)
  cat('"Config_Data.txt"', file = file_BaM , sep="\n", append = TRUE)
  cat('"Config_RemnantSigma.txt"', file = file_BaM , sep="\n", append = TRUE)                                       
  cat('"Config_MCMC.txt"', file = file_BaM , sep="\n", append = TRUE)                                            
  cat('"Config_Cooking.txt"', file = file_BaM , sep="\n", append = TRUE)
  cat('"Config_Summary.txt"', file = file_BaM , sep="\n", append = TRUE)
  cat('"Config_Residuals.txt"', file = file_BaM , sep="\n", append = TRUE)
  if( (predictionRC ==TRUE) | (predictionQt ==TRUE) |(predictionPrior ==TRUE) ) {
    cat('"Config_Pred_Master.txt"', file = file_BaM , sep="\n", append = TRUE)
  } else {
    cat('""', file = file_BaM , sep="\n", append = TRUE)
  }
  ###########################################################################################
  # PREDICTIONS :
  ###############
  file.Pred1 = paste(dir,"/Config_Pred_Master.txt",sep="")
  if ((predictionPrior == TRUE)&(predictionRC== TRUE)&(predictionQt == TRUE))  {
    cat(8, file =file.Pred1,sep="\n")
  } else if ((predictionPrior == TRUE) &(predictionRC == FALSE)&(predictionQt ==FALSE ))  {
    cat(2, file =file.Pred1,sep="\n")
  } else if ((predictionPrior == TRUE) & (predictionRC == TRUE)&(predictionPrior == FALSE)) {
    cat(5, file =file.Pred1,sep="\n")
  } else if ((predictionPrior == TRUE) & (predictionRC == FALSE) &(predictionQt == TRUE)) {
    cat(5, file =file.Pred1,sep="\n")
  } else if ((predictionPrior == FALSE)& (predictionRC == TRUE)&(predictionQt == TRUE)) {
    cat(6, file =file.Pred1,sep="\n")
  } else if ((predictionPrior == FALSE)& (predictionRC == FALSE)&(predictionQt == TRUE)) {
    cat(3, file =file.Pred1,sep="\n")
  } else if ((predictionPrior == FALSE)& (predictionRC == TRUE)&(predictionQt == FALSE)) {
    cat(3, file =file.Pred1,sep="\n")
  }
  
  ###########################
  if(predictionPrior==TRUE) {
    ###########################
    cat("'Config_Pred_Prior.txt'", file =file.Pred1, append = TRUE,sep="\n")
    cat("'Config_Pred_Prior_Qt.txt'", file =file.Pred1, append = TRUE,sep="\n")
    file.Pred11 = paste0(dir,"/Config_Pred_Prior.txt")
    cat("'BaM_BaRatin_2\\Hgrid.txt'", file =file.Pred11, sep="\n")
    cat(ngrid, file =file.Pred11, append = TRUE, sep="\n")
    cat("1", file = file.Pred11, append = TRUE,sep="\n")   #n of spaghetti
    cat(".true.", file = file.Pred11, append = TRUE,sep="\n")                              #!!! Propagate parametric uncertainty?
    cat(".false.", file = file.Pred11, append = TRUE,sep="\n")                             #!!! Propagate remnant uncertainty for each output variable? (size nY)
    cat(mcmc.prior, file = file.Pred11, append = TRUE,sep="\n")                          #!!! Nsim[prior]. If <=0: posterior sampling (nsim is given by mcmc sample); 
    #cat(-1, file = file.Pred11, append = TRUE,sep="\n")                                    #!!! Nsim[prior]. If <=0: posterior sampling (nsim is given by mcmc sample); 
    cat("'Qrc_Prior.spag'", file = file.Pred11, append = TRUE,sep="\n")                    # !!! Files containing spaghettis for each output variable (size nY)
    cat(".true.", file = file.Pred11, append = TRUE,sep="\n")                              #!!! Post-processing: transpose spag file (so that each column is a spaghetti)? 
    cat(".true.", file = file.Pred11, append = TRUE,sep="\n")                              #!!! Post-processing: create envelops? (size nY)
    cat("'Qrc_Prior.env'", file = file.Pred11, append = TRUE,sep="\n")                     # !!! Post-processing: name of envelop files (size nY)
    cat(".false.", file = file.Pred11, append = TRUE,sep="\n")                              #!!! Print progress in console during computations?
    cat(".false." , file = file.Pred11, append = TRUE,sep="\n")                           #!!! Do state prediction? (size nState)
    
    ###################################################################
    file.Pred21 = paste0(dir,"/Config_Pred_Prior_Qt.txt")
    cat("'BaM_BaRatin_2\\limni.txt'", file =file.Pred21,sep="\n")
    cat(nlimni, file =file.Pred21, append = TRUE, sep="\n")
    cat(nspag, file = file.Pred21, append = TRUE,sep="\n")   #n of spaghetti
    cat(".true.", file = file.Pred21, append = TRUE,sep="\n")                              #!!! Propagate parametric uncertainty?
    cat(".false.", file = file.Pred21, append = TRUE,sep="\n")                             #!!! Propagate remnant uncertainty for each output variable? (size nY)
    cat(mcmc.prior, file = file.Pred21, append = TRUE,sep="\n")                                #!!! Nsim[prior]. If <=0: posterior sampling (nsim is given by mcmc sample); 
    cat("'Qt_prior.spag'", file = file.Pred21, append = TRUE,sep="\n")                    # !!! Files containing spaghettis for each output variable (size nY)
    cat(".true.", file = file.Pred21, append = TRUE,sep="\n")                              #!!! Post-processing: transpose spag file (so that each column is a spaghetti)? 
    cat(".true.", file = file.Pred21, append = TRUE,sep="\n")                              #!!! Post-processing: create envelops? (size nY)
    cat("'Qt_prior.env'", file = file.Pred21, append = TRUE,sep="\n")                     # !!! Post-processing: name of envelop files (size nY)
    cat(".false.", file = file.Pred21, append = TRUE,sep="\n")                              #!!! Print progress in console during computations?
    cat(".false." , file = file.Pred21, append = TRUE,sep="\n")                           #!!! Do state prediction? (size nState)
  }
  
  
  ###################################################################
  if(predictionRC==TRUE) {  
    ###################################################################
    #cat("'Config_Pred_Prior.txt'", file = file.Pred1, append = TRUE,sep="\n")
    cat("'Config_Pred_RCMaxpost.txt'", file = file.Pred1, append = TRUE,sep="\n")
    cat("'Config_Pred_RCParamU.txt'", file =file.Pred1, append = TRUE,sep="\n")
    cat("'Config_Pred_RCTotalU.txt'", file =file.Pred1, append = TRUE,sep="\n")
    ##################################################################
    file.Pred3 = paste(dir,"/Config_Pred_RCMaxpost.txt",sep="")
    cat("'BaM_BaRatin_2\\Hgrid.txt'", file =file.Pred3, sep="\n")
    cat(ngrid, file =file.Pred3,sep="\n", append = TRUE)
    cat("1", file = file.Pred3, append = TRUE,sep="\n")   #n of spaghetti
    cat(".false.", file = file.Pred3, append = TRUE,sep="\n")                              #!!! Propagate parametric uncertainty?
    cat(".false.", file = file.Pred3, append = TRUE,sep="\n")                             #!!! Propagate remnant uncertainty for each output variable? (size nY)
    cat("-1", file = file.Pred3, append = TRUE,sep="\n")                                #!!! Nsim[prior]. If <=0: posterior sampling (nsim is given by mcmc sample); 
    cat("'Qrc_Maxpost.spag'", file = file.Pred3, append = TRUE,sep="\n")                    # !!! Files containing spaghettis for each output variable (size nY)
    cat(".true.", file = file.Pred3, append = TRUE,sep="\n")                              #!!! Post-processing: transpose spag file (so that each column is a spaghetti)? 
    cat(".true.", file = file.Pred3, append = TRUE,sep="\n")                              #!!! Post-processing: create envelops? (size nY)
    cat("'Qrc_Maxpost.env'", file = file.Pred3, append = TRUE,sep="\n")                     # !!! Post-processing: name of envelop files (size nY)
    cat(".false.", file = file.Pred3, append = TRUE,sep="\n")                              #!!! Print progress in console during computations?
    cat(".false." , file = file.Pred3, append = TRUE,sep="\n")                           #!!! Do state prediction? (size nState)
    ###################################################################
    file.Pred4 = paste(dir,"/Config_Pred_RCParamU.txt",sep="")
    cat("'BaM_BaRatin_2\\Hgrid.txt'", file =file.Pred4, sep="\n")
    cat(ngrid, file =file.Pred4,sep="\n", append = TRUE)
    cat("1", file = file.Pred4, append = TRUE,sep="\n")                                   #n of spaghetti
    cat(".true.", file = file.Pred4, append = TRUE,sep="\n")                              #!!! Propagate parametric uncertainty?
    cat(".false.", file = file.Pred4, append = TRUE,sep="\n")                             #!!! Propagate remnant uncertainty for each output variable? (size nY)
    cat("-1", file = file.Pred4, append = TRUE,sep="\n")                                  #!!! Nsim[prior]. If <=0: posterior sampling (nsim is given by mcmc sample); 
    cat("'Qrc_ParamU.spag'", file = file.Pred4, append = TRUE,sep="\n")                    # !!! Files containing spaghettis for each output variable (size nY)
    cat(".true.", file = file.Pred4, append = TRUE,sep="\n")                              #!!! Post-processing: transpose spag file (so that each column is a spaghetti)? 
    cat(".true.", file = file.Pred4, append = TRUE,sep="\n")                              #!!! Post-processing: create envelops? (size nY)
    cat("'Qrc_ParamU.env'", file = file.Pred4, append = TRUE,sep="\n")                     # !!! Post-processing: name of envelop files (size nY)
    cat(".false.", file = file.Pred4, append = TRUE,sep="\n")                              #!!! Print progress in console during computations?
    cat(".false." , file = file.Pred4, append = TRUE,sep="\n")                            #!!! Do state prediction? (size nState)
    ###################################################################
    file.Pred5 = paste(dir,"/Config_Pred_RCTotalU.txt",sep="")
    cat("'BaM_BaRatin_2\\Hgrid.txt'", file =file.Pred5,sep="\n")
    cat(ngrid, file =file.Pred5, append = TRUE, sep="\n")
    cat("1", file = file.Pred5, append = TRUE,sep="\n")                                  #n of spaghetti
    cat(".true.", file = file.Pred5, append = TRUE,sep="\n")                              #!!! Propagate parametric uncertainty?
    cat(".true.", file = file.Pred5, append = TRUE,sep="\n")                             #!!! Propagate remnant uncertainty for each output variable? (size nY)
    cat("-1", file = file.Pred5, append = TRUE,sep="\n")                                #!!! Nsim[prior]. If <=0: posterior sampling (nsim is given by mcmc sample); 
    cat("'Qrc_TotalU.spag'", file = file.Pred5, append = TRUE,sep="\n")                    # !!! Files containing spaghettis for each output variable (size nY)
    cat(".true.", file = file.Pred5, append = TRUE,sep="\n")                              #!!! Post-processing: transpose spag file (so that each column is a spaghetti)? 
    cat(".true.", file = file.Pred5, append = TRUE,sep="\n")                              #!!! Post-processing: create envelops? (size nY)
    cat("'Qrc_TotalU.env'", file = file.Pred5, append = TRUE,sep="\n")                     # !!! Post-processing: name of envelop files (size nY)
    cat(".false.", file = file.Pred5, append = TRUE,sep="\n")                              #!!! Print progress in console during computations?
    cat(".false." , file = file.Pred5, append = TRUE,sep="\n")                           #!!! Do state prediction? (size nState)
  }
  
  
  #############################################################################
  if(predictionQt==TRUE) {
    #############################################################################
    cat("'Config_Pred_Maxpost.txt'", file =file.Pred1, append = TRUE,sep="\n")
    cat("'Config_Pred_ParamU.txt'", file =file.Pred1, append = TRUE,sep="\n")
    cat("'Config_Pred_TotalU.txt'", file =file.Pred1, append = TRUE,sep="\n")
    
    file.Pred6 = paste0(dir,"/Config_Pred_Maxpost.txt")
    cat("'BaM_BaRatin_2\\limni.txt'", file =file.Pred6,sep="\n")
    cat(nlimni, file =file.Pred6, append = TRUE, sep="\n")
    cat(nspag, file = file.Pred6, append = TRUE,sep="\n")   #n of spaghetti
    cat(".false.", file = file.Pred6, append = TRUE,sep="\n")                              #!!! Propagate parametric uncertainty?
    cat(".false.", file = file.Pred6, append = TRUE,sep="\n")                             #!!! Propagate remnant uncertainty for each output variable? (size nY)
    cat("-1", file = file.Pred6, append = TRUE,sep="\n")                                #!!! Nsim[prior]. If <=0: posterior sampling (nsim is given by mcmc sample); 
    cat("'Qt_Maxpost.spag'", file = file.Pred6, append = TRUE,sep="\n")                    # !!! Files containing spaghettis for each output variable (size nY)
    cat(".true.", file = file.Pred6, append = TRUE,sep="\n")                              #!!! Post-processing: transpose spag file (so that each column is a spaghetti)? 
    cat(".true.", file = file.Pred6, append = TRUE,sep="\n")                              #!!! Post-processing: create envelops? (size nY) 
    cat("'Qt_Maxpost.env'", file = file.Pred6, append = TRUE,sep="\n")                     # !!! Post-processing: name of envelop files (size nY)
    cat(".false.", file = file.Pred6, append = TRUE,sep="\n")                              #!!! Print progress in console during computations?
    cat(".false." , file = file.Pred6, append = TRUE,sep="\n")                           #!!! Do state prediction? (size nState)
    
    
    file.Pred8 = paste0(dir,"/Config_Pred_ParamU.txt")
    cat("'BaM_BaRatin_2\\limni.txt'", file =file.Pred8,sep="\n")
    cat(nlimni, file =file.Pred8, append = TRUE, sep="\n")
    cat(nspag, file = file.Pred8, append = TRUE,sep="\n")   #n of spaghetti
    cat(".true.", file = file.Pred8, append = TRUE,sep="\n")                              #!!! Propagate parametric uncertainty?
    cat(".false.", file = file.Pred8, append = TRUE,sep="\n")                             #!!! Propagate remnant uncertainty for each output variable? (size nY)
    cat("-1", file = file.Pred8, append = TRUE,sep="\n")                                #!!! Nsim[prior]. If <=0: posterior sampling (nsim is given by mcmc sample); 
    cat("'Qt_ParamU.spag'", file = file.Pred8, append = TRUE,sep="\n")                    # !!! Files containing spaghettis for each output variable (size nY)
    cat(".true.", file = file.Pred8, append = TRUE,sep="\n")                              #!!! Post-processing: transpose spag file (so that each column is a spaghetti)? 
    cat(".true.", file = file.Pred8, append = TRUE,sep="\n")                              #!!! Post-processing: create envelops? (size nY)
    cat("'Qt_ParamU.env'", file = file.Pred8, append = TRUE,sep="\n")                     # !!! Post-processing: name of envelop files (size nY)
    cat(".false.", file = file.Pred8, append = TRUE,sep="\n")                              #!!! Print progress in console during computations?
    cat(".false." , file = file.Pred8, append = TRUE,sep="\n")                           #!!! Do state prediction? (size nState)
    
    
    file.Pred7 = paste(dir,"/Config_Pred_TotalU.txt",sep="")
    cat("'BaM_BaRatin_2\\limni.txt'", file =file.Pred7,sep="\n")
    cat(nlimni, file =file.Pred7, append = TRUE, sep="\n")
    cat(nspag, file = file.Pred7, append = TRUE,sep="\n")   #n of spaghetti
    cat(".true.", file = file.Pred7, append = TRUE,sep="\n")                              #!!! Propagate parametric uncertainty?
    cat(".true.", file = file.Pred7, append = TRUE,sep="\n")                             #!!! Propagate remnant uncertainty for each output variable? (size nY)
    cat("-1", file = file.Pred7, append = TRUE,sep="\n")                                #!!! Nsim[prior]. If <=0: posterior sampling (nsim is given by mcmc sample); 
    cat("'Qt_TotalU.spag'", file = file.Pred7, append = TRUE,sep="\n")                    # !!! Files containing spaghettis for each output variable (size nY)
    cat(".true.", file = file.Pred7, append = TRUE,sep="\n")                              #!!! Post-processing: transpose spag file (so that each column is a spaghetti)? 
    cat(".true.", file = file.Pred7, append = TRUE,sep="\n")                              #!!! Post-processing: create envelops? (size nY)
    cat("'Qt_TotalU.env'", file = file.Pred7, append = TRUE,sep="\n")                     # !!! Post-processing: name of envelop files (size nY)
    cat(".false.", file = file.Pred7, append = TRUE,sep="\n")                              #!!! Print progress in console during computations?
    cat(".false." , file = file.Pred7, append = TRUE,sep="\n")                           #!!! Do state prediction? (size nState)
  }
  
  
  
  ##########################################################################
  # RUN OPTIONS 
  ##########################################################################
  file.run = paste0(dir.exe,"/Config_RunOptions.txt")
  if (simMCMC == TRUE) {  
    cat(".true.", file =file.run,sep="\n")   # Do the MCMC simulations
  } else {
    cat(".false.", file =file.run,sep="\n")  # Do not do the MCMC simulation    
  }
  cat(".true.", file =file.run, append = TRUE, sep="\n")    # Do MCMC summary?
  cat(".true.", file =file.run, append = TRUE, sep="\n")    # Do Residual diagnostics?
  if ((predictionRC == TRUE) | (predictionQt == TRUE) | (predictionPrior ==TRUE)) {
    cat(".true.", file =file.run, append = TRUE, sep="\n")    # Do Predictions?
  } else {
    cat(".false.", file =file.run, append = TRUE, sep="\n")   # Do Predictions?
  }
  
  ###################################################################   RESIDUALS CONFIG
  file.residuals = paste(dir,"/Config_Residuals.txt",sep="")
  cat('"Results_Residuals.txt"' , file =file.residuals ,sep="\n")      # Result file
  
  ###################################################################   SUMMARY CONFIG
  file.summary = paste(dir,"/Config_Summary.txt",sep="")
  cat('"Results_Summary.txt"' , file =file.summary ,sep="\n")    #Summary stat results file name
  
  ###################################################################   COOKING CONFIG
  file.cooking = paste(dir,"/Config_Cooking.txt",sep="")
  cat('"Results_MCMC_Cooked.txt"' , file =file.cooking ,sep="\n")  # mcmc Results file name
  cat(0.5, file =file.cooking, append = TRUE, sep="\n")            # Burn factor
  cat(10, file =file.cooking, append = TRUE, sep="\n")             # Nslim
}




#########################################################################################################
BaRatin_configNspagPeriod <- function( dir, period, nperiod,
                                  nsim, 
                                  prior,
                                  ncontrol, 
                                  M,
                                  nobs, 
                                  Ncycles,
                                  ngrid,
                                  nlimni, 
                                  predictionRC, 
                                  predictionQt, 
                                  predictionPrior,
                                  simMCMC, 
                                  mcmc.prior,
                                  remnant.err.model, g1.prior, g2.prior, g1.distr.type, g2.distr.type,
                                  nspag,Hmax, dir_code, dir.bam.short) {
  
  #########################################################################################################
  writeConfigFiles(prior     = prior,        # Prior object produced by function 'fit'
                   start     = prior[[2]],        # starting vector produced by function propagate_XXX$start
                   ncontrol,                      # number of hydraulic controls
                   model     ='Config_Model.txt', # Model configuration file
                   names.bac = c('b','a','c'),    # base name of the 3 parameters for one control
                   dir)#,isVar = rep(F,6),nperiod = 1)
  
  npar = 3*ncontrol

  #**********************************************************************************************   DATA
  # file.name2 = paste0(dir_code,"/BaM_exe/BaM_BaRatin_2/Config_Data.txt")
  # cat("'BaM_BaRatin_2\\Gaugings_data.txt'", file =file.name2,sep="\n")
  file.name2 = paste0(dir,"/Config_Data.txt")
  cat(paste0("'",dir.bam.short,"Gaugings_data.txt'"), file =file.name2, sep="\n")
  cat(1,      file = file.name2, append = TRUE,sep="\n")
  cat(nobs,   file = file.name2, append = TRUE,sep="\n")      #nobs in the file  
  cat(4,      file = file.name2, append = TRUE,sep="\n")      #number of columns
  cat(1,      file = file.name2, append = TRUE,sep="\n")      #column with the input variable obs
  cat(0,      file = file.name2, append = TRUE,sep="\n")      
  cat(0,      file = file.name2, append = TRUE,sep="\n")
  cat(0,      file = file.name2, append = TRUE,sep="\n")
  cat(2,      file = file.name2, append = TRUE,sep="\n")      #column with output variable obs
  cat(3,      file = file.name2, append = TRUE,sep="\n")      #column with uQ (=uQrel/100*Q*0.5)
  cat(0,      file = file.name2, append = TRUE,sep="\n")
  cat(0,      file = file.name2, append = TRUE,sep="\n")
  ###################################################################   MCMC
  # file.mcmc = paste(dir_code,"/BaM_exe/BaM_BaRatin_2/Config_MCMC.txt",sep="")
  file.mcmc = paste(dir,"/Config_MCMC.txt",sep="")
  cat('"Results_MCMC.txt"', file =file.mcmc,sep="\n")
  cat(100, file = file.mcmc, append = TRUE,sep="\n")   #Nadapt
  cat(Ncycles, file = file.mcmc, append = TRUE,sep="\n")  #Ncycles
  cat(0.1, file =file.mcmc, append = TRUE,sep="\n")    #minMoveRate
  cat(0.5, file =file.mcmc, append = TRUE,sep="\n")    #maxMoveRate
  cat(0.9, file =file.mcmc, append = TRUE,sep="\n")    #DownMult
  cat(1.1, file =file.mcmc, append = TRUE,sep="\n")    #UpMult
  cat(0, file =file.mcmc, append = TRUE,sep="\n")      #mode for init jump distr
  cat("****", file =file.mcmc, append = TRUE,sep="\n") 
  cat(0.1, file =file.mcmc, append = TRUE,sep="\n")    #MultFact
  cat(0.1,file =file.mcmc, append = TRUE, sep=",")     #RC MultiFact
  cat(0.1,file =file.mcmc, append = TRUE, sep=",")     
  cat(0.1, file =file.mcmc, append = TRUE,sep="\n")
  cat(0.1,file =file.mcmc, append = TRUE, sep=",")      #Remnant MultiFact
  cat(0.1, file =file.mcmc, append = TRUE,sep="\n")
  ###################################################################   MATRIX
  # file.matrix = paste(dir_code,"/BaM_exe/BaM_BaRatin_2/Config_ControlMatrix.txt",sep="")
  file.matrix = paste(dir,"/Config_ControlMatrix.txt",sep="")
  write.table(M, file =file.matrix, row.names = FALSE, col.names = FALSE)
  cat(Hmax, file =file.matrix, append = TRUE,sep="\n")
  #---------------------------------------------------------------------- 
  # file.remnant = paste(dir_code,"/BaM_exe/BaM_BaRatin_2/Config_RemnantSigma.txt",sep="")
  file.remnant = paste(dir,"/Config_RemnantSigma.txt",sep="")
  if (remnant.err.model == "Linear") {
    cat("'Linear'", file = file.remnant, sep="\n")                         #! Function f used in sdev=f(Qrc) 
    cat(2, file = file.remnant, append = TRUE, sep="\n")                    #! Number of parameters gamma for f
    cat("gamma1", file = file.remnant, append = TRUE, sep="\n")             #! Parameter Name
    cat(g1.prior[3], file = file.remnant, append = TRUE, sep="\n")                   #! Initial Guess
    cat(g1.distr.type, file = file.remnant, append = TRUE, sep="\n")            #! Prior distribution
    cat(g1.prior[1],file =file.remnant, append = TRUE, sep=",")
    cat(",",file =file.remnant, append = TRUE, sep=",")
    cat(g1.prior[2],file =file.remnant, append = TRUE, sep="\n")
    
    cat("gamma2", file = file.remnant, append = TRUE, sep="\n")             #! Parameter Name
    cat(g2.prior[3], file = file.remnant, append = TRUE, sep="\n")                   #! Initial Guess
    cat(g2.distr.type, file = file.remnant, append = TRUE, sep="\n")             #! Initial Guess
    cat(g2.prior[1],file =file.remnant, append = TRUE, sep=",")
    cat(",",file =file.remnant, append = TRUE, sep=",")
    cat(g2.prior[2],file =file.remnant, append = TRUE, sep="\n")
    
  } else if (remnant.err.model == "Constant") {
    cat("'Constant'", file = file.remnant, sep="\n")                       #! Function f used in sdev=f(Qrc) 
    cat(1, file = file.remnant, append = TRUE, sep="\n")                    #! Number of parameters gamma for f
    cat("gamma1", file = file.remnant, append = TRUE, sep="\n")             #! Parameter Name
    cat(g1.prior[3], file = file.remnant, append = TRUE, sep="\n")                   #! Initial Guess
    cat(g1.distr.type, file = file.remnant, append = TRUE, sep="\n")            #! Prior distribution
    cat(g1.prior[1],file =file.remnant, append = TRUE, sep=",")
    cat(",",file =file.remnant, append = TRUE, sep=",")
    cat(g1.prior[2],file =file.remnant, append = TRUE, sep="\n")
  }
  ###################################################################  Config_BaM
  file_BaM <- paste0(dir_code,"/Config_BaM.txt")
  #creation of Config_BaM.txt
  cat(paste0("'",dir.bam.short,"'"), file =file_BaM , sep="\n", append = FALSE)
  cat('"Config_RunOptions.txt"', file = file_BaM , sep="\n", append = TRUE)    
  cat('"Config_Model.txt"', file = file_BaM , sep="\n", append = TRUE)
  cat('"Config_ControlMatrix.txt"', file = file_BaM , sep="\n", append = TRUE)
  cat('"Config_Data.txt"', file = file_BaM , sep="\n", append = TRUE)
  cat('"Config_RemnantSigma.txt"', file = file_BaM , sep="\n", append = TRUE)                                       
  cat('"Config_MCMC.txt"', file = file_BaM , sep="\n", append = TRUE)                                            
  cat('"Config_Cooking.txt"', file = file_BaM , sep="\n", append = TRUE)
  cat('"Config_Summary.txt"', file = file_BaM , sep="\n", append = TRUE)
  cat('"Config_Residuals.txt"', file = file_BaM , sep="\n", append = TRUE)
  if( (predictionRC ==TRUE) | (predictionQt ==TRUE) |(predictionPrior ==TRUE) ) {
    cat('"Config_Pred_Master.txt"', file = file_BaM , sep="\n", append = TRUE)
  } else {
    cat('""', file = file_BaM , sep="\n", append = TRUE)
  }
  ###########################################################################################
  # PREDICTIONS :
  ###############
  file.Pred1 = paste(dir,"/Config_Pred_Master.txt",sep="")
  if ((predictionPrior == TRUE)&(predictionRC== TRUE)&(predictionQt == TRUE))  {
    cat(8, file =file.Pred1,sep="\n")
  } else if ((predictionPrior == TRUE) &(predictionRC == FALSE)&(predictionQt ==FALSE ))  {
    cat(2, file =file.Pred1,sep="\n")
  } else if ((predictionPrior == TRUE) & (predictionRC == TRUE)&(predictionPrior == FALSE)) {
    cat(5, file =file.Pred1,sep="\n")
  } else if ((predictionPrior == TRUE) & (predictionRC == FALSE) &(predictionQt == TRUE)) {
    cat(5, file =file.Pred1,sep="\n")
  } else if ((predictionPrior == FALSE)& (predictionRC == TRUE)&(predictionQt == TRUE)) {
    cat(6, file =file.Pred1,sep="\n")
  } else if ((predictionPrior == FALSE)& (predictionRC == FALSE)&(predictionQt == TRUE)) {
    cat(3, file =file.Pred1,sep="\n")
  } else if ((predictionPrior == FALSE)& (predictionRC == TRUE)&(predictionQt == FALSE)) {
    cat(3, file =file.Pred1,sep="\n")
  }
  
  ###########################
  if(predictionPrior==TRUE) {
    ###########################
    cat("'Config_Pred_Prior.txt'", file =file.Pred1, append = TRUE,sep="\n")
    cat("'Config_Pred_Prior_Qt.txt'", file =file.Pred1, append = TRUE,sep="\n")
    file.Pred11 = paste0(dir,"/Config_Pred_Prior.txt")
    cat(paste0("'",dir.bam.short,"Hgrid.txt'"), file =file.Pred11, sep="\n")
    cat(ngrid, file =file.Pred11, append = TRUE, sep="\n")
    cat("1", file = file.Pred11, append = TRUE,sep="\n")   #n of spaghetti
    cat(".true.", file = file.Pred11, append = TRUE,sep="\n")                              #!!! Propagate parametric uncertainty?
    cat(".false.", file = file.Pred11, append = TRUE,sep="\n")                             #!!! Propagate remnant uncertainty for each output variable? (size nY)
    cat(mcmc.prior, file = file.Pred11, append = TRUE,sep="\n")                          #!!! Nsim[prior]. If <=0: posterior sampling (nsim is given by mcmc sample); 
    #cat(-1, file = file.Pred11, append = TRUE,sep="\n")                                    #!!! Nsim[prior]. If <=0: posterior sampling (nsim is given by mcmc sample); 
    cat("'Qrc_Prior.spag'", file = file.Pred11, append = TRUE,sep="\n")                    # !!! Files containing spaghettis for each output variable (size nY)
    cat(".true.", file = file.Pred11, append = TRUE,sep="\n")                              #!!! Post-processing: transpose spag file (so that each column is a spaghetti)? 
    cat(".true.", file = file.Pred11, append = TRUE,sep="\n")                              #!!! Post-processing: create envelops? (size nY)
    cat("'Qrc_Prior.env'", file = file.Pred11, append = TRUE,sep="\n")                     # !!! Post-processing: name of envelop files (size nY)
    cat(".false.", file = file.Pred11, append = TRUE,sep="\n")                              #!!! Print progress in console during computations?
    cat(".false." , file = file.Pred11, append = TRUE,sep="\n")                           #!!! Do state prediction? (size nState)
    
    ###################################################################
    file.Pred21 = paste0(dir,"/Config_Pred_Prior_Qt.txt")
    cat(paste0("'",dir.bam.short,"limni.txt'"), file =file.Pred21,sep="\n")
    cat(nlimni, file =file.Pred21, append = TRUE, sep="\n")
    cat(nspag, file = file.Pred21, append = TRUE,sep="\n")   #n of spaghetti
    cat(".true.", file = file.Pred21, append = TRUE,sep="\n")                              #!!! Propagate parametric uncertainty?
    cat(".false.", file = file.Pred21, append = TRUE,sep="\n")                             #!!! Propagate remnant uncertainty for each output variable? (size nY)
    cat(mcmc.prior, file = file.Pred21, append = TRUE,sep="\n")                                #!!! Nsim[prior]. If <=0: posterior sampling (nsim is given by mcmc sample); 
    cat("'Qt_prior.spag'", file = file.Pred21, append = TRUE,sep="\n")                    # !!! Files containing spaghettis for each output variable (size nY)
    cat(".true.", file = file.Pred21, append = TRUE,sep="\n")                              #!!! Post-processing: transpose spag file (so that each column is a spaghetti)? 
    cat(".true.", file = file.Pred21, append = TRUE,sep="\n")                              #!!! Post-processing: create envelops? (size nY)
    cat("'Qt_prior.env'", file = file.Pred21, append = TRUE,sep="\n")                     # !!! Post-processing: name of envelop files (size nY)
    cat(".false.", file = file.Pred21, append = TRUE,sep="\n")                              #!!! Print progress in console during computations?
    cat(".false." , file = file.Pred21, append = TRUE,sep="\n")                           #!!! Do state prediction? (size nState)
  }
  
  
  ###################################################################
  if(predictionRC==TRUE) {  
    ###################################################################
    #cat("'Config_Pred_Prior.txt'", file = file.Pred1, append = TRUE,sep="\n")
    cat("'Config_Pred_RCMaxpost.txt'", file = file.Pred1, append = TRUE,sep="\n")
    cat("'Config_Pred_RCParamU.txt'", file =file.Pred1, append = TRUE,sep="\n")
    cat("'Config_Pred_RCTotalU.txt'", file =file.Pred1, append = TRUE,sep="\n")
    ##################################################################
    file.Pred3 = paste(dir,"/Config_Pred_RCMaxpost.txt",sep="")
    cat(paste0("'",dir.bam.short,"Hgrid.txt'"), file =file.Pred3, sep="\n")
    cat(ngrid, file =file.Pred3,sep="\n", append = TRUE)
    cat("1", file = file.Pred3, append = TRUE,sep="\n")   #n of spaghetti
    cat(".false.", file = file.Pred3, append = TRUE,sep="\n")                              #!!! Propagate parametric uncertainty?
    cat(".false.", file = file.Pred3, append = TRUE,sep="\n")                             #!!! Propagate remnant uncertainty for each output variable? (size nY)
    cat("-1", file = file.Pred3, append = TRUE,sep="\n")                                #!!! Nsim[prior]. If <=0: posterior sampling (nsim is given by mcmc sample); 
    cat("'Qrc_Maxpost.spag'", file = file.Pred3, append = TRUE,sep="\n")                    # !!! Files containing spaghettis for each output variable (size nY)
    cat(".true.", file = file.Pred3, append = TRUE,sep="\n")                              #!!! Post-processing: transpose spag file (so that each column is a spaghetti)? 
    cat(".true.", file = file.Pred3, append = TRUE,sep="\n")                              #!!! Post-processing: create envelops? (size nY)
    cat("'Qrc_Maxpost.env'", file = file.Pred3, append = TRUE,sep="\n")                     # !!! Post-processing: name of envelop files (size nY)
    cat(".false.", file = file.Pred3, append = TRUE,sep="\n")                              #!!! Print progress in console during computations?
    cat(".false." , file = file.Pred3, append = TRUE,sep="\n")                           #!!! Do state prediction? (size nState)
    ###################################################################
    file.Pred4 = paste(dir,"/Config_Pred_RCParamU.txt",sep="")
    cat(paste0("'",dir.bam.short,"Hgrid.txt'"), file =file.Pred4, sep="\n")
    cat(ngrid, file =file.Pred4,sep="\n", append = TRUE)
    cat("1", file = file.Pred4, append = TRUE,sep="\n")                                   #n of spaghetti
    cat(".true.", file = file.Pred4, append = TRUE,sep="\n")                              #!!! Propagate parametric uncertainty?
    cat(".false.", file = file.Pred4, append = TRUE,sep="\n")                             #!!! Propagate remnant uncertainty for each output variable? (size nY)
    cat("-1", file = file.Pred4, append = TRUE,sep="\n")                                  #!!! Nsim[prior]. If <=0: posterior sampling (nsim is given by mcmc sample); 
    cat("'Qrc_ParamU.spag'", file = file.Pred4, append = TRUE,sep="\n")                    # !!! Files containing spaghettis for each output variable (size nY)
    cat(".true.", file = file.Pred4, append = TRUE,sep="\n")                              #!!! Post-processing: transpose spag file (so that each column is a spaghetti)? 
    cat(".true.", file = file.Pred4, append = TRUE,sep="\n")                              #!!! Post-processing: create envelops? (size nY)
    cat("'Qrc_ParamU.env'", file = file.Pred4, append = TRUE,sep="\n")                     # !!! Post-processing: name of envelop files (size nY)
    cat(".false.", file = file.Pred4, append = TRUE,sep="\n")                              #!!! Print progress in console during computations?
    cat(".false." , file = file.Pred4, append = TRUE,sep="\n")                            #!!! Do state prediction? (size nState)
    ###################################################################
    file.Pred5 = paste(dir,"/Config_Pred_RCTotalU.txt",sep="")
    cat(paste0("'",dir.bam.short,"Hgrid.txt'"), file =file.Pred5,sep="\n")
    cat(ngrid, file =file.Pred5, append = TRUE, sep="\n")
    cat("1", file = file.Pred5, append = TRUE,sep="\n")                                  #n of spaghetti
    cat(".true.", file = file.Pred5, append = TRUE,sep="\n")                              #!!! Propagate parametric uncertainty?
    cat(".true.", file = file.Pred5, append = TRUE,sep="\n")                             #!!! Propagate remnant uncertainty for each output variable? (size nY)
    cat("-1", file = file.Pred5, append = TRUE,sep="\n")                                #!!! Nsim[prior]. If <=0: posterior sampling (nsim is given by mcmc sample); 
    cat("'Qrc_TotalU.spag'", file = file.Pred5, append = TRUE,sep="\n")                    # !!! Files containing spaghettis for each output variable (size nY)
    cat(".true.", file = file.Pred5, append = TRUE,sep="\n")                              #!!! Post-processing: transpose spag file (so that each column is a spaghetti)? 
    cat(".true.", file = file.Pred5, append = TRUE,sep="\n")                              #!!! Post-processing: create envelops? (size nY)
    cat("'Qrc_TotalU.env'", file = file.Pred5, append = TRUE,sep="\n")                     # !!! Post-processing: name of envelop files (size nY)
    cat(".false.", file = file.Pred5, append = TRUE,sep="\n")                              #!!! Print progress in console during computations?
    cat(".false." , file = file.Pred5, append = TRUE,sep="\n")                           #!!! Do state prediction? (size nState)
  }
  
  
  #############################################################################
  if(predictionQt==TRUE) {
    #############################################################################
    cat("'Config_Pred_Maxpost.txt'", file =file.Pred1, append = TRUE,sep="\n")
    cat("'Config_Pred_ParamU.txt'", file =file.Pred1, append = TRUE,sep="\n")
    cat("'Config_Pred_TotalU.txt'", file =file.Pred1, append = TRUE,sep="\n")
    
    file.Pred6 = paste0(dir,"/Config_Pred_Maxpost.txt")
    cat(paste0("'",dir.bam.short,"limni.txt'"), file =file.Pred6,sep="\n")
    cat(nlimni, file =file.Pred6, append = TRUE, sep="\n")
    cat(nspag, file = file.Pred6, append = TRUE,sep="\n")   #n of spaghetti
    cat(".false.", file = file.Pred6, append = TRUE,sep="\n")                              #!!! Propagate parametric uncertainty?
    cat(".false.", file = file.Pred6, append = TRUE,sep="\n")                             #!!! Propagate remnant uncertainty for each output variable? (size nY)
    cat("-1", file = file.Pred6, append = TRUE,sep="\n")                                #!!! Nsim[prior]. If <=0: posterior sampling (nsim is given by mcmc sample); 
    cat("'Qt_Maxpost.spag'", file = file.Pred6, append = TRUE,sep="\n")                    # !!! Files containing spaghettis for each output variable (size nY)
    cat(".true.", file = file.Pred6, append = TRUE,sep="\n")                              #!!! Post-processing: transpose spag file (so that each column is a spaghetti)? 
    cat(".true.", file = file.Pred6, append = TRUE,sep="\n")                              #!!! Post-processing: create envelops? (size nY) 
    cat("'Qt_Maxpost.env'", file = file.Pred6, append = TRUE,sep="\n")                     # !!! Post-processing: name of envelop files (size nY)
    cat(".false.", file = file.Pred6, append = TRUE,sep="\n")                              #!!! Print progress in console during computations?
    cat(".false." , file = file.Pred6, append = TRUE,sep="\n")                           #!!! Do state prediction? (size nState)
    
    
    file.Pred8 = paste0(dir,"/Config_Pred_ParamU.txt")
    cat(paste0("'",dir.bam.short,"limni.txt'"), file =file.Pred8,sep="\n")
    cat(nlimni, file =file.Pred8, append = TRUE, sep="\n")
    cat(nspag, file = file.Pred8, append = TRUE,sep="\n")   #n of spaghetti
    cat(".true.", file = file.Pred8, append = TRUE,sep="\n")                              #!!! Propagate parametric uncertainty?
    cat(".false.", file = file.Pred8, append = TRUE,sep="\n")                             #!!! Propagate remnant uncertainty for each output variable? (size nY)
    cat("-1", file = file.Pred8, append = TRUE,sep="\n")                                #!!! Nsim[prior]. If <=0: posterior sampling (nsim is given by mcmc sample); 
    cat("'Qt_ParamU.spag'", file = file.Pred8, append = TRUE,sep="\n")                    # !!! Files containing spaghettis for each output variable (size nY)
    cat(".true.", file = file.Pred8, append = TRUE,sep="\n")                              #!!! Post-processing: transpose spag file (so that each column is a spaghetti)? 
    cat(".true.", file = file.Pred8, append = TRUE,sep="\n")                              #!!! Post-processing: create envelops? (size nY)
    cat("'Qt_ParamU.env'", file = file.Pred8, append = TRUE,sep="\n")                     # !!! Post-processing: name of envelop files (size nY)
    cat(".false.", file = file.Pred8, append = TRUE,sep="\n")                              #!!! Print progress in console during computations?
    cat(".false." , file = file.Pred8, append = TRUE,sep="\n")                           #!!! Do state prediction? (size nState)
    
    
    file.Pred7 = paste(dir,"/Config_Pred_TotalU.txt",sep="")
    cat(paste0("'",dir.bam.short,"limni.txt'"), file =file.Pred7,sep="\n")
    cat(nlimni, file =file.Pred7, append = TRUE, sep="\n")
    cat(nspag, file = file.Pred7, append = TRUE,sep="\n")   #n of spaghetti
    cat(".true.", file = file.Pred7, append = TRUE,sep="\n")                              #!!! Propagate parametric uncertainty?
    cat(".true.", file = file.Pred7, append = TRUE,sep="\n")                             #!!! Propagate remnant uncertainty for each output variable? (size nY)
    cat("-1", file = file.Pred7, append = TRUE,sep="\n")                                #!!! Nsim[prior]. If <=0: posterior sampling (nsim is given by mcmc sample); 
    cat("'Qt_TotalU.spag'", file = file.Pred7, append = TRUE,sep="\n")                    # !!! Files containing spaghettis for each output variable (size nY)
    cat(".true.", file = file.Pred7, append = TRUE,sep="\n")                              #!!! Post-processing: transpose spag file (so that each column is a spaghetti)? 
    cat(".true.", file = file.Pred7, append = TRUE,sep="\n")                              #!!! Post-processing: create envelops? (size nY)
    cat("'Qt_TotalU.env'", file = file.Pred7, append = TRUE,sep="\n")                     # !!! Post-processing: name of envelop files (size nY)
    cat(".false.", file = file.Pred7, append = TRUE,sep="\n")                              #!!! Print progress in console during computations?
    cat(".false." , file = file.Pred7, append = TRUE,sep="\n")                           #!!! Do state prediction? (size nState)
  }
  
  
  
  ##########################################################################
  # RUN OPTIONS 
  ##########################################################################
  file.run = paste0(dir,"/Config_RunOptions.txt")
  if (simMCMC == TRUE) {  
    cat(".true.", file =file.run,sep="\n")   # Do the MCMC simulations
  } else {
    cat(".false.", file =file.run,sep="\n")  # Do not do the MCMC simulation    
  }
  cat(".true.", file =file.run, append = TRUE, sep="\n")    # Do MCMC summary?
  cat(".true.", file =file.run, append = TRUE, sep="\n")    # Do Residual diagnostics?
  if ((predictionRC == TRUE) | (predictionQt == TRUE) | (predictionPrior ==TRUE)) {
    cat(".true.", file =file.run, append = TRUE, sep="\n")    # Do Predictions?
  } else {
    cat(".false.", file =file.run, append = TRUE, sep="\n")   # Do Predictions?
  }
  
  ###################################################################   RESIDUALS CONFIG
  file.residuals = paste(dir,"/Config_Residuals.txt",sep="")
  cat('"Results_Residuals.txt"' , file =file.residuals ,sep="\n")      # Result file
  
  ###################################################################   SUMMARY CONFIG
  file.summary = paste(dir,"/Config_Summary.txt",sep="")
  cat('"Results_Summary.txt"' , file =file.summary ,sep="\n")    #Summary stat results file name
  
  ###################################################################   COOKING CONFIG
  file.cooking = paste(dir,"/Config_Cooking.txt",sep="")
  cat('"Results_MCMC_Cooked.txt"' , file =file.cooking ,sep="\n")  # mcmc Results file name
  cat(0.5, file =file.cooking, append = TRUE, sep="\n")            # Burn factor
  cat(10, file =file.cooking, append = TRUE, sep="\n")             # Nslim
}






################################################################################################################
BaRatin.propagation <- function(dir, nsim, propagat, b.distr , a.distr, c.distr,  
                                a.prior, st_a.prior, c.prior, st_c.prior, b.prior, st_b.prior, #priors
                                Bw.prior, Cr.prior, g.prior, Bc.prior, KS.prior, S0.prior,
                                st_Bw.prior, st_Cr.prior, st_g.prior, st_Bc.prior, st_KS.prior, st_S0.prior,
                                ncontrol, M ) {
  ################################################################################################################  
  
  names.bac=c('b','a','c')
  margins.bac =c(b.distr ,a.distr, c.distr)
  # Define priors on "physical" parameters
  #***************************************
  b=NULL;c=NULL;Bw=NULL;Cr=NULL;g=NULL;Bc=NULL;KS=NULL;S0=NULL; a = NULL; # Priors are defined through Monte Carlo samples
  for (i in 1:ncontrol) {
    # To allow reproductivity, 
    # fixe a pseudo random generation for the rnorm() used later for the propagation of BaRatin parameters:           
    set.seed(b.prior[i])  
    b[[i]] <- rnorm(nsim, mean=b.prior[i], sd=st_b.prior[i]) # offset (stage for which Q=0 for that control)
    set.seed(c.prior[i]) 
    c[[i]] <- rnorm(nsim, mean=c.prior[i], sd=st_c.prior[i]) # exponent
    if (propagat == TRUE) {
      if (control.type[i] == "rect.weir") {
        set.seed(Bw.prior[i]) 
        Bw[[i]] <- rlnorm(nsim, meanlog=log(Bw.prior[i]), sdlog=st_b.prior[i])  # weir width
        set.seed(Cr.prior[i]) 
        Cr[[i]] <- rlnorm(nsim, meanlog=log(Cr.prior[i]), sdlog=st_Cr.prior[i]) # Discharge coefficient
        set.seed(g.prior[i]) 
        g[[i]]  <- rlnorm(nsim, meanlog=log(g.prior[i]),  sdlog=st_g.prior[i])  # gravity
      } else if (control.type[i] == "rect.channel") {
        set.seed(Bc.prior[i]) 
        Bc[[i]] <- rlnorm(nsim, meanlog=log(Bc.prior[i]), sdlog=st_Bc.prior[i]) # channel width
        set.seed(KS.prior[i]) 
        KS[[i]] <- rlnorm(nsim, meanlog=log(KS.prior[i]), sdlog=st_KS.prior[i]) # Strickler coefficient
        set.seed(S0.prior[i]) 
        S0[[i]] <- rlnorm(nsim, meanlog=log(S0.prior[i]), sdlog=st_S0.prior[i]) # slope
      }
    } else {
      set.seed(a.prior[i]) 
      a[[i]] <- rnorm(nsim, mean = a.prior[i], sd =st_a.prior[i]) # exponent
    }
  }  
  
  # Starting point for all parameters:
  #***********************************
  start <- list("test")
  for (i in 1:ncontrol) {
    start[[paste0("b",i)]] = b.prior[i] 
    start[[paste0("c",i)]] = c.prior[i]
    if (propagat == TRUE) {
      if (control.type[i] == "rect.weir") {
        start[[paste0("Bw",i)]] = Bw.prior[i]
        start[[paste0("Cr",i)]] = Cr.prior[i]
        start[[paste0("g",i)]] = g.prior[i]
      } else if (control.type[i] == "rect.channel") {
        start[[paste0("Bc",i)]] = Bc.prior[i]
        start[[paste0("KS",i)]] = KS.prior[i]
        start[[paste0("S0",i)]] = S0.prior[i]
      } else if (control.type[i] == "other.function") {
        start[[paste0("a",i)]]= a.prior[i]
      }
    } else {
      start[[paste0("a",i)]]= a.prior[i]
    }
  }
  
  if (propagat == TRUE) {
    # Perform Monte-Carlo propagation ( "Prior_Propagation.R")
    #**********************************************************
    # nsim=length(b[[1]]);  
    # initialise parameters that will be computed by propagation
    a = NULL
    for (i in 1:ncontrol) {
      a[[i]]= vector(mode='double',length=nsim); 
    }
    for(i in 1:nsim){
      # Compute a for each control
      for (n in 1:ncontrol) {
        if (control.type[n] == "rect.weir") {
          a[[n]][i]= Cr[[n]][i]*Bw[[n]][i]*sqrt(2*g[[n]][i])
        } else if (control.type[n] == "rect.channel") {
          a[[n]][i]=KS[[n]][i]*Bc[[n]][i]*sqrt(S0[[n]][i])
        }
      }
    }
  }
  
  #starting values:
  #----------------
  a.s = NULL
  if (propagat == TRUE) {
    for (j in 1:ncontrol) {
      if (control.type[j] == "rect.weir") {
        a.s[[j]] = start[[paste0("Cr",j)]] *start[[paste0("Bw",j)]]* sqrt(2*start[[paste0("g",j)]])
      } else if (control.type[j] == "rect.channel") {
        a.s[[j]] = start[[paste0("KS",j)]] *start[[paste0("Bc",j)]]* sqrt(start[[paste0("S0",j)]])
      }
    }
  }
  #final vectors:
  #-----------------------------------
  sim= data.frame(b[[1]], a[[1]], c[[1]])
  if (propagat == TRUE) {
    start2= c(start[paste0("b",1)], a.s[1], start[paste0("c",1)])
  } else {
    start2= c(start[paste0("b",1)], start[paste0("a",1)], start[paste0("c",1)])
  }
  if (ncontrol > 1) {
    for (j in 2:ncontrol) {
      sim = cbind(sim , b[[j]], a[[j]], c[[j]] )
      if (propagat==TRUE) {
        start2 = c(start2, start[paste0("b",j)], a.s[j], start[paste0("c",j)])
      } else {
        start2 = c(start2, start[paste0("b",j)], start[paste0("a",j)], start[paste0("c",j)])
      }
    }
  }
  MC = list(sim=sim, start=start2)
  
  
  # Fit prior distribution on MCMC samples
  #***************************************
  margins=c();names=c();k=0    # define marginal prior distributions and parameter names
  for(i in 1:ncontrol){
    for(j in 1:3){
      k=k+1
      margins=c(margins,margins.bac[j])
      names=c(names,paste(names.bac[j],i,sep=''))
    }
  }
  # fit multivariate prior distribution
  prior =fit(sim     = MC$sim, 
             margins = margins, 
             names   = names,
             ncol    = 6,
             rowmax  = 3,
             plots=F             ) #call function from "Prior_Propagation.R"
  
  return(list(prior, MC$start))
}



#######################################################################################################
writeConfigFiles<-function(prior, # prior object produced by function 'fit'
                           start, # starting vector produced by function propagate_XXX$start
                           ncontrol, # number of hydraulic controls
                           model     ='Config_Model.txt', # Model configuration file
                           names.bac = c('b','a','c'), # base name of the 3 parameters for one control
                           dir
){
  #######################################################################################################
  #^* PURPOSE: write configuration files used by BaM
  comment=c(
    'Model ID',
    'nX: number of input variables',
    'nY: number of output variables',
    'nPar: number of parameters theta'
  )
  val=list('BaRatinBAC',1,1,3*ncontrol)
  # specification for each parameter
  k=0;m=1
  for(i in 1:ncontrol){ # loop on each hydraulic control
    for(j in 1:length(names.bac)){ # loop on each b-a-c parameter of the control
      k=k+1
      # comments
      comment=c(comment,
                'Parameter name',
                'Initial guess',
                'Prior distribution',
                'Prior parameters')
      # parname
      pname=paste(names.bac[j],i,sep='')
      val=c(val,pname)
      # starting point
      val=c(val, start[m])
      # prior distribution
      val=c(val, prior$margins[m])
      # prior parameters
      val=c(val, list(prior$priorpar[[m]]))
      m=m+1                
    }
  }
  writeConfig(val,comment,dir,model)
}


##########################################################################################
RC_controls = function(theta, h, M, ncontrols, op=1) {
  ##########################################################################################
  Q        = 0*h
  mask     = list();
  #stop     = FALSE
  if  (ncontrols >1) {
    for (nc in 1:(ncontrols - 1)) {
      # if ((abs(theta[ncontrols*3 + 3 + nc] - theta[ncontrols*3 + 3 + nc+1]) 
      #      <= 0.01)|(any(h <= theta[ncontrols*3 + 3 +1]))) {
      #   
      #    break 
      #    stop= TRUE
      # } else {
      mask[[nc]] = which(((h >  theta[ncontrols*3 + 3 + nc]) +
                          (h <= theta[ncontrols*3 + 3 + nc+1])) == 2)
      # }
    }
  }
  
  #if (stop==FALSE){
  mask[[ncontrols]] = which(h > theta[ncontrols*3 + 3 + ncontrols])
  for (segm in 1:ncontrols){  # for each segment of stage h
    Q[mask[[segm]]] = 0 
    for (ccc in 1:ncontrols){  # for each control 
      if (M[segm, ccc] ==1){
        Q[mask[[segm]]] =   Q[mask[[segm]]]  +
          theta[3*ccc - 1]*(h[mask[[segm]]] - theta[3*ccc -2])^theta[3*ccc]
      }
    }
  }   
  # 
  if(op == 1) {
    resQ = sapply(Q, 
                  function(Q,theta){Q + rnorm(1, 
                                              mean = 0, 
                                              sd = theta[1]+ theta[2]*Q)}, 
                  theta = c(theta[3*ncontrols +1], theta[3*ncontrols +2]))
  } else {
    resQ = Q
  }
  # } else {
  #     resQ =NA
  # }
  return(resQ)
}


# Rating curve for maximum posterior:   
###########################################################################################
RC_controls_mp = function(theta, h, M, ncontrols, op=1){ 
  ###########################################################################################
  Q        = 0*h
  mask     = NULL
  if  (ncontrols >1) {
    for (nc in 1:(ncontrols - 1)) {
      mask[[nc]] = which(((h >  theta[ncontrols*3 + 2 + nc]) +
                            (h <= theta[ncontrols*3 + 2 + nc+1])) == 2)
    }
  }
  mask[[ncontrols]] = which(h> theta[ (ncontrols*3 + 2 + ncontrols)])
  for (segm in 1:ncontrols){  # for each segment of stage h
    Q[mask[[segm]]] = 0 
    for (ccc in 1:ncontrols){  # for each control 
      if (M[segm, ccc] ==1){
        Q[mask[[segm]]] =   Q[mask[[segm]]]  + 
          theta[3*ccc - 1]*(h[mask[[segm]]] - theta[3*ccc -2])^theta[3*ccc]
      }
    }
  }  
  
  return(Q)
}

  
  
#########
plot.RC = function(data.MCMC,data.MCMC.MaxPost,ncontrols,gaugings,M,hgrid,log=T,dir.save,per){
  
  nsample           = length(data.MCMC[,1])
  min.grid          = min(data.MCMC[,1]) 
  
  MCMC.save    =  matrix(NA, nrow = nsample, ncol = ncol(data.MCMC))  
  MaxPost.save =  matrix(NA, nrow = 1,  ncol = ncol(data.MCMC)-1)  
  MCMC.save = cbind(data.MCMC, rep(1, nsample))
  MaxPost.save = c(data.MCMC.MaxPost,1)
  
  RC.Post    =  apply(MCMC.save,MARGIN = 1,  RC_controls,h = hgrid,  M = M, ncontrols = ncontrols)
  RC.MaxPost =  RC_controls_mp(theta = MaxPost.save,h=hgrid,M=M,ncontrols = ncontrols)
  
  
  data.tmp            = apply(RC.Post, MARGIN=1, quantile, probs=c(0.025,0.975), na.rm=TRUE)
  data.tmp            = apply(data.tmp, MARGIN=c(1,2), function(x){ifelse(x<0,0,x)})
  List.RC.quants      = data.frame(cbind(hgrid, t(data.tmp),  RC.MaxPost))
  colnames(List.RC.quants) = c("h", "inf", "sup", "maxpost")
  
  data.RC = List.RC.quants
  
 gaugings$Period=as.factor(gaugings$Period)
  
  plot.RC = ggplot(data.RC)+
        geom_vline(xintercept = MaxPost.save[ncontrols*3 +2 + 1],
                   colour = "blue", size=0.6, linetype ="dashed")+
        geom_vline(xintercept = MaxPost.save[ ncontrols*3 +2 + 2],
                   colour = "blue", size=0.6)+
        # geom_vline(xintercept = MaxPost.save[ ncontrols*3 +2 + 3],
        #            colour = "blue", size=0.6)+
        geom_smooth(aes(x=h, y=maxpost, ymax=sup, ymin=inf), 
                  size=1, stat='identity', alpha = 0.3) +  #alpha=0.1
        geom_path(aes(x=h, y=maxpost), na.rm=TRUE, size=1)+
        ### Gaugings
        geom_linerange(aes(x=h, ymax=Q+2*uQ, ymin=Q-2*uQ, colour=Period), data=gaugings, size=0.5)+
        geom_point(aes(x=h,y=Q,colour=Period),data=gaugings,  shape=16, size=2)+
        #axis
        coord_cartesian(ylim=c(1,15000), xlim=c(-6,12))+
        scale_x_continuous(breaks=seq(-6,12,2), labels=seq(-6,12,2), expand=c(0,0))+
        scale_y_continuous(breaks=seq(1,15000,5000), labels=seq(1,15000,5000), expand=c(0,0))+
        #labs
        xlab(expression(paste("Stage h [m]",sep="")))+
        ylab(expression(paste("Discharge Q [",m^3,".",s^-1,"]",sep="")))+
        labs(colour = "Period")+
        theme_bw(base_size=10)+
        theme( axis.text        = element_text(size=15)
               ,axis.title       = element_text(size=20,face="bold")
               ,panel.grid.major = element_blank()  # element_line(size=1.2)
               ,panel.grid.minor = element_blank()  # element_line(size=0.8)
               ,legend.text      = element_text(size=20)
               ,legend.title     = element_text(size=30)
               ,plot.margin      = unit(c(0.5,0.5,0.5,0.5),"cm")
               ,legend.key.size  = unit(1.5, "cm")
               ,legend.position  = "none")

  pdf(paste0(dir.save,"/RC_P",per,".pdf"), 14, 8 ,useDingbats=F)
  print(plot.RC)
  dev.off()
  
  if(log==T){
  # RC plot Logarithmic scale:
  plot.RClog = plot.RC +
    scale_y_log10(breaks = c(100,1000, 10000) , labels =c(100,1000, 10000) , expand=c(0,0))+
    annotation_logticks(base = 10, sides = "l", scaled = TRUE, colour = "black",
                         size = 1, linetype = 1)
  
  pdf(paste0(dir.save,"/RClog_P",per,".pdf"), 14, 8 ,useDingbats=F)
  print(plot.RClog)
  dev.off()}
  
}
  

  
  
  
  
  