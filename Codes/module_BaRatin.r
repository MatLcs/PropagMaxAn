##############################################################################################################################
BaRatin_app <- function(dir,
                        t_Gaug, Q_Gaug, uQ_Gaug,h_Gaug,
                        nobs.gaug, nlimni,  #data of gaugings and limni
                        propagat,
                        b.distr ,a.distr, c.distr,
                        a.prior, st_a.prior,  #priors
                        c.prior, st_c.prior,
                        b.prior, st_b.prior,        #priors
                        Bw.prior, Cr.prior, g.prior, 
                        Bc.prior, KS.prior, S0.prior,
                        st_Bw.prior, st_Cr.prior, st_g.prior, 
                        st_Bc.prior, st_KS.prior, st_S0.prior,
                        ncontrol, M ,
                        remnant.err.model, g1.prior, g2.prior, 
                        g1.distr.type, g2.distr.type,  #remnant error model priors
                        predictionRC, 
                        predictionQt, 
                        predictionPrior,
                        simMCMC,   mcmc.prior,                                 #predictions choice
                        Ncycles,                                                              #mcmc oprtions
                        Hmin, Hmax,                                                           #grid limits for plots
                        iter                                                                  #Iteration index of the segm method
                        ) {
###############################################################################################################################
  Hgrid = seq(Hmin, Hmax, 0.01)
  ngrid = length(Hgrid)
  nsim  = 10000
  write.table(Hgrid, file =paste0(dir_code,"/BaM_exe/BaM_BaRatin_2/Hgrid.txt"), col.names = FALSE, row.names = FALSE)
  
  
  BaRatin_config(dir,
                 nsim = nsim, 
                 propagat, 
                 b.distr ,a.distr, c.distr,  
                 a.prior, st_a.prior, 
                 c.prior, st_c.prior,
                 b.prior, st_b.prior, #priors
                 Bw.prior, Cr.prior, g.prior, 
                 Bc.prior, KS.prior, S0.prior,
                 st_Bw.prior, st_Cr.prior, st_g.prior,
                 st_Bc.prior, st_KS.prior, st_S0.prior,
                 ncontrol, 
                 M,
                 nobs = nobs.gaug,
                 Ncycles, 
                 ngrid,
                 nlimni, 
                 predictionRC,
                 predictionQt,
                 predictionPrior,
                 simMCMC,
                 mcmc.prior,
                 remnant.err.model, 
                 g1.prior, g2.prior, g1.distr.type, g2.distr.type
                 )
  message("***************************************************************"); flush.console()
  message(c("iteration=",iter,"   Applying BaRatin !!!  Wait ... ")); flush.console()
  message("***************************************************************"); flush.console()
  system2(paste(dir_code,"/BaM_exe/BaM_2exp_pool2.exe",sep=""), stdout =NULL, stderr =NULL); 
}















#########################################################################################################
BaRatin_config <- function(dir, 
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
                           remnant.err.model, g1.prior, g2.prior, g1.distr.type, g2.distr.type) {
#########################################################################################################
  Hmax = 10  ##### CAREFUL !!!!!!!!!!!!!
  
  
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
  file.name2 = paste0(dir_code,"/BaM_exe/BaM_BaRatin_2/Config_Data.txt")
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
  file.mcmc = paste(dir_code,"/BaM_exe/BaM_BaRatin_2/Config_MCMC.txt",sep="")
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
  file.matrix = paste(dir_code,"/BaM_exe/BaM_BaRatin_2/Config_ControlMatrix.txt",sep="")
  write.table(M, file =file.matrix, row.names = FALSE, col.names = FALSE)
  cat(Hmax, file =file.matrix, append = TRUE,sep="\n")
  #---------------------------------------------------------------------- 
  file.remnant = paste(dir_code,"/BaM_exe/BaM_BaRatin_2/Config_RemnantSigma.txt",sep="")
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
  file_BaM <- paste(dir_code,"/BaM_exe/Config_BaM.txt",sep="")
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
  file.Pred1 = paste(dir_code,"/BaM_exe/BaM_BaRatin_2/Config_Pred_Master.txt",sep="")
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
    file.Pred11 = paste0(dir_code,"/BaM_exe/BaM_BaRatin_2/Config_Pred_Prior.txt")
    cat("'BaM_BaRatin_2\\Hgrid.txt'", file =file.Pred11, sep="\n")
    cat(ngrid, file =file.Pred11, append = TRUE, sep="\n")
    cat("1", file = file.Pred11, append = TRUE,sep="\n")   #n of spaghetti
    cat(".true.", file = file.Pred11, append = TRUE,sep="\n")                              #!!! Propagate parametric uncertainty?
    cat(".true.", file = file.Pred11, append = TRUE,sep="\n")                             #!!! Propagate remnant uncertainty for each output variable? (size nY)
    cat(mcmc.prior, file = file.Pred11, append = TRUE,sep="\n")                          #!!! Nsim[prior]. If <=0: posterior sampling (nsim is given by mcmc sample); 
    #cat(-1, file = file.Pred11, append = TRUE,sep="\n")                                    #!!! Nsim[prior]. If <=0: posterior sampling (nsim is given by mcmc sample); 
    cat("'Qrc_Prior.spag'", file = file.Pred11, append = TRUE,sep="\n")                    # !!! Files containing spaghettis for each output variable (size nY)
    cat(".true.", file = file.Pred11, append = TRUE,sep="\n")                              #!!! Post-processing: transpose spag file (so that each column is a spaghetti)? 
    cat(".true.", file = file.Pred11, append = TRUE,sep="\n")                              #!!! Post-processing: create envelops? (size nY)
    cat("'Qrc_Prior.env'", file = file.Pred11, append = TRUE,sep="\n")                     # !!! Post-processing: name of envelop files (size nY)
    cat(".true.", file = file.Pred11, append = TRUE,sep="\n")                              #!!! Print progress in console during computations?
    cat(".false." , file = file.Pred11, append = TRUE,sep="\n")                           #!!! Do state prediction? (size nState)

    ###################################################################
    file.Pred21 = paste0(dir_code,"/BaM_exe/BaM_BaRatin_2/Config_Pred_Prior_Qt.txt")
    cat("'BaM_BaRatin_2\\limni.txt'", file =file.Pred21,sep="\n")
    cat(nlimni, file =file.Pred21, append = TRUE, sep="\n")
    cat("1", file = file.Pred21, append = TRUE,sep="\n")   #n of spaghetti
    cat(".true.", file = file.Pred21, append = TRUE,sep="\n")                              #!!! Propagate parametric uncertainty?
    cat(".true.", file = file.Pred21, append = TRUE,sep="\n")                             #!!! Propagate remnant uncertainty for each output variable? (size nY)
    cat(mcmc.prior, file = file.Pred21, append = TRUE,sep="\n")                                #!!! Nsim[prior]. If <=0: posterior sampling (nsim is given by mcmc sample); 
    cat("'Qt_prior.spag'", file = file.Pred21, append = TRUE,sep="\n")                    # !!! Files containing spaghettis for each output variable (size nY)
    cat(".true.", file = file.Pred21, append = TRUE,sep="\n")                              #!!! Post-processing: transpose spag file (so that each column is a spaghetti)? 
    cat(".true.", file = file.Pred21, append = TRUE,sep="\n")                              #!!! Post-processing: create envelops? (size nY)
    cat("'Qt_prior.env'", file = file.Pred21, append = TRUE,sep="\n")                     # !!! Post-processing: name of envelop files (size nY)
    cat(".true.", file = file.Pred21, append = TRUE,sep="\n")                              #!!! Print progress in console during computations?
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
    file.Pred3 = paste(dir_code,"/BaM_exe/BaM_BaRatin_2/Config_Pred_RCMaxpost.txt",sep="")
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
    file.Pred4 = paste(dir_code,"/BaM_exe/BaM_BaRatin_2/Config_Pred_RCParamU.txt",sep="")
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
    file.Pred5 = paste(dir_code,"/BaM_exe/BaM_BaRatin_2/Config_Pred_RCTotalU.txt",sep="")
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

    file.Pred6 = paste0(dir_code,"/BaM_exe/BaM_BaRatin_2/Config_Pred_Maxpost.txt")
    cat("'BaM_BaRatin_2\\limni.txt'", file =file.Pred6,sep="\n")
    cat(nlimni, file =file.Pred6, append = TRUE, sep="\n")
    cat("1", file = file.Pred6, append = TRUE,sep="\n")   #n of spaghetti
    cat(".false.", file = file.Pred6, append = TRUE,sep="\n")                              #!!! Propagate parametric uncertainty?
    cat(".false.", file = file.Pred6, append = TRUE,sep="\n")                             #!!! Propagate remnant uncertainty for each output variable? (size nY)
    cat("-1", file = file.Pred6, append = TRUE,sep="\n")                                #!!! Nsim[prior]. If <=0: posterior sampling (nsim is given by mcmc sample); 
    cat("'Qt_maxpost.spag'", file = file.Pred6, append = TRUE,sep="\n")                    # !!! Files containing spaghettis for each output variable (size nY)
    cat(".true.", file = file.Pred6, append = TRUE,sep="\n")                              #!!! Post-processing: transpose spag file (so that each column is a spaghetti)? 
    cat(".true.", file = file.Pred6, append = TRUE,sep="\n")                              #!!! Post-processing: create envelops? (size nY) 
    cat("'Qt_Maxpost.env'", file = file.Pred6, append = TRUE,sep="\n")                     # !!! Post-processing: name of envelop files (size nY)
    cat(".false.", file = file.Pred6, append = TRUE,sep="\n")                              #!!! Print progress in console during computations?
    cat(".false." , file = file.Pred6, append = TRUE,sep="\n")                           #!!! Do state prediction? (size nState)

    
    file.Pred8 = paste0(dir_code,"/BaM_exe/BaM_BaRatin_2/Config_Pred_ParamU.txt")
    cat("'BaM_BaRatin_2\\limni.txt'", file =file.Pred8,sep="\n")
    cat(nlimni, file =file.Pred8, append = TRUE, sep="\n")
    cat("1", file = file.Pred8, append = TRUE,sep="\n")   #n of spaghetti
    cat(".true.", file = file.Pred8, append = TRUE,sep="\n")                              #!!! Propagate parametric uncertainty?
    cat(".false.", file = file.Pred8, append = TRUE,sep="\n")                             #!!! Propagate remnant uncertainty for each output variable? (size nY)
    cat("-1", file = file.Pred8, append = TRUE,sep="\n")                                #!!! Nsim[prior]. If <=0: posterior sampling (nsim is given by mcmc sample); 
    cat("'Qt_ParamU.spag'", file = file.Pred8, append = TRUE,sep="\n")                    # !!! Files containing spaghettis for each output variable (size nY)
    cat(".true.", file = file.Pred8, append = TRUE,sep="\n")                              #!!! Post-processing: transpose spag file (so that each column is a spaghetti)? 
    cat(".true.", file = file.Pred8, append = TRUE,sep="\n")                              #!!! Post-processing: create envelops? (size nY)
    cat("'Qt_ParamU.env'", file = file.Pred8, append = TRUE,sep="\n")                     # !!! Post-processing: name of envelop files (size nY)
    cat(".false.", file = file.Pred8, append = TRUE,sep="\n")                              #!!! Print progress in console during computations?
    cat(".false." , file = file.Pred8, append = TRUE,sep="\n")                           #!!! Do state prediction? (size nState)

    
    file.Pred7 = paste(dir_code,"/BaM_exe/BaM_BaRatin_2/Config_Pred_TotalU.txt",sep="")
    cat("'BaM_BaRatin_2\\limni.txt'", file =file.Pred7,sep="\n")
    cat(nlimni, file =file.Pred7, append = TRUE, sep="\n")
    cat("1", file = file.Pred7, append = TRUE,sep="\n")   #n of spaghetti
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
  file.run = paste0(dir_code,"/BaM_exe/BaM_BaRatin_2/Config_RunOptions.txt")
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
  file.residuals = paste(dir_code,"/BaM_exe/BaM_BaRatin_2/Config_Residuals.txt",sep="")
  cat('"Results_Residuals.txt"' , file =file.residuals ,sep="\n")      # Result file
  
  ###################################################################   SUMMARY CONFIG
  file.summary = paste(dir_code,"/BaM_exe/BaM_BaRatin_2/Config_Summary.txt",sep="")
  cat('"Results_Summary.txt"' , file =file.summary ,sep="\n")    #Summary stat results file name
  
  ###################################################################   COOKING CONFIG
  file.cooking = paste(dir_code,"/BaM_exe/BaM_BaRatin_2/Config_Cooking.txt",sep="")
  cat('"Results_MCMC_Cooked.txt"' , file =file.cooking ,sep="\n")  # mcmc Results file name
  cat(0.5, file =file.cooking, append = TRUE, sep="\n")            # Burn factor
  cat(10, file =file.cooking, append = TRUE, sep="\n")             # Nslim
}















#########################################################################################
BaM_BaRatin_noPred <- function(nobs,ngrid) {
#########################################################################################
  file_BaM_nopred <- paste(dir_code,"/BaM_exe/Config_BaM.txt",sep="")
  #creation of Config_BaM.txt
  cat('"BaM_BaRatin_2/"', file =file_BaM_nopred , sep="\n", append = FALSE)
  cat('"Config_RunOptions.txt"', file = file_BaM_nopred , sep="\n", append = TRUE)    
  cat('"Config_Model.txt"', file = file_BaM_nopred , sep="\n", append = TRUE)
  cat('"Config_ControlMatrix.txt"', file = file_BaM_nopred , sep="\n", append = TRUE)
  cat('"Config_Data.txt"', file = file_BaM_nopred , sep="\n", append = TRUE)
  cat('"Config_RemnantSigma.txt"', file = file_BaM_nopred , sep="\n", append = TRUE)                                       
  cat('"Config_MCMC.txt"', file = file_BaM_nopred , sep="\n", append = TRUE)                                            
  cat('"Config_Cooking.txt"', file = file_BaM_nopred , sep="\n", append = TRUE)
  cat('"Config_Summary.txt"', file = file_BaM_nopred , sep="\n", append = TRUE)
  cat('"Config_Residuals.txt"', file = file_BaM_nopred , sep="\n", append = TRUE)
  cat('"Config_Pred_Master.txt"', file = file_BaM_nopred , sep="\n", append = TRUE)
  
  file.Pred1 = paste(dir_code,"/BaM_exe/BaM_BaRatin_2/Config_Pred_Master.txt",sep="")
  cat('4', file =file.Pred1,sep="\n")
  cat("'Config_Pred_Prior.txt'", file = file.Pred1, append = TRUE,sep="\n")
  cat("'Config_Pred_RCMaxpost.txt'", file = file.Pred1, append = TRUE,sep="\n")
  cat("'Config_Pred_RCParamU.txt'", file =file.Pred1, append = TRUE,sep="\n")
  cat("'Config_Pred_RCTotalU.txt'", file =file.Pred1, append = TRUE,sep="\n")
  
  file.Pred2 = paste(dir_code,"/BaM_exe/BaM_BaRatin_2/Config_Pred_Prior.txt",sep="")
  cat("'BaM_BaRatin_2\\Hgrid.txt'", file =file.Pred2, sep="\n")
  cat(ngrid, file =file.Pred2, append = TRUE, sep="\n")
  cat("1", file = file.Pred2, append = TRUE,sep="\n")   #n of spaghetti
  cat(".true.", file = file.Pred2, append = TRUE,sep="\n")                              #!!! Propagate parametric uncertainty?
  cat(".false.", file = file.Pred2, append = TRUE,sep="\n")                             #!!! Propagate remnant uncertainty for each output variable? (size nY)
  cat("1000", file = file.Pred2, append = TRUE,sep="\n")                                #!!! Nsim[prior]. If <=0: posterior sampling (nsim is given by mcmc sample); 
  cat("'Qrc_Prior.spag'", file = file.Pred2, append = TRUE,sep="\n")                    # !!! Files containing spaghettis for each output variable (size nY)
  cat(".true.", file = file.Pred2, append = TRUE,sep="\n")                              #!!! Post-processing: transpose spag file (so that each column is a spaghetti)? 
  cat(".true.", file = file.Pred2, append = TRUE,sep="\n")                              #!!! Post-processing: create envelops? (size nY)
  cat("'Qrc_Prior.env'", file = file.Pred2, append = TRUE,sep="\n")                     # !!! Post-processing: name of envelop files (size nY)
  cat(".true.", file = file.Pred2, append = TRUE,sep="\n")                              #!!! Print progress in console during computations?
  cat(".false." , file = file.Pred2, append = TRUE,sep="\n")                           #!!! Do state prediction? (size nState)
  
  file.Pred3 = paste(dir_code,"/BaM_exe/BaM_BaRatin_2/Config_Pred_RCMaxpost.txt",sep="")
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
  cat(".true.", file = file.Pred3, append = TRUE,sep="\n")                              #!!! Print progress in console during computations?
  cat(".false." , file = file.Pred3, append = TRUE,sep="\n")                           #!!! Do state prediction? (size nState)
  
  file.Pred4 = paste(dir_code,"/BaM_exe/BaM_BaRatin_2/Config_Pred_RCParamU.txt",sep="")
  cat("'BaM_BaRatin_2\\Hgrid.txt'", file =file.Pred4, sep="\n")
  cat(ngrid, file =file.Pred4,sep="\n", append = TRUE)
  cat("1", file = file.Pred4, append = TRUE,sep="\n")   #n of spaghetti
  cat(".true.", file = file.Pred4, append = TRUE,sep="\n")                              #!!! Propagate parametric uncertainty?
  cat(".false.", file = file.Pred4, append = TRUE,sep="\n")                             #!!! Propagate remnant uncertainty for each output variable? (size nY)
  cat("-1", file = file.Pred4, append = TRUE,sep="\n")                                #!!! Nsim[prior]. If <=0: posterior sampling (nsim is given by mcmc sample); 
  cat("'Qrc_ParamU.spag'", file = file.Pred4, append = TRUE,sep="\n")                    # !!! Files containing spaghettis for each output variable (size nY)
  cat(".true.", file = file.Pred4, append = TRUE,sep="\n")                              #!!! Post-processing: transpose spag file (so that each column is a spaghetti)? 
  cat(".true.", file = file.Pred4, append = TRUE,sep="\n")                              #!!! Post-processing: create envelops? (size nY)
  cat("'Qrc_ParamU.env'", file = file.Pred4, append = TRUE,sep="\n")                     # !!! Post-processing: name of envelop files (size nY)
  cat(".true.", file = file.Pred4, append = TRUE,sep="\n")                              #!!! Print progress in console during computations?
  cat(".false." , file = file.Pred4, append = TRUE,sep="\n")                           #!!! Do state prediction? (size nState)
  
  file.Pred5 = paste(dir_code,"/BaM_exe/BaM_BaRatin_2/Config_Pred_RCTotalU.txt",sep="")
  cat("'BaM_BaRatin_2\\Hgrid.txt'", file =file.Pred5,sep="\n")
  cat(ngrid, file =file.Pred5, append = TRUE, sep="\n")
  cat("1", file = file.Pred5, append = TRUE,sep="\n")   #n of spaghetti
  cat(".true.", file = file.Pred5, append = TRUE,sep="\n")                              #!!! Propagate parametric uncertainty?
  cat(".true.", file = file.Pred5, append = TRUE,sep="\n")                             #!!! Propagate remnant uncertainty for each output variable? (size nY)
  cat("-1", file = file.Pred5, append = TRUE,sep="\n")                                #!!! Nsim[prior]. If <=0: posterior sampling (nsim is given by mcmc sample); 
  cat("'Qrc_TotalU.spag'", file = file.Pred5, append = TRUE,sep="\n")                    # !!! Files containing spaghettis for each output variable (size nY)
  cat(".true.", file = file.Pred5, append = TRUE,sep="\n")                              #!!! Post-processing: transpose spag file (so that each column is a spaghetti)? 
  cat(".true.", file = file.Pred5, append = TRUE,sep="\n")                              #!!! Post-processing: create envelops? (size nY)
  cat("'Qrc_TotalU.env'", file = file.Pred5, append = TRUE,sep="\n")                     # !!! Post-processing: name of envelop files (size nY)
  cat(".true.", file = file.Pred5, append = TRUE,sep="\n")                              #!!! Print progress in console during computations?
  cat(".false." , file = file.Pred5, append = TRUE,sep="\n")                           #!!! Do state prediction? (size nState)
}






###############################################################################
Config_predictions <- function(ncontrol,a,k,c,st_a,st_k,st_c,nobs) {
###############################################################################
  file.name1 = paste(dir_code,"/BaM_exe/BaM_BaRatin_2/Config_Model.txt",sep="")
  cat('"BaRatin"', file =file.name1,sep="\n")
  cat(1, file = file.name1, append = TRUE,sep="\n")
  cat(1, file =file.name1, append = TRUE,sep="\n")
  cat(npar, file =file.name1, append = TRUE,sep="\n")
  
  cat('"k1"', file =file.name1, append = TRUE,sep="\n")
  cat(k[1], file =file.name1, append = TRUE,sep="\n")
  cat("'Gaussian'", file =file.name1, append = TRUE,sep="\n")
  cat(k[1],file =file.name1, append = TRUE, sep=",")
  cat(",",file =file.name1, append = TRUE, sep=",")
  cat(st_k[1],file =file.name1, append = TRUE, sep="\n")  
  
  cat('"a1"', file =file.name1, append = TRUE,sep="\n")
  cat(a[1], file =file.name1, append = TRUE,sep="\n")
  cat("'Gaussian'", file =file.name1, append = TRUE,sep="\n")
  cat(a[1],file =file.name1, append = TRUE, sep=",")
  cat(",",file =file.name1, append = TRUE, sep=",")
  cat(st_a[1],file =file.name1, append = TRUE, sep="\n")
  
  cat('"c1"', file =file.name1, append = TRUE,sep="\n")
  cat(c[1], file =file.name1, append = TRUE,sep="\n")
  cat("'Gaussian'", file =file.name1, append = TRUE,sep="\n")
  cat(c[1],file =file.name1, append = TRUE, sep=",")
  cat(",",file =file.name1, append = TRUE, sep=",")
  cat(st_c[1],file =file.name1, append = TRUE, sep="\n")
  
  cat('"k2"', file =file.name1, append = TRUE,sep="\n")
  cat(k[2], file =file.name1, append = TRUE,sep="\n")
  cat("'Gaussian'", file =file.name1, append = TRUE,sep="\n")
  cat(k[2],file =file.name1, append = TRUE, sep=",")
  cat(",",file =file.name1, append = TRUE, sep=",")
  cat(st_k[2],file =file.name1, append = TRUE, sep="\n")
  
  cat('"a2"', file =file.name1, append = TRUE,sep="\n")
  cat(a[2], file =file.name1, append = TRUE,sep="\n")
  cat("'Gaussian'", file =file.name1, append = TRUE,sep="\n")
  cat(a[2],file =file.name1, append = TRUE, sep=",")
  cat(",",file =file.name1, append = TRUE, sep=",")
  cat(st_a[2],file =file.name1, append = TRUE, sep="\n")
  
  cat('"c2"', file =file.name1, append = TRUE,sep="\n")
  cat(c[2], file =file.name1, append = TRUE,sep="\n")
  cat("'Gaussian'", file =file.name1, append = TRUE,sep="\n")
  cat(c[2],file =file.name1, append = TRUE, sep=",")
  cat(",",file =file.name1, append = TRUE, sep=",")
  cat(st_c[2],file =file.name1, append = TRUE, sep="\n")
  
  cat('"k3"', file =file.name1, append = TRUE,sep="\n")
  cat(k[3], file =file.name1, append = TRUE,sep="\n")
  cat("'Gaussian'", file =file.name1, append = TRUE,sep="\n")
  cat(k[3],file =file.name1, append = TRUE, sep=",")
  cat(",",file =file.name1, append = TRUE, sep=",")
  cat(st_k[3],file =file.name1, append = TRUE, sep="\n")
  
  cat('"a3"', file =file.name1, append = TRUE,sep="\n")
  cat(a[3], file =file.name1, append = TRUE,sep="\n")
  cat("'Gaussian'", file =file.name1, append = TRUE, sep="\n")
  cat(a[3],file =file.name1, append = TRUE, sep=",")
  cat(",",file =file.name1, append = TRUE, sep=",")
  cat(st_a[3],file =file.name1, append = TRUE, sep="\n")
  
  cat('"c3"', file =file.name1, append = TRUE,sep="\n")
  cat(c[3], file =file.name1, append = TRUE,sep="\n")
  cat("'Gaussian'", file =file.name1, append = TRUE,sep="\n")
  cat(c[3],file =file.name1, append = TRUE, sep=",")
  cat(",",file =file.name1, append = TRUE, sep=",")
  cat(st_c[3],file =file.name1, append = TRUE, sep="\n")
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
           file.save = dir) #call function from "Prior_Propagation.R"

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

















#***************************************************************************************
Qtimeseries <- function(time,Qt,m){
  #***************************************************************************************
  plot(time,Qt, type = "l", ylim = c(0,150))
  #AVERAGING:
  Qmean = averageh(Qt,m)
  tmean = averaget(time,m)
  plot(tmean[650:1000],Qmean[650:1000], ylim =c(0,150),type ="l")
  return(Qt)
}




  
  
  
  
  
  












############################################################################################################
#BaRatin - SPD (Mansanarez et al., 2019)                   
BaRatin_SPD.bac_app <- function(dir_code, 
                                dir.BaM, 
                                dir.SPD.config, 
                                dir.SPD.results,  
                                save.all.results   = TRUE,
                                plot.results.only  = TRUE
) {
  #############################################################################################################
  # Directories:
  dir.create(paste0(dir.case_study,"/Results/segmentation_gaugings/SPD"))
  dir.SPD.exe <- paste0(dir_code,"/BaM_exe/BaRatin_SPD")
  dir.SPD.segm.gaug.results <- paste0(dir.case_study,"/Results/segmentation_gaugings/SPD")
  write.table(g.s.results.1[,1:4], paste0(dir.SPD.segm.gaug.results,"/Gaugings_data_SPD.txt"), 
              sep ="\t", row.names=FALSE, col.names = c("h","Q","uQ","Period"))
  write.table(g.s.results.1[,1:4], paste0(dir.SPD.exe,"/Gaugings_data_SPD.txt"), 
              sep ="\t", row.names=FALSE, col.names = c("h","Q","uQ","Period")
  )
  #b.priorSPD= c(3, 7.6), st_b.priorSPD= c(1, 1) #Maxau
  #b.priorSPD=  b.prior # Wairau  # c(0, 0, 1.2); Meyras
  #st_b.priorSPD = st_b.prior#c(1, 1,0.2) #Meyras
  b.priorSPD= c(0, 0, 1.2); st_b.priorSPD= c(1, 1,0.2) #Meyras !!!!!!!!!!!!!!!!!
  
  # Read the gaugings data:
  data4BaRatinSPD <- read.table(paste0(dir.SPD.results,"/Gaugings_data_SPD.txt"), header= TRUE) 
  #write.table(data4BaRatinSPD, paste0(dir.SPD.config,"/Gaugings_data_SPD.txt"),sep="\t",row.names=FALSE)
  names.bac=c('b','a','c')
  # Define priors on "physical" parameters:
  b=NULL;c=NULL;Bw=NULL;Cr=NULL;g=NULL;Bc=NULL;KS=NULL;S0=NULL;  # Priors are defined through Monte Carlo samples
  for (i in 1:ncontrol) {
    b[[i]] <- rnorm(nsim, mean=b.prior[i], sd=st_b.prior[i]) # offset (stage for which Q=0 for that control)
    c[[i]] <- rnorm(nsim, mean=c.prior[i], sd=st_c.prior[i]) # exponent
    # a[[i]] <- rnorm(nsim, mean=a.prior[i], sd=st_a.prior[i]) # parameter
    #IN CASE OF BARATIN :
    if (control.type[i] == "rect.weir") {
      Bw[[i]] <- rlnorm(nsim, meanlog=log(Bw.prior[i]), sdlog=st_b.prior[i])  # weir width
      Cr[[i]] <- rlnorm(nsim, meanlog=log(Cr.prior[i]), sdlog=st_Cr.prior[i]) # Discharge coefficient
      g[[i]]  <- rlnorm(nsim, meanlog=log(g.prior[i]),  sdlog=st_g.prior[i])  # gravity
    } else if (control.type[i] == "rect.channel") {
      Bc[[i]] <- rlnorm(nsim, meanlog=log(Bc.prior[i]), sdlog=st_Bc.prior[i]) # channel width
      KS[[i]] <- rlnorm(nsim, meanlog=log(KS.prior[i]), sdlog=st_KS.prior[i]) # Strickler coefficient
      S0[[i]] <- rlnorm(nsim, meanlog=log(S0.prior[i]), sdlog=st_S0.prior[i]) # slope
    }
  }
  # Priors for incremental global changes:
  if (global.change == TRUE) {
    d.g <- matrix(rnorm(nsim*(nperiod-1), mean=dg.prior, sd=st_dg.prior), nrow=nsim, ncol=nperiod-1)
  } else {
    d.g = NULL
  }
  if (local.change == TRUE) {
    d.l <- matrix(rnorm(nsim*(nperiod-1), mean=dl.prior, sd=st_dl.prior), nrow=nsim, ncol=nperiod-1)
  } else {
    d.l =NULL
  }
  # Priors for incremental global changes:
  if (width.change == TRUE) {
    d.Bc1 <- matrix(rlnorm(nsim*(nperiod-1), meanlog=d.Bc1.prior, sd=st_dBc1.prior), nrow=nsim, ncol=nperiod-1)
  } else {
    d.Bc1 = NULL
  }
  # d.b1 <- matrix(rnorm(nsim*(nperiod-1), mean=d.b1.prior, sd=st_db1.prior), nrow=nsim, ncol=nperiod-1)
  # d.a1 <- matrix(rnorm(nsim*(nperiod-1), mean=d.a1.prior, sd=st_d.a1.prior), nrow=nsim, ncol=nperiod-1)
  
  
  # Starting point for all parameters:
  start <- list()
  for (i in 1:ncontrol) {
    if (control.type[i] == "rect.weir") {
      start[[paste0("b",i)]] = b.prior[i] 
      start[[paste0("c",i)]] = c.prior[i]
      start[[paste0("Bw",i)]] = Bw.prior[i]
      start[[paste0("Cr",i)]] = Cr.prior[i]
      start[[paste0("g",i)]] = g.prior[i]
    } else if (control.type[i] == "rect.channel") {
      start[[paste0("b",i)]] =  b.prior[i] 
      start[[paste0("c",i)]]= c.prior[i]
      start[[paste0("Bc",i)]] = Bc.prior[i]
      start[[paste0("KS",i)]] = KS.prior[i]
      start[[paste0("S0",i)]] = S0.prior[i]
    } else if (control.type[i] == "other.function") {
      start[[paste0("b",i)]] =  b.prior[i] 
      start[[paste0("c",i)]]= c.prior[i]
      start[[paste0("a",i)]]= a.prior[i]
    }
  }
  #changes parameters:
  for (i in 1:(nperiod-1)) {
    if (global.change == TRUE) {
      start[[paste0("dg", i)]] = dg.prior[1]
    } 
    if (local.change == TRUE) {
      start[[paste0("dl", i)]] = dl.prior[1]
    }
    if (width.change == TRUE) {
      start[paste0("dBc1", i)] = d.Bc1.prior
    }
  }
  
  # Perform Monte-Carlo propagation   (call function from "Prior_Propagation.R")
  MC = propagation(b, c,
                   Bw,Cr,g,  # rectangular weir
                   Bc,KS,S0, # large rectangular channel
                   #d.b1, d.Bc1, #Wairau
                   d.g, d.l,   # Meyras: global and local incremental changes
                   d.Bc1,
                   start,    # starting vector
                   changes.method, #choice between "cumsum" and "cumprod" methods ofr incremental changes   
                   n.parvar = n.parvar,
                   global.change,
                   local.change,
                   width.change,
                   case_study_name)
  
  # Fit prior distribution on MCMC samples:
  margins=c();names=c();k=0    # define marginal prior distributions and parameter names
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
  # fit multivariate prior distribution :  call function from "Prior_Propagation.R"
  prior=fit(sim=MC$sim, margins=margins, names=names, ncol=6, rowmax=3, file.save = dir.SPD.results)  
  
  # write config files from "Prior_Propagation.R" module
  colperiod=4     # column where period will be stored in BaM data file
  writeConfigFilesSPD(prior =prior, 
                      start=MC$start, 
                      ncontrol=ncontrol,
                      nperiod= nperiod, 
                      colperiod=colperiod, 
                      isVar =isVar, 
                      #isVar=rep(F,ncontrol*nperiod), 
                      model='Config_Model.txt', 
                      correl='PriorCorrelation.txt',
                      names.bac=c('b','a','c'), 
                      directory = dir.SPD.config)
  nobs <- length(data4BaRatinSPD$Q)
  BaRatin_SPD_config(dir.BaM, dir.SPD.config , 
                     pred, nobs, M, 
                     remnant, g1.prior, g2.prior, g1.distr.type, g2.distr.type,
                     Ncycles)
  # Exection of BaM.exe:
  setwd(dir.BaM)
  message("**********************************************************"); flush.console()
  message("Applying BaRatin-SPD (Mansanarez et al., 2019) !!!  Wait ... "); flush.console()
  system2(paste0(dir.BaM,"/BaM_2exp_pool2.exe"), stdout =NULL, stderr =NULL)
  
  
  #####################################################################################################          
  # Plotting results of BaRatin-SPD in terms of rating curves :
  plot.SPD(dir.SPD.exe, 
           dir.SPD.results=dir.SPD.segm.gaug.results, 
           dir.SPD.config = dir.SPD.exe,
           nperiod = nperiod.from.segm.gaugings.1, 
           df.limni,
           FinalColors =  colo , #  c( "red" , "blue" , "chartreuse", "orange" , "darkviolet", "cyan", "saddlebrown" ), #colo, 
           ylim.wind =grid_RC.ylim, xlim.wind =grid_RC.xlim, 
           breaks.lin.x = seq(grid_RC.xlim[1],grid_RC.xlim[2],grid_RC.xstep), 
           labels.lin.x = seq(grid_RC.xlim[1],grid_RC.xlim[2],grid_RC.xstep), 
           breaks.lin.y = seq(grid_RC.ylim[1], grid_RC.ylim[2], grid_RC.ystep), 
           labels.lin.y = seq(grid_RC.ylim[1], grid_RC.ylim[2], grid_RC.ystep),
           ylim.log.wind =c(grid_RC.ylim.log[1], grid_RC.ylim.log[2]),  
           breaks.log= ticks_RC.y.log,   
           labels.log= ticks_RC.y.log,
           case_study_name) 
  
  
  # Boxplots of b1 and b2:
  SPD.summary = read.table(file=paste0(dir.SPD.segm.gaug.results, "/Results_Summary.txt"))
  SPD.mcmc.cooked = read.table(file=paste0(dir.SPD.segm.gaug.results, "/Results_MCMC_Cooked.txt"), header=TRUE)
  # b1:
  SPD.mcmc.cooked.b1 = SPD.mcmc.cooked[,1:nperiod.from.segm.gaugings.1]; 
  names(SPD.mcmc.cooked.b1) = seq(1,nperiod.from.segm.gaugings.1,1)
  SPD.mcmc.cooked.b1 = SPD.mcmc.cooked.b1 %>% gather(period, b1, na.rm = FALSE, convert = FALSE)
  # b2:
  #SPD.mcmc.cooked.b2 = data.frame(SPD.mcmc.cooked[,(nperiod.from.segm.gaugings.1+3)])
  #names(SPD.mcmc.cooked.b2) = 1
  #SPD.mcmc.cooked.b2 = cbind(SPD.mcmc.cooked.b2, replicate(nperiod.from.segm.gaugings.1, SPD.mcmc.cooked.b2$'1'))
  #
  SPD.mcmc.cooked.b2 = SPD.mcmc.cooked[,(nperiod.from.segm.gaugings.1+3 ):(2*nperiod.from.segm.gaugings.1+2)]; 
  names(SPD.mcmc.cooked.b2) = seq(1,nperiod.from.segm.gaugings.1,1)
  SPD.mcmc.cooked.b2 = SPD.mcmc.cooked.b2 %>% gather(period, b2, na.rm = FALSE, convert = FALSE)
  
  SPD.bt = ggplot()+
    geom_boxplot(data=SPD.mcmc.cooked.b1, aes(x= period, y= b1, color="b1"), outlier.size = 0.1, lwd =0.2)+
    geom_boxplot(data=SPD.mcmc.cooked.b2, aes(x= period, y= b2, color="b2"), outlier.size = 0.1, lwd =0.2)+
    stat_summary(fun.y = mean, geom="point", shape=16, size=2)+
    scale_y_continuous(name = TeX("$b_1, b_2$"))+
    scale_color_manual(name = "Legend", labels=c("b1", "b2") , values = c("blue", "red"))+
    scale_x_discrete(name = "Period", limits = seq(1,nperiod.from.segm.gaugings.1, 1) )+
    theme_bw()
  ggsave(SPD.bt, filename = paste0(dir.SPD.segm.gaug.results,"/bt_boxplots.png"), 
         bg = "transparent", device = "png", width = 6, height =4, dpi = 400)
  
  
  
  #df$asDate = as.Date(dtm[,1], "%d.%m.%Y")   #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! DATE FORMAT !!!
  #plotting river bed b(t) from BaRatin-SPD:
  bt.SPD(nperiod         = nperiod.from.segm.gaugings.1,  
         df.limni        = df.limni,  
         dir.SPD.results = dir.SPD.segm.gaug.results,
         t.shift.for.b   = data.annotate.gaug.adjust.1,
         h_G             = g.s.results.1$X.h. , 
         t_G             = g.s.results.1$X.tP., 
         color_G         = g.s.results.1$color, 
         times.uncert    = TRUE, 
         officialShiftsTime = officialShiftsTime ,
         ylimits         = c(-1, 4))#grid_RC.xlim)
  
}  
























###########################################################################################
BaRatin_SPD_config <- function(dir.BaM, dir.SPD.config, pred, nobs, M,                    #
                               remnant, g1.prior, g2.prior,                               #
                               g1.distr.type, g2.distr.type,                              # 
                               Ncycles) {                                                 #
###########################################################################################
  # creation of Config_BaM.txt
  file.BaM <- paste(dir.BaM,"/Config_BaM.txt",sep="")
  cat('"BaRatin_SPD/"', file = file.BaM, sep="\n", append = FALSE)
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
  cat("'BaRatin_SPD\\Gaugings_data_SPD.txt'", file =file.name2, sep="\n")
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
  cat("10.", file = file.matrix,  append = TRUE, sep="\n")  #hmax
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
  cat(".true.", file =file.run, sep="\n")                             # Do MCMC?
  cat(".true.", file =file.run, append = TRUE, sep="\n")              # Do MCMC summary?
  cat(".true.", file =file.run, append = TRUE, sep="\n")              # Do Residual diagnostics?
  if (pred == TRUE) {
    cat(".true.", file =file.run, append = TRUE, sep="\n")         # Do Predictions?
  } else {
    cat(".false.", file =file.run, append = TRUE, sep="\n")        # Do Predictions?
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

































############################################################################################################
# Plotting multiperiod rating curves in linear and log scales:                                             #
plot.SPD <- function(dir.SPD.exe, dir.SPD.results,  dir.SPD.config, nperiod, df.limni, FinalColors,        #
                     ylim.wind, xlim.wind, breaks.lin.x, labels.lin.x, breaks.lin.y, labels.lin.y,         #
                     ylim.log.wind, breaks.log, labels.log, case_study_name) {                             #
############################################################################################################
  # Settings :
  require(ggplot2);require(reshape2);require(scales);require(grid);require(RColorBrewer)
  op_cooked=TRUE
  nburn=0.5; 
  nslim=10
  err.per=0.1; 
  color.inter="#F4A582"; color.interBorder="#F4A582"; color.gaugings="#0571B0"; color.RC="#CA0020";
  color.measurement="#92C5DE"; color.param.inter="#FEE0B6"; color.param.interBorder="#FEE0B6";
  color.param.inter = "#B2ABD2"; color.param.interBorder = "#B2ABD2"; face.wind="bold";
  na.rm=TRUE; shape.gaugings=18; size.gaugings=4; alpha.inter=0.5; 
  size.RC=2; size.wind=20; size.axis=20; shape.meas=18; size.meas=4;
  width.errorbar=0.005; alpha.param.inter=0.5; alpha.inter=0.5; width.wind=30; height.wind=30
  hydr.op=1; time.unit="h"; blank=TRUE; param.plot = TRUE; error.bar.op=NULL;  inter.hydr=NULL;
  blank=FALSE;  param.plot=FALSE; hydr.op=NULL; na.rm.meas=TRUE; time.unit="s"; 
  # Data loading:
  if(op_cooked){
    data.MCMC.cooked=as.matrix(read.table(paste0(dir.SPD.exe,"/Results_MCMC_Cooked.txt")  #### Cooked MCMC
                                          ,header=TRUE,dec=".", sep=""))
  }else{
    data.MCMC.raw=as.matrix(read.table(paste0(dir.SPD.exe,"/Results_MCMC.txt") #### Raw MCMC
                                       ,header=TRUE,dec=".", sep=""))
    nsim.tmp=length(data.MCMC.raw[,1])
    data.MCMC.cooked=data.MCMC.raw[seq(nburn*nsim.tmp, nsim.tmp, nslim),]
    rm("data.MCMC.raw")
  }
  logPost.ncol=which(colnames(data.MCMC.cooked)=="LogPost")
  data.MCMC.MaxPost=as.numeric(read.table(paste(dir.SPD.exe,"/Results_Summary.txt",sep="")
                                          ,row.names=1,dec=".",sep="", skip = 16))
  nsample=length(data.MCMC.cooked[,1])
  
  # Propagation: Convert MCMC into rating curve results
  hgrid=seq(xlim.wind[1], xlim.wind[2]+1, 0.01) 
  
  
  if ((case_study_name == "Meyras") | (case_study_name == "Illinois_Tinley_Creek")) {
    RC_3controls = function(theta, h, op=1){ # Rating curve with discharge error propagation:
      Q = 0*h
      mask1 = which(((h>theta[1])+(h<=theta[13]))==2)
      mask2 = which(((h>theta[13])+(h<=theta[7]))==2)
      mask3 = which(h>theta[7])
      Q[mask1] = theta[2]*(h[mask1]-theta[1])^theta[3] # section control
      Q[mask2] = theta[5]*(h[mask2]-theta[4])^theta[6] # channel control
      Q[mask3] = theta[5]*(h[mask3]-theta[4])^theta[6] +
        theta[8]*(h[mask3]-theta[7])^theta[9] # channel + floodway controls
      #print(c(Q[mask1],Q[mask2],Q[mask3]))
      if(op==1) {
        resQ = sapply(Q,function(Q,theta){Q + rnorm(1, mean=0, sd=theta[1]+theta[2]*Q)}, theta=c(theta[10],theta[11]))
      }else{
        resQ = Q
      }
      return(resQ)
    }
    RC_3controls_mp=function(theta,h){  # Rating curve for maximum posterior:
      Q=0*h
      mask1 = which(((h>theta[1])+(h<=theta[12]))==2)
      mask2 = which(((h>theta[12])+(h<=theta[7]))==2)
      mask3 = which(h>theta[7])
      Q[mask1] = theta[2]*(h[mask1]-theta[1])^theta[3] # section control
      Q[mask2] = theta[5]*(h[mask2]-theta[4])^theta[6] # channel control
      Q[mask3] = theta[5]*(h[mask3]-theta[4])^theta[6] +
        theta[8]*(h[mask3]-theta[7])^theta[9] # channel + floodway controls
      return(Q)
    }
    #Initialisation:
    MCMC.save = matrix(NA, nrow=nperiod*nsample, ncol=14) ### b1 a1 c1 b2 a2 c2 b3 a3 c3 gamma1 gamma2 logpost k2   (nper)
    # 1  2  3  4  5  6  7  8  9    10     11      12   13    (14)
    #please notice that in summary file there is not the logpost column: !!!!
    MaxPost.save = matrix(NA, nrow=nperiod, ncol=13)      ### b1 a1 c1 b2 a2 c2 b3 a3 c3 gamma1 gamma2 k2 (nper)
    # 1  2  3  4  5  6  7  8  9    10     11   12
    for(i in 1:nperiod){
      #for all mcmc:
      col.num = c( i, nperiod+1, nperiod+2,   #b1, a1, c1
                   nperiod+2+i, 2*nperiod+2+1, 2*nperiod+2+2, #b2, a2, c2
                   2*nperiod+4+1, 2*nperiod+4+2, 2*nperiod+4+3, #b3, a3, c3
                   2*nperiod+4+4, 2*nperiod+4+5 , #gamma1, gamma2
                   2*nperiod+4+6 ,  #logpost
                   2*nperiod+4+6+ 6*(i-1)+2 ) #k2
      MCMC.save[(nsample*(i-1)+1):(i*nsample),] = cbind(data.MCMC.cooked[,col.num], rep(i, nsample) )
      #for maxpost:
      col.num.mp = c( i, nperiod+1, nperiod+2, #b1, a1, c1
                      nperiod+2+i, 2*nperiod+2+1, 2*nperiod+2+2, #b2, a2, c2
                      2*nperiod+4+1, 2*nperiod+4+2, 2*nperiod+4+3, #b3, a3, c3
                      2*nperiod+4+4, 2*nperiod+4+5 , #gamma1, gamma2
                      2*nperiod+4+5 + 6*(i-1) +2 ) #k2
      MaxPost.save[i,] = c(data.MCMC.MaxPost[col.num.mp],i)
    }
    RC.Post <-  apply(MCMC.save, MARGIN=1, RC_3controls, h=hgrid) #apply RC for hgrid with structural error:
    RC.MaxPost= apply(MaxPost.save, MARGIN=1, RC_3controls_mp, h=hgrid) # Maximum posterior computing
    data.gaug=read.table(paste(dir.SPD.exe,"/Gaugings_data_SPD.txt",sep=""),header=TRUE,dec=".",sep="") # Gaugings loading
    #data.gaug=data.gaug[-which(data.gaug$Q==-9999),]
  }  
  
  if ((case_study_name == "Wairau")) {
    RC_2controls = function(theta, h, op=1){ # Rating curve with discharge error propagation:
      Q = 0*h
      mask1 = which(((h>theta[1])+(h<=theta[4]))==2)
      mask2 = which(h>theta[4])
      
      Q[mask1] = theta[2]*(h[mask1]-theta[1])^theta[3] # channel control
      Q[mask2] = theta[2]*(h[mask2]-theta[1])^theta[3] + theta[5]*(h[mask2]-theta[4])^theta[6] # channel control + floodways
      
      #print(mask2)
      if(op==1) {
        resQ = sapply(Q,function(Q,theta){Q + rnorm(1, mean=0, sd=theta[1]+theta[2]*Q)}, theta=c(theta[7],theta[8]))
      }else{
        resQ = Q
      }
      return(resQ)
    }
    RC_2controls_mp=function(theta,h){  # Rating curve for maximum posterior:
      Q=0*h
      mask1 = which(((h> theta[1])+(h <= theta[4]))==2)
      mask2 = which(h> theta[4])
      
      Q[mask1] = theta[2]*(h[mask1]-theta[1])^theta[3]
      Q[mask2] = theta[2]*(h[mask2]-theta[1])^theta[3] + theta[5]*(h[mask2]-theta[4])^theta[6] # channel control
      return(Q)
    }
    MCMC.save = matrix(NA, nrow=nperiod*nsample, ncol=10) ### b1 a1 c1 b2 a2 c2 gamma1 gamma2 logpost (nper)
    MaxPost.save = matrix(NA, nrow=nperiod, ncol=9)      ### b1 a1 c1 b2 a2 c2 gamma1 gamma2 (nper)
    for(i in 1:nperiod){
      #for all mcmc:
      col.num = c( i, nperiod+i, 2*nperiod+1,
                   2*nperiod+2, 2*nperiod+3, 2*nperiod+4,
                   2*nperiod+5, 2*nperiod+6 , 2*nperiod+7)
      
      MCMC.save[(nsample*(i-1)+1):(i*nsample),] = cbind(data.MCMC.cooked[,col.num], rep(i, nsample) )
      #for maxpost:
      col.num.mp = c( i, nperiod+i, 2*nperiod +1,
                      2*nperiod+2, 2*nperiod+3, 2*nperiod+4,
                      2*nperiod+5, 2*nperiod+6)
      
      MaxPost.save[i,] = c(data.MCMC.MaxPost[col.num.mp],i)
    }
    RC.Post <-  apply(MCMC.save, MARGIN=1, RC_3controls, h=hgrid) #apply RC for hgrid with structural error:
    RC.MaxPost= apply(MaxPost.save, MARGIN=1, RC_3controls_mp, h=hgrid) # Maximum posterior computing
    data.gaug=read.table(paste(dir.SPD.exe,"/Gaugings_data_SPD.txt",sep=""),header=TRUE,dec=".",sep="") # Gaugings loading
    #data.gaug=data.gaug[-which(data.gaug$Q==-9999),]
  }
  
  if ((case_study_name == "Maxau")) {
    RC_2controls = function(theta, h, op=1){ # Rating curve with discharge error propagation:
      Q = 0*h
      mask1 = which(((h>theta[1])+(h<=theta[4]))==2)
      mask2 = which(h>theta[4])
      
      Q[mask1] = theta[2]*(h[mask1]-theta[1])^theta[3] # channel control
      Q[mask2] = theta[2]*(h[mask2]-theta[1])^theta[3] + theta[5]*(h[mask2]-theta[4])^theta[6] # channel control + floodways
      
      #print(mask2)
      if(op==1) {
        resQ = sapply(Q,function(Q,theta){Q + rnorm(1, mean=0, sd=theta[1]+theta[2]*Q)}, theta=c(theta[7],theta[8]))
      }else{
        resQ = Q
      }
      return(resQ)
    }
    #---------------------------------
    RC_2controls_mp=function(theta,h){  # Rating curve for maximum posterior:
      Q=0*h
      mask1 = which(((h> theta[1])+(h <= theta[4]))==2)
      mask2 = which(h> theta[4])
      
      Q[mask1] = theta[2]*(h[mask1]-theta[1])^theta[3]
      Q[mask2] = theta[2]*(h[mask2]-theta[1])^theta[3] + theta[5]*(h[mask2]-theta[4])^theta[6] # channel control
      return(Q)
    }
    #---------------------------------
    MCMC.save = matrix(NA, nrow=nperiod*nsample, ncol=10) ### b1 a1 c1 b2 a2 c2 gamma1 gamma2 logpost (nper)
    MaxPost.save = matrix(NA, nrow=nperiod, ncol=9)      ### b1 a1 c1 b2 a2 c2 gamma1 gamma2 (nper)
    for(i in 1:nperiod){
      #for all mcmc:
      col.num = c( i, nperiod+1, nperiod+2,
                   nperiod+3, nperiod+4, nperiod+5,
                   nperiod+6, nperiod+7 , nperiod+8)
      
      MCMC.save[(nsample*(i-1)+1):(i*nsample),] = cbind(data.MCMC.cooked[,col.num], rep(i, nsample) )
      #for maxpost:
      col.num.mp = c( i, nperiod+1, nperiod +2,
                      nperiod+3, nperiod+4, nperiod+5,
                      nperiod+6, nperiod+7)
      
      MaxPost.save[i,] = c(data.MCMC.MaxPost[col.num.mp],i)
    }
    RC.Post <-  apply(MCMC.save, MARGIN=1, RC_2controls, h=hgrid) #apply RC for hgrid with structural error:
    RC.MaxPost= apply(MaxPost.save, MARGIN=1, RC_2controls_mp, h=hgrid) # Maximum posterior computing
    data.gaug=read.table(paste(dir.SPD.exe,"/Gaugings_data_SPD.txt",sep=""),header=TRUE,dec=".",sep="") # Gaugings loading
    #data.gaug=data.gaug[-which(data.gaug$Q==-9999),]
  }
  
  #Escalier_Mat_section_channel: 2 controls
  # # #***************************************
  # RC_2controls = function(theta, h, op=1){ # Rating curve with discharge error propagation:
  #   #***************************************
  #   Q = 0*h
  #   mask1 = which(((h>theta[1])+(h<=theta[9]))==2)
  #   mask2 = which(h>theta[9])
  #   Q[mask1] = theta[2]*(h[mask1]-theta[1])^theta[3] # section control
  #   Q[mask2] = theta[5]*(h[mask2]-theta[4])^theta[6] # channel control
  #   if(op==1) {
  #     resQ = sapply(Q,function(Q,theta){Q + rnorm(1, mean=0, sd=theta[1]+theta[2]*Q)}, theta=c(theta[7],theta[8]))
  #   }else{
  #     resQ = Q
  #   }
  #   return(resQ)
  # }
  # 
  # #***********************************
  # RC_2controls_mp=function(theta,h){  # Rating curve for maximum posterior:
  #   #************************************
  #   Q=0*h
  #   mask1 = which(((h>theta[1])+(h<=theta[9]))==2)
  #   mask2 = which(h>theta[9])
  #   Q[mask1] = theta[2]*(h[mask1]-theta[1])^theta[3] # section control
  #   Q[mask2] = theta[5]*(h[mask2]-theta[4])^theta[6] # channel control
  #   return(Q)
  # }
  # 
  # #Initialisation:
  # #******************
  # #MCMC.save = matrix(NA, nrow=nperiod*nsample, ncol=10) ### b1 a1 c1 b2 a2 c2 gamma1 gamma2  k2 (nper)
  # MCMC.save = matrix(NA, nrow=nperiod*nsample, ncol=10)  ### b1 a1 c1 b2 a2 c2 gamma1 gamma2  k2 (nper)
  # # 1  2  3  4  5  6  7  8  9    10     11      12   13    (14)
  # #MaxPost.save = matrix(NA, nrow=nperiod, ncol=10)      ### b1 a1 c1 b2 a2 c2 gamma1 gamma2 k2 (nper)
  # MaxPost.save = matrix(NA, nrow=nperiod, ncol=10)      ### b1 a1 c1 b2 a2 c2 gamma1 gamma2 k2 (nper)
  # # 1  2  3  4  5  6  7  8  9    10     11   12
  # 
  # for(i in 1:nperiod){
  #   #***********************************************************
  #   #for all mcmc:
  #   #***********************************************************
  #   # col.num = c( i, nperiod+1, nperiod+2,
  #   #              nperiod+2+i, 2*nperiod+3, 2*nperiod+4,
  #   #              2*nperiod+5, 2*nperiod+6,
  #   #              2*nperiod+6+4*(i-1)+2 )
  #   col.num = c( i, nperiod+1, nperiod+2,
  #                nperiod+2+i, nperiod+2+nperiod+1, nperiod+2+nperiod+2,
  #                2*nperiod+2+2+1, 2*nperiod+2+2+2,
  #                2*nperiod+2+2+3+4*(i-1)+2)
  # 
  #   MCMC.save[(nsample*(i-1)+1):(i*nsample),] = cbind(data.MCMC.cooked[,col.num], rep(i, nsample) )
  # 
  #   #***********************************************************
  #   #for maxpost:
  #   #***********************************************************
  #   # col.num.mp =  c( i, nperiod+1, nperiod+2,
  #   #                  nperiod+2+i, 2*nperiod+3, 2*nperiod+4,
  #   #                  2*nperiod+5, 2*nperiod+6,
  #   #                  2*nperiod+6+4*(i-1)+2)
  #   col.num.mp = c( i, nperiod+1, nperiod+2,
  #                   nperiod+2+i, nperiod+2+nperiod+1, nperiod+2+nperiod+2,
  #                   2*nperiod+2+2+1, 2*nperiod+2+2+2,
  #                   2*nperiod+2+2+3+4*(i-1)+1 )
  #   MaxPost.save[i,] = c(data.MCMC.MaxPost[col.num.mp],i)
  # }
  # RC.Post <-  apply(MCMC.save, MARGIN=1, RC_2controls, h=hgrid) #apply RC for hgrid with structural error:
  # RC.MaxPost= apply(MaxPost.save, MARGIN=1, RC_2controls_mp, h=hgrid) # Maximum posterior computing
  # data.gaug=read.table(paste(dir.SPD.exe,"/Gaugings_data_SPD.txt",sep=""),header=TRUE,dec=".",sep="") # Gaugings loading
  # #data.gaug=data.gaug[-which(data.gaug$Q==-9999),]
  
  if ((case_study_name == "Escalier_Mat_section_channel")) {
    #****************************************************************************************
    # #Escalier Mat: 2 controls
    RC_2controls = function(theta, h, op=1){ # Rating curve with discharge error propagation:
      Q = 0*h
      mask1 = which(((h>theta[1])+(h<=theta[9]))==2)
      mask2 = which(h>theta[9])
      Q[mask1] = theta[2]*(h[mask1]-theta[1])^theta[3] # channel control
      Q[mask2] = theta[5]*(h[mask2]-theta[4])^theta[6] # channel control
      if(op==1) {
        resQ = sapply(Q,function(Q,theta){Q + rnorm(1, mean=0, sd=theta[1]+theta[2]*Q)}, theta=c(theta[7],theta[8]))
      }else{
        resQ = Q
      }
      return(resQ)
    }
    RC_2controls_mp=function(theta,h){  # Rating curve for maximum posterior:
      Q=0*h
      mask1 = which(((h>theta[1])+(h<=theta[9]))==2)
      mask2 = which(h>theta[9])
      Q[mask1] = theta[2]*(h[mask1]-theta[1])^theta[3] # channel control
      Q[mask2] = theta[5]*(h[mask2]-theta[4])^theta[6] # channel control
      return(Q)
    }
    #MCMC.save = matrix(NA, nrow=nperiod*nsample, ncol=10) ### b1 a1 c1 b2 a2 c2 gamma1 gamma2  k2 (nper)
    MCMC.save = matrix(NA, nrow=nperiod*nsample, ncol=10)  ### b1 a1 c1 b2 a2 c2 gamma1 gamma2  k2 (nper)
    # 1  2  3  4  5  6  7  8  9    10     11      12   13    (14)
    #MaxPost.save = matrix(NA, nrow=nperiod, ncol=10)      ### b1 a1 c1 b2 a2 c2 gamma1 gamma2 k2 (nper)
    MaxPost.save = matrix(NA, nrow=nperiod, ncol=10)      ### b1 a1 c1 b2 a2 c2 gamma1 gamma2 k2 (nper)
    # 1  2  3  4  5  6  7  8  9    10     11   12
    for(i in 1:nperiod){
      #for all mcmc:
      # col.num = c( i, nperiod+1, nperiod+2,
      #              nperiod+2+i, 2*nperiod+3, 2*nperiod+4,
      #              2*nperiod+5, 2*nperiod+6,
      #              2*nperiod+6+4*(i-1)+2 )
      col.num = c( i, nperiod+1, nperiod+2,
                   nperiod+3, nperiod+4, nperiod+5,
                   nperiod+6, nperiod+7,
                   nperiod+8+4*(i-1)+2)
      # col.num = c( i, nperiod+1, nperiod+2,
      #                nperiod+2+i, nperiod+2+nperiod+1, nperiod+2+nperiod+2,
      #                2*nperiod+2+2+1, 2*nperiod+2+2+2,
      #                2*nperiod+2+2+3+4*(i-1)+2)
      
      MCMC.save[(nsample*(i-1)+1):(i*nsample),] = cbind(data.MCMC.cooked[,col.num], rep(i, nsample) )
      #for maxpost:
      # col.num.mp =  c( i, nperiod+1, nperiod+2,
      #                  nperiod+2+i, 2*nperiod+3, 2*nperiod+4,
      #                  2*nperiod+5, 2*nperiod+6,
      #                  2*nperiod+6+4*(i-1)+2)
      col.num.mp = c( i, nperiod+1, nperiod+2,
                      nperiod+3, nperiod+4, nperiod+5,
                      nperiod+6, nperiod+7,
                      nperiod+7+4*(i-1)+2)
      # col.num.mp = c( i, nperiod +1, nperiod+2,
      #                 nperiod+2+i, nperiod+2+nperiod+1, nperiod+2+nperiod+2,
      #                 2*nperiod+2+2+1, 2*nperiod+2+2+2,
      #                 2*nperiod+2+2+3+4*(i-1)+1 )
      MaxPost.save[i,] = c(data.MCMC.MaxPost[col.num.mp],i)
    }
    RC.Post <-  apply(MCMC.save, MARGIN=1, RC_2controls, h=hgrid) #apply RC for hgrid with structural error:
    RC.MaxPost= apply(MaxPost.save, MARGIN=1, RC_2controls_mp, h=hgrid) # Maximum posterior computing
    data.gaug=read.table(paste(dir.SPD.exe,"/Gaugings_data_SPD.txt",sep=""),header=TRUE,dec=".",sep="") # Gaugings loading
    #data.gaug=data.gaug[-which(data.gaug$Q==-9999),]
  }
  
  # Quantiles for Figure RC: 
  message("Plotting BaRatin-SPD !!!  Wait ... "); flush.console()
  List.RC.quants = list(NULL)
  for(i in 1:nperiod){
    data.tmp = apply(RC.Post[,(nsample*(i-1)+1):(nsample*i)], MARGIN=1, quantile, probs=c(0.025,0.975) ) #, na.rm=TRUE)
    data.tmp = apply(data.tmp, MARGIN=c(1,2), function(x){ifelse(x<0,0,x)})
    List.RC.quants[[i]] = data.frame(cbind(hgrid, t(data.tmp), RC.MaxPost[,i]))
    colnames(List.RC.quants[[i]]) = c("h", "inf", "sup", "maxpost")
  }
  inter.per=seq(1,nperiod,1)   # inter.per=c(25:35)
  nRC=length(inter.per)
  # Palette for plot:
  #palette.per= grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)] #brewer.pal(nRC,"Spectral")
  #palette.per= colorRampPalette(c("red","orange","yellow","green","blue","grey","purple")) 
  # palette.per = distinctColorPalette(k = nperiod, altCol = FALSE, runTsne = FALSE)
  # palette.per = c("red","blue","green","orange","brown3","lightblue","black","pink", "yellow","chartreuse", "brown", "violet", "grey","blue4",
  #                 "cyan", "blueviolet", "cyan4", "darkgreen")
  palette.per = FinalColors
  palette.per = palette.per[1:nperiod]
  data.plot.RC = data.frame(do.call("rbind", List.RC.quants[inter.per]), Period=rep(inter.per,each=length(hgrid)))
  inter.null=which(data.plot.RC$maxpost==0)
  # Initialise the plot:
  err.per=0.1; color.inter="#F4A582"; color.interBorder="#F4A582"; alpha.inter=0.5; color.gaugings="#0571B0"
  na.rm=TRUE; shape.gaugings=18; size.gaugings=4; color.RC="#CA0020"; size.RC=2; size.wind=20; size.axis=20
  face.wind="bold"; hydr.op=NULL; color.measurement="#92C5DE"; na.rm.meas=TRUE; shape.meas=18; size.meas=4
  error.bar.op=NULL; width.errorbar=0.005; time.unit="s"; inter.hydr=NULL;blank=FALSE; width.wind=20
  height.wind=20; param.plot=FALSE; color.param.inter="#FEE0B6"; color.param.interBorder="#FEE0B6"
  alpha.param.inter=0.5; hydr.op=1;time.unit="h"; blank=TRUE; alpha.inter=0.5; width.wind=40; height.wind=20
  param.plot = TRUE; color.param.inter = "#B2ABD2"; color.param.interBorder = "#B2ABD2";
  vect.break=seq(1,nperiod,1) #seq(1,9,1)
  vect.emp=rep("",22)
  breaks.wind=c(vect.break*0.01,vect.break*0.1,vect.break,vect.break*10,vect.break*100,vect.break*1000)
  labels.wind=c(0.01,vect.emp,0.1,vect.emp,1,vect.emp,10,vect.emp,100,vect.emp,1000,vect.emp)
  pos.num=function(x.int){inter=which(x.int==data.gaug$period);return(inter)}
  inter.gaug=unlist(sapply(inter.per,pos.num), recursive = TRUE, use.names = TRUE)
  # data for the plot:
  data.obs=data.gaug[inter.gaug,]
  data.obs$period=factor(data.obs$period)
  data.RC =data.plot.RC #[-inter.null,]
  data.RC$Period=factor(data.RC$Period)
  data.gaug=read.table(paste(dir.SPD.exe,"/Gaugings_data_SPD.txt",sep=""),header=TRUE,dec=".",sep="") # Gaugings loading
  #data.gaug=read.table(paste(dir.SPD.exe,"/Gaugings_data_SPD_official.txt",sep=""),header=TRUE,dec=".",sep="") # Gaugings loading
  data.gaug$Period=as.factor(data.gaug$Period)
  # save results of the RCs for each stable period in one file with Maxpost, Q2.5, Q97.5 total
  write.table( data.RC, paste0(dir.SPD.results,"/RC_SPD_env.txt"), sep ="\t", row.names=FALSE)
  
  # RC ggplot::
  plot.RC= ggplot(data.RC)
  if ((case_study_name == "Meyras") | (case_study_name == "Illinois_Tinley_Creek")) {
    # plot.RC = plot.RC + 
    # geom_vline(xintercept = MaxPost.save[inter.per,9],colour=palette.per,size=1)+
    # geom_vline(xintercept = MaxPost.save[inter.per,1], colour=palette.per, size=1, linetype="longdash")
  } else if ( case_study_name == "Escalier") {
    plot.RC = plot.RC + 
      geom_vline(xintercept = MaxPost.save[inter.per,9],colour=palette.per,size=0.6)+
      geom_vline(xintercept = MaxPost.save[inter.per,1], colour=palette.per, size=0.6, linetype="longdash")
  }
  plot.RC = plot.RC + 
    geom_smooth(aes(x=h, y=maxpost, ymax=sup, ymin=inf, group=Period, fill=Period), 
                size=1, stat='identity', alpha=0.2) +  #alpha=0.1
    geom_path(aes(x=h, y=maxpost, group=Period, colour=Period), na.rm=TRUE, size=1)+
    ### Gaugings
    geom_linerange(aes(x=h, ymax=Q+2*uQ, ymin=Q-2*uQ, colour=Period), data=data.gaug, size=0.8)+
    geom_point(aes(x=h,y=Q,colour=Period),data=data.gaug, na.rm=na.rm, shape=16, size=size.gaugings)+
    ### Labels
    xlab(expression(paste("Stage h [m]",sep="")))+
    ylab(expression(paste("Discharge Q [",m^3,".",s^-1,"]",sep="")))+
    labs(colour = "Period")+
    scale_colour_manual(values = palette.per)+
    scale_fill_manual(values = palette.per)+
    #scale_y_continuous(breaks=breaks.wind,labels=labels.wind)+
    coord_cartesian(ylim=ylim.wind, xlim=xlim.wind)+
    scale_x_continuous(breaks=breaks.lin.x, labels=labels.lin.x, expand=c(0,0))+
    scale_y_continuous(breaks=breaks.lin.y, labels=labels.lin.y, expand=c(0,0))+
    ### Theme
    theme_bw(base_size=20)+
    theme( axis.text=element_text(size=15)
           ,axis.title=element_text(size=20,face="bold")
           ,panel.grid.major=element_line(size=1.2)
           ,panel.grid.minor=element_line(size=0.8)
           ,legend.text=element_text(size=20)
           ,legend.title=element_text(size=30)
           ,plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm")
           ,legend.key.size=unit(1.5, "cm")
           ,legend.position="none")
  ggsave(plot.RC, filename=paste0(dir.SPD.results,"/RC_SPD.png"), 
         device = "png", width = 10, height =8, dpi = 400, units = "in")
  # RC plot Logarithmic scale:
  plot.RClog = plot.RC + 
    coord_cartesian(ylim=ylim.log.wind, xlim=xlim.wind) +
    scale_y_log10(breaks=breaks.log, labels=labels.log, expand=c(0,0))+
    annotation_logticks(base = 10, sides = "l", scaled = TRUE,colour = "black", 
                        size = 1, linetype = 1)
  # ggsave(plot.RClog, filename=paste0(dir.SPD.results,"/RClog_SPD.png"), 
  # device = "png", width = 14, height =8, dpi = 400, units = "in")
  pdf(paste0(dir.SPD.results,"/RClog_SPD.pdf"), 14, 8 ,useDingbats=F)
  print(plot.RClog)
  dev.off()
  
  #save results:
  list.files.pool <- c( paste0(dir.SPD.exe,"/Results_MCMC_Cooked.txt"),
                        paste0(dir.SPD.exe,"/Results_Residuals.txt"),
                        paste0(dir.SPD.exe,"/Results_Summary.txt"),
                        paste0(dir.SPD.exe,"/Config_Model.txt"),
                        paste0(dir.SPD.exe,"/Gaugings_data_SPD.txt"))
  #paste(dir.BaM.rec,"/tgrid.txt", sep=""),paste(dir.BaM.rec,"/ht_Maxpost.spag", sep=""),
  #paste(dir.BaM.rec,"/ht_ParamU.spag", sep=""),paste(dir.BaM.rec,"/ht_ParamU.env", sep=""),
  #paste(dir.BaM.rec,"/ht_TotalU.env", sep=""),paste(dir.BaM.rec,"/ht_TotalU.spag", sep=""))
  for (i in 1:length(list.files.pool)) {
    file.copy(list.files.pool[i], dir.SPD.results ,overwrite = TRUE)
  }
  setwd(dir.SPD.results)
  Results_summary.SPD <- read.table(paste(dir.SPD.results,"/Results_Summary.txt",sep=""), header =TRUE)
  Results_mcmc.SPD <- read.table(paste(dir.SPD.results,"/Results_MCMC_Cooked.txt",sep=""), header =TRUE)
  Results_residuals.SPD <- read.table(paste(dir.SPD.results,"/Results_Residuals.txt",sep=""), header =TRUE)
  
  #Save MCMC plots:
  message("Plotting mcmc of SPD !!!  Wait ... "); flush.console()
  setwd(dir.SPD.results)
  mcmc=MCMCplot(doLogPost=F,doDPar=T, MCMCfile= "Results_MCMC_Cooked.txt" , type="trace")
  ggsave(mcmc, filename =paste0(dir.SPD.results,"/mcmc_SPD.png"), 
         bg = "transparent", width = 15, height =25, dpi = 200)
} 























#########################################################################################################
bt.SPD <- function(nperiods, df.limni, dir.SPD.results,                                                 #
                   t.shift.for.b, 
                   h_G, t_G, color_G,                                       
                   times.uncert,                                                                        #
                   officialShiftsTime, 
                   ylimits) {                                                                           #
  #########################################################################################################
  
  
  #Reading and plotting b1 and b2:
  #*******************************
  # Boxplots of b1 and b2:
  SPD.summary = read.table(file=paste0(dir.SPD.results, "/Results_Summary.txt"))
  SPD.mcmc.cooked = read.table(file=paste0(dir.SPD.results, "/Results_MCMC_Cooked.txt"), header=TRUE)
  ##################
  # b1:
  SPD.mcmc.cooked.b1 = SPD.mcmc.cooked[,1:nperiods]; 
  names(SPD.mcmc.cooked.b1) = seq(1,nperiods,1)
  SPD.mcmc.cooked.b1 = SPD.mcmc.cooked.b1 %>% gather(period, b1, na.rm = FALSE, convert = FALSE)
  
  bt1.mcmc = SPD.mcmc.cooked[,1:nperiods]
  bt1.summary = SPD.summary[,1:nperiods]
  bt1.df = data.frame(  t(c(quantile(bt1.mcmc[,1], p = c(0.025, 0.5, 0.975)),
                            mean=    mean(bt1.mcmc[,1]), 
                            stdev =   std(bt1.mcmc[,1]), 
                            maxpost = bt1.summary[16,1])))
  names(bt1.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
  for (i in 2:nperiods){
    bt1.df = rbind(bt1.df, t(c(quantile(bt1.mcmc[,i], p = c(0.025, 0.5, 0.975)), 
                               mean = mean(bt1.mcmc[,i]), 
                               stdev =std(bt1.mcmc[,i]), 
                               maxpost = bt1.summary[16,i]) ))
  }
  write.table(bt1.df, paste0(dir.SPD.results,"/bt1_df.txt"), 
              sep ="\t", row.names=FALSE)
  ##################
  SPD.mcmc.cooked.b2 = SPD.mcmc.cooked[,(nperiods+3 ):(2*nperiods+2)]; 
  names(SPD.mcmc.cooked.b2) = seq(1,nperiods,1)
  SPD.mcmc.cooked.b2 = SPD.mcmc.cooked.b2 %>% gather(period, b2, na.rm = FALSE, convert = FALSE)
  
  bt2.mcmc = SPD.mcmc.cooked[,(nperiods+3 ):(2*nperiods+2)]; 
  bt2.summary = SPD.summary[,(nperiods+3 ):(2*nperiods+2)]; 
  bt2.df = data.frame(  t(c(quantile(bt2.mcmc[,1], p = c(0.025, 0.5, 0.975)),
                            mean =    mean(bt2.mcmc[,1]), 
                            stdev =   std(bt2.mcmc[,1]), 
                            maxpost = bt2.summary[16,1])))
  names(bt2.df) = c("2.5%", "50%", "97.5%", "mean", "stdev" , "maxpost")
  for (i in 2:nperiods){
    bt2.df = rbind(bt2.df, t(c(quantile(bt2.mcmc[,i], p = c(0.025, 0.5, 0.975)), 
                               mean = mean(bt2.mcmc[,i]), 
                               stdev =std(bt2.mcmc[,i]), 
                               maxpost = bt2.summary[16,i]) ))
  }
  write.table(bt2.df, paste0(dir.SPD.results,"/bt2_df.txt"), 
              sep ="\t", row.names=FALSE)
  #boxplots of b1 and b2 parameters for each period:
  bt.boxplot = ggplot()+
    geom_boxplot(data=SPD.mcmc.cooked.b1, aes(x= period, y= b1, color="b1"), outlier.size = 0.1, lwd =0.2)+
    geom_boxplot(data=SPD.mcmc.cooked.b2, aes(x= period, y= b2, color="b2"), outlier.size = 0.1, lwd =0.2)+
    stat_summary(fun.y = mean, geom="point", shape=16, size=2)+
    scale_y_continuous(name = TeX("$b_1, b_2$"))+
    scale_color_manual(name = "Legend", labels=c("b1", "b2") , values = c("blue", "red"))+
    scale_x_discrete(name = "Period", limits = seq(1,nperiods, 1) )+
    theme_bw()
  ggsave(bt.boxplot, filename = paste0(dir.SPD.results,"/bt_boxplots.png"), 
         bg = "transparent", device = "png", width = 6, height =4, dpi = 400)
  # df$asDate = as.Date(dtm[,1], "%d.%m.%Y")   #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! DATE FORMAT !!!
  
  # Read shift times:
  shifts.all = t.shift.for.b$t.adj
  
  if (times.uncert ==TRUE) {   # /!\ this has to be changed !!!!!!!!!!!!!!
    t.shifts.before = c(0,t.shift.for.b$t.adj)
    t.shifts.plus = c(t.shift.for.b$t.adj, tail(df.limni$t_limni,1))
    bt2.df =cbind(bt2.df, 
                  t.shifts.before = t.shifts.before , 
                  t.shifts.plus=t.shifts.plus)
    
    bt1.df =cbind(bt1.df, 
                  t.shifts.before = t.shifts.before , 
                  t.shifts.plus=    t.shifts.plus)
  } else {
    t.shifts = shifts[-length(shifts)]
    t.shifts.before = c(0,t.shifts)
    t.shifts.plus = sort(c(t.shifts, tail(df.limni$t_limni,1)))
  }
  
  # index.hmax=0
  # for (aaaa in 1:length(t.shift.for.b$t.adj)) {
  #   for (oooo in 2:length(df.limni$h_limni)) {
  #     if ((df.limni$t_limni[oooo-1] <= t.shift.for.b$t.adj) & (df.limni$t_limni[oooo] >= t.shift.for.b$t.adj)) {
  #         index.hmax[aaaa] = oooo
  #     }
  #   }
  # }
  # 
  # PLot of b(t) with limni
  #########################
  bt.plot = ggplot()
  if (times.uncert ==TRUE) {
    #geom_line(data = df.bt, aes(x=t.shifts.before, y =bt.MAP), color = "red", size = 1 )+
    bt.plot = bt.plot +
      # geom_rect(mapping= aes(xmin= shifts$t2, xmax=shifts$t97 ,ymin=-Inf, ymax=Inf), 
      #   fill="blue", alpha=0.1) +
      #geom_vline(aes(xintercept = shifts$tMAP), color = "blue", lwd =0.5, linetype = "solid")+
      
      #geom_vline(aes(xintercept = bt2.df$t.shifts.before[-1]), color = "red", lwd =0.3, linetype = "dashed")+
      # geom_segment(data = bt2.df, mapping= aes(x =t.shifts.before , 
      #                                          y = maxpost, 
      #                                          xend = t.shifts.plus, 
      #                                          yend = maxpost),color = "red", size = 1) +
      # geom_rect(data = bt2.df , mapping = aes(xmin= t.shifts.before, 
      #                                         xmax=t.shifts.plus,
    #                                         ymin= `2.5%`, 
    #                                         ymax= `97.5%`), 
    #           fill="red", alpha=0.3) +
    geom_segment(data = bt1.df, mapping= aes(x =t.shifts.before , 
                                             y = maxpost, 
                                             xend = t.shifts.plus, 
                                             yend = maxpost), 
                 color = unique(color_G), 
                 size = 1) +
      # geom_rect(data = bt1.df , mapping = aes(xmin= t.shifts.before, 
      #                                           xmax=t.shifts.plus,
      #                                           ymin= `2.5%`, 
      #                                           ymax= `97.5%`), 
      #           fill=color_G, alpha=0.3) +
      geom_line(data = df.limni, aes(x = t_limni, y= h_limni), color="gray40", size =0.3)+
      geom_point(aes(x=t_G, y = h_G), color= color_G, size=4)+
      scale_y_continuous(name = expression("Stage h [m]"), limits =ylimits, expand = c(0,0)) + 
      scale_x_continuous(name = expression("Time [days]"), expand = c(0,0))+
      theme_bw(base_size=20)+
      theme(axis.text=element_text(size=15),
            axis.title=element_text(size=20,face="bold")
            ,legend.text=element_text(size=20),legend.title=element_text(size=30)
            ,legend.key.size=unit(1.5, "cm"),legend.position="none"
            ,panel.grid = element_blank())+
      geom_segment(aes(x =t.shift.for.b$t.adj ,   
                       y = bt1.df$maxpost[- length(bt1.df$maxpost)] ,
                       xend = t.shift.for.b$t.adj, 
                       yend = bt1.df$maxpost[-1]),
                   arrow = arrow(length = unit(0.2, "cm"), ends = "both"), 
                   size =0.8,
                   color= "black")
    # geom_point(aes(x = t.shift.for.b$t.adj, 
    #                y = df.limni$h_limni[index.hmax]),
    #              size = 10, color="red", pch=21, fill="NA")
  } else {
    bt.plot = bt.plot +
      geom_vline(aes(xintercept = shifts.all), color = "blue", lwd =0.7, linetype = "solid") +
      geom_segment(mapping = aes(x =t.shifts.before , y = bt.MAP, 
                                 xend = t.shifts.plus, yend = bt.MAP), color = "red", size = 1) +
      geom_rect(mapping = aes(xmin= t.shifts.before, xmax=t.shifts.plus, ymin=bt.q10,
                              ymax= bt.q90), fill="red", alpha=0.3) +
      geom_line(data = df.limni, aes(x = t_limni, y= h_limni), color="black", size =0.3)+
      scale_y_continuous(name = expression("Stage h [m]"), limits =ylimits, expand = c(0,0)) + 
      scale_x_continuous(name = expression("Time [days]"), expand = c(0,0))+
      geom_point(aes(x=t_G, y = h_G), color=color_G, size=4)+
      theme_bw(base_size=20)+
      theme(axis.text=element_text(size=15),axis.title=element_text(size=20,face="bold")
            ,legend.text=element_text(size=20),legend.title=element_text(size=30)
            ,legend.key.size=unit(1.5, "cm"),legend.position="none"
            ,panel.grid = element_blank())
  }
  # ggsave(bt.plot, filename = paste0(dir.SPD.results,"/bt.png"), 
  #        width = 16, height =9, dpi = 400)
  pdf(paste0(dir.SPD.results,"/bt.pdf"), 16, 9 ,useDingbats=F)
  print(bt.plot)
  dev.off()
  
  # return the parameter b statistics
  return(list(bt1.df=bt1.df, 
              bt2.df=bt2.df,
              t.shift.for.b=t.shift.for.b))
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






#########################################################################################################
BaRatin_configNspag <- function(dir, 
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
                           nspag) {
  #########################################################################################################
  Hmax = 10  ##### CAREFUL !!!!!!!!!!!!!
  
  
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
  file.name2 = paste0(dir_code,"/BaM_exe/BaM_BaRatin_2/Config_Data.txt")
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
  file.mcmc = paste(dir_code,"/BaM_exe/BaM_BaRatin_2/Config_MCMC.txt",sep="")
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
  file.matrix = paste(dir_code,"/BaM_exe/BaM_BaRatin_2/Config_ControlMatrix.txt",sep="")
  write.table(M, file =file.matrix, row.names = FALSE, col.names = FALSE)
  cat(Hmax, file =file.matrix, append = TRUE,sep="\n")
  #---------------------------------------------------------------------- 
  file.remnant = paste(dir_code,"/BaM_exe/BaM_BaRatin_2/Config_RemnantSigma.txt",sep="")
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
  file_BaM <- paste(dir_code,"/BaM_exe/Config_BaM.txt",sep="")
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
  file.Pred1 = paste(dir_code,"/BaM_exe/BaM_BaRatin_2/Config_Pred_Master.txt",sep="")
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
    file.Pred11 = paste0(dir_code,"/BaM_exe/BaM_BaRatin_2/Config_Pred_Prior.txt")
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
    file.Pred21 = paste0(dir_code,"/BaM_exe/BaM_BaRatin_2/Config_Pred_Prior_Qt.txt")
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
    file.Pred3 = paste(dir_code,"/BaM_exe/BaM_BaRatin_2/Config_Pred_RCMaxpost.txt",sep="")
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
    file.Pred4 = paste(dir_code,"/BaM_exe/BaM_BaRatin_2/Config_Pred_RCParamU.txt",sep="")
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
    file.Pred5 = paste(dir_code,"/BaM_exe/BaM_BaRatin_2/Config_Pred_RCTotalU.txt",sep="")
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
    
    file.Pred6 = paste0(dir_code,"/BaM_exe/BaM_BaRatin_2/Config_Pred_Maxpost.txt")
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
    
    
    file.Pred8 = paste0(dir_code,"/BaM_exe/BaM_BaRatin_2/Config_Pred_ParamU.txt")
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
    
    
    file.Pred7 = paste(dir_code,"/BaM_exe/BaM_BaRatin_2/Config_Pred_TotalU.txt",sep="")
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
  file.run = paste0(dir_code,"/BaM_exe/BaM_BaRatin_2/Config_RunOptions.txt")
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
  file.residuals = paste(dir_code,"/BaM_exe/BaM_BaRatin_2/Config_Residuals.txt",sep="")
  cat('"Results_Residuals.txt"' , file =file.residuals ,sep="\n")      # Result file
  
  ###################################################################   SUMMARY CONFIG
  file.summary = paste(dir_code,"/BaM_exe/BaM_BaRatin_2/Config_Summary.txt",sep="")
  cat('"Results_Summary.txt"' , file =file.summary ,sep="\n")    #Summary stat results file name
  
  ###################################################################   COOKING CONFIG
  file.cooking = paste(dir_code,"/BaM_exe/BaM_BaRatin_2/Config_Cooking.txt",sep="")
  cat('"Results_MCMC_Cooked.txt"' , file =file.cooking ,sep="\n")  # mcmc Results file name
  cat(0.5, file =file.cooking, append = TRUE, sep="\n")            # Burn factor
  cat(10, file =file.cooking, append = TRUE, sep="\n")             # Nslim
}









