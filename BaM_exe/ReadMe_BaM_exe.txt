The uncertainty quantification at mostly all steps of the tools in BayDERS is performed through BaM (BAyesian Modelling), developed by Benjamin Renard (INRAE) with the use of a few Fortran libraries developed by Dmitri Kavetski (University of Adelaide).

BaM represents a probabilistic Bayesian approach and is based on the Monte Carlo Markov Chain algorithm (MCMC) to explore the posterior distribution of inferred parameters. For a better description of this algorithm, please see Renard et al. [2006]. 
https://github.com/BaM-tools/BaM


#############
IMPORTANT !!!
#############                                     

Please, notice that you will need to include the following libraries .dll (freely available in the web) 
in this folder "BaM_exe": 

- libgcc_s_dw2-1.dll
- libgcc_s_sjlj-1.dll
- libgfortran-3.dll
- libquadmath-0.dll
- libwinpthread-1.dll 


In some cases, you may also need to upload the following libraries:
libatomic-1.dll 
libcharset-1.dll
libexpat-1.dll
libgettextlib-0-18-3.dll
libgettextpo-0.dll 
libgettextsrc-0-18-3.dll
libgmp-10.dll
libgomp-1.dll
libiconv-2.dll
libintl-8.dll
libisl-15.dll
libltdl-7.dll
libmingwex-0.dll
libmpc-3.dll
libmpfr-4.dll
libssp-0.dll
libstdc++-6.dll
