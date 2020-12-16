# ----------------------------------------------------------------
# Code script presented in the publication:
#
# "Impact of Equatorial Atlantic Variability on ENSO Predictive Skill"
# by Eleftheria Exarchou, Pablo Ortega, Belen Rodrıguez de Fonseca,
# Teresa Losada, Irene Polo, and Chloe Prodhomme 
# in Nature Communications
#
# Data from EC-Earth predictions are available on request from 
# Eleftheria Exarchou. 
#
# HadISSTv1.1 and ERSSTv4 are publicly available
# here:https://www.metoffice.gov.uk/hadobs/hadisst/data/ and 
# https://www1.ncdc.noaa.gov/pub/data/.
#
# The data have been formatted to follow a particular data
# structure that allows the function "Load" by the s2dverification
# package to read the data. All data have been loaded with "Load" 
# and saved as an R object.      
#
# Additional information:
# Period: 1981 - 2018
# Ensemble members: 15 
#
# Date: 15 Dec 2020
# Author: Eleftheria Exarchou (eleftheria.exarchou@gmail.com)
# ----------------------------------------------------------------

rm(list=ls())
gc()

library(s2dverification)
library(boot)
library(ncdf4)
library(RColorBrewer)
library(abind)
library(multiApply)
library(SpecsVerification)

#-- Load DATA 
load("EC-Earth-indices.Rdata")

#-- Significance for correlation 
corr.significance  <- function(exp1,exp2,obs) {
	# Seasonal means 
	tmp.exp1 = Smoothing (exp1, runmeanlen = 3, numdimt = 4 )
	tmp.exp2 = Smoothing (exp2, runmeanlen = 3, numdimt = 4 )
	tmp.obs  = Smoothing (obs, runmeanlen = 3, numdimt = 4 )
	# Remove NA values
	tmp.exp1 [ is.na (tmp.exp1 )] = 0
	tmp.exp2 [ is.na (tmp.exp2 )] = 0
	tmp.obs  [ is.na (tmp.obs  )] = 0
	# Name dimensions 
	names(dim(tmp.exp1)) <- c ('dat1' , 'dat2', 'sdates' , 'ltimes')
	names(dim(tmp.exp2)) <- c ('dat1' , 'dat2', 'sdates' , 'ltimes')
	names(dim(tmp.obs )) <- c ('dat1' , 'dat2', 'sdates' , 'ltimes')
	# Call the function
	tmp   <- Apply(list(tmp.exp1,
						tmp.exp2,
						tmp.obs), 
				 list('sdates', 'sdates', 'sdates'), 
				 CorrDiff,conf.level = 0.95,handle.na='only.complete.triplets')$output1
    return(tmp)
} 

#-- Calculate anomalies for a period of 1981:cutoff_date
anomalies.obs <- function(data,cutoff_sdate) {
    tmp1 = InsertDim(InsertDim (data$obs[1,1,1:cutoff_sdate,],1,1),1,1)
	tmp =  Ano_CrossValid(tmp1 , tmp1)$ano_obs
  return(tmp)
}
anomalies.exp <- function(data,cutoff_sdate) {
    tmp1 = InsertDim( data$mod[1,,1:cutoff_sdate,],1,1)
	tmp =  InsertDim (Mean1Dim (Ano_CrossValid(tmp1 , tmp1, memb=F)$ano_exp ,2),1,1)
  return(tmp)
}

 atl3.obs.ano      <- anomalies.obs (atl3.obs , 1,38 ) 
 nino3.obs.ano     <- anomalies.obs (nino3.obs , 1,38 ) 
 nino34.obs.ano    <- anomalies.obs (nino34.obs , 1,38 ) 
 
 atl3.ctr.ano      <- anomalies.exp(atl3.ctr , 1,38 ) 
 nino3.ctr.ano     <- anomalies.exp(nino3.ctr , 1,38 ) 
 nino34.ctr.ano    <- anomalies.exp(nino34.ctr , 1,38 ) 
 
 atl3.nud1.ano      <- anomalies.exp(atl3.nud1 , 1,38 ) 
 nino3.nud1.ano     <- anomalies.exp(nino3.nud1 , 1,38 ) 
 nino34.nud1.ano    <- anomalies.exp(nino34.nud1 , 1,38 ) 
 
 atl3.nud2.ano      <- anomalies.exp(atl3.nud2 , 1,38 ) 
 nino3.nud2.ano     <- anomalies.exp(nino3.nud2 , 1,38 ) 
 nino34.nud2.ano    <- anomalies.exp(nino34.nud2 , 1,38 ) 
 

#-- Prediction Skill (Anonaly Correlation Coefficient)
Corr  <-  s2dverification::Corr

atl3.ctr.skill.seas <- Corr( Smoothing (atl3.ctr.ano, runmeanlen = 3, numdimt = 4), 	                        
                              Smoothing (atl3.obs.ano,runmeanlen = 3,numdimt = 4),
                              posloop=1,poscor=3)

atl3.nud1.skill.seas <- Corr( Smoothing (atl3.nud1.ano, runmeanlen = 3, numdimt = 4), 	                        
                              Smoothing (atl3.obs.ano,runmeanlen = 3,numdimt = 4),
                              posloop=1,poscor=3)

atl3.nud2.skill.seas <- Corr( Smoothing (atl3.nud2.ano, runmeanlen = 3, numdimt = 4), 	                        
                              Smoothing (atl3.obs.ano,runmeanlen = 3,numdimt = 4),
                              posloop=1,poscor=3)


nino3.ctr.skill.seas <- Corr( Smoothing (nino3.ctr.ano, runmeanlen = 3, numdimt = 4), 	                        
                              Smoothing (nino3.obs.ano,runmeanlen = 3,numdimt = 4),
                              posloop=1,poscor=3)

nino3.nud1.skill.seas <- Corr( Smoothing (nino3.nud1.ano, runmeanlen = 3, numdimt = 4), 	                        
                              Smoothing (nino3.obs.ano,runmeanlen = 3,numdimt = 4),
                              posloop=1,poscor=3)

nino3.nud2.skill.seas <- Corr( Smoothing (nino3.nud2.ano, runmeanlen = 3, numdimt = 4), 	                        
                              Smoothing (nino3.obs.ano,runmeanlen = 3,numdimt = 4),
                              posloop=1,poscor=3)



nino34.ctr.skill.seas <- Corr( Smoothing (nino34.ctr.ano, runmeanlen = 3, numdimt = 4), 	                        
                              Smoothing (nino34.obs.ano,runmeanlen = 3,numdimt = 4),
                              posloop=1,poscor=3)

nino34.nud1.skill.seas <- Corr( Smoothing (nino34.nud1.ano, runmeanlen = 3, numdimt = 4), 	                        
                              Smoothing (nino34.obs.ano,runmeanlen = 3,numdimt = 4),
                              posloop=1,poscor=3)

nino34.nud2.skill.seas <- Corr( Smoothing (nino34.nud2.ano, runmeanlen = 3, numdimt = 4), 	                        
                              Smoothing (nino34.obs.ano,runmeanlen = 3,numdimt = 4),
                              posloop=1,poscor=3)


# -- Correlations differences between ctr and SST-nudged
#   forecasts, assess the uncertainty using 
#   the function CorrDiff of the SpecsVerification R-package 
atl3.corr.sign.ctr.nud1    <- corr.significance(atl3.ctr.ano,  atl3.nud1.ano,    atl3.obs.ano)
nino3.corr.sign.ctr.nud1   <- corr.significance(nino3.ctr.ano, nino3.nud1.ano,   nino3.obs.ano)
nino34.corr.sign.ctr.nud1  <- corr.significance(nino34.ctr.ano,nino34.nud1.ano,  nino34.obs.ano)
atl3.corr.sign.ctr.nud2    <- corr.significance(atl3.ctr.ano,  atl3.nud2.ano,    atl3.obs.ano)
nino3.corr.sign.ctr.nud2   <- corr.significance(nino3.ctr.ano, nino3.nud2.ano,   nino3.obs.ano)
nino34.corr.sign.ctr.nud2  <- corr.significance(nino34.ctr.ano,nino34.nud2.ano,  nino34.obs.ano)

# -- RMSS skill 
atl3.ctr.rmsss   <- RMSSS(Smoothing(atl3.ctr.ano, runmeanlen = 3, numdimt = 4), 
                          Smoothing(atl3.obs.ano, runmeanlen = 3, numdimt = 4) ,
                          posloop = 1, posRMS = 3, pval = TRUE)
nino3.ctr.rmsss   <- RMSSS(Smoothing(nino3.ctr.ano, runmeanlen = 3, numdimt = 4), 
                          Smoothing(nino3.obs.ano, runmeanlen = 3, numdimt = 4) ,
                          posloop = 1, posRMS = 3, pval = TRUE)
nino34.ctr.rmsss   <- RMSSS(Smoothing(nino34.ctr.ano, runmeanlen = 3, numdimt = 4), 
                          Smoothing(nino34.obs.ano, runmeanlen = 3, numdimt = 4) ,
                          posloop = 1, posRMS = 3, pval = TRUE)
atl3.nud1.rmsss   <- RMSSS(Smoothing(atl3.nud1.ano, runmeanlen = 3, numdimt = 4), 
                          Smoothing(atl3.obs.ano, runmeanlen = 3, numdimt = 4) ,
                          posloop = 1, posRMS = 3, pval = TRUE)
nino3.nud1.rmsss   <- RMSSS(Smoothing(nino3.nud1.ano, runmeanlen = 3, numdimt = 4), 
                          Smoothing(nino3.obs.ano, runmeanlen = 3, numdimt = 4) ,
                          posloop = 1, posRMS = 3, pval = TRUE)
nino34.nud1.rmsss   <- RMSSS(Smoothing(nino34.nud1.ano, runmeanlen = 3, numdimt = 4), 
                          Smoothing(nino34.obs.ano, runmeanlen = 3, numdimt = 4) ,
                          posloop = 1, posRMS = 3, pval = TRUE)
atl3.nud2.rmsss   <- RMSSS(Smoothing(atl3.nud2.ano, runmeanlen = 3, numdimt = 4), 
                          Smoothing(atl3.obs.ano, runmeanlen = 3, numdimt = 4) ,
                          posloop = 1, posRMS = 3, pval = TRUE)
nino3.nud2.rmsss   <- RMSSS(Smoothing(nino3.nud2.ano, runmeanlen = 3, numdimt = 4), 
                          Smoothing(nino3.obs.ano, runmeanlen = 3, numdimt = 4) ,
                          posloop = 1, posRMS = 3, pval = TRUE)
nino34.nud2.rmsss   <- RMSSS(Smoothing(nino34.nud2.ano, runmeanlen = 3, numdimt = 4), 
                          Smoothing(nino34.obs.ano, runmeanlen = 3, numdimt = 4) ,
                          posloop = 1, posRMS = 3, pval = TRUE)
# -- Plot ACC skill 
atl3.exp1   = atl3.nud1.skill.seas   
atl3.exp2   = atl3.ctr.skill.seas   
atl3.exp3   = atl3.nud2.skill.seas   
nino3.exp1  = nino3.nud1.skill.seas   
nino3.exp2  = nino3.cre.skill.seas   
nino3.exp3  = nino3.nud2.skill.seas   
nino34.exp1 = nino34.nud1.skill.seas   
nino34.exp2 = nino34.ctr.skill.seas   
nino34.exp3 = nino34.nud2.skill.seas   

  atl3.seas        <-  array(NA, dim=c(3,6)) 
  nino3.seas       <-  array(NA, dim=c(3,6)) 
  nino34.seas      <-  array(NA, dim=c(3,6)) 
  atl3.seas.sign   <-  array(NA, dim=c(3,6)) 
  nino3.seas.sign  <-  array(NA, dim=c(3,6)) 
  nino34.seas.sign <-  array(NA, dim=c(3,6)) 

atl3.seas[1,]    <-  atl3.exp1   [ 1,1,2,1,2:7]
atl3.seas[2,]    <-  atl3.exp2   [ 1,1,2,1,2:7]
atl3.seas[3,]    <-  atl3.exp3   [ 1,1,2,1,2:7]
nino3.seas[1,]   <-  nino3.exp1  [ 1,1,2,1,2:7]
nino3.seas[2,]   <-  nino3.exp2  [ 1,1,2,1,2:7]
nino3.seas[3,]   <-  nino3.exp3  [ 1,1,2,1,2:7]
nino34.seas[1,]  <-  nino34.exp1 [ 1,1,2,1,2:7]
nino34.seas[2,]  <-  nino34.exp2 [ 1,1,2,1,2:7]
nino34.seas[3,]  <-  nino34.exp3 [ 1,1,2,1,2:7]

atl3.seas.sign[1,]      <-  atl3.exp1  [ 1,1,1,1,2:7]   * atl3.exp1  [ 1,1,3,1,2:7]
atl3.seas.sign[2,]      <-  atl3.exp2  [ 1,1,1,1,2:7]   * atl3.exp2  [ 1,1,3,1,2:7]
atl3.seas.sign[3,]      <-  atl3.exp3  [ 1,1,1,1,2:7]   * atl3.exp3  [ 1,1,3,1,2:7]
nino3.seas.sign[1,]     <-  nino3.exp1 [ 1,1,1,1,2:7]   * nino3.exp1  [ 1,1,3,1,2:7]
nino3.seas.sign[2,]     <-  nino3.exp2 [ 1,1,1,1,2:7]   * nino3.exp2  [ 1,1,3,1,2:7]
nino3.seas.sign[3,]     <-  nino3.exp3 [ 1,1,1,1,2:7]   * nino3.exp3 [ 1,1,3,1,2:7]
nino34.seas.sign[1,]    <-  nino34.exp1 [ 1,1,1,1,2:7]  * nino34.exp1 [ 1,1,3,1,2:7]
nino34.seas.sign[2,]    <-  nino34.exp2 [ 1,1,1,1,2:7]  * nino34.exp2 [ 1,1,3,1,2:7]
nino34.seas.sign[3,]    <-  nino34.exp3 [ 1,1,1,1,2:7]  * nino34.exp3 [ 1,1,3,1,2:7]

legend1="NUD-var"
legend2="CTR"
legend3="NUD-JJAS"

cols<- brewer.pal(4,"Dark2")
postscript( "Figure3-ACC.ps" ,width=7,height=5)
main.title1="Prediction skill in SST in ATL3" 
main.title2="Prediction skill in SST in Nino3.4"  
main.title3="Prediction skill in SST in Nino3"  
sign.threshold = 0.05

# -------- 1st page 
time=seq(1,6)
# First page 
plot(time,atl3.seas[1,],type="l",ylab="Correlation
	 coefficient",xlab="",ylim=c(0.,1),lwd=2,col=cols[1],main=main.title1,font.main=1,xaxt='n')
points(time [atl3.seas.sign[1,]   >0  ] ,   atl3.seas[1,][atl3.seas.sign[1,]   >0 ] ,lwd=2,col=cols[1])

lines(time,atl3.seas[2,],lwd=2,lty=1,col=cols[2])
points(time [atl3.seas.sign[2,]   >0  ] ,   atl3.seas[2,][atl3.seas.sign[2,]   >0 ] ,lwd=2,col=cols[2])

lines(time,atl3.seas[3,],lwd=2,lty=1,col=cols[3])
points(time [atl3.seas.sign[3,]   >0  ] ,   atl3.seas[3,][atl3.seas.sign[3,]   >0 ] ,lwd=2,col=cols[3])

### Here if p < 0.05  (2nd value) , correlation differences are significant,
points(time [atl3.seas.sign [1,] >0 & atl3.corr.sign.ctr.nud1 [2,1,1,2:7]   < sign.threshold   ],
	   atl3.seas[1,][atl3.seas.sign[1,] >0 & atl3.corr.sign.ctr.nud1 [2,1,1,2:7] < sign.threshold ],pch=19, cex=1.5, col=cols[1])

points(time [atl3.seas.sign [3,] >0 & atl3.corr.sign.ctr.nud2 [2,1,1,2:7]   < sign.threshold   ],
	   atl3.seas[3,][atl3.seas.sign[3,] >0 & atl3.corr.sign.ctr.nud2 [2,1,1,2:7] < sign.threshold ],pch=19, cex=1.5, col=cols[3])

legend(4,0.99,c(legend1,legend2,legend3 )
      ,lty=c(1),cex=c(1.)
      ,lwd=c(2),col=c(cols[1], cols[2], cols[3]  )
      ,seg.len=c(2.5),  bty = "n")

axis(1, at=time,labels=c('JJA','JAS','ASO', 'SON', 'OND', 'NDJ'))

# ------ Second page ;;  nino 3.4   
plot(time,nino34.seas[1,],type="l",ylab="Correlation coefficient",xlab="",ylim=c(0.75,1),lwd=2,col=cols[1],main=main.title2,font.main=1,xaxt='n')
points(time [nino34.seas.sign[1,]   >0  ] ,   nino34.seas[1,][nino34.seas.sign[1,]   >0 ] ,lwd=2,col=cols[1])

lines(time,nino34.seas[2,],lwd=2,lty=1,col=cols[2])
points(time [nino34.seas.sign[2,]   >0  ] ,   nino34.seas[2,][nino34.seas.sign[2,]   >0 ] ,lwd=2,col=cols[2])

lines(time,nino34.seas[3,],lwd=2,lty=1,col=cols[3])
points(time [nino34.seas.sign[3,]   >0  ] ,  nino34.seas[3,][nino34.seas.sign[3,]   >0 ] ,lwd=2,col=cols[3])

### Here if p < 0.05  (2nd value) I consider it significant
points(time [nino34.seas.sign [1,] >0 & nino34.corr.sign.ctr.nud1 [2,1,1,2:7]   < sign.threshold   ],
	   nino34.seas[1,][nino34.seas.sign[1,] >0 & nino34.corr.sign.ctr.nud1 [2,1,1,2:7] < sign.threshold ],pch=19, cex=1.5, col=cols[1])

points(time [nino34.seas.sign [3,] >0 & nino34.corr.sign.ctr.nud2 [2,1,1,2:7]   < sign.threshold   ],
	   nino34.seas[3,][nino34.seas.sign[3,] >0 & nino34.corr.sign.ctr.nud2 [2,1,1,2:7] < sign.threshold ],pch=19, cex=1.5, col=cols[3])


#############

legend(3.,0.99,c(legend1, legend2, legend3)
      ,lty=c(1),cex=c(1.)
      ,lwd=c(2),col=c(cols[1], cols[2], cols[3]  )
      ,seg.len=c(2.5),  bty = "n")

axis(1, at=time,labels=c('JJA','JAS','ASO', 'SON', 'OND', 'NDJ'))


# -------- Third page ;; nino 3  
plot(time,nino3.seas[1,],type="l",ylab="Correlation coefficient",xlab="",ylim=c(0.75,1),lwd=2,col=cols[1],main=main.title3,font.main=1,xaxt='n')
points(time [nino3.seas.sign[1,]   >0  ] ,   nino3.seas[1,][nino3.seas.sign[1,]   >0 ] ,lwd=2,col=cols[1])

lines(time,nino3.seas[2,],lwd=2,lty=1,col=cols[2])
points(time [nino3.seas.sign[2,]   >0  ] ,   nino3.seas[2,][nino3.seas.sign[2,]   >0 ] ,lwd=2,col=cols[2])

lines(time,nino3.seas[3,],lwd=2,lty=1,col=cols[3])
points(time [nino3.seas.sign[3,]   >0  ] ,  nino3.seas[3,][nino3.seas.sign[3,]   >0 ] ,lwd=2,col=cols[3])

### Here if p < 0.05  (2nd value) I consider it significant,
points(time [nino3.seas.sign [1,] >0 & nino3.corr.sign.ctr.nud1 [2,1,1,2:7]   < sign.threshold   ],
	   nino3.seas[1,][nino3.seas.sign[1,] >0 & nino3.corr.sign.ctr.nud1 [2,1,1,2:7] < sign.threshold ],pch=19, cex=1.5, col=cols[1])


points(time [nino3.seas.sign [3,] >0 & nino3.corr.sign.ctr.nud2 [2,1,1,2:7]   < sign.threshold   ],
	   nino3.seas[3,][nino3.seas.sign[3,] >0 & nino3.corr.sign.ctr.nud2 [2,1,1,2:7] < sign.threshold ],pch=19, cex=1.5, col=cols[3])


legend(3.,0.99,c(legend1, legend2, legend3)
      ,lty=c(1),cex=c(1.)
      ,lwd=c(2),col=c(cols[1], cols[2] , cols[3]  )
      ,seg.len=c(2.5),  bty = "n")

axis(1, at=time,labels=c('JJA','JAS','ASO', 'SON', 'OND', 'NDJ'))


dev.off()


# ------- Plot RMSSS  
atl3.exp1   =  atl3.nud1.rmsss    
atl3.exp2   =  atl3.ctr.rmsss 
atl3.exp3   =  atl3.nud2.rmsss 
nino3.exp1  =  nino3.nud1.rmsss
nino3.exp2  =  nino3.ctr.rmsss
nino3.exp3  =  nino3.nud2.rmsss
nino34.exp1 =  nino34.nud1.rmsss 
nino34.exp2 =  nino34.ctr.rmsss 
nino34.exp3 =  nino34.nud2.rmsss 

  atl3.seas        <-  array(NA, dim=c(3,6)) 
  nino3.seas       <-  array(NA, dim=c(3,6)) 
  nino34.seas      <-  array(NA, dim=c(3,6)) 
  atl3.seas.sign   <-  array(NA, dim=c(3,6)) 
  nino3.seas.sign  <-  array(NA, dim=c(3,6)) 
  nino34.seas.sign <-  array(NA, dim=c(3,6)) 

atl3.seas[1,]    <-  atl3.exp1   [ 1,1,1,1,2:7]
atl3.seas[2,]    <-  atl3.exp2   [ 1,1,1,1,2:7]
atl3.seas[3,]    <-  atl3.exp3   [ 1,1,1,1,2:7]
nino3.seas[1,]   <-  nino3.exp1  [ 1,1,1,1,2:7]
nino3.seas[2,]   <-  nino3.exp2  [ 1,1,1,1,2:7]
nino3.seas[3,]   <-  nino3.exp3  [ 1,1,1,1,2:7]
nino34.seas[1,]  <-  nino34.exp1 [ 1,1,1,1,2:7]
nino34.seas[2,]  <-  nino34.exp2 [ 1,1,1,1,2:7]
nino34.seas[3,]  <-  nino34.exp3 [ 1,1,1,1,2:7]

atl3.seas.sign[1,]      <-   atl3.exp1   [ 1,1,2,1,2:7]
atl3.seas.sign[2,]      <-   atl3.exp2   [ 1,1,2,1,2:7]
atl3.seas.sign[3,]      <-   atl3.exp3   [ 1,1,2,1,2:7]
nino3.seas.sign[1,]     <-   nino3.exp1  [ 1,1,2,1,2:7]
nino3.seas.sign[2,]     <-   nino3.exp2  [ 1,1,2,1,2:7]
nino3.seas.sign[3,]     <-   nino3.exp3  [ 1,1,2,1,2:7]
nino34.seas.sign[1,]    <-   nino34.exp1 [ 1,1,2,1,2:7]
nino34.seas.sign[2,]    <-   nino34.exp2 [ 1,1,2,1,2:7]
nino34.seas.sign[3,]    <-   nino34.exp3 [ 1,1,2,1,2:7]

legend1="NUD-var"
legend2="CTR"
legend3="NUD-JJAS"

cols<- brewer.pal(4,"Dark2")
postscript( "Figure3-RMSSS.ps" ,width=7,height=5)
main.title1="RMSS skill in SST in ATL3" 
main.title2="RMSS skill in SST in Nino3.4"  
main.title3="RMSS skill in SST in Nino3"  

# -------- 1st page 
time=seq(1,6)
# First page 
plot(time,atl3.seas[1,],type="l",ylab="Correlation
	 coefficient",xlab="",ylim=c(-0.5,1),lwd=2,col=cols[1],main=main.title1,font.main=1,xaxt='n')
points(time [atl3.seas.sign[1,]   < 0.05  ] ,   atl3.seas[1,][atl3.seas.sign[1,]  < 0.05 ] ,lwd=2,col=cols[1])

lines(time,atl3.seas[2,],lwd=2,lty=1,col=cols[2])
points(time [atl3.seas.sign[2,]   <0.05  ] ,   atl3.seas[2,][atl3.seas.sign[2,]   <0.05 ] ,lwd=2,col=cols[2])

lines(time,atl3.seas[3,],lwd=2,lty=1,col=cols[3])
points(time [atl3.seas.sign[3,]   <0.05  ] ,   atl3.seas[3,][atl3.seas.sign[3,]   <0.05 ] ,lwd=2,col=cols[3])


legend(4,0.3,c(legend1,legend2,legend3 )
      ,lty=c(1),cex=c(1.)
      ,lwd=c(2),col=c(cols[1], cols[2], cols[3]  )
      ,seg.len=c(2.5),  bty = "n")

axis(1, at=time,labels=c('JJA','JAS','ASO', 'SON', 'OND', 'NDJ'))

# ------ Second page ;;  nino 3.4   
plot(time,nino34.seas[1,],type="l",ylab="Correlation coefficient",xlab="",ylim=c(0.2,0.5),lwd=2,col=cols[1],main=main.title2,font.main=1,xaxt='n')
points(time [nino34.seas.sign[1,]  <0.05  ] ,   nino34.seas[1,][nino34.seas.sign[1,]   <0.05] ,lwd=2,col=cols[1])

lines(time,nino34.seas[2,],lwd=2,lty=1,col=cols[2])
points(time [nino34.seas.sign[2,]   <0.05 ] ,   nino34.seas[2,][nino34.seas.sign[2,]   <0.05] ,lwd=2,col=cols[2])

lines(time,nino34.seas[3,],lwd=2,lty=1,col=cols[3])
points(time [nino34.seas.sign[3,]   <0.05 ] ,  nino34.seas[3,][nino34.seas.sign[3,]   <0.05] ,lwd=2,col=cols[3])

#############

legend(4.,0.3,c(legend1, legend2, legend3)
      ,lty=c(1),cex=c(1.)
      ,lwd=c(2),col=c(cols[1], cols[2], cols[3]  )
      ,seg.len=c(2.5),  bty = "n")

axis(1, at=time,labels=c('JJA','JAS','ASO', 'SON', 'OND', 'NDJ'))


# -------- Third page ;; nino 3  
plot(time,nino3.seas[1,],type="l",ylab="Correlation coefficient",xlab="",ylim=c(0.3,.5),lwd=2,col=cols[1],main=main.title3,font.main=1,xaxt='n')
points(time [nino3.seas.sign[1,]  <0.05 ] ,   nino3.seas[1,][nino3.seas.sign[1,]   <0.05] ,lwd=2,col=cols[1])

lines(time,nino3.seas[2,],lwd=2,lty=1,col=cols[2])
points(time [nino3.seas.sign[2,]  <0.05  ] ,   nino3.seas[2,][nino3.seas.sign[2,]  <0.05 ] ,lwd=2,col=cols[2])

lines(time,nino3.seas[3,],lwd=2,lty=1,col=cols[3])
points(time [nino3.seas.sign[3,]  <0.05  ] ,  nino3.seas[3,][nino3.seas.sign[3,]   <0.05] ,lwd=2,col=cols[3])

legend(3.,0.99,c(legend1, legend2, legend3)
      ,lty=c(1),cex=c(1.)
      ,lwd=c(2),col=c(cols[1], cols[2] , cols[3]  )
      ,seg.len=c(2.5),  bty = "n")

axis(1, at=time,labels=c('JJA','JAS','ASO', 'SON', 'OND', 'NDJ'))


dev.off()
 
