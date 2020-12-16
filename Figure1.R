# ----------------------------------------------------------------
# Code script presented in the publication:
#
# "Impact of Equatorial Atlantic Variability on ENSO Predictive Skill"
# by Eleftheria Exarchou, Pablo Ortega, Belen Rodrıguez de Fonseca,
# Teresa Losada, Irene Polo, and Chloe Prodhomme 
# in Nature Communications
#
# NMME data have been downloaded from the IRI database
# (http://iridl.ldeo.columbia.edu/SOURCES/.Models/.NMME/). EUROSIP and
# ERA-interim data have been downloaded from the MARS computing facility
# of ECMWF (https://www.ecmwf.int/en/forecasts/). 
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
# Period: 1981 - 2011
# Ensemble members: all available (information on Supplementary Table 1)
#
# Date: 15 Dec 2020
# Author: Eleftheria Exarchou (eleftheria.exarchou@gmail.com)
# ----------------------------------------------------------------

rm(list=ls())
gc()

library(s2dverification)
library(RColorBrewer)
library(ncdf4)

# --- Load data: atl3, nino3 and nino34 area-mean SST
load ("NMME-Eurosip-ECEarth.Rdata")

# --- Anomalies 

atl3.mod.ano    <-  Ano_CrossValid (atl3$mod,atl3$obs)$ano_exp
nino3.mod.ano   <-  Ano_CrossValid (nino3$mod,nino3$obs)$ano_exp
nino34.mod.ano  <-  Ano_CrossValid (nino34$mod,nino34$obs)$ano_exp

atl3.obs.ano    <-  Ano_CrossValid (atl3$obs,atl3$obs)$ano_obs
nino3.obs.ano   <-  Ano_CrossValid (nino3$obs,nino3$obs)$ano_obs
nino34.obs.ano  <-  Ano_CrossValid (nino34$obs,nino34$obs)$ano_obs

# --- Correlations  
atl3.corr   <- Corr(atl3.mod.ano,   atl3.obs.ano    , posloop = 1, poscor = 3)
nino3.corr  <- Corr(nino3.mod.ano,  nino3.obs.ano   , posloop = 1, poscor = 3)
nino34.corr <- Corr(nino34.mod.ano, nino34.obs.ano  , posloop = 1, poscor = 3)

# --- PLOT 

exps      <- c ( 
                 'cfs_v2'
              ,  'ecmwf/system4' 
              ,  'cancm4'
              ,  'cancm3'
              ,  'cm2p5-flor-a06'
              ,  'cm2p5-flor-b01'
              ,  'rsmas-ccsm4'
              ,  'rsmas-ccsm3'
              ,  'cm2p1'
              ,  'cm2p1-aer04'  
              ,  'echam4p5'   
              ,  'gmao-062012'  
              ,  'ecmwf/system5' 
              ,  'ncep'  
              ,  'EC-Earth' 
 )

time=seq(1,10)

color1=c("blue4", "red", "yellow", "magenta",  "green")
pal <- colorRampPalette(color1)
cols=pal(16)
 
postscript(paste("Figure1.ps",sep=""),width=7,height=5)

# First page Nino3 
plot(time,nino3.corr[15,1,2,1,],type="l",ylab="Anomaly correlation coefficient", 
	 xlab="",ylim=c(0.,1),lwd=2,col=cols[15], 
	 main="Prediction skill Nino3",font.main=1.5,xaxt='n',
	 cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)

points(time [nino3.corr[15,1,1,,1] * nino3.corr[15,1,3,1,] >0  ] , nino3.corr[15,1,2,1,] [
	   nino3.corr[15,1,1,1,] * nino3.corr[15,1,3,1,]  >0  ] ,lwd=2,col=cols[15])

for ( ens in 1:14 )
{ 
lines(time, nino3.corr[ens,1,2,1,],lwd=2,lty=1,col=cols[ens])
points(time [nino3.corr[ens,1,1,,1] * nino3.corr[ens,1,3,1,] >0  ] , nino3.corr[ens,1,2,1,] [
	   nino3.corr[ens,1,1,1,] * nino3.corr[ens,1,3,1,]  >0  ] ,lwd=2,col=cols[ens])
} 

# Persistence 
lines(time, corr.per.nino3  ,lwd=4,lty=1,col='black')

legend(0.7,0.33,c(exps[1:5] )
      ,lty=c(1),cex=c(0.8)
      ,lwd=c(2),col=c(cols[1:5])
      ,seg.len=c(2),  bty = "n")

legend(3.7,0.33,c(exps[6:10] )
      ,lty=c(1),cex=c(0.8)
      ,lwd=c(2),col=c(cols[6:10])
      ,seg.len=c(2),  bty = "n")

legend(7.7,0.33,c(exps[11:15] , 'persistence')
      ,lty=c(1),cex=c(0.8)
      ,lwd=c(2),col=c(cols[11:15],'black')
      ,seg.len=c(2),  bty = "n")

axis(1, at=time,cex.axis=1.3, labels=c('Jun','Jul','Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb', 'Mar'  ))

# Second page Nino3.4 
plot(time,nino34.corr[15,1,2,1,],type="l",ylab="Anomaly correlation coefficient", 
	 xlab="",ylim=c(0.,1),lwd=2,col=cols[15], 
	 main="Prediction skill Nino3.4",font.main=1.5,xaxt='n',
	 cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)

points(time [nino34.corr[15,1,1,,1] * nino34.corr[15,1,3,1,] >0  ] , nino34.corr[15,1,2,1,] [
	   nino34.corr[15,1,1,1,] * nino34.corr[15,1,3,1,]  >0  ] ,lwd=2,col=cols[15])

for ( ens in 1:14 )
{ 
lines(time, nino34.corr[ens,1,2,1,],lwd=2,lty=1,col=cols[ens])
points(time [nino34.corr[ens,1,1,,1] * nino34.corr[ens,1,3,1,] >0  ] , nino34.corr[ens,1,2,1,] [
	   nino34.corr[ens,1,1,1,] * nino34.corr[ens,1,3,1,]  >0  ] ,lwd=2,col=cols[ens])
} 

# Persistence 
lines(time, corr.per.nino34  ,lwd=4,lty=1,col='black')

legend(0.7,0.33,c(exps[1:5] )
      ,lty=c(1),cex=c(0.8)
      ,lwd=c(2),col=c(cols[1:5])
      ,seg.len=c(2),  bty = "n")

legend(3.7,0.33,c(exps[6:10] )
      ,lty=c(1),cex=c(0.8)
      ,lwd=c(2),col=c(cols[6:10])
      ,seg.len=c(2),  bty = "n")

legend(7.7,0.33,c(exps[11:15] , 'persistence')
      ,lty=c(1),cex=c(0.8)
      ,lwd=c(2),col=c(cols[11:15],'black')
      ,seg.len=c(2),  bty = "n")

axis(1, at=time,cex.axis=1.3, labels=c('Jun','Jul','Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb', 'Mar'  ))

# Third page Atl3
plot(time,atl3.corr[15,1,2,1,],type="l",ylab="Anomaly correlation coefficient", 
	 xlab="",ylim=c(0.,1),lwd=2,col=cols[15], 
	 main="Prediction skill ATL3",font.main=1.5,xaxt='n',
	 cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)

points(time [atl3.corr[15,1,1,,1] * atl3.corr[15,1,3,1,] >0  ] , atl3.corr[15,1,2,1,] [
	   atl3.corr[15,1,1,1,] * atl3.corr[15,1,3,1,]  >0  ] ,lwd=2,col=cols[15])

for ( ens in 1:14 )
{ 
lines(time, atl3.corr[ens,1,2,1,],lwd=2,lty=1,col=cols[ens])
points(time [atl3.corr[ens,1,1,,1] * atl3.corr[ens,1,3,1,] >0  ] , atl3.corr[ens,1,2,1,] [
	   atl3.corr[ens,1,1,1,] * atl3.corr[ens,1,3,1,]  >0  ] ,lwd=2,col=cols[ens])
} 

# Persistence 
lines(time, corr.per.atl3  ,lwd=4,lty=1,col='black')
axis(1, at=time,cex.axis=1.3, labels=c('Jun','Jul','Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb', 'Mar'  ))
dev.off()

