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
# The script below calculates the correlations between the JJA ATL3 index
# and the SST/uas/vas fields and saves the results in netCDF files
#
# Date: 15 Dec 2020
# Author: Eleftheria Exarchou (eleftheria.exarchou@gmail.com)
# ----------------------------------------------------------------

rm(list=ls())
gc()

library(s2dverification)
library(boot)
library(RColorBrewer)
library(abind)
library(ncdf4)
library(startR)
library(multiApply)
library(SpecsVerification)
library(multiApply)
library(plyr)
library(parallel)
library(doParallel)
library(foreach)
source ("Apply.R")
source ("zzz.R")
#--- Load data 
load ("Spatial-data.Rdata")

#--- Concatenate the ensemble members onto a single array 
dim1  <-  dim(tos.ctr.ano)[1] 
dim2  <-  dim(tos.ctr.ano)[2] 
dim3  <-  dim(tos.ctr.ano)[3] 
dim4  <-  dim(tos.ctr.ano)[4] 
dim5  <-  dim(tos.ctr.ano)[5] 
dim6  <-  dim(tos.ctr.ano)[6] 

concat  <- function(exp) {

atl3.concat = array(NA, dim=c( 1, 1,  dim2*dim3, dim4  ))
tos.concat = array(NA, dim=c( 1, 1,  dim2*dim3, dim4, dim5,  dim6   ))
uas.concat = array(NA, dim=c( 1, 1,  dim2*dim3, dim4, dim5,  dim6   ))
vas.concat = array(NA, dim=c( 1, 1,  dim2*dim3, dim4, dim5,  dim6   ))

atl3.ano= switch (exp,tos.ctr.atl3.ano,tos.nud1.atl3.ano,tos.nud2.atl3.ano)
tos.ano = switch (exp,tos.ctr.ano,     tos.nud1.ano,     tos.nud2.ano    )
uas.ano = switch (exp,uas.ctr.ano,     uas.nud1.ano,     uas.nud2.ano    )
vas.ano = switch (exp,vas.ctr.ano,     vas.nud1.ano,     vas.nud2.ano    )

for ( i  in 1:dim4   )
{             
atl3.concat [ 1, 1, , i ]  =   c (atl3.ano[1,,, i ])
} 


for ( i  in 1:dim4  )
{             
for ( j  in 1:dim5  )
{             
for ( k  in 1:dim6  )
{ 
tos.concat[ 1, 1, , i, j, k  ]  =   c (tos.ano[1, ,, i, j,k ])
uas.concat[ 1, 1, , i, j, k  ]  =   c (uas.ano[1, ,, i, j,k ])
vas.concat[ 1, 1, , i, j, k  ]  =   c (vas.ano[1, ,, i, j,k ])
} 
} 
} 

#--- Seasonal means 
#--- tos 
 atl3.concat.jja   <- Season(  atl3.concat, posdim = 4, 6, 6, 8)
 tos.concat.jja    <- Season(  tos.concat, posdim = 4, 6, 6, 8)
 tos.concat.jas    <- Season(  tos.concat, posdim = 4, 6, 7, 9)
 tos.concat.aso    <- Season(  tos.concat, posdim = 4, 6, 8, 10)
 tos.concat.son    <- Season(  tos.concat, posdim = 4, 6, 9, 11)
 tos.concat.ond    <- Season(  tos.concat, posdim = 4, 6, 10, 12)
 tos.concat.ndj    <- Season(  tos.concat, posdim = 4, 6, 11, 1)

 #--- uas 
 uas.concat.jja    <- Season(  uas.concat, posdim = 4, 6, 6, 8)
 uas.concat.jas    <- Season(  uas.concat, posdim = 4, 6, 7, 9)
 uas.concat.aso    <- Season(  uas.concat, posdim = 4, 6, 8, 10)
 uas.concat.son    <- Season(  uas.concat, posdim = 4, 6, 9, 11)
 uas.concat.ond    <- Season(  uas.concat, posdim = 4, 6, 10, 12)
 uas.concat.ndj    <- Season(  uas.concat, posdim = 4, 6, 11, 1)
 #--- uas 
 vas.concat.jja    <- Season(  vas.concat, posdim = 4, 6, 6, 8)
 vas.concat.jas    <- Season(  vas.concat, posdim = 4, 6, 7, 9)
 vas.concat.aso    <- Season(  vas.concat, posdim = 4, 6, 8, 10)
 vas.concat.son    <- Season(  vas.concat, posdim = 4, 6, 9, 11)
 vas.concat.ond    <- Season(  vas.concat, posdim = 4, 6, 10, 12)
 vas.concat.ndj    <- Season(  vas.concat, posdim = 4, 6, 11, 1)

    results <- list()
    results$atl3     <- atl3.concat.jja
    results$tos.jja  <- tos.concat.jja
    results$tos.jas  <- tos.concat.jas
    results$tos.aso  <- tos.concat.aso
    results$tos.son  <- tos.concat.son
    results$tos.ond  <- tos.concat.ond
    results$tos.ndj  <- tos.concat.ndj

    results$uas.jja  <- uas.concat.jja
    results$uas.jas  <- uas.concat.jas
    results$uas.aso  <- uas.concat.aso
    results$uas.son  <- uas.concat.son
    results$uas.ond  <- uas.concat.ond
    results$uas.ndj  <- uas.concat.ndj
    
    results$vas.jja  <- vas.concat.jja
    results$vas.jas  <- vas.concat.jas
    results$vas.aso  <- vas.concat.aso
    results$vas.son  <- vas.concat.son
    results$vas.ond  <- vas.concat.ond
    results$vas.ndj  <- vas.concat.ndj

    return (results)
} 

tmp=concat(1)
tos.ctr.atl3.ano.concat.jja = tmp$atl3
tos.ctr.ano.concat.jja = tmp$tos.jja  
tos.ctr.ano.concat.jas = tmp$tos.jas  
tos.ctr.ano.concat.aso = tmp$tos.aso  
tos.ctr.ano.concat.son = tmp$tos.son  

uas.ctr.ano.concat.jja = tmp$uas.jja  
uas.ctr.ano.concat.jas = tmp$uas.jas  
uas.ctr.ano.concat.aso = tmp$uas.aso  
uas.ctr.ano.concat.son = tmp$uas.son  

vas.ctr.ano.concat.jja = tmp$vas.jja  
vas.ctr.ano.concat.jas = tmp$vas.jas  
vas.ctr.ano.concat.aso = tmp$vas.aso  
vas.ctr.ano.concat.son = tmp$vas.son  

tmp=concat(2)
tos.nud1.atl3.ano.concat.jja = tmp$atl3
tos.nud1.ano.concat.jja = tmp$tos.jja  
tos.nud1.ano.concat.jas = tmp$tos.jas  
tos.nud1.ano.concat.aso = tmp$tos.aso  
tos.nud1.ano.concat.son = tmp$tos.son  

uas.nud1.ano.concat.jja = tmp$uas.jja  
uas.nud1.ano.concat.jas = tmp$uas.jas  
uas.nud1.ano.concat.aso = tmp$uas.aso  
uas.nud1.ano.concat.son = tmp$uas.son  

vas.nud1.ano.concat.jja = tmp$vas.jja  
vas.nud1.ano.concat.jas = tmp$vas.jas  
vas.nud1.ano.concat.aso = tmp$vas.aso  
vas.nud1.ano.concat.son = tmp$vas.son  

#--- Calculate Correlation between ATL3 JJA and the spatial SST/uas/vas
#    fields 
tos.cor.ctr.jja  <-  array (0, dim = c(4, dim5,dim6 ) )
tos.cor.ctr.jas  <-  array (0, dim = c(4, dim5,dim6 ) )
tos.cor.ctr.aso  <-  array (0, dim = c(4, dim5,dim6 ) )
tos.cor.ctr.son  <-  array (0, dim = c(4, dim5,dim6 ) )
tos.cor.ctr.ond  <-  array (0, dim = c(4, dim5,dim6 ) )
tos.cor.ctr.ndj  <-  array (0, dim = c(4, dim5,dim6 ) )

uas.cor.ctr.jja  <-  array (0, dim = c(4, dim5,dim6 ) )
uas.cor.ctr.jas  <-  array (0, dim = c(4, dim5,dim6 ) )
uas.cor.ctr.aso  <-  array (0, dim = c(4, dim5,dim6 ) )
uas.cor.ctr.son  <-  array (0, dim = c(4, dim5,dim6 ) )
uas.cor.ctr.ond  <-  array (0, dim = c(4, dim5,dim6 ) )
uas.cor.ctr.ndj  <-  array (0, dim = c(4, dim5,dim6 ) )

vas.cor.ctr.jja  <-  array (0, dim = c(4, dim5,dim6 ) )
vas.cor.ctr.jas  <-  array (0, dim = c(4, dim5,dim6 ) )
vas.cor.ctr.aso  <-  array (0, dim = c(4, dim5,dim6 ) )
vas.cor.ctr.son  <-  array (0, dim = c(4, dim5,dim6 ) )
vas.cor.ctr.ond  <-  array (0, dim = c(4, dim5,dim6 ) )
vas.cor.ctr.ndj  <-  array (0, dim = c(4, dim5,dim6 ) )
###
tos.nud1.ctr.jja  <-  array (0, dim = c(4, dim5,dim6 ) )
tos.nud1.ctr.jas  <-  array (0, dim = c(4, dim5,dim6 ) )
tos.nud1.ctr.aso  <-  array (0, dim = c(4, dim5,dim6 ) )
tos.nud1.ctr.son  <-  array (0, dim = c(4, dim5,dim6 ) )
tos.nud1.ctr.ond  <-  array (0, dim = c(4, dim5,dim6 ) )
tos.nud1.ctr.ndj  <-  array (0, dim = c(4, dim5,dim6 ) )

uas.nud1.ctr.jja  <-  array (0, dim = c(4, dim5,dim6 ) )
uas.nud1.ctr.jas  <-  array (0, dim = c(4, dim5,dim6 ) )
uas.nud1.ctr.aso  <-  array (0, dim = c(4, dim5,dim6 ) )
uas.nud1.ctr.son  <-  array (0, dim = c(4, dim5,dim6 ) )
uas.nud1.ctr.ond  <-  array (0, dim = c(4, dim5,dim6 ) )
uas.nud1.ctr.ndj  <-  array (0, dim = c(4, dim5,dim6 ) )

vas.nud1.ctr.jja  <-  array (0, dim = c(4, dim5,dim6 ) )
vas.nud1.ctr.jas  <-  array (0, dim = c(4, dim5,dim6 ) )
vas.nud1.ctr.aso  <-  array (0, dim = c(4, dim5,dim6 ) )
vas.nud1.ctr.son  <-  array (0, dim = c(4, dim5,dim6 ) )
vas.nud1.ctr.ond  <-  array (0, dim = c(4, dim5,dim6 ) )
vas.nud1.ctr.ndj  <-  array (0, dim = c(4, dim5,dim6 ) )
Corr = s2dverification::Corr 

for ( i  in 1:dim5 )
{ 
for ( j  in 1:dim6 )
{ 
tos.cor.ctr.jja[,i,j]<-Corr(InsertDim(tos.ctr.ano.concat.jja[1,1,,1,i,j],1,1),
                             InsertDim(tos.ctr.atl3.ano.concat.jja[1,1,,1],1,1), posloop=1,poscor=2) 
tos.cor.ctr.jas[,i,j]<-Corr(InsertDim(tos.ctr.ano.concat.jas[1,1,,1,i,j],1,1),
                             InsertDim(tos.ctr.atl3.ano.concat.jja[1,1,,1],1,1), posloop=1,poscor=2 ) 
tos.cor.ctr.aso[,i,j]<-Corr(InsertDim(tos.ctr.ano.concat.aso[1,1,,1,i,j],1,1),
                             InsertDim(tos.ctr.atl3.ano.concat.jja[1,1,,1],1,1), posloop=1,poscor=2 ) 
tos.cor.ctr.son[,i,j]<-Corr(InsertDim(tos.ctr.ano.concat.son[1,1,,1,i,j],1,1),
                             InsertDim(tos.ctr.atl3.ano.concat.jja[1,1,,1],1,1), posloop=1,poscor=2 ) 

uas.cor.ctr.jja[,i,j]<-Corr(InsertDim(uas.ctr.ano.concat.jja[1,1,,1,i,j],1,1),
                             InsertDim(tos.ctr.atl3.ano.concat.jja[1,1,,1],1,1), posloop=1,poscor=2 ) 
uas.cor.ctr.jas[,i,j]<-Corr(InsertDim(uas.ctr.ano.concat.jas[1,1,,1,i,j],1,1),
                             InsertDim(tos.ctr.atl3.ano.concat.jja[1,1,,1],1,1), posloop=1,poscor=2 ) 
uas.cor.ctr.aso[,i,j]<-Corr(InsertDim(uas.ctr.ano.concat.aso[1,1,,1,i,j],1,1),
                             InsertDim(tos.ctr.atl3.ano.concat.jja[1,1,,1],1,1), posloop=1,poscor=2 ) 
uas.cor.ctr.son[,i,j]<-Corr(InsertDim(uas.ctr.ano.concat.son[1,1,,1,i,j],1,1),
                             InsertDim(tos.ctr.atl3.ano.concat.jja[1,1,,1],1,1), posloop=1,poscor=2 ) 

vas.cor.ctr.jja[,i,j]<-Corr(InsertDim(vas.ctr.ano.concat.jja[1,1,,1,i,j],1,1),
                             InsertDim(tos.ctr.atl3.ano.concat.jja[1,1,,1],1,1), posloop=1,poscor=2 ) 
vas.cor.ctr.jas[,i,j]<-Corr(InsertDim(vas.ctr.ano.concat.jas[1,1,,1,i,j],1,1),
                             InsertDim(tos.ctr.atl3.ano.concat.jja[1,1,,1],1,1), posloop=1,poscor=2 ) 
vas.cor.ctr.aso[,i,j]<-Corr(InsertDim(vas.ctr.ano.concat.aso[1,1,,1,i,j],1,1),
                             InsertDim(tos.ctr.atl3.ano.concat.jja[1,1,,1],1,1), posloop=1,poscor=2 ) 
vas.cor.ctr.son[,i,j]<-Corr(InsertDim(vas.ctr.ano.concat.son[1,1,,1,i,j],1,1),
                             InsertDim(tos.ctr.atl3.ano.concat.jja[1,1,,1],1,1), posloop=1,poscor=2 ) 
} 
} 

for ( i  in 1:dim5 )
{ 
for ( j  in 1:dim6 )
{ 
tos.cor.nud1.jja[,i,j]<-Corr(InsertDim(tos.nud1.ano.concat.jja[1,1,,1,i,j],1,1),
                             InsertDim(tos.nud1.atl3.ano.concat.jja[1,1,,1],1,1), posloop=1,poscor=2) 
tos.cor.nud1.jas[,i,j]<-Corr(InsertDim(tos.nud1.ano.concat.jas[1,1,,1,i,j],1,1),
                             InsertDim(tos.nud1.atl3.ano.concat.jja[1,1,,1],1,1), posloop=1,poscor=2 ) 
tos.cor.nud1.aso[,i,j]<-Corr(InsertDim(tos.nud1.ano.concat.aso[1,1,,1,i,j],1,1),
                             InsertDim(tos.nud1.atl3.ano.concat.jja[1,1,,1],1,1), posloop=1,poscor=2 ) 
tos.cor.nud1.son[,i,j]<-Corr(InsertDim(tos.nud1.ano.concat.son[1,1,,1,i,j],1,1),
                             InsertDim(tos.nud1.atl3.ano.concat.jja[1,1,,1],1,1), posloop=1,poscor=2 ) 

uas.cor.nud1.jja[,i,j]<-Corr(InsertDim(uas.nud1.ano.concat.jja[1,1,,1,i,j],1,1),
                             InsertDim(tos.nud1.atl3.ano.concat.jja[1,1,,1],1,1), posloop=1,poscor=2 ) 
uas.cor.nud1.jas[,i,j]<-Corr(InsertDim(uas.nud1.ano.concat.jas[1,1,,1,i,j],1,1),
                             InsertDim(tos.nud1.atl3.ano.concat.jja[1,1,,1],1,1), posloop=1,poscor=2 ) 
uas.cor.nud1.aso[,i,j]<-Corr(InsertDim(uas.nud1.ano.concat.aso[1,1,,1,i,j],1,1),
                             InsertDim(tos.nud1.atl3.ano.concat.jja[1,1,,1],1,1), posloop=1,poscor=2 ) 
uas.cor.nud1.son[,i,j]<-Corr(InsertDim(uas.nud1.ano.concat.son[1,1,,1,i,j],1,1),
                             InsertDim(tos.nud1.atl3.ano.concat.jja[1,1,,1],1,1), posloop=1,poscor=2 ) 

vas.cor.nud1.jja[,i,j]<-Corr(InsertDim(vas.nud1.ano.concat.jja[1,1,,1,i,j],1,1),
                             InsertDim(tos.nud1.atl3.ano.concat.jja[1,1,,1],1,1), posloop=1,poscor=2 ) 
vas.cor.nud1.jas[,i,j]<-Corr(InsertDim(vas.nud1.ano.concat.jas[1,1,,1,i,j],1,1),
                             InsertDim(tos.nud1.atl3.ano.concat.jja[1,1,,1],1,1), posloop=1,poscor=2 ) 
vas.cor.nud1.aso[,i,j]<-Corr(InsertDim(vas.nud1.ano.concat.aso[1,1,,1,i,j],1,1),
                             InsertDim(tos.nud1.atl3.ano.concat.jja[1,1,,1],1,1), posloop=1,poscor=2 ) 
vas.cor.nud1.son[,i,j]<-Corr(InsertDim(vas.nud1.ano.concat.son[1,1,,1,i,j],1,1),
                             InsertDim(tos.nud1.atl3.ano.concat.jja[1,1,,1],1,1), posloop=1,poscor=2 ) 
} 
} 

#--- Bootstrapping 

Rpar=1000
tmp1.ctr = tos.ctr.atl3.ano.concat.jja[1,1,,1] 
tmp2.ctr = array(0,c(Rpar,dim2*dim3)) 
tmp1.nud1 = tos.nud1.atl3.ano.concat.jja[1,1,,1] 
tmp2.nud1 = array(0,c(Rpar,dim2*dim3)) 

 for (j in 1:Rpar) 
 {
	 for (i in 1:113) 
 {
 tmp2.ctr[j,(5*i-4):(5*i)] <- sample(tmp1.ctr[(5*i-4):(5*i)])
 tmp2.nud1[j,(5*i-4):(5*i)] <- sample(tmp1.nud1[(5*i-4):(5*i)])
 }
 }


cor.ctr.jas <- array (0, dim = c(Rpar, 4, dim5,dim6 ) )
cor.ctr.aso <- array (0, dim = c(Rpar, 4, dim5,dim6 ) )
cor.ctr.son <- array (0, dim = c(Rpar, 4, dim5,dim6 ) )

cor.ctr.jas.uas <- array (0, dim = c(Rpar, 4, dim5,dim6 ) )
cor.ctr.aso.uas <- array (0, dim = c(Rpar, 4, dim5,dim6 ) )
cor.ctr.son.uas <- array (0, dim = c(Rpar, 4, dim5,dim6 ) )

cor.ctr.jas.vas <- array (0, dim = c(Rpar, 4, dim5,dim6 ) )
cor.ctr.aso.vas <- array (0, dim = c(Rpar, 4, dim5,dim6 ) )
cor.ctr.son.vas <- array (0, dim = c(Rpar, 4, dim5,dim6 ) )

cor.nud1.jas <- array (0, dim = c(Rpar, 4, dim5,dim6 ) )
cor.nud1.aso <- array (0, dim = c(Rpar, 4, dim5,dim6 ) )
cor.nud1.son <- array (0, dim = c(Rpar, 4, dim5,dim6 ) )

cor.nud1.jas.uas <- array (0, dim = c(Rpar, 4, dim5,dim6 ) )
cor.nud1.aso.uas <- array (0, dim = c(Rpar, 4, dim5,dim6 ) )
cor.nud1.son.uas <- array (0, dim = c(Rpar, 4, dim5,dim6 ) )

cor.nud1.jas.vas <- array (0, dim = c(Rpar, 4, dim5,dim6 ) )
cor.nud1.aso.vas <- array (0, dim = c(Rpar, 4, dim5,dim6 ) )
cor.nud1.son.vas <- array (0, dim = c(Rpar, 4, dim5,dim6 ) )


#--- I need to name the dimensions 
names(dim(tos.ctr.ano.concat.jja))[3] <- 'rep'
names(dim(tos.ctr.ano.concat.jas))[3] <- 'rep'
names(dim(tos.ctr.ano.concat.aso))[3] <- 'rep'
names(dim(tos.ctr.ano.concat.son))[3] <- 'rep'
#
names(dim(tos.ctr.ano.concat.jja))[4] <- 'ftime'
names(dim(tos.ctr.ano.concat.jas))[4] <- 'ftime'
names(dim(tos.ctr.ano.concat.aso))[4] <- 'ftime'
names(dim(tos.ctr.ano.concat.son))[4] <- 'ftime'
#
names(dim(tos.ctr.ano.concat.jja))[5] <- 'lat'
names(dim(tos.ctr.ano.concat.jas))[5] <- 'lat'
names(dim(tos.ctr.ano.concat.aso))[5] <- 'lat'
names(dim(tos.ctr.ano.concat.son))[5] <- 'lat'
#
names(dim(tos.ctr.ano.concat.jja))[6] <- 'lon'
names(dim(tos.ctr.ano.concat.jas))[6] <- 'lon'
names(dim(tos.ctr.ano.concat.aso))[6] <- 'lon'
names(dim(tos.ctr.ano.concat.son))[6] <- 'lon'

names(dim(uas.ctr.ano.concat.jja))[3] <- 'rep'
names(dim(uas.ctr.ano.concat.jas))[3] <- 'rep'
names(dim(uas.ctr.ano.concat.aso))[3] <- 'rep'
names(dim(uas.ctr.ano.concat.son))[3] <- 'rep'
#
names(dim(uas.ctr.ano.concat.jja))[4] <- 'ftime'
names(dim(uas.ctr.ano.concat.jas))[4] <- 'ftime'
names(dim(uas.ctr.ano.concat.aso))[4] <- 'ftime'
names(dim(uas.ctr.ano.concat.son))[4] <- 'ftime'
#
names(dim(uas.ctr.ano.concat.jja))[5] <- 'lat'
names(dim(uas.ctr.ano.concat.jas))[5] <- 'lat'
names(dim(uas.ctr.ano.concat.aso))[5] <- 'lat'
names(dim(uas.ctr.ano.concat.son))[5] <- 'lat'
#
names(dim(uas.ctr.ano.concat.jja))[6] <- 'lon'
names(dim(uas.ctr.ano.concat.jas))[6] <- 'lon'
names(dim(uas.ctr.ano.concat.aso))[6] <- 'lon'
names(dim(uas.ctr.ano.concat.son))[6] <- 'lon'
##################
names(dim(vas.ctr.ano.concat.jja))[3] <- 'rep'
names(dim(vas.ctr.ano.concat.jas))[3] <- 'rep'
names(dim(vas.ctr.ano.concat.aso))[3] <- 'rep'
names(dim(vas.ctr.ano.concat.son))[3] <- 'rep'
#
names(dim(vas.ctr.ano.concat.jja))[4] <- 'ftime'
names(dim(vas.ctr.ano.concat.jas))[4] <- 'ftime'
names(dim(vas.ctr.ano.concat.aso))[4] <- 'ftime'
names(dim(vas.ctr.ano.concat.son))[4] <- 'ftime'
#
names(dim(vas.ctr.ano.concat.jja))[5] <- 'lat'
names(dim(vas.ctr.ano.concat.jas))[5] <- 'lat'
names(dim(vas.ctr.ano.concat.aso))[5] <- 'lat'
names(dim(vas.ctr.ano.concat.son))[5] <- 'lat'
#
names(dim(vas.ctr.ano.concat.jja))[6] <- 'lon'
names(dim(vas.ctr.ano.concat.jas))[6] <- 'lon'
names(dim(vas.ctr.ano.concat.aso))[6] <- 'lon'
names(dim(vas.ctr.ano.concat.son))[6] <- 'lon'
##################

names (dim(tmp2.ctr))[1] <- 'tmp'
names (dim(tmp2.ctr))[2] <- 'member'

names(dim(tos.nud1.ano.concat.jja))[3] <- 'rep'
names(dim(tos.nud1.ano.concat.jas))[3] <- 'rep'
names(dim(tos.nud1.ano.concat.aso))[3] <- 'rep'
names(dim(tos.nud1.ano.concat.son))[3] <- 'rep'
#
names(dim(tos.nud1.ano.concat.jja))[4] <- 'ftime'
names(dim(tos.nud1.ano.concat.jas))[4] <- 'ftime'
names(dim(tos.nud1.ano.concat.aso))[4] <- 'ftime'
names(dim(tos.nud1.ano.concat.son))[4] <- 'ftime'
#
names(dim(tos.nud1.ano.concat.jja))[5] <- 'lat'
names(dim(tos.nud1.ano.concat.jas))[5] <- 'lat'
names(dim(tos.nud1.ano.concat.aso))[5] <- 'lat'
names(dim(tos.nud1.ano.concat.son))[5] <- 'lat'
#
names(dim(tos.nud1.ano.concat.jja))[6] <- 'lon'
names(dim(tos.nud1.ano.concat.jas))[6] <- 'lon'
names(dim(tos.nud1.ano.concat.aso))[6] <- 'lon'
names(dim(tos.nud1.ano.concat.son))[6] <- 'lon'

names(dim(uas.nud1.ano.concat.jja))[3] <- 'rep'
names(dim(uas.nud1.ano.concat.jas))[3] <- 'rep'
names(dim(uas.nud1.ano.concat.aso))[3] <- 'rep'
names(dim(uas.nud1.ano.concat.son))[3] <- 'rep'
#
names(dim(uas.nud1.ano.concat.jja))[4] <- 'ftime'
names(dim(uas.nud1.ano.concat.jas))[4] <- 'ftime'
names(dim(uas.nud1.ano.concat.aso))[4] <- 'ftime'
names(dim(uas.nud1.ano.concat.son))[4] <- 'ftime'
#
names(dim(uas.nud1.ano.concat.jja))[5] <- 'lat'
names(dim(uas.nud1.ano.concat.jas))[5] <- 'lat'
names(dim(uas.nud1.ano.concat.aso))[5] <- 'lat'
names(dim(uas.nud1.ano.concat.son))[5] <- 'lat'
#
names(dim(uas.nud1.ano.concat.jja))[6] <- 'lon'
names(dim(uas.nud1.ano.concat.jas))[6] <- 'lon'
names(dim(uas.nud1.ano.concat.aso))[6] <- 'lon'
names(dim(uas.nud1.ano.concat.son))[6] <- 'lon'
##################
names(dim(vas.nud1.ano.concat.jja))[3] <- 'rep'
names(dim(vas.nud1.ano.concat.jas))[3] <- 'rep'
names(dim(vas.nud1.ano.concat.aso))[3] <- 'rep'
names(dim(vas.nud1.ano.concat.son))[3] <- 'rep'
#
names(dim(vas.nud1.ano.concat.jja))[4] <- 'ftime'
names(dim(vas.nud1.ano.concat.jas))[4] <- 'ftime'
names(dim(vas.nud1.ano.concat.aso))[4] <- 'ftime'
names(dim(vas.nud1.ano.concat.son))[4] <- 'ftime'
#
names(dim(vas.nud1.ano.concat.jja))[5] <- 'lat'
names(dim(vas.nud1.ano.concat.jas))[5] <- 'lat'
names(dim(vas.nud1.ano.concat.aso))[5] <- 'lat'
names(dim(vas.nud1.ano.concat.son))[5] <- 'lat'
#
names(dim(vas.nud1.ano.concat.jja))[6] <- 'lon'
names(dim(vas.nud1.ano.concat.jas))[6] <- 'lon'
names(dim(vas.nud1.ano.concat.aso))[6] <- 'lon'
names(dim(vas.nud1.ano.concat.son))[6] <- 'lon'
##################

names (dim(tmp2.nud1))[1] <- 'tmp'
names (dim(tmp2.nud1))[2] <- 'member'

#--- Correlation function 

cor.fun   <- function(tmpx, tmpy) 
{
	#if (length(tmpx) > 0 && length(tmpy) > 0 && !all(is.na(tmpx)) && !all(is.na(tmpy))) {
	if (any(is.na(tmpx)) || any(is.na(tmpy))) {
		NA
	} else {
	    cor(tmpy,  tmpx , use="pairwise.complete.obs" ) 
	}
} 

# --- Actual correlations for 1000 iterations 
tmp      = tos.ctr.ano.concat.jja
tmp [is.na(tos.ctr.ano.concat.jja)]=0.
cor.ctr.jja <- Apply(list(tmp2.ctr, tmp[,1,,1,, ] ),
			 target_dims = list(2, 1),
			 cor.fun,
			 ncores = 20)
tmp      = tos.ctr.ano.concat.jas
tmp [is.na(tos.ctr.ano.concat.jas)]=0.
cor.ctr.jas <- Apply(list(tmp2.ctr, tmp[,1,,1,, ] ),
			 target_dims = list(2, 1),
			 cor.fun,
			 ncores = 20)
tmp      = tos.ctr.ano.concat.aso
tmp [is.na(tos.ctr.ano.concat.aso)]=0.
cor.ctr.aso <- Apply(list(tmp2.ctr, tmp[,1,,1,, ] ),
			 target_dims = list(2, 1),
			 cor.fun,
			 ncores = 20)
tmp      = tos.ctr.ano.concat.son
tmp [is.na(tos.ctr.ano.concat.son)]=0.
cor.ctr.son <- Apply(list(tmp2.ctr, tmp[,1,,1,, ] ),
			 target_dims = list(2, 1),
			 cor.fun,
			 ncores = 20)
##################
tmp = uas.ctr.ano.concat.jja
tmp [is.na(uas.ctr.ano.concat.jja)]=0.
#cor.ctr.jja.uas <- Apply(list(tmp2.ctr, uas.ctr.ano.concat.jja[,1,,1,, ] ),
cor.ctr.jja.uas <- Apply(list(tmp2.ctr, tmp[,1,,1,, ] ),
			 target_dims = list(2, 1),
			 cor.fun,
			 ncores = 20)
tmp      = uas.ctr.ano.concat.jas
tmp [is.na(uas.ctr.ano.concat.jas)]=0.
cor.ctr.jas.uas <- Apply(list(tmp2.ctr, tmp[,1,,1,, ] ),
			 target_dims = list(2, 1),
			 cor.fun,
			 ncores = 20)
tmp      = uas.ctr.ano.concat.aso
tmp [is.na(uas.ctr.ano.concat.aso)]=0.
cor.ctr.aso.uas <- Apply(list(tmp2.ctr, tmp[,1,,1,, ] ),
			 target_dims = list(2, 1),
			 cor.fun,
			 ncores = 20)
tmp      = uas.ctr.ano.concat.son
tmp [is.na(uas.ctr.ano.concat.son)]=0.
cor.ctr.son.uas <- Apply(list(tmp2.ctr, tmp[,1,,1,, ] ),
			 target_dims = list(2, 1),
			 cor.fun,
			 ncores = 20)

##############
tmp      = vas.ctr.ano.concat.jja
tmp [is.na(vas.ctr.ano.concat.jja)]=0.
cor.ctr.jja.vas <- Apply(list(tmp2.ctr, tmp[,1,,1,, ] ),
			 target_dims = list(2, 1),
			 cor.fun,
			 ncores = 20)
tmp      = vas.ctr.ano.concat.jas
tmp [is.na(vas.ctr.ano.concat.jas)]=0.
cor.ctr.jas.vas <- Apply(list(tmp2.ctr, tmp[,1,,1,, ] ),
			 target_dims = list(2, 1),
			 cor.fun,
			 ncores = 20)
tmp      = vas.ctr.ano.concat.aso
tmp [is.na(vas.ctr.ano.concat.aso)]=0.
cor.ctr.aso.vas <- Apply(list(tmp2.ctr, tmp[,1,,1,, ] ),
			 target_dims = list(2, 1),
			 cor.fun,
			 ncores = 20)
tmp      = vas.ctr.ano.concat.son
tmp [is.na(vas.ctr.ano.concat.son)]=0.
cor.ctr.son.vas <- Apply(list(tmp2.ctr, tmp[,1,,1,, ] ),
			 target_dims = list(2, 1),
			 cor.fun,
			 ncores = 20)
##################
tmp      = tos.nud1.ano.concat.jja
tmp [is.na(tos.nud1.ano.concat.jja)]=0.
cor.nud1.jja <- Apply(list(tmp2.nud1, tmp[,1,,1,, ] ),
			 target_dims = list(2, 1),
			 cor.fun,
			 ncores = 20)
tmp      = tos.nud1.ano.concat.jas
tmp [is.na(tos.nud1.ano.concat.jas)]=0.
cor.nud1.jas <- Apply(list(tmp2.nud1, tmp[,1,,1,, ] ),
			 target_dims = list(2, 1),
			 cor.fun,
			 ncores = 20)
tmp      = tos.nud1.ano.concat.aso
tmp [is.na(tos.nud1.ano.concat.aso)]=0.
cor.nud1.aso <- Apply(list(tmp2.nud1, tmp[,1,,1,, ] ),
			 target_dims = list(2, 1),
			 cor.fun,
			 ncores = 20)
tmp      = tos.nud1.ano.concat.son
tmp [is.na(tos.nud1.ano.concat.son)]=0.
cor.nud1.son <- Apply(list(tmp2.nud1, tmp[,1,,1,, ] ),
			 target_dims = list(2, 1),
			 cor.fun,
			 ncores = 20)
##################
tmp = uas.nud1.ano.concat.jja
tmp [is.na(uas.nud1.ano.concat.jja)]=0.
#cor.nud1.jja.uas <- Apply(list(tmp2.nud1, uas.nud1.ano.concat.jja[,1,,1,, ] ),
cor.nud1.jja.uas <- Apply(list(tmp2.nud1, tmp[,1,,1,, ] ),
			 target_dims = list(2, 1),
			 cor.fun,
			 ncores = 20)
tmp      = uas.nud1.ano.concat.jas
tmp [is.na(uas.nud1.ano.concat.jas)]=0.
cor.nud1.jas.uas <- Apply(list(tmp2.nud1, tmp[,1,,1,, ] ),
			 target_dims = list(2, 1),
			 cor.fun,
			 ncores = 20)
tmp      = uas.nud1.ano.concat.aso
tmp [is.na(uas.nud1.ano.concat.aso)]=0.
cor.nud1.aso.uas <- Apply(list(tmp2.nud1, tmp[,1,,1,, ] ),
			 target_dims = list(2, 1),
			 cor.fun,
			 ncores = 20)
tmp      = uas.nud1.ano.concat.son
tmp [is.na(uas.nud1.ano.concat.son)]=0.
cor.nud1.son.uas <- Apply(list(tmp2.nud1, tmp[,1,,1,, ] ),
			 target_dims = list(2, 1),
			 cor.fun,
			 ncores = 20)

##############
tmp      = vas.nud1.ano.concat.jja
tmp [is.na(vas.nud1.ano.concat.jja)]=0.
cor.nud1.jja.vas <- Apply(list(tmp2.nud1, tmp[,1,,1,, ] ),
			 target_dims = list(2, 1),
			 cor.fun,
			 ncores = 20)
tmp      = vas.nud1.ano.concat.jas
tmp [is.na(vas.nud1.ano.concat.jas)]=0.
cor.nud1.jas.vas <- Apply(list(tmp2.nud1, tmp[,1,,1,, ] ),
			 target_dims = list(2, 1),
			 cor.fun,
			 ncores = 20)
tmp      = vas.nud1.ano.concat.aso
tmp [is.na(vas.nud1.ano.concat.aso)]=0.
cor.nud1.aso.vas <- Apply(list(tmp2.nud1, tmp[,1,,1,, ] ),
			 target_dims = list(2, 1),
			 cor.fun,
			 ncores = 20)
tmp      = vas.nud1.ano.concat.son
tmp [is.na(vas.nud1.ano.concat.son)]=0.
cor.nud1.son.vas <- Apply(list(tmp2.nud1, tmp[,1,,1,, ] ),
			 target_dims = list(2, 1),
			 cor.fun,
			 ncores = 20)
##################


## -- Do the ranking --- ###############

signif.jja <- array (0, dim = c(dim5,dim6 ) )
signif.jas <- array (0, dim = c(dim5,dim6 ) )
signif.aso <- array (0, dim = c(dim5,dim6 ) )
signif.son <- array (0, dim = c(dim5,dim6 ) )

#
diff.jja <- cor.nud1.jja$output1 - cor.ctr.jja$output1
diff.jas <- cor.nud1.jas$output1 - cor.ctr.jas$output1
diff.aso <- cor.nud1.aso$output1 - cor.ctr.aso$output1
diff.son <- cor.nud1.son$output1 - cor.ctr.son$output1


signif.jja <- apply(diff.jja[,,],c(2,3),function(x)(quantile(x,0.05,na.rm=TRUE)>=0)||(quantile(x,0.95,na.rm=TRUE)<=0))
signif.jas <- apply(diff.jas[,,],c(2,3),function(x)(quantile(x,0.05,na.rm=TRUE)>=0)||(quantile(x,0.95,na.rm=TRUE)<=0))
signif.aso <- apply(diff.aso[,,],c(2,3),function(x)(quantile(x,0.05,na.rm=TRUE)>=0)||(quantile(x,0.95,na.rm=TRUE)<=0))
signif.son <- apply(diff.son[,,],c(2,3),function(x)(quantile(x,0.05,na.rm=TRUE)>=0)||(quantile(x,0.95,na.rm=TRUE)<=0))

signif.jja.uas <- array (0, dim = c(dim5,dim6 ) )
signif.jas.uas <- array (0, dim = c(dim5,dim6 ) )
signif.aso.uas <- array (0, dim = c(dim5,dim6 ) )
signif.son.uas <- array (0, dim = c(dim5,dim6 ) )

signif.jja.vas <- array (0, dim = c(dim5,dim6 ) )
signif.jas.vas <- array (0, dim = c(dim5,dim6 ) )
signif.aso.vas <- array (0, dim = c(dim5,dim6 ) )
signif.son.vas <- array (0, dim = c(dim5,dim6 ) )

diff.jja.uas <- cor.nud1.jja.uas$output1 - cor.ctr.jja.uas$output1
diff.jas.uas <- cor.nud1.jas.uas$output1 - cor.ctr.jas.uas$output1
diff.aso.uas <- cor.nud1.aso.uas$output1 - cor.ctr.aso.uas$output1
diff.son.uas <- cor.nud1.son.uas$output1 - cor.ctr.son.uas$output1
            
diff.jja.vas <- cor.nud1.jja.vas$output1 - cor.ctr.jja.vas$output1
diff.jas.vas <- cor.nud1.jas.vas$output1 - cor.ctr.jas.vas$output1
diff.aso.vas <- cor.nud1.aso.vas$output1 - cor.ctr.aso.vas$output1
diff.son.vas <- cor.nud1.son.vas$output1 - cor.ctr.son.vas$output1

signif.jja.uas <- apply(diff.jja.uas[,,],c(2,3),function(x)(quantile(x,0.05,na.rm=TRUE)>=0)||(quantile(x,0.95,na.rm=TRUE)<=0))
signif.jas.uas <- apply(diff.jas.uas[,,],c(2,3),function(x)(quantile(x,0.05,na.rm=TRUE)>=0)||(quantile(x,0.95,na.rm=TRUE)<=0))
signif.aso.uas <- apply(diff.aso.uas[,,],c(2,3),function(x)(quantile(x,0.05,na.rm=TRUE)>=0)||(quantile(x,0.95,na.rm=TRUE)<=0))
signif.son.uas <- apply(diff.son.uas[,,],c(2,3),function(x)(quantile(x,0.05,na.rm=TRUE)>=0)||(quantile(x,0.95,na.rm=TRUE)<=0))

signif.jja.vas <- apply(diff.jja.vas[,,],c(2,3),function(x)(quantile(x,0.05,na.rm=TRUE)>=0)||(quantile(x,0.95,na.rm=TRUE)<=0))
signif.jas.vas <- apply(diff.jas.vas[,,],c(2,3),function(x)(quantile(x,0.05,na.rm=TRUE)>=0)||(quantile(x,0.95,na.rm=TRUE)<=0))
signif.aso.vas <- apply(diff.aso.vas[,,],c(2,3),function(x)(quantile(x,0.05,na.rm=TRUE)>=0)||(quantile(x,0.95,na.rm=TRUE)<=0))
signif.son.vas <- apply(diff.son.vas[,,],c(2,3),function(x)(quantile(x,0.05,na.rm=TRUE)>=0)||(quantile(x,0.95,na.rm=TRUE)<=0))


#---  Saving the correlations on netCDF files for plotting with a
#     separate NCL script 

## ----------- Common part for all files ---------------------------- ###############
# define dimensions
timenc =seq (1,1)
londim     <-   ncdim_def("lon","degrees_east",as.double(tos.lon )) 
latdim     <-   ncdim_def("lat","degrees_north",as.double(tos.lat)) 
timedim    <-   ncdim_def("time","days since 1950-01-01 00:00:00",timenc,unlim=TRUE)

# define variables
fillvalue=1.e+36
tosname <- "Sea Surf Temperature"
uasname <- "Surface zonal wind"
vasname <- "Surface meridional wind"
dotname <- "Significance (1 for significant, 0 for not )"

tmp.def.tos <- ncvar_def("tos","deg_C",list(londim,latdim,timedim),fillvalue,tosname,prec="single")
tmp.def.uas <- ncvar_def("uas","m/s",list(londim,latdim,timedim),fillvalue,uasname,prec="single")
tmp.def.vas <- ncvar_def("vas","m/s",list(londim,latdim,timedim),fillvalue,vasname,prec="single")
#

save.cor.file  <- function(exp) {

 #1        CTR 
 #2        NUD1
 #3        NUD2

expname  = switch( exp, "ctr", "nud1", "nud2") 
### JJA  #############################################
tos.tmp2 = switch( exp,  tos.cor.ctr.jja, tos.cor.nud1.jja, tos.cor.nud2.jja ) 
uas.tmp2 = switch( exp,  uas.cor.ctr.jja, uas.cor.nud1.jja, uas.cor.nud2.jja ) 
vas.tmp2 = switch( exp,  vas.cor.ctr.jja, vas.cor.nud1.jja, vas.cor.nud2.jja ) 

tos.tmp = tos.tmp2[2,, ]
tos.tmp [ tos.tmp2[1,,] * tos.tmp2[3,,] <=0 ]= NA 
tos.tmp [which ( is.na ( tos.ctr$mod[1, 1, 1, 1, ,] ) )] = NA  

uas.tmp = uas.tmp2[2,,]
vas.tmp = vas.tmp2[2,,]
uas.tmp[,] = NA
vas.tmp[,] = NA

uas.tmp3 = uas.tmp2
vas.tmp3 = vas.tmp2
uas.tmp3 [which (is.na(uas.tmp2))] = 0 
vas.tmp3 [which (is.na(vas.tmp2))] = 0 
for ( i  in 1:dim5   )
{ 
for ( j  in 1:dim6   )
{ 
if (         (uas.tmp3[1,i,j] * uas.tmp3[3,i,j] >= 0 ) 
	 	   &&  (vas.tmp3[ 1 , i, j ] * vas.tmp3[3,i,j] >= 0 ) )  
		  { 
   uas.tmp [i,j] = uas.tmp2[2,i,j]
   vas.tmp [i,j] = vas.tmp2[2,i,j]
		  } 
} 
} 
# create netCDF file and put arrays
ncfname.tos <- paste0("correlation_tos_",expname,"_jja.nc")
ncfname.uas <- paste0("correlation_uas_",expname,"_jja.nc")
ncfname.vas <- paste0("correlation_vas_",expname,"_jja.nc")

ncout.tos <- nc_create(ncfname.tos,list(tmp.def.tos))
ncout.uas <- nc_create(ncfname.uas,list(tmp.def.uas))
ncout.vas <- nc_create(ncfname.vas,list(tmp.def.vas))

# put variables
ncvar_put(ncout.tos,tmp.def.tos,t(tos.tmp)) ## transpose of tmp 
ncvar_put(ncout.uas,tmp.def.uas,t(uas.tmp))
ncvar_put(ncout.vas,tmp.def.vas,t(vas.tmp))

ncatt_put(ncout.tos,"lon","axis","X") #,verbose=FALSE) #,definemode=FALSE)
ncatt_put(ncout.tos,"lat","axis","Y")
ncatt_put(ncout.tos,"time","axis","T")
##
ncatt_put(ncout.uas,"lon","axis","X") #,verbose=FALSE) #,definemode=FALSE)
ncatt_put(ncout.uas,"lat","axis","Y")
ncatt_put(ncout.uas,"time","axis","T")
##
ncatt_put(ncout.vas,"lon","axis","X") #,verbose=FALSE) #,definemode=FALSE)
ncatt_put(ncout.vas,"lat","axis","Y")
ncatt_put(ncout.vas,"time","axis","T")
##

# close the file, writing data to disk
nc_close(ncout.tos)
nc_close(ncout.uas)
nc_close(ncout.vas)

### JAS  #############################################
tos.tmp2 = switch( exp,  tos.cor.ctr.jas, tos.cor.nud1.jas, tos.cor.nud2.jas ) 
uas.tmp2 = switch( exp,  uas.cor.ctr.jas, uas.cor.nud1.jas, uas.cor.nud2.jas ) 
vas.tmp2 = switch( exp,  vas.cor.ctr.jas, vas.cor.nud1.jas, vas.cor.nud2.jas ) 

tos.tmp = tos.tmp2[2,, ]
tos.tmp [ tos.tmp2[1,,] * tos.tmp2[3,,] <=0 ]= NA 
tos.tmp [which ( is.na ( tos.ctr$mod[1, 1, 1, 1, ,] ) )] = NA  

uas.tmp = uas.tmp2[2,,]
vas.tmp = vas.tmp2[2,,]
uas.tmp[,] = NA
vas.tmp[,] = NA

uas.tmp3 = uas.tmp2
vas.tmp3 = vas.tmp2
uas.tmp3 [which (is.na(uas.tmp2))] = 0 
vas.tmp3 [which (is.na(vas.tmp2))] = 0 

for ( i  in 1:dim5   )
{ 
for ( j  in 1:dim6   )
{ 
if (         (uas.tmp3[1,i,j] * uas.tmp3[3,i,j] >= 0 ) 
	 	   &&  (vas.tmp3[ 1 , i, j ] * vas.tmp3[3,i,j] >= 0 ) )  
		  { 
   uas.tmp [i,j] = uas.tmp2[2,i,j]
   vas.tmp [i,j] = vas.tmp2[2,i,j]
		  } 
} 
} 
# create netCDF file and put arrays
ncfname.tos <- paste0("correlation_tos_",expname,"_jas.nc")
ncfname.uas <- paste0("correlation_uas_",expname,"_jas.nc")
ncfname.vas <- paste0("correlation_vas_",expname,"_jas.nc")
ncout.tos <- nc_create(ncfname.tos,list(tmp.def.tos))
ncout.uas <- nc_create(ncfname.uas,list(tmp.def.uas))
ncout.vas <- nc_create(ncfname.vas,list(tmp.def.vas))

# put variables
ncvar_put(ncout.tos,tmp.def.tos,t(tos.tmp)) ## transpose of tmp 
ncvar_put(ncout.uas,tmp.def.uas,t(uas.tmp))
ncvar_put(ncout.vas,tmp.def.vas,t(vas.tmp))
# put additional attributes into dimension and data variables
ncatt_put(ncout.tos,"lon","axis","X") #,verbose=FALSE) #,definemode=FALSE)
ncatt_put(ncout.tos,"lat","axis","Y")
ncatt_put(ncout.tos,"time","axis","T")
##
ncatt_put(ncout.uas,"lon","axis","X") #,verbose=FALSE) #,definemode=FALSE)
ncatt_put(ncout.uas,"lat","axis","Y")
ncatt_put(ncout.uas,"time","axis","T")
##
ncatt_put(ncout.vas,"lon","axis","X") #,verbose=FALSE) #,definemode=FALSE)
ncatt_put(ncout.vas,"lat","axis","Y")
ncatt_put(ncout.vas,"time","axis","T")
##

# close the file, writing data to disk
nc_close(ncout.tos)
nc_close(ncout.uas)
nc_close(ncout.vas)
#nc_close(ncout.dot)

### ASO  #############################################
tos.tmp2 = switch( exp,  tos.cor.ctr.aso, tos.cor.nud1.aso, tos.cor.nud2.aso ) 
uas.tmp2 = switch( exp,  uas.cor.ctr.aso, uas.cor.nud1.aso, uas.cor.nud2.aso ) 
vas.tmp2 = switch( exp,  vas.cor.ctr.aso, vas.cor.nud1.aso, vas.cor.nud2.aso ) 

tos.tmp = tos.tmp2[2,,]
tos.tmp [ tos.tmp2[1,,] * tos.tmp2[3,,] <=0 ]= NA 
tos.tmp [which ( is.na ( tos.ctr$mod[1, 1, 1, 1, ,] ) )] = NA  

uas.tmp = uas.tmp2[2,,]
vas.tmp = vas.tmp2[2,,]
uas.tmp[,] = NA
vas.tmp[,] = NA

uas.tmp3 = uas.tmp2
vas.tmp3 = vas.tmp2
uas.tmp3 [which (is.na(uas.tmp2))] = 0 
vas.tmp3 [which (is.na(vas.tmp2))] = 0 


for ( i  in 1:dim5   )
{ 
for ( j  in 1:dim6   )
{ 
if (         (uas.tmp3[1,i,j] * uas.tmp3[3,i,j] >= 0 ) 
	 	   &&  (vas.tmp3[ 1 , i, j ] * vas.tmp3[3,i,j] >= 0 ) )  
		  { 
   uas.tmp [i,j] = uas.tmp2[2,i,j]
   vas.tmp [i,j] = vas.tmp2[2,i,j]
		  } 
} 
} 
# create netCDF file and put arrays
ncfname.tos <- paste0("correlation_tos_",expname,"_aso.nc")
ncfname.uas <- paste0("correlation_uas_",expname,"_aso.nc")
ncfname.vas <- paste0("correlation_vas_",expname,"_aso.nc")
ncout.tos <- nc_create(ncfname.tos,list(tmp.def.tos))
ncout.uas <- nc_create(ncfname.uas,list(tmp.def.uas))
ncout.vas <- nc_create(ncfname.vas,list(tmp.def.vas))

# put variables
ncvar_put(ncout.tos,tmp.def.tos,t(tos.tmp)) ## transpose of tmp 
ncvar_put(ncout.uas,tmp.def.uas,t(uas.tmp))
ncvar_put(ncout.vas,tmp.def.vas,t(vas.tmp))
# put additional attributes into dimension and data variables
ncatt_put(ncout.tos,"lon","axis","X") #,verbose=FALSE) #,definemode=FALSE)
ncatt_put(ncout.tos,"lat","axis","Y")
ncatt_put(ncout.tos,"time","axis","T")
##
ncatt_put(ncout.uas,"lon","axis","X") #,verbose=FALSE) #,definemode=FALSE)
ncatt_put(ncout.uas,"lat","axis","Y")
ncatt_put(ncout.uas,"time","axis","T")
##
ncatt_put(ncout.vas,"lon","axis","X") #,verbose=FALSE) #,definemode=FALSE)
ncatt_put(ncout.vas,"lat","axis","Y")
ncatt_put(ncout.vas,"time","axis","T")
##

# close the file, writing data to disk
nc_close(ncout.tos)
nc_close(ncout.uas)
nc_close(ncout.vas)
#nc_close(ncout.dot)

### SON  #############################################
tos.tmp2 = switch( exp, tos.cor.ctr.son, tos.cor.nud1.son, tos.cor.nud2.son ) 
uas.tmp2 = switch( exp, uas.cor.ctr.son, uas.cor.nud1.son, uas.cor.nud2.son ) 
vas.tmp2 = switch( exp, vas.cor.ctr.son, vas.cor.nud1.son, vas.cor.nud2.son ) 


tos.tmp = tos.tmp2[2,,]
tos.tmp [ tos.tmp2[1,,] * tos.tmp2[3,,] <=0 ]= NA 
tos.tmp [which ( is.na ( tos.ctr$mod[1, 1, 1, 1, ,] ) )] = NA  

uas.tmp = uas.tmp2[2,,]
vas.tmp = vas.tmp2[2,,]
uas.tmp[,] = NA
vas.tmp[,] = NA

uas.tmp3 = uas.tmp2
vas.tmp3 = vas.tmp2
uas.tmp3 [which (is.na(uas.tmp2))] = 0 
vas.tmp3 [which (is.na(vas.tmp2))] = 0 

for ( i  in 1:dim5   )
{ 
for ( j  in 1:dim6   )
{ 
if (         (uas.tmp3[1,i,j] * uas.tmp3[3,i,j] >= 0 ) 
	 	  &&  (vas.tmp3[ 1 , i, j ] * vas.tmp3[3,i,j] >= 0 ) )  
		  { 
   uas.tmp [i,j] = uas.tmp2[2,i,j]
   vas.tmp [i,j] = vas.tmp2[2,i,j]
		  } 
} 
} 
# create netCDF file and put arrays
ncfname.tos <- paste0("../Ncfiles/regression_correlation_ATL3/correlation_tos_",expname,"_son.nc")
ncfname.uas <- paste0("../Ncfiles/regression_correlation_ATL3/correlation_uas_",expname,"_son.nc")
ncfname.vas <- paste0("../Ncfiles/regression_correlation_ATL3/correlation_vas_",expname,"_son.nc")
ncout.tos <- nc_create(ncfname.tos,list(tmp.def.tos))
ncout.uas <- nc_create(ncfname.uas,list(tmp.def.uas))
ncout.vas <- nc_create(ncfname.vas,list(tmp.def.vas))

# put variables
ncvar_put(ncout.tos,tmp.def.tos,t(tos.tmp)) ## transpose of tmp 
ncvar_put(ncout.uas,tmp.def.uas,t(uas.tmp))
ncvar_put(ncout.vas,tmp.def.vas,t(vas.tmp))
# put additional attributes into dimension and data variables
ncatt_put(ncout.tos,"lon","axis","X") #,verbose=FALSE) #,definemode=FALSE)
ncatt_put(ncout.tos,"lat","axis","Y")
ncatt_put(ncout.tos,"time","axis","T")
##
ncatt_put(ncout.uas,"lon","axis","X") #,verbose=FALSE) #,definemode=FALSE)
ncatt_put(ncout.uas,"lat","axis","Y")
ncatt_put(ncout.uas,"time","axis","T")
##
ncatt_put(ncout.vas,"lon","axis","X") #,verbose=FALSE) #,definemode=FALSE)
ncatt_put(ncout.vas,"lat","axis","Y")
ncatt_put(ncout.vas,"time","axis","T")
##

# close the file, writing data to disk
nc_close(ncout.tos)
nc_close(ncout.uas)
nc_close(ncout.vas)
#nc_close(ncout.dot)

} # finish the function


#--- Saving the significance assessment bassed on the bootstrap, on separate netCDF files 

## ----------- Common part for all files ---------------------------- 

## start function 
save.dot.file  <- function(exp) {

 #1        NUD1
 #2        NUD2

# define dimensions

## define variables
fillvalue=1.e+36
expname  = switch( exp, "nud1","nud2") 
#

### JJA  #############################################
tos.tmp1  = switch( exp, signif.jja     , signif.jja2 )
uas.tmp1  = switch( exp, signif.jja.uas , signif.jja.uas2 )
vas.tmp1  = switch( exp, signif.jja.vas , signif.jja.vas2 )

tos.tmp2 = switch(exp,tos.cor.nud1.jja[1,,]*tos.cor.nud1.jja[3,,],tos.cor.nud2.jja[1,,]*tos.cor.nud2.jja[3,,])
uas.tmp2 = switch(exp,uas.cor.nud1.jja[1,,]*uas.cor.nud1.jja[3,,],uas.cor.nud2.jja[1,,]*uas.cor.nud2.jja[3,,])
vas.tmp2 = switch(exp,vas.cor.nud1.jja[1,,]*vas.cor.nud1.jja[3,,],vas.cor.nud2.jja[1,,]*vas.cor.nud2.jja[3,,])

tos.tmp3 = switch(exp,tos.cor.nud1.jja[2,,]-tos.cor.ctr.jja[2,,],tos.cor.nud2.jja[2,,]-tos.cor.ctr.jja[2,,] ) 
uas.tmp3 = switch(exp,uas.cor.nud1.jja[2,,]-uas.cor.ctr.jja[2,,],uas.cor.nud2.jja[2,,]-uas.cor.ctr.jja[2,,] ) 
vas.tmp3 = switch(exp,vas.cor.nud1.jja[2,,]-vas.cor.ctr.jja[2,,],vas.cor.nud2.jja[2,,]-vas.cor.ctr.jja[2,,] ) 

tos.tmp3[which ( is.na ( tos.nud1$mod[1, 1, 1, 1, ,] ) )] = NA
tos.tmp3[tos.tmp2 <= 0.] = NA
tos.tmp3[tos.tmp1 == F] = NA

uas.tmp3[uas.tmp2 <= 0. ] = NA
uas.tmp3[uas.tmp1 == F  ] = NA
vas.tmp3[uas.tmp2 <= 0.] = NA
vas.tmp3[uas.tmp1 == F] = NA
uas.tmp3[vas.tmp2 <= 0. ] = NA
uas.tmp3[vas.tmp1 == F  ] = NA
vas.tmp3[vas.tmp2 <= 0.] = NA
vas.tmp3[vas.tmp1 == F] = NA

# create netCDF file and put arrays
ncfname.tos.diff <- paste0("correlation_diff_tos_",expname,"_jja.nc")
ncfname.uas.diff <- paste0("correlation_diff_uas_",expname,"_jja.nc")
ncfname.vas.diff <- paste0("correlation_diff_vas_",expname,"_jja.nc")

ncout.tos.diff <- nc_create(ncfname.tos.diff,list(tmp.def.tos))
ncout.uas.diff <- nc_create(ncfname.uas.diff,list(tmp.def.uas))
ncout.vas.diff <- nc_create(ncfname.vas.diff,list(tmp.def.vas))

# put variables
ncvar_put(ncout.tos.diff, tmp.def.tos, t(tos.tmp3))
ncvar_put(ncout.uas.diff, tmp.def.uas, t(uas.tmp3))
ncvar_put(ncout.vas.diff, tmp.def.vas, t(vas.tmp3))

# put additional attributes into dimension and data variables
ncatt_put(ncout.tos.diff,"lon","axis","X") #,verbose=FALSE) #,definemode=FALSE)
ncatt_put(ncout.tos.diff,"lat","axis","Y")

ncatt_put(ncout.uas.diff,"lon","axis","X") #,verbose=FALSE) #,definemode=FALSE)
ncatt_put(ncout.uas.diff,"lat","axis","Y")

ncatt_put(ncout.vas.diff,"lon","axis","X") #,verbose=FALSE) #,definemode=FALSE)
ncatt_put(ncout.vas.diff,"lat","axis","Y")
##

# close the file, writing data to disk
nc_close(ncout.tos.diff)
nc_close(ncout.uas.diff)
nc_close(ncout.vas.diff)


### JAS  #############################################
tos.tmp1  = switch( exp, signif.jas     , signif.jas2 )
uas.tmp1  = switch( exp, signif.jas.uas , signif.jas.uas2 )
vas.tmp1  = switch( exp, signif.jas.vas , signif.jas.vas2 )

tos.tmp2 = switch(exp,tos.cor.nud1.jas[1,,]*tos.cor.nud1.jas[3,,],tos.cor.nud2.jas[1,,]*tos.cor.nud2.jas[3,,])
uas.tmp2 = switch(exp,uas.cor.nud1.jas[1,,]*uas.cor.nud1.jas[3,,],uas.cor.nud2.jas[1,,]*uas.cor.nud2.jas[3,,])
vas.tmp3 = switch(exp,vas.cor.nud1.jas[1,,]*vas.cor.nud1.jas[3,,],vas.cor.nud2.jas[1,,]*vas.cor.nud2.jas[3,,])

tos.tmp3 = switch(exp,tos.cor.nud1.jas[2,,]-tos.cor.ctr.jas[2,,],tos.cor.nud2.jas[2,,]-tos.cor.ctr.jas[2,,] ) 
uas.tmp3 = switch(exp,uas.cor.nud1.jas[2,,]-uas.cor.ctr.jas[2,,],uas.cor.nud2.jas[2,,]-uas.cor.ctr.jas[2,,] ) 
vas.tmp3 = switch(exp,vas.cor.nud1.jas[2,,]-vas.cor.ctr.jas[2,,],vas.cor.nud2.jas[2,,]-vas.cor.ctr.jas[2,,] ) 

tos.tmp3[which ( is.na ( tos.nud1$mod[1, 1, 1, 1, ,] ) )] = NA
tos.tmp3[tos.tmp2 <= 0.] = NA
tos.tmp3[tos.tmp1 == F] = NA

uas.tmp3[uas.tmp2 <= 0. ] = NA
uas.tmp3[uas.tmp1 == F  ] = NA
vas.tmp3[uas.tmp2 <= 0.] = NA
vas.tmp3[uas.tmp1 == F] = NA
uas.tmp3[vas.tmp2 <= 0. ] = NA
uas.tmp3[vas.tmp1 == F  ] = NA
vas.tmp3[vas.tmp2 <= 0.] = NA
vas.tmp3[vas.tmp1 == F] = NA


# create netCDF file and put arrays
ncfname.tos.diff <- paste0("correlation_diff_tos_",expname,"_jas.nc")
ncfname.uas.diff <- paste0("correlation_diff_uas_",expname,"_jas.nc")
ncfname.vas.diff <- paste0("correlation_diff_vas_",expname,"_jas.nc")

ncout.tos.diff <- nc_create(ncfname.tos.diff,list(tmp.def.tos))
ncout.uas.diff <- nc_create(ncfname.uas.diff,list(tmp.def.uas))
ncout.vas.diff <- nc_create(ncfname.vas.diff,list(tmp.def.vas))

# put variables
ncvar_put(ncout.tos.diff, tmp.def.tos, t(tos.tmp3))
ncvar_put(ncout.uas.diff, tmp.def.uas, t(uas.tmp3))
ncvar_put(ncout.vas.diff, tmp.def.vas, t(vas.tmp3))

# put additional attributes into dimension and data variables
ncatt_put(ncout.tos.diff,"lon","axis","X") #,verbose=FALSE) #,definemode=FALSE)
ncatt_put(ncout.tos.diff,"lat","axis","Y")

ncatt_put(ncout.uas.diff,"lon","axis","X") #,verbose=FALSE) #,definemode=FALSE)
ncatt_put(ncout.uas.diff,"lat","axis","Y")

ncatt_put(ncout.vas.diff,"lon","axis","X") #,verbose=FALSE) #,definemode=FALSE)
ncatt_put(ncout.vas.diff,"lat","axis","Y")
##

# close the file, writing data to disk
nc_close(ncout.tos.diff)
nc_close(ncout.uas.diff)
nc_close(ncout.vas.diff)


### ASO  #############################################
tos.tmp1  = switch( exp, signif.aso     , signif.aso2 )
uas.tmp1  = switch( exp, signif.aso.uas , signif.aso.uas2 )
vas.tmp1  = switch( exp, signif.aso.vas , signif.aso.vas2 )

tos.tmp2 = switch(exp,tos.cor.nud1.aso[1,,]*tos.cor.nud1.aso[3,,],tos.cor.nud2.aso[1,,]*tos.cor.nud2.aso[3,,])
uas.tmp2 = switch(exp,uas.cor.nud1.aso[1,,]*uas.cor.nud1.aso[3,,],uas.cor.nud2.aso[1,,]*uas.cor.nud2.aso[3,,])
vas.tmp3 = switch(exp,vas.cor.nud1.aso[1,,]*vas.cor.nud1.aso[3,,],vas.cor.nud2.aso[1,,]*vas.cor.nud2.aso[3,,])

tos.tmp3 = switch(exp,tos.cor.nud1.aso[2,,]-tos.cor.ctr.aso[2,,],tos.cor.nud2.aso[2,,]-tos.cor.ctr.aso[2,,] ) 
uas.tmp3 = switch(exp,uas.cor.nud1.aso[2,,]-uas.cor.ctr.aso[2,,],uas.cor.nud2.aso[2,,]-uas.cor.ctr.aso[2,,] ) 
vas.tmp3 = switch(exp,vas.cor.nud1.aso[2,,]-vas.cor.ctr.aso[2,,],vas.cor.nud2.aso[2,,]-vas.cor.ctr.aso[2,,] ) 

tos.tmp3[which ( is.na ( tos.nud1$mod[1, 1, 1, 1, ,] ) )] = NA
tos.tmp3[tos.tmp2 <= 0.] = NA
tos.tmp3[tos.tmp1 == F] = NA


uas.tmp3[uas.tmp2 <= 0. ] = NA
uas.tmp3[uas.tmp1 == F  ] = NA
vas.tmp3[uas.tmp2 <= 0.] = NA
vas.tmp3[uas.tmp1 == F] = NA
uas.tmp3[vas.tmp2 <= 0. ] = NA
uas.tmp3[vas.tmp1 == F  ] = NA
vas.tmp3[vas.tmp2 <= 0.] = NA
vas.tmp3[vas.tmp1 == F] = NA



# create netCDF file and put arrays
ncfname.tos.diff <- paste0("correlation_diff_tos_",expname,"_aso.nc")
ncfname.uas.diff <- paste0("correlation_diff_uas_",expname,"_aso.nc")
ncfname.vas.diff <- paste0("correlation_diff_vas_",expname,"_aso.nc")

ncout.tos.diff <- nc_create(ncfname.tos.diff,list(tmp.def.tos))
ncout.uas.diff <- nc_create(ncfname.uas.diff,list(tmp.def.uas))
ncout.vas.diff <- nc_create(ncfname.vas.diff,list(tmp.def.vas))

# put variables
ncvar_put(ncout.tos.diff, tmp.def.tos, t(tos.tmp3))
ncvar_put(ncout.uas.diff, tmp.def.uas, t(uas.tmp3))
ncvar_put(ncout.vas.diff, tmp.def.vas, t(vas.tmp3))

# put additional attributes into dimension and data variables
ncatt_put(ncout.tos.diff,"lon","axis","X") #,verbose=FALSE) #,definemode=FALSE)
ncatt_put(ncout.tos.diff,"lat","axis","Y")

ncatt_put(ncout.uas.diff,"lon","axis","X") #,verbose=FALSE) #,definemode=FALSE)
ncatt_put(ncout.uas.diff,"lat","axis","Y")

ncatt_put(ncout.vas.diff,"lon","axis","X") #,verbose=FALSE) #,definemode=FALSE)
ncatt_put(ncout.vas.diff,"lat","axis","Y")
##

# close the file, writing data to disk
nc_close(ncout.tos.diff)
nc_close(ncout.uas.diff)
nc_close(ncout.vas.diff)


### SON  #############################################
tos.tmp1  = switch( exp, signif.son     , signif.son2 )
uas.tmp1  = switch( exp, signif.son.uas , signif.son.uas2 )
vas.tmp1  = switch( exp, signif.son.vas , signif.son.vas2 )

tos.tmp2 = switch(exp,tos.cor.nud1.son[1,,]*tos.cor.nud1.son[3,,],tos.cor.nud2.son[1,,]*tos.cor.nud2.son[3,,])
uas.tmp2 = switch(exp,uas.cor.nud1.son[1,,]*uas.cor.nud1.son[3,,],uas.cor.nud2.son[1,,]*uas.cor.nud2.son[3,,])
vas.tmp3 = switch(exp,vas.cor.nud1.son[1,,]*vas.cor.nud1.son[3,,],vas.cor.nud2.son[1,,]*vas.cor.nud2.son[3,,])

tos.tmp3 = switch(exp,tos.cor.nud1.son[2,,]-tos.cor.ctr.son[2,,],tos.cor.nud2.son[2,,]-tos.cor.ctr.son[2,,] ) 
uas.tmp3 = switch(exp,uas.cor.nud1.son[2,,]-uas.cor.ctr.son[2,,],uas.cor.nud2.son[2,,]-uas.cor.ctr.son[2,,] ) 
vas.tmp3 = switch(exp,vas.cor.nud1.son[2,,]-vas.cor.ctr.son[2,,],vas.cor.nud2.son[2,,]-vas.cor.ctr.son[2,,] ) 

tos.tmp3[which ( is.na ( tos.nud1$mod[1, 1, 1, 1, ,] ) )] = NA
tos.tmp3[tos.tmp2 <= 0.] = NA
tos.tmp3[tos.tmp1 == F] = NA


uas.tmp3[uas.tmp2 <= 0. ] = NA
uas.tmp3[uas.tmp1 == F  ] = NA
vas.tmp3[uas.tmp2 <= 0.] = NA
vas.tmp3[uas.tmp1 == F] = NA
uas.tmp3[vas.tmp2 <= 0. ] = NA
uas.tmp3[vas.tmp1 == F  ] = NA
vas.tmp3[vas.tmp2 <= 0.] = NA
vas.tmp3[vas.tmp1 == F] = NA

# create netCDF file and put arrays
ncfname.tos.diff <- paste0("correlation_diff_tos_",expname,"_son.nc")
ncfname.uas.diff <- paste0("correlation_diff_uas_",expname,"_son.nc")
ncfname.vas.diff <- paste0("correlation_diff_vas_",expname,"_son.nc")

ncout.tos.diff <- nc_create(ncfname.tos.diff,list(tmp.def.tos))
ncout.uas.diff <- nc_create(ncfname.uas.diff,list(tmp.def.uas))
ncout.vas.diff <- nc_create(ncfname.vas.diff,list(tmp.def.vas))

# put variables
ncvar_put(ncout.tos.diff, tmp.def.tos, t(tos.tmp3))
ncvar_put(ncout.uas.diff, tmp.def.uas, t(uas.tmp3))
ncvar_put(ncout.vas.diff, tmp.def.vas, t(vas.tmp3))

# put additional attributes into dimension and data variables
ncatt_put(ncout.tos.diff,"lon","axis","X") #,verbose=FALSE) #,definemode=FALSE)
ncatt_put(ncout.tos.diff,"lat","axis","Y")

ncatt_put(ncout.uas.diff,"lon","axis","X") #,verbose=FALSE) #,definemode=FALSE)
ncatt_put(ncout.uas.diff,"lat","axis","Y")

ncatt_put(ncout.vas.diff,"lon","axis","X") #,verbose=FALSE) #,definemode=FALSE)
ncatt_put(ncout.vas.diff,"lat","axis","Y")
##

# close the file, writing data to disk
nc_close(ncout.tos.diff)
nc_close(ncout.uas.diff)
nc_close(ncout.vas.diff)

#
} # finish the function


#--- Save data 
save.cor.file(1)
save.cor.file(2)
save.cor.file(3)

save.dot.file(1)
save.dot.file(2)

