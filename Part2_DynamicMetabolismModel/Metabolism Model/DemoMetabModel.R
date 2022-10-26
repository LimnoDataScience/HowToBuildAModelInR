## GSTP Metabolism Model 2013
## Edited October 2022 by PC Hanson
## Model simulation of epilimnetic lake metabolism for Lake Mendota
## Observational data from NTL LTER

# O2 Model: dO2/dt = NPP - Rtot + Fatm
# O2(t) = (O2(t-1) + NPP(t-1) - Rtot(t-1) + Fatm(t-1)) * dt

# Phytoplankton Model: dPhyto/dt = NPP - RPhyto - Settling
# Phyto(t) = (Phyto(t-1) + NPP(t-1) - RPhyto(t-1) - Settling(t-1)) * dt

# NPP = f(PAR,phosphorus)
# Rtot = (RDOC + RPhyto) * dt; RDOC = f(DOC); RPhyto = f(Phyto)
# Settling = f(zMix,Phyto)

##################################################
# (1) Setup constants and parameters

# These three constants you might normally get from observational data
phosphorus = 0.050 # Constant value of phosphorus (ug/L)
DOC = 5.0  # Constant value for dissolved organic carbon (mg/L)
zMix = 10  # Constant depth of the lake's upper mixed layer (m)

# Model parameters
# Fatm = k * (DO deviation from saturation) / mixed layer
kDO = 1.0   # Parameter, gas exchange piston velocity (meters/day), often calculated from, e.g., wind speed

# NPP = PAR * P * pNPP * Theta^Temperature
pNPP = 0.08  # Parameter, converts light and phosphorus to NPP (mgC/unitP/unitLight)
thetaNPP = 1.2 # Parameter, Arrhenius coefficient for temperature adjustment for NPP

# RPhyto = Phyto * phytoR * Theta^Temperature
phytoR = 0.8   # Parameter, 1st order decay of phytoplankton
thetaR = 1.08  # Parameter, Arrhenius coefficient for temperature adjustment for respiration

# SettlingPhyto = Phyto * settlingPhyto / zMix
settlingPhyto = 1.0 # Parmaeter, Settling in m/d

# RDOC = DOC * docR
docR = 0.02  # Parameter, 1st order decay of DOC
 
# Other necessary constants for the model
CtoO2 = 32/12 # Convert C values to DO values
dt = 1/48     # Model time step (day fraction); 1/48 = 30 minutes; for now, has to be synchronized with 30 min met data
nSteps = 10500   # Known apriori (number of records in the series)
LoadData = TRUE 

##################################################
# (2) Pre-simulation setup

# If you're using Rstudio, change to current directory
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# Only load the data if it hasn't already been loaded, otherwise the code runs slowly
if (exists('doObs')){
  LoadData = FALSE
}

#Load the data
if (LoadData){
  doObs = read.table("Data/Mendota.doobs", #This is the file name
                     header=T, #This says that there *is* a header in the file
                     sep="\t", #This says the column separator is a tab
                     colClasses=c(dateTime="POSIXct")) #This says the column "dateTime" is in a POSIXct format, which is date/time in R
  wtr = read.table("Data/Mendota.wtr",header=T,sep="\t",colClasses=c(dateTime="POSIXct"))
  PAR = read.table("Data/Mendota.par",header=T,sep="\t",colClasses=c(dateTime="POSIXct"))
  # Plot DO as a function of temperature
  # plot(wTemp,doObs$DO,ylab='DO (mg/L)',xlab='Water temp (c)',main="",cex=1.3)
  ##Now lets plot the data quick
  par(mfrow=c(3,1),lend=2,mai = c(0.35,0.8, 0.08, 0.05),oma = c(2,1,0.2,0.2), cex = 0.9)
  plot(doObs$dateTime,doObs[,2],type="l",ylab='DO (mg/L)')
  # title("doObs")
  plot(wtr$dateTime,wtr[,3],type="l",ylab='Water temp (C)')
  # title("Water Temp")
  plot(PAR$dateTime,PAR[,2],type="l",ylab='PAR')
  # title("PAR")
  # Plot DO as a function of T, just to show it
  par(mfrow=c(1,1),lend=2,mai = c(0.8,0.8, 0.08, 0.05),oma = c(0.5,0.5,0.2,0.2), cex = 0.9)
  myLinearModel = lm(doObs[,2] ~ wtr[,3])
  plot(wtr[,3],doObs[,2], xlab = 'Water T (C)',ylab = 'DO (mg/L)',pch='.')
  abline(a=myLinearModel$coefficients[1],b=myLinearModel$coefficients[2],lty=2,lwd=3)
  
  # Simplify the vectors
  wTemp = wtr[,3]
  PAR = PAR[,2]
}

# Setup variables to track states and fluxes
# State variables
doPredic = vector(mode="double",length=nSteps)
phyto = vector(mode="double",length=nSteps)
# Fluxes
fatm = vector(mode="double",length=nSteps)
npp = vector(mode="double",length=nSteps)
RDOC = vector(mode="double",length=nSteps)
Rphyto = vector(mode="double",length=nSteps)
Rtot = vector(mode="double",length=nSteps)
# DO at saturation
dosat = vector(mode="double",length=nSteps)

# Initialize values for state variables and fluxes
doPredic[1] = doObs[1,2] # Set the first value to the first observed value
phyto[1] = 0.1 # g/m3
fatm[1] = NA # Initialize atmospheric exchange
npp[1] = NA # Initialize primary production
RDOC[1] = NA  # Initialize respiration
Rphyto[1] = NA  # Initialize respiration
Rtot[1] = NA  # Initialize respiration
dosat[1] = NA # Initialize DO at saturation

##################################################
# (3) Run the model

# Cycle through the model with 48 time steps per day (i.e., every 30 min)
for(i in 2:nSteps){
  # if water temp is missing, use previous value
  if(is.na(wTemp[i])){
    wTemp[i] = wTemp[i-1]
  }
  # if PAR is missing, use previous value
  if(is.na(PAR[i])){
    PAR[i] = PAR[i-1]
  }
  
  # Calculate the saturation value of dissolved oxygen as a function of temperature 
  dosat[i] = -0.00006 * wTemp[i]^3 + 0.0069 * wTemp[i]^2 - 0.3906 * wTemp[i] + 14.578
  
  # Atmospheric exchange
  # Fatm = k * (DO deviation from saturation) / mixed layer
  fatm[i] = dt * kDO * (dosat[i] - doPredic[i-1]) / zMix # units: m/d * g/m3 * 1/m should be mg/L/30min 
  
  # Primary production
  # NPP = PAR * P * pNPP * Theta^Temperature
  npp[i] = dt * PAR[i] * phosphorus * pNPP * thetaNPP^(wTemp[i] - 20)
  
  # Respiration terms
  # RDOC = [DOC] * docR
  RDOC[i] = dt * DOC * docR
  # RPhyto = Phyto * phytoR * Theta^Temperature
  Rphyto[i] = dt * phyto[i-1] * phytoR * thetaR^(wTemp[i] - 20)
  Rtot[i] = RDOC[i] + Rphyto[i]
  
  # SettlingPhyto = Phyto * settlingPhyto/zMix
  Settling = dt * phyto[i-1] * settlingPhyto/zMix
  
  # Mass balance equation for dissolved oxygen
  doPredic[i] = doPredic[i-1] + fatm[i] + npp[i]*CtoO2 - Rtot[i]*CtoO2
  # Mass balance equaation for phytos
  phyto[i] = phyto[i-1] + npp[i] - Rphyto[i] - Settling
}

##################################################
# (4) Plot results

par(mfrow=c(3,1),lend=2,mai = c(0.35,0.8, 0.08, 0.05),oma = c(2,1,0.2,0.2), cex = 0.9)

# Plot the predictions
plot(doObs[1:nSteps,1],doObs[1:nSteps,2],ylim=c(0,25),xlab="",ylab="DO (mg/L)", type="l")
lines(doObs[1:nSteps,1],dosat,col='grey',lwd=2)
lines(doObs[1:nSteps,1],doPredic,col="red")
lines(doObs[1:nSteps,1],phyto,col='green')
legend('topright',legend=c('Observed','Modeled','Saturation','Phyto (mgC/L)'),lty=c(1,1),col=c('black','red','grey','green'))

# Plot the fluxes
myYLim = c(min(c(fatm/dt,-Rtot/dt),na.rm=TRUE),max(npp/dt,na.rm=TRUE))
plot(doObs[1:nSteps,1],npp/dt,ylim=myYLim,type='l',col='green',xlab="",ylab="Rates (mgC/L/d)",main="")
lines(doObs[1:nSteps,1],-Rtot/dt,type='l',col='red',lwd=2)
lines(doObs[1:nSteps,1],fatm/dt,type='l',col='black')
abline(h=0,lty=2)
legend('topright',legend=c('NPP','Rtot','Fatm'),lty=c(1,1,1),col=c('green','red','black'))

# Plot the cumulative fluxes
# Calculate the y limits
CumSumMax = max(c(abs(cumsum(npp[2:nSteps])), abs(cumsum(Rtot[2:nSteps])), abs(cumsum(fatm[2:nSteps]))))
myYLim = c(-CumSumMax,CumSumMax)
plot(doObs[2:nSteps,1],cumsum(npp[2:nSteps]),ylim=myYLim,type='l',col='green',xlab="",ylab="Cumulative rates (mgC)",main="")
lines(doObs[2:nSteps,1],cumsum(-Rtot[2:nSteps]),type='l',col='red',lwd=2)
lines(doObs[2:nSteps,1],cumsum(fatm[2:nSteps]),type='l',col='black',lwd=2)
lines(doObs[2:nSteps,1],cumsum(-Rtot[2:nSteps]+npp[2:nSteps]),type='l',col='blue',lwd=2)
abline(h=0,lty=2)
legend('topleft',legend=c('NPP','Rtot','Fatm','NEP'),lty=c(1,1,1,1),col=c('green','red','black','blue'))



