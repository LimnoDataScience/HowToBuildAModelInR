# Load libraries
library(tidyverse)
library(patchwork)
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
  doObs = read_tsv("Data/Mendota.doobs") 
  wtr = read_tsv("Data/Mendota.wtr")
  PAR = read_tsv("Data/Mendota.par") |> mutate(PAR = as.numeric(PAR)) #dataset has non-numeric characters
  
  ##Now lets plot the data quickly
  p1 = ggplot(doObs) + geom_line(aes(x = dateTime, y = DO)) + ylab('DO (mg/L)')
  p2 = ggplot(wtr) + geom_line(aes(x = dateTime, y = temp0.5)) + ylab('Water temp 0.5 m (°C)')
  p3 = ggplot(PAR) + geom_line(aes(x = dateTime, y = PAR)) + ylab('PAR')
  # patchwork plots
  print(p1 / p2 / p3 & theme_bw())
  
  # Let's join the data, so we have one dataframe
  in.data = doObs |> left_join(PAR) |> left_join(wtr)
  
  # Plot DO as a function of T, just to show it
  p4 = ggplot(in.data, aes(x = temp0.5, y = DO)) + 
    geom_point(size = 1) +
    geom_smooth(method = 'lm') + 
    xlab('Water temp 0.5 m (°C)') + ylab('DO (mg/L)') +
    theme_bw()
  print(p4) # You have to print since these are in a loop 
  
  #What data do we actually want to use in the model
  # Create vectors
  wTemp = wtr$temp0.5
  PAR = PAR$PAR
}

# Setup variables to track states and fluxes
# Initialize values for state variables
doPredic = doObs$DO[1] # Set the first value to the first observed value
phyto = 0.1 # g/m3

# Initialize values for fluxes
fatm = NA_real_ # Initialize atmospheric exchange
npp = NA_real_ # Initialize primary production
RDOC = NA_real_  # Initialize respiration
Rphyto = NA_real_  # Initialize respiration
Rtot = NA_real_  # Initialize respiration
dosat = NA_real_ # Initialize DO at saturation

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

# Setup variables to track fluxes
out.model = doObs |> slice(1:nSteps) |> 
  mutate(doPredic = doPredic,
         phyto = phyto,
         fatm = fatm, 
         npp = npp, 
         RDOC = RDOC, 
         Rphyto = Rphyto, 
         Rtot = Rtot, 
         dosat = dosat)


##################################################
# (4) Plot results
# There are many ways to do this in ggplot. Ideally you'd want to convert the
# dataframes from wide to long, but the code below uses the wide format

# Plot the predictions
p5 = ggplot(out.model) +
  geom_line(aes(x = dateTime, y = DO, color = 'black')) +
  geom_line(aes(x = dateTime, y = dosat, color = 'grey50')) +
  geom_line(aes(x = dateTime, y = doPredic, color = 'red4')) +
  geom_line(aes(x = dateTime, y = phyto, color = 'darkgreen')) +
  scale_color_identity(guide = "legend",
                       name = 'Predictions',
                       labels = c("Observed", "Phyto", "Saturation","Modeled")) +
  ylab('DO (mg/L)') +
  scale_x_datetime(date_labels = '%b', breaks = 'month') +
  theme_bw(base_size = 9) +
  theme(axis.title.x = element_blank())

# Plot the fluxes
p6 = ggplot(out.model) +
  geom_line(aes(x = dateTime, y = npp/dt, color = 'green4')) +
  geom_line(aes(x = dateTime, y = -Rtot/dt, color = 'red4')) +
  geom_line(aes(x = dateTime, y = fatm/dt, color = 'black')) +
  scale_color_identity(guide = "legend",
                       name = 'Fluxes',
                       labels = c("Fatm", "NPP", "Rtot")) +
  ylab('DO (mg/L)') +
  scale_x_datetime(date_labels = '%b', breaks = 'month') +
  theme_bw(base_size = 9) +
  theme(axis.title.x = element_blank())

# Plot the cumulative fluxes
p7 = ggplot(out.model |> slice(-1)) +
  geom_line(aes(x = dateTime, y = cumsum(npp), color = 'green4')) +
  geom_line(aes(x = dateTime, y = cumsum(-Rtot), color = 'red4')) +
  geom_line(aes(x = dateTime, y = cumsum(fatm), color = 'black')) +
  geom_line(aes(x = dateTime, y = cumsum(-Rtot + npp), color = 'darkslateblue')) +
  scale_color_identity(guide = "legend",
                       name = 'Cum. Fluxes',
                       labels = c("Fatm", "NEP", "NPP", "Rtot")) +
  ylab('DO (mg/L)') +
  scale_x_datetime(date_labels = '%b', breaks = 'month') +
  theme_bw(base_size = 9) +
  theme(axis.title.x = element_blank())

# Patchwork plots together
p5 / p6 / p7 + 
  plot_annotation(tag_levels = 'a', tag_suffix = ')') & 
  theme(plot.tag = element_text(size = 8))


