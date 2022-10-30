library(hydroGOF) #for RMSE function

metabFunction <- function(pars, doObs, wTemp, PAR) {
  
  phytoR = pars[1]
  kDO = pars[2]
  # (1) Setup constants and parameters
  # These three constants you might normally get from observational data
  phosphorus = 0.050 # Constant value of phosphorus (ug/L)
  DOC = 5.0  # Constant value for dissolved organic carbon (mg/L)
  zMix = 10  # Constant depth of the lake's upper mixed layer (m)
  
  # Model parameters
  # Fatm = k * (DO deviation from saturation) / mixed layer
  # kDO = 1.0   # Parameter, gas exchange piston velocity (meters/day), often calculated from, e.g., wind speed
  
  # NPP = PAR * P * pNPP * Theta^Temperature
  pNPP = 0.08  # Parameter, converts light and phosphorus to NPP (mgC/unitP/unitLight)
  thetaNPP = 1.2 # Parameter, Arrhenius coefficient for temperature adjustment for NPP
  
  # RPhyto = Phyto * phytoR * Theta^Temperature
  # phytoR = 0.8   # Parameter, 1st order decay of phytoplankton
  thetaR = 1.08  # Parameter, Arrhenius coefficient for temperature adjustment for respiration
  
  # SettlingPhyto = Phyto * settlingPhyto / zMix
  settlingPhyto = 1.0 # Parmaeter, Settling in m/d
  
  # RDOC = DOC * docR
  docR = 0.02  # Parameter, 1st order decay of DOC
  
  # Other necessary constants for the model
  CtoO2 = 32/12 # Convert C values to DO values
  dt = 1/48     # Model time step (day fraction); 1/48 = 30 minutes; for now, has to be synchronized with 30 min met data
  nSteps = 10500
  
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
  
  # Make dataframe of output vectors 
  out.model = doObs %>% slice(1:nSteps) %>% 
    mutate(doPredic = doPredic,
           phyto = phyto,
           fatm = fatm, 
           npp = npp, 
           RDOC = RDOC, 
           Rphyto = Rphyto, 
           Rtot = Rtot, 
           dosat = dosat)
  
  out.rmse = rmse(sim = out.model$doPredic, obs = out.model$DO)
  return(out.rmse)
}

metabFunction(pars = c(0.9,1), doObs = doObs, 
              wTemp = wTemp, PAR = PAR)

optim(par = c(0.9,1), fn = metabFunction, doObs = doObs, 
      wTemp = wTemp, PAR = PAR, method = 'BFGS')$par

