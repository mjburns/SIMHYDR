

##### SIMHYD in R
##### MJB, 2014

##### INPUT DATA
# Daily rainfall
# Daily potential evapotranspiration
# Model parameters include:

# insc - Interception store capacity (mm)
# coeff - Maximum infiltration loss (mm)
# sq - Infiltration loss exponent
# smsc - Soil moisture store capacity (mm)
# sub - Constant of proportionality in interflow equation
# crak - Constant of proportionality in groundwater recharge equation
# k - Baseflow linear recession parameter

##### DEVELOPMENT NOTES
# Model was developed based on the diagram from: "/Users/mjburns/Dropbox/UNGAUGED/useful_docs/57449_1.pdf"
# Value ranges from "hessd-8-8701-2011-print.pdf"

##### DEFINE THE FUNCTION
simhyd <- function(rainfall_input, et_input, insc, coeff, sq, smsc, sub, crak, k, store) {
  
  # The input to SIMHYD includes daily rainfall, daily et, and the model parameters	
  # Give the user warnings if any paramaters are outside safe ranges
  #if(insc < 0.5 | insc > 5) cat("insc", "is outside the range of typical values", "\n")
  #if(coeff < 50 | coeff > 400) cat("coeff", "is outside the range of typical values", "\n")
  #if(sq < 0 | sq > 6) cat("sq", "is outside the range of typical values", "\n")
  #if(smsc < 50 | smsc > 500) cat("smsc", "is outside the range of typical values", "\n")
  #if(sub < 0 | sub > 1) cat("sub", "is outside the range of typical values", "\n")
  #if(crak < 0 | crak > 1) cat("crak", "is outside the range of typical values", "\n")
  #if(k < 0.003 | k > 0.3) cat("k", "is outside the range of typical values", "\n")
  
  #Predefine the size of vector, this is based on the input data
  size_vector <- NROW(rainfall_input)
  
  #Allocate stores
  soil_moisture_store <- rep(NA, size_vector)   #Soil moisture store
  ground_store <- rep(NA, size_vector)     #Groundwater store
  total_flow <- rep(NA, size_vector) #Runoff store
  
  # Start the rainfall-runoff model
  
  for(i in 1:size_vector) {
    # Need to deal with the initial time-step
    
    if(i == 1) {
      #Make the soil moisture store initial value half of capacity (value given)
      #this is similar to the chiew spreadsheet
      soil_initial <- 0.5 * smsc
      
      #Ground initial (in looking at Chiew spreadsheet, this is the value they use)
      ground_initial <- 5
      
      #Calculate imax (maximum interception)
      imax <- min(insc, et_input[i])
      
      #Calculate int (interception)
      int <- min(imax, rainfall_input[i])
      
      #Calculate inr (rainfall not lost in the interception store)
      inr <- rainfall_input[i] - int
      
      #Calculate rmo (water infiltrated)
      rmo <- min((coeff*exp(-sq*(soil_initial/smsc))), inr)
      
      #Calculate irun (infiltration excess runoff)
      irun <- inr - rmo
      
      #Calculate srun (saturation excess runoff and interflow)
      srun <- sub * (soil_initial/smsc) * rmo
      
      #Calculate rec (water going to the groundwater store)
      rec <- crak * (soil_initial/smsc) * (rmo-srun) 
      
      #Calculate water going into moisture store
      smf <- rmo - srun - rec
      
      #Add water going to moisture store (smf) to the store
      #Need to account for overflow
      soil_moisture_store[i] <- ifelse(((soil_initial + smf) > smsc), smsc, soil_initial + smf)
      any_overflow <- ifelse(((soil_initial + smf) > smsc), ((soil_initial + smf) - smsc), 0)
      
      #Calculate pot (account for water lost due to interception)
      pot <- et_input[i] - int
      
      #Calculate et (assume that et is taken out of the store after any overflow)
      et <- min((10*(soil_initial/smsc)), pot)
      
      #Reduce the soil moisture store by ET
      soil_moisture_store[i] <- soil_moisture_store[i] - et
      
      #Add recharge and any overflow from the soil moisture store to the groundwater store
      ground_store[i] <- ground_initial + rec + any_overflow
      
      #Calculate baseflow
      base <- k * ground_store[i]
      
      #Take baseflow out from store
      ground_store[i] <- ground_store[i] - base
      
      #Add all the flow compon together
      total_flow[i] <- base + srun + irun    	
      
    }	
    
    else {	
      
      #Calculate imax (maximum interception)
      imax <- min(insc, et_input[i])
      
      #Calculate int (interception)
      int <- min(imax, rainfall_input[i])
      
      #Calculate inr (rainfall not lost in the interception store)
      inr <- rainfall_input[i] - int
      
      #Calculate rmo (water infiltrated)
      rmo <- min((coeff*exp(-sq*(soil_moisture_store[i-1]/smsc))), inr)
      
      #Calculate irun (infiltration excess runoff)
      irun <- inr - rmo
      
      #Calculate srun (saturation excess runoff and interflow)
      srun <- sub * (soil_moisture_store[i-1]/smsc) * rmo
      
      #Calculate rec (water going to the groundwater store)
      rec <- crak * (soil_moisture_store[i-1]/smsc) * (rmo-srun) 
      
      #Calculate water going into moisture store
      smf <- rmo - srun - rec
      
      #Add water going to moisture store (smf) to the store
      #Need to account for overflow
      soil_moisture_store[i] <- ifelse(((soil_moisture_store[i-1] + smf) > smsc), smsc, soil_moisture_store[i-1] + smf)
      any_overflow <- ifelse(((soil_moisture_store[i-1] + smf) > smsc), ((soil_moisture_store[i-1] + smf) - smsc), 0)
      
      #Calculate pot (account for water lost due to interception)
      pot <- et_input[i] - int
      
      #Calculate et (assume that et is taken out of the store after any overflow)
      et <- min((10*(soil_moisture_store[i]/smsc)), pot)
      
      #Reduce the soil moisture store by ET
      soil_moisture_store[i] <- soil_moisture_store[i] - et
      
      #Add recharge and any overflow from the soil moisture store to the groundwater store
      ground_store[i] <- ground_store[i-1] + rec + any_overflow
      
      #Calculate baseflow
      base <- k * ground_store[i]
      
      #Take baseflow out from store
      ground_store[i] <- ground_store[i] - base
      
      #Add all the flow compon together
      total_flow[i] <- base + srun + irun
      
      # End the non first time-step
    } 
    # now end total loop
  }
  
  #End the function now and return total flow
  #Try returning the soil moisture store
  if(store == 0) {
    return(total_flow) }
  else {
    if(store == 1) {
      return(soil_moisture_store) } else {}} 
  #End the function
}

