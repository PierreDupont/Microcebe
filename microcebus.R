##--------------------------------------------------------------------------------------------------------------------
## Microcebus CJS model with unequal time intervals between sessions
modelCode <- nimbleCode({
  ##---------------------
  ## DEMOGRAPHIC PROCESS 
  for(a in 1:n.age){
    for(s in 1:2){
      for(w in 1:n.seasons){
        # Age, sex & season-specific survival hazard rates
        PHI[a,s,w] ~ dunif(0,1)                         
        for(t in 1:(n.sessions-1)){
          # Session-specific survival probabilities
          phi[a,s,w,t] <- pow(PHI[a,s,w], dt[t])  # allows for unequal session intervals (dt[t])
        }# time
      }# weather
    }# sex
  }# age
      
  for(i in 1:N){
    # State at first detection 
    z[i,f[i]] ~ dbern(1)    # (= 1 because no augmentation here ; model conditional on first capture)                  
    for(t in (f[i]+1):n.sessions){
      # Subsequent states
      z[i,t] ~ dbern(z[i,t-1] * phi[age[i,t], sex[i], season[t], t-1])  
    }#t
  }#i
  
  ##------------------
  ## DETECTION PROCESS
  for(a in 1:n.age){
    for(s in 1:2){
      for(f in 1:n.fragments){
        lambda[a,s,f] ~ dgamma(0.01,0.01)                   # Detection Hazard rate
        for(t in 1:n.sessions){
          p[a,s,f,t] <- 1-exp(-lambda[a,s,f] * n.days[t]) # allows for unequal session duration (n.days[t])
        }# session
      }# fragment
    }# sex
  }# age
  
  for(i in 1:N){
    for(t in (f[i]+1):n.sessions){
      y[i,t] ~ dbern(p[age[i,t], sex[i], fragment[i], t] * z[i,t])
    }#t
  }#i
  
  ##------------------------------------------------------------------------------------------------------------------
})

##--------------------------------------------------------------------------------------------------------------------
## Microcebus JSSA model with unequal time intervals between sessions
modelCode <- nimbleCode({
  ##---------------------
  ## DEMOGRAPHIC PROCESS 
  for(w in 1:n.season){
    for(a in 1:n.age){
      for(s in 1:2){
        # Age, sex & season-specific survival hazard rates
        PHI[a,s,w] ~ dunif(0,1)                           
      }# weather
    }# sex
  }# age
  
  for(w in 1:n.season){
    for(t in 1:(n.sessions-1)){
      # season and session-specific recruitment hazard rates
      GAMMA[w,t] ~ dunif(0,1)                           
    }# weather
  }# time
  
  for(i in 1:M){
    z[i,1] ~ dbern(gamma[1])                  # probability of inclusion in year 1. 
    a[i, 1] <- (1 - z[i,1])
    for(t in 2:n.sessions){
      phi[i,t-1] <- pow(PHI[age[i,t-1],sex[i],season[t-1]], dt[t-1])
      gamma[i,t] <- 1-exp(-GAMMA[a,s,f] * n.days[t])
      z[i,t] ~ dbern((phi[i,t-1] * z[i,t-1]) + (gamma[i,t] * a[i,t-1]))
    }#t
  }#i
  
  ##------------------
  ## DETECTION PROCESS
  for(a in 1:n.age){
    for(s in 1:2){
      for(f in 1:n.fragments){
        lambda[a,s,f] ~ dgamma(0.01,0.01)                 # Detection Hazard rates
        for(t in 1:n.sessions){
          p[a,s,f,t] <- 1-exp(-lambda[a,s,f] * n.days[t]) # allows for unequal session duration (n.days[t])
        }# session
      }# fragment
    }# sex
  }# age
  
  for(i in 1:M){
    for(t in 1:n.sessions){
      y[i,t] ~ dbern(p[age[i,t], sex[i], fragment[i], t] * z[i,t])
    }#t
  }#i
  
  ##------------------------------------------------------------------------------------------------------------------
})

##--------------------------------------------------------------------------------------------------------------------
## Microcebus Multi-Site model with unequal time intervals between sessions
modelCode <- nimbleCode({
  ##---------------------
  ## DEMOGRAPHIC PROCESS 
  for(a in 1:2){
    for(s in 1:2){
      for(w in 1:2){
        # Age, sex & season-specific survival hazard rates
        PHI[a,s,w] ~ dunif(0,1)                           
        
        # Session-specific survival probabilities
        for(t in 1:(n.sessions-1)){
          phi[a,s,w,t] <- pow(PHI[a,s,w], dt[t])  # allows for unequal session intervals (dt[t])
          
          ## DO THE SAME FOR TRANSITION PROBABILITIES 
          psi12[a,s,w,t] <- pow(PSI[a,s,w], dt[t])   # allows for unequal session intervals (dt[t])
          psi13[a,s,w,t] <- pow(PSI[a,s,w], dt[t])   # allows for unequal session intervals (dt[t])
          psi14[a,s,w,t] <- pow(PSI[a,s,w], dt[t])   # allows for unequal session intervals (dt[t])
          
          psi21[a,s,w,t] <- pow(PSI[a,s,w], dt[t])   # allows for unequal session intervals (dt[t])
          psi23[a,s,w,t] <- pow(PSI[a,s,w], dt[t])   # allows for unequal session intervals (dt[t])
          psi24[a,s,w,t] <- pow(PSI[a,s,w], dt[t])   # allows for unequal session intervals (dt[t])
          
          psi31[a,s,w,t] <- pow(PSI[a,s,w], dt[t])   # allows for unequal session intervals (dt[t])
          psi32[a,s,w,t] <- pow(PSI[a,s,w], dt[t])   # allows for unequal session intervals (dt[t])
          psi34[a,s,w,t] <- pow(PSI[a,s,w], dt[t])   # allows for unequal session intervals (dt[t])
          
          psi41[a,s,w,t] <- pow(PSI[a,s,w], dt[t])   # allows for unequal session intervals (dt[t])
          psi42[a,s,w,t] <- pow(PSI[a,s,w], dt[t])   # allows for unequal session intervals (dt[t])
          psi43[a,s,w,t] <- pow(PSI[a,s,w], dt[t])   # allows for unequal session intervals (dt[t])
          
          ## WRITE TRANSITION MATRIX (OMEGA)
          OMEGA[1,1:4,a,s,w,t] <- c(,,,,)
          OMEGA[1,1:4,a,s,w,t] <- c(,,,,)
          OMEGA[1,1:4,a,s,w,t] <- c(,,,,)
          OMEGA[1,1:4,a,s,w,t] <- c(,,,,)
          
          }# time
      }# weather
    }# sex
  }# age
  
  for(i in 1:N){
    # State at first detection 
    z[i,f[i]] ~ dcat()    # (= 1 because no augmentation here ; model conditional on first capture)                  
    for(t in (f[i]+1):n.sessions){
      
      # Subsequent states
      z[i,t] ~ dcat(OMEGA[z[i,t-1], 1:n.sites, age[i,t], sex[i], w[t-1], t-1])  
    }#t
  }#i
  
  ##------------------
  ## DETECTION PROCESS
  for(a in 1:n.age){
    for(s in 1:2){
      for(f in 1:n.fragments){
        lambda[a,s,f] ~ dgamma(0.01,0.01)                   # Detection Hazard rate
        for(t in 1:n.sessions){
          p[a,s,f,t] <- 1-exp(-lambda[a,s,f] * n.days[t]) # allows for unequal session duration (n.days[t])
        }# session
      }# fragment
    }# sex
  }# age
  
  for(i in 1:N){
    for(t in (f[i]+1):n.sessions){
      y[i,t] ~ dbern(p[age[i,t], sex[i], fragment[i], t] * z[i,t])
    }#t
  }#i
  
  ##------------------------------------------------------------------------------------------------------------------
})

##--------------------------------------------------------------------------------------------------------------------