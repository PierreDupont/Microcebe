## This script simulates a Microcebus-like capture-recapture dataset
## and analyze it using a custom capture-recapture model in NIMBLE 
rm(list=ls())

## Load Libraries 
library(nimble)
library(basicMCMCplots)
library(coda)


## Load personal working directories
source("workingDirectories.R")


WD <- file.path(simDir, "simulation1")

## -----------------------------------------------------------------------------
## ----- I. Simulate Data -----
## Here, we simulate a population of Microcebus
## based (approx.) on the life-history of Mandena's Microcebes. 
N0 <- 50                               ## Initial population size       
n.years <- 6                           ## Number of years simulated
season <- rep(c(1,1,1,1,1,1,1,1,2,2,2,2), n.years) ## Wet/dry seasons
phi0 <- c(0.98,0.8)                    ## Time- and season-specific survival
beta <- c(-0.02,-0.0)
r <- 0.08                              ## Reproduction probability
litterSize <- c(0.5,0.35,0.15)         ## Vector of litter size probabilities
meanDuration <- 4
lambdaDet <- 0.45
propSession <- 0.46

for(rep in 1:30){
## -----  1. Simulate population dynamics -----
## We start by initializing simulation objects.
POP <- list()                          ## List of population composition 
POP[[1]] <- rep(1, N0)                 ## Initial population composition
N <- vector()                          ## List of population size
N[1] <- N0                             ## Initial population size

## Then, we simulate individual survival, recruitment, and litter size.
## From that we derive session-specific population sizes.
n.months <- 12*n.years
phi <- NULL 
start.int <- end.int <- NULL 
for(t in 2:n.months){
  phi[t-1] <- 1/(1+exp(-(logit(phi0[season[t]]) + beta[season[t]] * t)))
  
  Z <- B <- rep(0,length(POP[[t-1]]))                ## Initialize vectors of individuals alive, reproduction indicator and number of newborns         
  for(i in 1:length(POP[[t-1]])){       
    if(POP[[t-1]][i] == 1){                          ## Sample states only if individual is alive
      Z[i] <- rbinom(1, 1, phi[t-1])                 ## Sample individual survival
      repro <- rbinom(1, 1, r*Z[i])                  ## Sample individual repro (only if individual survived)
      if(repro==1){B[i] <- sample( x = 1:3,
                                   size = 1,
                                   prob = litterSize)}## Sample individual litter size
    }#if
  }#i
  if(sum(B[]) > 0){Z <- c(Z, rep(1,sum(B[])))}       ## Add new recruits to the population vector
  N[t] <- sum(Z)                                     ## New population size
  POP[[t]] <- Z                                      ## New population composition
  if(length(POP[[t]]) > 200000)stop("pop size waayyyyyy too big...something is wrong!!!") 
}#t

## Plot of the population trajectory
plot(1:n.months, N, type = "l", axes = F,
     ylim = c(0, max(N)), lwd = 2,
     ylab = "Population size", xlab = "years")
axis(side = 1, at = (0:n.years)*12, labels = 0:n.years, xlab = "years")
axis(side = 2)


## Finally, we re-arrange individual life-histories in a (BIG) matrix.
z <- matrix(0, length(POP[[n.months]]), n.months)
for (t in 1:n.months){
  z[1:length(POP[[t]]), t] <- POP[[t]]
}#t



## -----  2. Simulate population monitoring -----
### Generate CMR dataset 
## Here, we randomly sample in which months the trapping sessions occured.
## (First and last sessions must be capture sessions for the model to make sense)
sessionIndex <- c( 1,
                   sample( x = 2:(n.months-1),
                           size = round(propSession*(n.months-2))),
                   n.months)
sessionIndex <- sessionIndex[order(sessionIndex)]
n.sessions <- length(sessionIndex)

## Then, we randomly sample the length of each capture session 
## and calculate the corresponding detection probability
sessionDuration <- rpois(n = length(sessionIndex), lambda = meanDuration)+1
p <- 1 - exp(-lambdaDet * sessionDuration)

## Finally, we sample individual observations only for capture sessions
y <- matrix(0, dim(z)[1], n.sessions)
for(i in 1:dim(z)[1]){
  for (t in 1:n.sessions){
    y[i,t] <- rbinom(size = 1, n = 1, prob = p[t] * z[i,sessionIndex[t]])
  }#t
}#i

## Keep only detected individuals
detected <- apply(y, 1, function(x)any(x == 1))
y.detected <- y[detected, ]

## Extract the first detection occasion per individual
f <- apply(y.detected, 1, function(x)min(which(x == 1)))

## Remove individuals detected for the first time on the last session
y.detected <- y.detected[which(f != dim(y.detected)[2]), ]
f <- f[which(f != dim(y.detected)[2])]



## -----  3. Extract monitoring information -----
start.int <- end.int <- dt1 <- dt2 <- NULL
for(t in 1:(length(sessionIndex)-1)){
  start.int[t] <- sessionIndex[t]
  end.int[t] <- sessionIndex[t+1]-1
  dt1[t] <- sum(season[start.int[t]:end.int[t]] == 1)
  dt2[t] <- sum(season[start.int[t]:end.int[t]] == 2)
}#t

## -----------------------------------------------------------------------------
## ----- II. Fit the model with NIMBLE -----
## -----  1. NIMBLE model code -----
modelCode <- nimbleCode({
  ##---------------------
  ## DEMOGRAPHIC PROCESS 
  phi0[1] <- ilogit(logit.phi0[1]) #~ dunif(0,1)
  phi0[2] <- ilogit(logit.phi0[2]) #~ dunif(0,1)
  logit.phi0[1] ~ dnorm(0,0.01)    #<- logit(phi0[1])
  logit.phi0[2] ~ dnorm(0,0.01)    #<- logit(phi0[2])
  beta[1] ~ dnorm(0,0.01)
  beta[2] ~ dnorm(0,0.01)
  for(m in 1:n.months){
    # PHI[m] <- phi0[season[m]] + beta[season[m]] * m
    logit(PHI[m]) <- logit.phi0[season[m]] + beta[season[m]] * m
  }# months
  
  for(t in 1:n.intervals){
    phi[t] <- prod(PHI[start.int[t]:end.int[t]])
  }# time
  
  for(i in 1:n.individuals){
    z[i,f[i]] ~ dbern(1)  
    for(t in f[i]:n.intervals){
      z[i,t+1] ~ dbern(z[i,t] * phi[t])  
    }#t
  }#i
  
  ##------------------
  ## DETECTION PROCESS
  ## Detection Hazard rate
  lambda ~ dgamma(0.01,0.01)                   
  for(t in 1:n.sessions){
    ## Allows for unequal sampling sessions (sessionDuration[t])
    p[t] <- 1-exp(-lambda * sessionDuration[t]) 
  }# session
  
  for(i in 1:n.individuals){
    for(t in (f[i]+1):n.sessions){
      y[i,t] ~ dbern(p[t] * z[i,t])
    }#t
  }#i
  
  ##----------------------------------------------------------------------------
})


getPhi <- nimbleFunction(
  run = function(phiVec = double(1)){
    returnType(double(0))
    lengthPhi <- length(PhiVec)
    out <- 1
    for(i in 1:lengthPhi){
      out <- out*phiVec[i]
    }
    return(out) 
  })

## -----  2. Format the data for NIMBLE -----
# Then, we organize the simulated data in a format useable by NIMBLE
# i.e. we create objects `constants`, `data`, and `inits` for later 
# use in the function `nimbleModel`. 
nimData <- list( y = y.detected,
                 sessionDuration = sessionDuration )

nimConstants <- list( n.individuals = dim(y.detected)[1],
                      n.months = n.months,
                      n.intervals = dim(y.detected)[2]-1,
                      n.sessions = dim(y.detected)[2],
                      season = season,
                      f = f,
                      start.int = start.int,
                      end.int = end.int)


l <- apply(nimData$y, 1, function(x)max(which(x == 1)))
z.init <- matrix(data =  NA,
                 nrow = nrow(nimData$y),
                 ncol = ncol(nimData$y))
for(i in 1:dim(z.init)[1]){
  z.init[i,f[i]:l[i]] <- 1
  if(l[i]<dim(z.init)[2]){
    z.init[i,(l[i]+1):(dim(z.init)[2])] <- 0
  }
}#i
nimInits <- list( z = z.init,
                  phi0 = c(0.85,0.85),
                  beta = c(0,0),
                  lambda = 0.5)


## -----  3. Create a NIMBLE model object -----
## Now, we can create the NIMBLE model object, using the model structure
## defined in modelCode, and the constants, data, and initial values.
Rmodel <- nimbleModel( code = modelCode,
                       constants = nimConstants,
                       data = nimData,
                       inits = nimInits,
                       calculate = F)
Rmodel$calculate()

### Configure and Build MCMC objects 
## We configure an MCMC algorithm to the `Rmodel` model object.
## We assign MCMC monitors to $/phi$, $/gamma$, $lambda$, and $/psi$.
conf <- configureMCMC(Rmodel, monitors = c("phi0", "beta", "lambda"), print = FALSE)
Rmcmc <- buildMCMC(conf)

### Compile and Run MCMC 
## Finally, we compile both the model and MCMC objects and
## execute the compiled MCMC for 10 000 iterations.
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
MCMC_runtime <- system.time(
  samples <- runMCMC( Cmcmc,
                      niter = 4000,
                      nburnin = 0,
                      nchains = 2,
                      samplesAsCodaMCMC = T))

plot(samples)

save(samples, file = file.path(simDir, "simulation1",
                               paste("sim_",rep,".RData",sep="")))

}#rep


## -----------------------------------------------------------------------------
## ----- III. Plot Relative Bias -----
parm_true <- c(beta, lambdaDet, phi0)
plot(1, type = "n", axes = F, ylab = "", xlab = "", xlim = c(0.5,5.5), ylim = c(-2,2))
axis(side = 1, at = 1:5, labels = c("beta[1]","beta[2]","lambda","phi0[1]","phi0[2]"))
axis(side = 2, at = seq(-2,2,1), labels = seq(-2,2,1))
abline(h = 0, lty = 2)
RB <- matrix(NA, nrow = 30, ncol = 5)
for(rep in 1:30){
  if(! rep %in% c(24,25)){
    load(file.path(WD,paste("sim_", rep, ".RData", sep = "")))
    temp <- as.mcmc.list(lapply(samples, function(x)as.mcmc(x[2000:4000, ])))
    #plot(temp)
    RB[rep, ] <- ((summary(temp)[[1]][ ,1]+0.01) - (parm_true + 0.01))/(parm_true+0.01)
    points(x = 1:5, y = RB[rep, ], pch = 19, col = "navyblue")
  }
}#rep



## -----------------------------------------------------------------------------