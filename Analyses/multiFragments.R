#######################################################
##### ------------ MANDENA MICROCEBE ------------ #####
##### -- PRELIMINARY ANALYSIS for FRAGMENT M13 -- #####
#######################################################
rm(list = ls())

## ------ LIBRARIES ------
library(dplyr)
library(tidyr)
library(magrittr)
library(purrr)
library(plyr)
library(coda)
library(nimble)


## ------ WORKING DIRECTORIES ------
source("workingDirectories.R")


## -----------------------------------------------------------------------------
## ------ I. LOAD & CLEAN DATA ------
## ------   1. Individual captures ------
##-- Load individual captures data 
capture_data <- read.csv(file.path(dataDir, "cmr_mhc_input_data_microcebus.csv"),
                         row.names = 1, header = TRUE)

##-- Remove duplicates from the capture data
capture_data <- capture_data[!duplicated(capture_data[ ]), ] %>% droplevels()

##-- Remove individuals without ID
capture_data <- capture_data[!is.na(capture_data$transponder), ] %>% droplevels()
capture_data <- capture_data[!(capture_data$transponder == ""), ] %>% droplevels()
capture_data <- capture_data[!(capture_data$transponder == "//"), ] %>% droplevels()
capture_data <- capture_data[!(capture_data$transponder == "/?/"), ] %>% droplevels()

##-- Format dates
capture_data$date <- as.POSIXct(strptime(capture_data$date, "%m/%d/%Y"))



## ------   2. Capture sessions ------
##-- Load capture sessions data
capture_sessions <- read.csv(file.path(dataDir, "sessions_dates_sites.csv"),h=T)
names(capture_sessions) <- c("start.date", "end.date", "site")

##-- Format dates
capture_sessions$end.date <- as.POSIXct(strptime(capture_sessions$end.date, "%m/%d/%Y"))
capture_sessions$start.date <- as.POSIXct(strptime(capture_sessions$start.date, "%m/%d/%Y"))

##-- Identify start months and years
capture_sessions$start.year <- as.numeric(format(capture_sessions$start.date,"%Y"))
capture_sessions$start.month <- as.numeric(format(capture_sessions$start.date,"%m"))

##-- Calculate the duration of each capture session
capture_sessions$duration <- difftime(time1 = capture_sessions$end.date,
                            time2 = capture_sessions$start.date,
                            units = "days") 
capture_sessions$duration <- as.numeric(capture_sessions$duration + 1) ## because at least one day of capture

##-- Aggregate capture sessions that happened in the same month
startAggSessions <- aggregate(start.date ~ start.month + start.year + site, 
                              data = capture_sessions, FUN = min)
endAggSessions <- aggregate(end.date ~ start.month + start.year + site,
                            data = capture_sessions, FUN = max)
durationAggSessions <- aggregate(duration ~ start.month + start.year + site,
                                 data = capture_sessions, FUN = sum)
aggSessions <- merge( startAggSessions, endAggSessions,
                 by = c("start.month", "start.year", "site"))
aggSessions <- merge( aggSessions, durationAggSessions,
                 by = c("start.month", "start.year", "site"))


## ------   3. Clean data for each fragment ------
fragments <- unique(capture_data$site)
minYear <- min(capture_sessions$start.year)
minMonth <- min(capture_sessions$start.month[capture_sessions$start.year == minYear])

data_list <- list()
for(f in 1:length(fragments)){
  data_list[[f]] <- list()
  ##-- Subset to fragment "s"
  print(fragments[f])
  thisData <- capture_data[capture_data$site == fragments[f], ]
  thisSessions <- capture_sessions[capture_sessions$site == fragments[f], ]
  
  ##-- Ignore fragments with less than 2 sessions 
  if(dim(thisSessions)[1] < 2){next}
  
  ##-- Ensure sessions are ordered by start.date
  thisSessions <- thisSessions[order(thisSessions$start.date), ]
  
  ##-- Give an index to each capture session
  thisSessions$index <- 1:nrow(thisSessions)
  
  ##-- Identify end months and years of each aggregated capture session
  thisSessions$end.year <- as.numeric(format(thisSessions$end.date, "%Y"))
  thisSessions$end.month <- as.numeric(format(thisSessions$end.date, "%m"))
  
  ##-- Store number of sessions for this fragment
  n.sessions <- nrow(thisSessions)
  data_list[[f]]$n.sessions <- n.sessions
  
  ##-- Calculate months index of each capture session
  for(s in 1:n.sessions){
    thisSessions$start.month.index[s] <- (thisSessions$start.year[s]-minYear)*12 +
      thisSessions$start.month[s] - minMonth  + 1
    thisSessions$end.month.index[s] <- (thisSessions$end.year[s]-minYear)*12 +
      thisSessions$end.month[s] - minMonth + 1
  }#s
  
  ##-- Calculate range of years covered by the study in this fragment
  years <- minYear:max(thisSessions$start.year)
  data_list[[f]]$years <- years
  data_list[[f]]$n.years <- length(years)
  
  ##-- Calculate range of months covered by the study in this fragment
  data_list[[f]]$months <- 1:max(thisSessions$start.month.index)
  data_list[[f]]$n.months <- length(months)-1

  
  
  ## ---------------------------------------------------------------------------
  ## ------ II. CREATE CAPTURE HISTORY ------
  ##-- Identify in which session each individual was captured
  for(c in 1:nrow(thisData)){
    temp <- thisSessions$index[thisSessions$start.date <= thisData$date[c] &
                                 thisSessions$end.date >= thisData$date[c]]
    thisData$session[c] <- ifelse(length(temp) == 0, NA, temp)
  }#c
  
  ##-- Create a dummy dataset
  dummy <- data.frame( transponder  = "dummy",
                       session = thisSessions$index)
  
  ##-- Combine real and dummy datasets
  thisData.dummy <- rbind.fill(thisData, dummy)
  
  ##-- Create the matrix of capture history
  thisCH <- table(thisData.dummy$transponder, thisData.dummy$session)
  
  ##-- Remove dummy individual
  thisCH <- thisCH[-which(dimnames(thisCH)[[1]] == "dummy"), ] 
  
  ##-- Extract the first detection session for each individual
  first <- apply(thisCH, 1, function(x)min(which(x >= 1)))
  
  ##-- Remove individuals detected for the first time on the last session
  #thisCH <- thisCH[which(first != dim(thisCH)[2]), ]
  #data_list[[f]]$f <- first[which(first != dim(thisCH)[2])]
  
  ##-- Reorder and turn into a matrix
  thisCH <- thisCH[order(dimnames(thisCH)[[1]]), ]
  thisCH <- as.matrix(thisCH)
  thisCH[thisCH > 0] <- 1
  data_list[[f]]$CH <- thisCH
  dim(thisCH)
  
  ##-- List individuals
  ids <- dimnames(thisCH)[[1]] 
  data_list[[f]]$ids <- ids
  data_list[[f]]$n.individuals <- length(ids)
  
  ##-- Sex
  sex <- unique(thisData[thisData$transponder %in% ids, c("transponder", "Sexe")])
  sex <- sex[order(sex$transponder), ] 
  all(dimnames(thisCH)[[1]] == sex$transponder)
  data_list[[f]]$sex <- ifelse(sex$Sexe == "f", 1, 2) 
  
  ##-- Age 
  # age <- matrix(data = NA, nrow = n.individuals, ncol = n.sessions)
  # for(i in ids){
  #   for(s in thisSessions$index){
  #     temp <- thisData[thisData$transponder == i & thisData$session == s, ]
  #     if(length(temp) != 0){ age[i,s] <- unique(temp$age_estimation_2) }
  #   }
  # }
  
  ##-- Interval lengths between capture sessions
  data_list[[f]]$start.int <- data_list[[f]]$end.int <- NULL
  for(i in 1:(n.sessions-1)){
    data_list[[f]]$start.int[i] <- thisSessions$start.month.index[i]
    data_list[[f]]$end.int[i] <- thisSessions$start.month.index[i+1]-1
  }#t
  
  ##-- Check the data with plots
  hist(rowSums(thisCH))                       ## num. of detections per id
  plot(thisSessions$duration, colSums(thisCH))## num. of ids detected per session duration
}#f



##-- Identify seasons for each month of the study
# Wet season from October until June!
season <- rep(c(1,1,1,1,1,1,2,2,2,1,1,1), n.years)
season <- season[minMonth:(n.months+minMonth-1)]

data_list[[6]]
## -----------------------------------------------------------------------------
## ------ III. NIMBLE MODEL ,------
##-- Write the CR model in NIMBLE
nimModel <- nimbleCode({
  ## DEMOGRAPHIC PROCESS 
  phi0[1] ~ dunif(0,1)
  phi0[2] ~ dunif(0,1)
  logit.phi0[1] <- logit(phi0[1])
  logit.phi0[2] <- logit(phi0[2])
  beta ~ dnorm(0,0.01)
  for(m in 1:n.months){
    logit(PHI[m]) <- logit.phi0[season[m]] + beta * m
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
  
  ## DETECTION PROCESS
  ## Detection Hazard rate
  lambda ~ dunif(0,5)                   
  for(t in 1:n.sessions){
    ## Allows for unequal sampling sessions (sessionDuration[t])
    p[t] <- 1-exp(-lambda * sessionDuration[t]) 
  }# session
  
  for(i in 1:n.individuals){
    for(t in (f[i]+1):n.sessions){
      y[i,t] ~ dbern(p[t] * z[i,t])
    }#t
  }#i
  
})

##-- Format the data for NIMBLE
nimData <- list( y = thisCH,
                 sessionDuration = sess13$duration)

nimConstants <- list( n.individuals = dim(thisCH)[1],
                      n.intervals = dim(thisCH)[2]-1,
                      n.sessions = dim(thisCH)[2],
                      n.months = n.months,
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
                  beta = c(0),
                  lambda = 0.5)


##-- Create a NIMBLE model object
Rmodel <- nimbleModel( code = nimModel,
                       constants = nimConstants,
                       data = nimData,
                       inits = nimInits,
                       calculate = F)
Rmodel$calculate()

##-- Configure and Build MCMC objects
conf <- configureMCMC(Rmodel, monitors = c("phi0", "beta", "lambda"), print = FALSE)
Rmcmc <- buildMCMC(conf)
23,,
##-- Compile and Run MCMC
## Finally, we compile both the model and MCMC objects and execute the compiled 
## MCMC for 50 000 iterations and 3 chains.
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
MCMC_runtime <- system.time(
  nimOutput <- runMCMC( Cmcmc,
                        niter = 100000,
                        nburnin = 50000,
                        nchains = 3,
                        thin = 5,
                        samplesAsCodaMCMC = T)
)
plot(nimOutput)
save(nimOutput, MCMC_runtime,
     file = file.path(analysisDir, "m16.RData"))

## -----------------------------------------------------------------------------
## For recruitment, take the approach of N.Hostetter
## use monthly betas, with sum(betas[1:n.months]) = 1
## then derive gammas from betas using:
## gamma[t] <- 1-sum(beta[1:t]) ... or something similar
## check in the age model!!!!