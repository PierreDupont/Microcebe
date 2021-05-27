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

##-- Clean-up fragments names
unique(capture_data$site)
capture_data$site[capture_data$site == "M13R"] <- "M13"
capture_data$site[capture_data$site == "M15A"] <- "M15a"
capture_data$site[capture_data$site == "M15D"] <- "M15d"
capture_data$site[capture_data$site == "M15E"] <- "M15e"
capture_data$site[capture_data$site == "M16 Vao"] <- "M16"
capture_data$site[capture_data$site == "Bush clearing"] <- "Bush"
capture_data$site[capture_data$site == "Bush cl"] <- "Bush"

fragments <- c("M13","M15a","M15b","M15c","M15d","M15e","M16")
# fragments <- c("M3","M4","M5","M6","M7","M9","M11","M12","M13",
#                "M15a","M15b","M15c","M15d","M15e","M16","M16_Bakuli","M20",
#                "Bush clearing","Bush cl","Cage","Mine","Corr  R.","Corr  E.",
#                "Tokon'ala")
capture_data <- capture_data[capture_data$site %in% fragments, ] 



## ------   2. Capture sessions ------
##-- Load capture sessions data
capture_sessions <- read.csv(file.path(dataDir, "sessions_dates_sites.csv"),h=T)
names(capture_sessions) <- c("start.date", "end.date", "site")

##-- Clean-up fragments names
unique(capture_sessions$site)
# [1] "M16"                "M3"                 "M15a"               "M4"                
# [5] "M5"                 "M20"                "M16_Bakuli"         "Corridor"          
# [9] "M6"                 "M13"                "M1"                 "M15b"              
# [13] "Bush clearing"      "Cage"               "Mine"               "Corr"              
# [17] "Corr R."            "M13R"               "Corr R. et M13R"    "M13R et M15E"      
# [21] "Tokon'ala"          "Corr E."            "Corr M16"           "M16 Vao"           
# [25] "M15c"               "M15d"               "M7"                 "M12"               
# [29] "M11"                "M9"                 "Herbarium"          "Salle des graines "

# capture_sessions$site[capture_sessions$site == "M13R et M15E" ] <- "???"

capture_sessions[capture_sessions$site %in% fragments, ]

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
minYear <- min(capture_sessions$start.year)
minMonth <- min(capture_sessions$start.month[capture_sessions$start.year == minYear])

data_list <- list()
for(f in 1:length(fragments)){
  data_list[[f]] <- list()
  ##-- Subset to fragment "s"
  print(fragments[f])
  thisData <- capture_data[capture_data$site == fragments[f], ]
  thisSessions <- capture_sessions[capture_sessions$site == fragments[f], ]
  
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
  years <- min(thisSessions$start.year):max(thisSessions$start.year)
  data_list[[f]]$years <- years
  data_list[[f]]$n.years <- length(years)
  
  ##-- Calculate range of months covered by the study in this fragment
  months <- min(thisSessions$start.month.index):max(thisSessions$start.month.index)
  data_list[[f]]$months <- months
  data_list[[f]]$n.months <- length(months)-1

  
  
  ## ---------------------------------------------------------------------------
  ## ------ II. CREATE CAPTURE HISTORY ------
  ##-- Identify in which session each individual was captured
  for(c in 1:nrow(thisData)){
    temp <- thisSessions$index[thisSessions$start.date <= thisData$date[c] &
                                 thisSessions$end.date >= thisData$date[c]]
    ## [PD]: WARNING!!!
    ## SOME CAPTURES WERE RECORDED OUTSIDE THE OFFICIAL CAPTURE SESSIONS
    if(length(temp) == 0)stop()
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
  thisCH <- thisCH[which(first != dim(thisCH)[2]), ]
  data_list[[f]]$first <- first[which(first != dim(thisCH)[2])]
  
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
  sex <- unique(thisData[thisData$transponder %in% ids, c("transponder","Sexe")])
  sex <- sex[order(sex$transponder), ] 
  all(dimnames(thisCH)[[1]] == sex$transponder)
  data_list[[f]]$sex <- ifelse(sex$Sexe == "f",1,2) 
  
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
n.months <- max(unlist(lapply(data_list, function(x)length(x$months))))
season <- rep(c(1,1,1,1,1,1,2,2,2,1,1,1), length(years))
season <- season[minMonth:(n.months + minMonth - 1)]



## -----------------------------------------------------------------------------
## ------ II. NIMBLE MODEL ------
##-- Write the CR model in NIMBLE
nimModel <- nimbleCode({

  ## DEMOGRAPHIC PROCESS 
  phi0[1,1,1] ~ dunif(0,1)  ## Baseline survival female/protected/wet
  phi0[1,1,2] ~ dunif(0,1)  ## Baseline survival female/protected/dry
  phi0[1,2,1] ~ dunif(0,1)  ## Baseline survival female/disturbed/wet
  phi0[1,2,2] ~ dunif(0,1)  ## Baseline survival female/disturbed/dry
  phi0[2,1,1] ~ dunif(0,1)  ## Baseline survival male/protected/wet
  phi0[2,1,2] ~ dunif(0,1)  ## Baseline survival male/protected/dry
  phi0[2,2,1] ~ dunif(0,1)  ## Baseline survival male/disturbed/wet
  phi0[2,2,2] ~ dunif(0,1)  ## Baseline survival male/disturbed/dry
  
  beta[1,1] ~ dnorm(0,0.01) ## Temporal effect female/protected
  beta[1,2] ~ dnorm(0,0.01) ## Temporal effect female/disturbed
  beta[2,1] ~ dnorm(0,0.01) ## Temporal effect male/protected
  beta[2,2] ~ dnorm(0,0.01) ## Temporal effect male/disturbed
  
  for(m in 1:n.months){
    logit(PHI[1,1,m]) <- logit(phi0[1,1,season[m]]) + beta[1,1] * m
    logit(PHI[1,2,m]) <- logit(phi0[1,2,season[m]]) + beta[1,2] * m
    logit(PHI[2,1,m]) <- logit(phi0[2,1,season[m]]) + beta[2,1] * m
    logit(PHI[2,2,m]) <- logit(phi0[2,2,season[m]]) + beta[2,2] * m
  }# months
  
  ## Multi-sites model
  for(f in 1:n.fragments){
    for(t in 1:n.intervals[f]){
      phi[1,t,f] <- prod(PHI[1,status[f],start.int[t,f]:end.int[t,f]])
      phi[2,t,f] <- prod(PHI[2,status[f],start.int[t,f]:end.int[t,f]])
    }# time
    
    for(i in 1:n.individuals[f]){
      z[i,first[i,f],f] ~ dbern(1)  
      for(t in first[i,f]:n.intervals[f]){
        z[i,t+1,f] ~ dbern(z[i,t,f] * phi[sex[i],t,f])  
      }#t
    }#i
  }#f
  
  
  ## DETECTION PROCESS
  lambda[1,1] ~ dunif(0,5) ## Detection Hazard rate female/protected
  lambda[1,2] ~ dunif(0,5) ## Detection Hazard rate female/disturbed
  lambda[2,1] ~ dunif(0,5) ## Detection Hazard rate male/protected
  lambda[2,2] ~ dunif(0,5) ## Detection Hazard rate male/disturbed
  
  ## Multi-sites model
  for(f in 1:n.fragments){
    for(t in 1:n.sessions[f]){
      ## Allows for unequal sampling sessions (sessionDuration[t])
      p[1,t,f] <- 1-exp(-lambda[1,status[f]] * sessionDuration[t,f]) 
      p[2,t,f] <- 1-exp(-lambda[2,status[f]] * sessionDuration[t,f]) 
    }# session
    
    for(i in 1:n.individuals[f]){
      for(t in (first[i]+1):n.sessions[f]){
        y[i,t,f] ~ dbern(p[sex[i],t,f] * z[i,t,f])
      }#t
    }#i
  }#f
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

