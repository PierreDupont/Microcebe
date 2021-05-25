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

## ------ 

## -----------------------------------------------------------------------------
## ------ I. LOAD & CLEAN DATA ------
## ------   1. Individual captures ------
##-- Load individual captures data 
capture_data <- read.csv(file.path(dataDir, "cmr_mhc_input_data_microcebus.csv"),
                         row.names = 1, header = TRUE)

##-- Remove duplicates from the capture data
capture_data <- capture_data[!duplicated(capture_data[ ]),] %>% droplevels()

##-- Subset to fragment M13
m13 <- capture_data[capture_data$site == "M16", ]

##-- Format dates
m13$date <- as.POSIXct(strptime(m13$date, "%m/%d/%Y"))


## ------   2. Capture sessions ------
##-- Load capture sessions data
capture_sessions <- read.csv(file.path(dataDir, "sessions_dates_sites.csv"),h=T)
names(capture_sessions) <- c("start.date", "end.date", "site")

##-- Subset to fragment M13
sess13 <- capture_sessions[capture_sessions$site == "M16", ]

##-- Format dates
sess13$end.date <- as.POSIXct(strptime(sess13$end.date, "%m/%d/%Y"))
sess13$start.date <- as.POSIXct(strptime(sess13$start.date, "%m/%d/%Y"))

##-- Identify start months and years
sess13$start.year <- as.numeric(format(sess13$start.date,"%Y"))
sess13$start.month <- as.numeric(format(sess13$start.date,"%m"))

##-- Calculate the duration of each capture session
sess13$duration <- difftime(time1 = sess13$end.date,
                            time2 = sess13$start.date,
                            units = "days") 
sess13$duration <- as.numeric(sess13$duration + 1) ## because at least one day of capture

##-- Aggregate capture sessions that happened in the same month
startAggSessions <- aggregate(start.date ~ start.month + start.year, data = sess13, FUN = min)
endAggSessions <- aggregate(end.date ~ start.month + start.year, data = sess13, FUN = max)
durationAggSessions <- aggregate(duration ~ start.month + start.year, data = sess13, FUN = sum)
sess13 <- merge(startAggSessions, endAggSessions, by = c("start.month", "start.year"))
sess13 <- merge(sess13, durationAggSessions, by = c("start.month", "start.year"))

##-- Ensure sessions are ordered by start.date
sess13 <- sess13[order(sess13$start.date), ]

##-- Give an index to each capture session
sess13$index <- 1:nrow(sess13)

##-- Identify end months and years of each aggregated capture session
sess13$end.year <- as.numeric(format(sess13$end.date,"%Y"))
sess13$end.month <- as.numeric(format(sess13$end.date,"%m"))
n.sessions <- nrow(sess13)

##-- Calculate months index of each capture session
minYear <- min(sess13$start.year)
minMonth <- min(sess13$start.month[sess13$start.year == minYear])
for(s in 1:n.sessions){
  sess13$start.month.index[s] <- (sess13$start.year[s]-minYear)*12 +
    sess13$start.month[s] - minMonth  + 1
  sess13$end.month.index[s] <- (sess13$end.year[s]-minYear)*12 +
    sess13$end.month[s] - minMonth + 1
}#s

##-- Calculate range of years covered by the study
years <- minYear:max(sess13$start.year)
n.years <- length(years)

##-- Calculate range of months covered by the study 
months <- 1:max(sess13$start.month.index)
n.months <- length(months)-1

##-- Identify seasons for each month of the study
# Wet season from October until June!
season <- rep(c(1,1,1,1,1,1,2,2,2,1,1,1), n.years)
season <- season[minMonth:(n.months+minMonth-1)]


## -----------------------------------------------------------------------------
## ------ II. CREATE CAPTURE HISTORY ------
##-- Identify in which session each individual was captured
for(c in 1:nrow(m13)){
  m13$session[c] <- sess13$index[sess13$start.date <= m13$date[c] & sess13$end.date >= m13$date[c]]
}#c

##-- Create a dummy dataset
dummy <- data.frame( transponder  = "dummy",
                     session = sess13$index)

##-- Combine real and dummy datasets
m13.dummy <- rbind.fill(m13, dummy)

##-- Create the matrix of capture history
ch13 <- table(m13.dummy$transponder, m13.dummy$session)

##-- Remove dummy individual
ch13 <- ch13[-which(dimnames(ch13)[[1]] == "dummy"), ] 

##-- Extract the first detection session for each individual
f <- apply(ch13, 1, function(x)min(which(x >= 1)))

##-- Remove individuals detected for the first time on the last session
ch13 <- ch13[which(f != dim(ch13)[2]), ]
f <- f[which(f != dim(ch13)[2])]

##-- Reorder and turn into a matrix
ch13 <- ch13[order(dimnames(ch13)[[1]]), ]
ch13 <- as.matrix(ch13)
ch13[ch13 > 0] <- 1
dim(ch13)

##-- List individuals
ids <- dimnames(ch13)[[1]] 
n.individuals <- length(ids)

##-- Sex
sex <- unique(m13[m13$transponder %in% ids, c("transponder", "Sexe")])
sex <- sex[order(sex$transponder), ] 
all(dimnames(ch13)[[1]] == sex$transponder)
sex <- ifelse(sex$Sexe == "f", 1, 2) 


# ##-- Age 
# age <- matrix(data = NA, nrow = n.individuals, ncol = n.sessions)
# for(i in ids){
#   for(s in sess13$index){
#     temp <- m13[m13$transponder == i & m13$session == s, ]
#     if(length(temp) != 0){
#       age[i,s] <- unique(temp$age_estimation_2)
#     }
#   }
# }

##-- Interval lengths between capture sessions
start.int <- end.int <- NULL
for(i in 1:(n.sessions-1)){
  start.int[i] <- sess13$start.month.index[i]
  end.int[i] <- sess13$start.month.index[i+1]-1
}#t


##-- Check the data with plots
hist(rowSums(ch13))                  ## num. of detections per individual
plot(sess13$duration, colSums(ch13)) ## num. of ids detected per session duration



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
nimData <- list( y = ch13,
                 sessionDuration = sess13$duration)

nimConstants <- list( n.individuals = dim(ch13)[1],
                      n.intervals = dim(ch13)[2]-1,
                      n.sessions = dim(ch13)[2],
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