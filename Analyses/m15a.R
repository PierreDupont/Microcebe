#######################################################
##### ------------ MANDENA MICROCEBE ------------ #####
##### -- PRELIMINARY ANALYSIS for FRAGMENT M15a -- #####
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
capture_data <- capture_data[!duplicated(capture_data[ ]),] %>% droplevels()

##-- Subset to fragment M15a
m15a <- capture_data[capture_data$site == "M15a", ]

##-- Remove non marked individuals and problem individuals with no session assigned)
m15a <- m15a[!(m15a$transponder %in% c("//", "/16A90A2/", "/16AABA8/")),]

##-- Format dates
m15a$date <- as.POSIXct(strptime(m15a$date, "%m/%d/%Y"))


## ------   2. Capture sessions ------
##-- Load capture sessions data
capture_sessions <- read.csv(file.path(dataDir, "sessions_dates_sites.csv"),h=T)
names(capture_sessions) <- c("start.date", "end.date", "site")

##-- Subset to fragment M15a
sess15a <- capture_sessions[capture_sessions$site == "M15a", ]

##-- Format dates
sess15a$end.date <- as.POSIXct(strptime(sess15a$end.date, "%m/%d/%Y"))
sess15a$start.date <- as.POSIXct(strptime(sess15a$start.date, "%m/%d/%Y"))

##-- Identify start months and years
sess15a$start.year <- as.numeric(format(sess15a$start.date,"%Y"))
sess15a$start.month <- as.numeric(format(sess15a$start.date,"%m"))

##-- Calculate the duration of each capture session
sess15a$duration <- difftime(time1 = sess15a$end.date,
                            time2 = sess15a$start.date,
                            units = "days") 
sess15a$duration <- as.numeric(sess15a$duration + 1) ## because at least one day of capture

##-- Aggregate capture sessions that happened in the same month
startAggSessions <- aggregate(start.date ~ start.month + start.year, data = sess15a, FUN = min)
endAggSessions <- aggregate(end.date ~ start.month + start.year, data = sess15a, FUN = max)
durationAggSessions <- aggregate(duration ~ start.month + start.year, data = sess15a, FUN = sum)
sess15a <- merge(startAggSessions, endAggSessions, by = c("start.month", "start.year"))
sess15a <- merge(sess15a, durationAggSessions, by = c("start.month", "start.year"))

##-- Ensure sessions are ordered by start.date
sess15a <- sess15a[order(sess15a$start.date), ]

##-- Give an index to each capture session
sess15a$index <- 1:nrow(sess15a)

##-- Identify end months and years of each aggregated capture session
sess15a$end.year <- as.numeric(format(sess15a$end.date,"%Y"))
sess15a$end.month <- as.numeric(format(sess15a$end.date,"%m"))
n.sessions <- nrow(sess15a)

##-- Calculate months index of each capture session
minYear <- min(sess15a$start.year)
minMonth <- min(sess15a$start.month[sess15a$start.year == minYear])
for(s in 1:n.sessions){
  sess15a$start.month.index[s] <- (sess15a$start.year[s]-minYear)*12 +
    sess15a$start.month[s] - minMonth  + 1
  sess15a$end.month.index[s] <- (sess15a$end.year[s]-minYear)*12 +
    sess15a$end.month[s] - minMonth + 1
}#s

##-- Calculate range of years covered by the study
years <- minYear:max(sess15a$start.year)
n.years <- length(years)

##-- Calculate range of months covered by the study 
months <- 1:max(sess15a$start.month.index)
n.months <- length(months)-1

##-- Identify seasons for each month of the study
# Wet season from October until June!
season <- rep(c(1,1,1,1,1,1,2,2,2,1,1,1), n.years)
season <- season[minMonth:(n.months+minMonth-1)]


## -----------------------------------------------------------------------------
## ------ II. CREATE CAPTURE HISTORY ------
##-- Identify in which session each individual was captured
for(c in 1:nrow(m15a)){
  m15a$session[c] <- sess15a$index[sess15a$start.date <= m15a$date[c] & sess15a$end.date >= m15a$date[c]]
}#c

##-- Create a dummy dataset
dummy <- data.frame( transponder  = "dummy",
                     session = sess15a$index)

##-- Combine real and dummy datasets
m15a.dummy <- rbind.fill(m15a, dummy)

##-- Create the matrix of capture history
ch15a <- table(m15a.dummy$transponder, m15a.dummy$session)

##-- Remove dummy individual
ch15a <- ch15a[-which(dimnames(ch15a)[[1]] == "dummy"), ] 

##-- Extract the first detection session for each individual
f <- apply(ch15a, 1, function(x)min(which(x >= 1)))

##-- Remove individuals detected for the first time on the last session
ch15a <- ch15a[which(f != dim(ch15a)[2]), ]
f <- f[which(f != dim(ch15a)[2])]

##-- Reorder and turn into a matrix
ch15a <- ch15a[order(dimnames(ch15a)[[1]]), ]
ch15a <- as.matrix(ch15a)
ch15a[ch15a > 0] <- 1
dim(ch15a)

##-- List individuals
ids <- dimnames(ch15a)[[1]] 
n.individuals <- length(ids)

##-- Sex
sex <- unique(m15a[m15a$transponder %in% ids, c("transponder", "Sexe")])
sex <- sex[order(sex$transponder), ] 
all(dimnames(ch15a)[[1]] == sex$transponder)
sex <- ifelse(sex$Sexe == "f", 1, 2) 


# ##-- Age 
# age <- matrix(data = NA, nrow = n.individuals, ncol = n.sessions)
# for(i in ids){
#   for(s in sess15a$index){
#     temp <- m15a[m15a$transponder == i & m15a$session == s, ]
#     if(length(temp) != 0){
#       age[i,s] <- unique(temp$age_estimation_2)
#     }
#   }
# }

##-- Interval lengths between capture sessions
start.int <- end.int <- NULL
for(i in 1:(n.sessions-1)){
  start.int[i] <- sess15a$start.month.index[i]
  end.int[i] <- sess15a$start.month.index[i+1]-1
}#t


##-- Check the data with plots
hist(rowSums(ch15a))                  ## num. of detections per individual
plot(sess15a$duration, colSums(ch15a)) ## num. of ids detected per session duration



## -----------------------------------------------------------------------------
## ------ III. NIMBLE MODEL ------
##-- Write the CR model in NIMBLE
nimModel <- nimbleCode({
  ## DEMOGRAPHIC PROCESS 
  phi0[1,1] ~ dunif(0,1)
  phi0[2,1] ~ dunif(0,1)
  phi0[1,2] ~ dunif(0,1)
  phi0[2,2] ~ dunif(0,1)
  
  beta[1] ~ dnorm(0,0.01)
  beta[2] ~ dnorm(0,0.01)
  
  for(m in 1:n.months){
    logit(PHI[m,1]) <- logit(phi0[season[m],1]) + beta[1] * m
    logit(PHI[m,2]) <- logit(phi0[season[m],2]) + beta[2] * m
  }# months
  
  for(t in 1:n.intervals){
    phi[t,1] <- prod(PHI[start.int[t]:end.int[t],1])
    phi[t,2] <- prod(PHI[start.int[t]:end.int[t],2])
  }# time
  
  for(i in 1:n.individuals){
    z[i,f[i]] ~ dbern(1)  
    for(t in f[i]:n.intervals){
      z[i,t+1] ~ dbern(z[i,t] * phi[t,sex[i]])  
    }#t
  }#i
  
  ## DETECTION PROCESS
  ## Detection Hazard rate
  lambda[1] ~ dunif(0,5)  
  lambda[2] ~ dunif(0,5)                   
  for(t in 1:n.sessions){
    ## Allows for unequal sampling sessions (sessionDuration[t])
    p[t,1] <- 1-exp(-lambda[1] * sessionDuration[t]) 
    p[t,2] <- 1-exp(-lambda[2] * sessionDuration[t]) 
  }# session
  
  for(i in 1:n.individuals){
    for(t in (f[i]+1):n.sessions){
      y[i,t] ~ dbern(p[t,sex[i]] * z[i,t])
    }#t
  }#i
  
})

##-- Format the data for NIMBLE
nimData <- list( y = ch15a,
                 sessionDuration = sess15a$duration)

nimConstants <- list( n.individuals = dim(ch15a)[1],
                      n.intervals = dim(ch15a)[2]-1,
                      n.sessions = dim(ch15a)[2],
                      n.months = n.months,
                      season = season,
                      sex = sex,
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
                  phi0 = matrix(rep(0.85,4), ncol = 2),
                  beta = c(0,0),
                  lambda = c(0.2,0.2))


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


MCMCtrace(nimOutput, 
          pdf = TRUE, 
          open_pdf = TRUE, 
          filename = 'M15a',
          wd="C:/Users/anvargas/Dropbox/Mouse lemur CMR data/06_Results/01_Model fragment")

## -----------------------------------------------------------------------------