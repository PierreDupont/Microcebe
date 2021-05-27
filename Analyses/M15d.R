#######################################################
##### ------------ MANDENA MICROCEBE ------------ #####
##### -- PRELIMINARY ANALYSIS for FRAGMENT M15d -- #####
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
library(MCMCvis)


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

##-- Subset to fragment M15d
m15d <- capture_data[capture_data$site == "M15d", ]
unique(m15d$transponder)
#-- Remove non marked individuals too
m15d <- m15d[!(m15d$transponder %in% "//"),]


##-- Format dates
m15d$date <- as.POSIXct(strptime(m15d$date, "%m/%d/%Y"))


## ------   2. Capture sessions ------
##-- Load capture sessions data
capture_sessions <- read.csv(file.path(dataDir, "sessions_dates_sites.csv"),h=T)
names(capture_sessions) <- c("start.date", "end.date", "site")

##-- Subset to fragment M15d
sess15d <- capture_sessions[capture_sessions$site == "M15d", ]

##-- Format dates
sess15d$end.date <- as.POSIXct(strptime(sess15d$end.date, "%m/%d/%Y"))
sess15d$start.date <- as.POSIXct(strptime(sess15d$start.date, "%m/%d/%Y"))

##-- Identify start months and years
sess15d$start.year <- as.numeric(format(sess15d$start.date,"%Y"))
sess15d$start.month <- as.numeric(format(sess15d$start.date,"%m"))

##-- Calculate the duration of each capture session
sess15d$duration <- difftime(time1 = sess15d$end.date,
                            time2 = sess15d$start.date,
                            units = "days") 
sess15d$duration <- as.numeric(sess15d$duration + 1) ## because at least one day of capture

length(sess15d$duration)

##-- Aggregate capture sessions that happened in the same month
startAggSessions <- aggregate(start.date ~ start.month + start.year, data = sess15d, FUN = min)
endAggSessions <- aggregate(end.date ~ start.month + start.year, data = sess15d, FUN = max)
durationAggSessions <- aggregate(duration ~ start.month + start.year, data = sess15d, FUN = sum)
sess15d <- merge(startAggSessions, endAggSessions, by = c("start.month", "start.year"))
sess15d <- merge(sess15d, durationAggSessions, by = c("start.month", "start.year"))

##-- Ensure sessions are ordered by start.date
sess15d <- sess15d[order(sess15d$start.date), ]

##-- Give an index to each capture session
sess15d$index <- 1:nrow(sess15d)

##-- Identify end months and years of each aggregated capture session
sess15d$end.year <- as.numeric(format(sess15d$end.date,"%Y"))
sess15d$end.month <- as.numeric(format(sess15d$end.date,"%m"))
n.sessions <- nrow(sess15d)

##-- Calculate months index of each capture session
minYear <- min(sess15d$start.year)
minMonth <- min(sess15d$start.month[sess15d$start.year == minYear])
for(s in 1:n.sessions){
  sess15d$start.month.index[s] <- (sess15d$start.year[s]-minYear)*12 +
    sess15d$start.month[s] - minMonth  + 1
  sess15d$end.month.index[s] <- (sess15d$end.year[s]-minYear)*12 +
    sess15d$end.month[s] - minMonth + 1
}#s

##-- Calculate range of years covered by the study
years <- minYear:max(sess15d$start.year)
n.years <- length(years)

##-- Calculate range of months covered by the study 
months <- 1:max(sess15d$start.month.index)
n.months <- length(months)-1

##-- Identify seasons for each month of the study
# Wet season from October until June!
season <- rep(c(1,1,1,1,1,1,2,2,2,1,1,1), n.years)
season <- season[minMonth:(n.months+minMonth-1)]


## -----------------------------------------------------------------------------
## ------ II. CREATE CAPTURE HISTORY ------
##-- Identify in which session each individual was captured
c=1
for(c in 1:nrow(m15d)){
  m15d$session[c] <- sess15d$index[sess15d$start.date <= m15d$date[c] & sess15d$end.date >= m15d$date[c]]
}#c

##-- Create a dummy dataset
dummy <- data.frame( transponder  = "dummy",
                     session = sess15d$index)

##-- Combine real and dummy datasets
m15d.dummy <- rbind.fill(m15d, dummy)

##-- Create the matrix of capture history
ch15d <- table(m15d.dummy$transponder, m15d.dummy$session)

##-- Remove dummy individual
ch15d <- ch15d[-which(dimnames(ch15d)[[1]] == "dummy"), ] 

##-- Extract the first detection session for each individual
f <- apply(ch15d, 1, function(x)min(which(x >= 1)))

##-- Remove individuals detected for the first time on the last session
ch15d <- ch15d[which(f != dim(ch15d)[2]), ]
f <- f[which(f != dim(ch15d)[2])]

##-- Reorder and turn into a matrix
ch15d <- ch15d[order(dimnames(ch15d)[[1]]), ]
ch15d <- as.matrix(ch15d)
ch15d[ch15d > 0] <- 1
dim(ch15d)

##-- List individuals
ids <- dimnames(ch15d)[[1]] 
n.individuals <- length(ids)

##-- Sex
sex <- unique(m15d[m15d$transponder %in% ids, c("transponder", "Sexe")])
sex <- sex[order(sex$transponder), ] 
all(dimnames(ch15d)[[1]] == sex$transponder)
sex <- ifelse(sex$Sexe == "f", 1, 2) 


# ##-- Age 
# age <- matrix(data = NA, nrow = n.individuals, ncol = n.sessions)
# for(i in ids){
#   for(s in sess15d$index){
#     temp <- m15d[m15d$transponder == i & m15d$session == s, ]
#     if(length(temp) != 0){
#       age[i,s] <- unique(temp$age_estimation_2)
#     }
#   }
# }

##-- Interval lengths between capture sessions
start.int <- end.int <- NULL
for(i in 1:(n.sessions-1)){
  start.int[i] <- sess15d$start.month.index[i]
  end.int[i] <- sess15d$start.month.index[i+1]-1
}#t


##-- Check the data with plots
hist(rowSums(ch15d))                  ## num. of detections per individual
plot(sess15d$duration, colSums(ch15d)) ## num. of ids detected per session duration



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
nimData <- list( y = ch15d,
                 sessionDuration = sess15d$duration)

nimConstants <- list( n.individuals = dim(ch15d)[1],
                      n.intervals = dim(ch15d)[2]-1,
                      n.sessions = dim(ch15d)[2],
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
          filename = 'M15d',
          wd="C:/Users/anvargas/Dropbox/Mouse lemur CMR data/06_Results/01_Model fragment")
## -----------------------------------------------------------------------------