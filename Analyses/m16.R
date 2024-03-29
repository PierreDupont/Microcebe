#######################################################
##### ------------ MANDENA MICROCEBE ------------ #####
##### -- PRELIMINARY ANALYSIS for FRAGMENT m16 -- #####
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

nrow(capture_data)
##-- Remove duplicates from the capture data
capture_data <- capture_data[!duplicated(capture_data[ ]),] %>% droplevels()

##-- Subset to fragment m16
m16 <- capture_data[capture_data$site %in% c("M16", "M16 Vao"), ]

#-- Remove non marked individuals
m16 <- m16[!(m16$transponder %in% "//"),]

##-- Format dates
m16$date <- as.POSIXct(strptime(m16$date, "%m/%d/%Y"))


## ------   2. Capture sessions ------
##-- Load capture sessions data
capture_sessions <- read.csv(file.path(dataDir, "sessions_dates_sites.csv"),h=T)
names(capture_sessions) <- c("start.date", "end.date", "site")

##-- Subset to fragment m16
sess16 <- capture_sessions[capture_sessions$site%in% c("M16", "M16 Vao"), ]

##-- Format dates
sess16$end.date <- as.POSIXct(strptime(sess16$end.date, "%m/%d/%Y"))
sess16$start.date <- as.POSIXct(strptime(sess16$start.date, "%m/%d/%Y"))

##-- Identify start months and years
sess16$start.year <- as.numeric(format(sess16$start.date,"%Y"))
sess16$start.month <- as.numeric(format(sess16$start.date,"%m"))

##-- Calculate the duration of each capture session
sess16$duration <- difftime(time1 = sess16$end.date,
                            time2 = sess16$start.date,
                            units = "days") 
sess16$duration <- as.numeric(sess16$duration + 1) ## because at least one day of capture

##-- Aggregate capture sessions that happened in the same month
startAggSessions <- aggregate(start.date ~ start.month + start.year, data = sess16, FUN = min)
endAggSessions <- aggregate(end.date ~ start.month + start.year, data = sess16, FUN = max)
durationAggSessions <- aggregate(duration ~ start.month + start.year, data = sess16, FUN = sum)
sess16 <- merge(startAggSessions, endAggSessions, by = c("start.month", "start.year"))
sess16 <- merge(sess16, durationAggSessions, by = c("start.month", "start.year"))

##-- Ensure sessions are ordered by start.date
sess16 <- sess16[order(sess16$start.date), ]

##-- Give an index to each capture session
sess16$index <- 1:nrow(sess16)

##-- Identify end months and years of each aggregated capture session
sess16$end.year <- as.numeric(format(sess16$end.date,"%Y"))
sess16$end.month <- as.numeric(format(sess16$end.date,"%m"))
n.sessions <- nrow(sess16)

##-- Calculate months index of each capture session
minYear <- min(sess16$start.year)
minMonth <- min(sess16$start.month[sess16$start.year == minYear])
for(s in 1:n.sessions){
  sess16$start.month.index[s] <- (sess16$start.year[s]-minYear)*12 +
    sess16$start.month[s] - minMonth  + 1
  sess16$end.month.index[s] <- (sess16$end.year[s]-minYear)*12 +
    sess16$end.month[s] - minMonth + 1
}#s

##-- Calculate range of years covered by the study
years <- minYear:max(sess16$start.year)
n.years <- length(years)

##-- Calculate range of months covered by the study 
months <- 1:max(sess16$start.month.index)
n.months <- length(months)-1

##-- Identify seasons for each month of the study
# Wet season from October until June!
season <- rep(c(1,1,1,1,1,1,2,2,2,1,1,1), n.years)
season <- season[minMonth:(n.months+minMonth-1)]


## -----------------------------------------------------------------------------
## ------ II. CREATE CAPTURE HISTORY ------
##-- Identify in which session each individual was captured
for(c in 1:nrow(m16)){
  m16$session[c] <- sess16$index[sess16$start.date <= m16$date[c] & sess16$end.date >= m16$date[c]]
}#c

##-- Create a dummy dataset
dummy <- data.frame( transponder  = "dummy",
                     session = sess16$index)

##-- Combine real and dummy datasets
m16.dummy <- rbind.fill(m16, dummy)

##-- Create the matrix of capture history
ch16 <- table(m16.dummy$transponder, m16.dummy$session)

##-- Remove dummy individual
ch16 <- ch16[-which(dimnames(ch16)[[1]] == "dummy"), ] 

##-- Extract the first detection session for each individual
f <- apply(ch16, 1, function(x)min(which(x >= 1)))

##-- Remove individuals detected for the first time on the last session
ch16 <- ch16[which(f != dim(ch16)[2]), ]
f <- f[which(f != dim(ch16)[2])]

##-- Reorder and turn into a matrix
ch16 <- ch16[order(dimnames(ch16)[[1]]), ]
ch16 <- as.matrix(ch16)
ch16[ch16 > 0] <- 1
dim(ch16)

##-- List individuals
ids <- dimnames(ch16)[[1]] 
n.individuals <- length(ids)

##-- Sex
sex <- unique(m16[m16$transponder %in% ids, c("transponder", "Sexe")])
sex <- sex[order(sex$transponder), ] 
all(dimnames(ch16)[[1]] == sex$transponder)
sex <- ifelse(sex$Sexe == "f", 1, 2) 


# ##-- Age 
# age <- matrix(data = NA, nrow = n.individuals, ncol = n.sessions)
# for(i in ids){
#   for(s in sess16$index){
#     temp <- m16[m16$transponder == i & m16$session == s, ]
#     if(length(temp) != 0){
#       age[i,s] <- unique(temp$age_estimation_2)
#     }
#   }
# }

##-- Interval lengths between capture sessions
start.int <- end.int <- NULL
for(i in 1:(n.sessions-1)){
  start.int[i] <- sess16$start.month.index[i]
  end.int[i] <- sess16$start.month.index[i+1]-1
}#t


##-- Check the data with plots
hist(rowSums(ch16))                  ## num. of detections per individual
plot(sess16$duration, colSums(ch16)) ## num. of ids detected per session duration



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
nimData <- list( y = ch16,
                 sessionDuration = sess16$duration)

nimConstants <- list( n.individuals = dim(ch16)[1],
                      n.intervals = dim(ch16)[2]-1,
                      n.sessions = dim(ch16)[2],
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
          filename = 'M16',
          wd="C:/Users/anvargas/Dropbox/Mouse lemur CMR data/06_Results/01_Model fragment")

## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
## For recruitment, take the approach of N.Hostetter
## use monthly betas, with sum(betas[1:n.months]) = 1
## then derive gammas from betas using:
## gamma[t] <- 1-sum(beta[1:t]) ... or something similar
## check in the age model!!!!