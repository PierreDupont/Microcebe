########################################################
##### ------------ MANDENA MICROCEBE ------------- #####
##### SIMULTANEOUS ANALYSIS for MULTIPLE FRAGMENTS #####
########################################################
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
modelName <- "microcebe_temp_time_transloc"
source("C:/My_documents/RovQuant/Temp/PD/FUNCTIONS/FunctionScripts/wildMap.R")
myCols <- wildMap(4)


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
capture_data$month <- as.numeric(format(capture_data$date,"%m"))

##-- Clean-up fragments names
capture_data$site <- toupper(capture_data$site)
unique(capture_data$site)
capture_data$site[capture_data$site == "M16 VAO"] <- "M16"

## Fragments to remove:
# "M3" "M4" "M6" "M7" "M9" "M11" "M12" "M13R" "M15E" "M16_Bakuli"   
# "Bush clearing" "Bush cl" "Cage" "Mine" "Corr  R." "Corr  E." "Tokon'ala"             
fragments <- c("M5","M13","M15A","M15B","M15C","M15D","M16","M20")
capture_data <- capture_data[capture_data$site %in% fragments, ] 



## ------   2. Capture sessions ------
##-- Load capture sessions data
capture_sessions <- read.csv(file.path(dataDir, "sessions_dates_sites.csv"),h=T)
names(capture_sessions) <- c("start.date", "end.date", "site")

##-- Clean-up fragments names
capture_sessions$site <- toupper(capture_sessions$site)
unique(capture_sessions$site)
capture_sessions$site[capture_sessions$site == "M16 VAO"] <- "M16"
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
capture_sessions <- merge( aggSessions, durationAggSessions,
                           by = c("start.month", "start.year", "site"))

##-- Extract season for each capture session
capture_sessions$season <- ifelse(capture_sessions$start.month %in% c(6,7,8,9), 2, 1)



## ------   3. Weather data -----
weather <- read.csv(file.path(dataDir, "climate_fort_dauphin.csv"),h=T)
weather$MONTH <- rep(1:12,22)
names(weather)

weather <- weather[ ,c("YEAR", "MONTH", "T", "TM", "Tm", "PP", "DAYS_WITH_DATA")]
plot(weather$MONTH,weather$T)

weather$meanPP <- weather$PP/weather$DAYS_WITH_DATA
plot(weather$MONTH,weather$meanPP)

##-- Extract mean temperature for each capture session
capture_sessions <- merge( capture_sessions, weather,
                           by.x = c("start.month", "start.year"),
                           by.y = c("MONTH", "YEAR"),
                           all.x = T, all.y = F)
names(capture_sessions)



## ------   4. Clean data for each fragment ------
## Identify first year of capture
minYear <- min(capture_sessions$start.year)
## Identify last year of capture
maxYear <- max(capture_sessions$start.year)
## Identify first month of capture in the first year
minMonth <- min(capture_sessions$start.month[capture_sessions$start.year == minYear])
firstMonth <- min(capture_sessions$start.month[capture_sessions$start.year == minYear])
## Identify last month of capture in the last year
lastMonth <- max(capture_sessions$start.month[capture_sessions$start.year == maxYear])

## Store data for each fragment in a list
data_list <- list()
n.fragments <- length(fragments)

for(f in 1:n.fragments){
  data_list[[f]] <- list()
  
  ##-- Subset data & capture sessions to fragment "s"
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
  data_list[[f]]$start.month.index <- thisSessions$start.month.index
  data_list[[f]]$end.month.index <- thisSessions$end.month.index
  data_list[[f]]$duration <- thisSessions$duration 
  data_list[[f]]$seas <- thisSessions$season 
  data_list[[f]]$temp2 <- thisSessions$T
  data_list[[f]]$precip <- thisSessions$meanPP
  
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
    #if(length(temp) == 0)stop()
    thisData$session[c] <- ifelse(length(temp) == 0,NA,temp)
  }#c
  
  ## Remove detections outside capture sessions
  if(any(is.na(thisData$session))){print(thisData[is.na(thisData$session), ])}
  thisData <- thisData[!is.na(thisData$session), ]
  
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
  sex <- NULL
  for(i in 1:length(ids)){
    sex[i] <- unique(thisData$Sexe[thisData$transponder == ids[i]])
    if(length(sex[i]) != 1)stop(paste("sex issue with id:", ids[i]))
  }#i
  data_list[[f]]$sex <- ifelse(sex == "f", 1, 2) 
  
  ##-- Translocation
  transloc <- NULL
  for(i in 1:length(ids)){
    transloc[i] <- unique(thisData$residency[thisData$transponder == ids[i]])
    if(length(transloc[i]) != 1)stop("translocation issue")
  }#i
  data_list[[f]]$transloc <- ifelse(transloc == "resident", 1, 2) 
  
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
names(data_list) <- fragments



## ------   5. Format data in arrays ------
##-- Extract the number of ids detected per fragment
n.individuals <- unlist(lapply(data_list, function(x)x$n.individuals))
maxIDs <- max(n.individuals)

##-- Extract the number of capture sessions per fragment
n.sessions <- unlist(lapply(data_list, function(x)x$n.sessions))
maxSessions <- max(n.sessions)

##-- Extract the number of months
months <- unlist(lapply(data_list, function(x)x$months))
maxMonth <- max(months)

##-- Put into ragged arrays
y <- array(NA, c(maxIDs, maxSessions, n.fragments))
sex <- first <- transloc <- matrix(NA, nrow = maxIDs, ncol = length(data_list))
start.int <- end.int <- duration <- seas <- temp2 <- matrix( NA,
                                                             nrow = maxSessions,
                                                             ncol = length(data_list))
for(f in 1:length(data_list)){
  y[1:n.individuals[f],1:n.sessions[f],f] <- data_list[[f]]$CH
  sex[1:n.individuals[f],f] <- data_list[[f]]$sex
  transloc[1:n.individuals[f],f] <- data_list[[f]]$transloc
  first[1:n.individuals[f],f] <- data_list[[f]]$first
  start.int[1:n.sessions[f],f] <- data_list[[f]]$start.month.index
  end.int[1:n.sessions[f],f] <- data_list[[f]]$end.month.index
  duration[1:n.sessions[f],f] <- data_list[[f]]$duration
  seas[1:n.sessions[f],f] <- data_list[[f]]$seas
  temp2[1:n.sessions[f],f] <- data_list[[f]]$temp2
}#f

##-- Identify seasons for each month of the study
years <- unique(unlist(lapply(data_list, function(x)x$years)))
years <- min(years):max(years)

##-- Wet season from October until May!
season <- rep(c(1,1,1,1,1,2,2,2,2,1,1,1), length(years)) 
season <- season[minMonth:(maxMonth + minMonth - 1)]

##-- Temperature & precipitation
weather$index <- (weather$YEAR-minYear)*12 + weather$MONTH- minMonth  + 1
temp <- weather$T[weather$index >= 1 & weather$index <= maxMonth]
precip <- weather$meanPP[weather$index >= 1 & weather$index <= maxMonth]

##-- Identify protection status for each fragment
status <- c(1,1,2,2,2,2,2,1)


## -----------------------------------------------------------------------------
## ------ II. NIMBLE MODEL ------
## ------   1. MODEL ------
nimModel <- nimbleCode({
  
  ## DEMOGRAPHIC PROCESS 
  logit.phi0[1,1] ~ dnorm(0,0.01) ## Baseline survival female/disturbed
  logit.phi0[2,1] ~ dnorm(0,0.01) ## Baseline survival male/disturbed
  logit.phi0[1,2] ~ dnorm(0,0.01) ## Baseline survival female/protected
  logit.phi0[2,2] ~ dnorm(0,0.01) ## Baseline survival male/protected
  
  beta.time[1,1] ~ dnorm(0,0.01)  ## Temporal effect female/disturbed
  beta.time[2,1] ~ dnorm(0,0.01)  ## Temporal effect male/disturbed
  beta.time[1,2] ~ dnorm(0,0.01)  ## Temporal effect female/protected
  beta.time[2,2] ~ dnorm(0,0.01)  ## Temporal effect male/protected
  
  beta.temp[1,1] ~ dnorm(0,0.01)  ## Temperature effect female/disturbed
  beta.temp[2,1] ~ dnorm(0,0.01)  ## Temperature effect male/disturbed
  beta.temp[1,2] ~ dnorm(0,0.01)  ## Temperature effect female/protected
  beta.temp[2,2] ~ dnorm(0,0.01)  ## Temperature effect male/protected
  
  beta.transloc ~ dnorm(0,0.01)
  
  ## Monthly survival probabilities
  for(m in 1:n.months){
    logit(PHI[1,1,1,m]) <- logit.phi0[1,1] + beta.temp[1,1] * temp[m] + beta.time[1,1] * m
    logit(PHI[2,1,1,m]) <- logit.phi0[2,1] + beta.temp[2,1] * temp[m] + beta.time[2,1] * m
    logit(PHI[1,2,1,m]) <- logit.phi0[1,2] + beta.temp[1,2] * temp[m] + beta.time[1,2] * m
    logit(PHI[2,2,1,m]) <- logit.phi0[2,2] + beta.temp[2,2] * temp[m] + beta.time[2,2] * m
    logit(PHI[1,1,2,m]) <- logit(PHI[1,1,1,m]) + beta.transloc
    logit(PHI[2,1,2,m]) <- logit(PHI[2,1,1,m]) + beta.transloc
    logit(PHI[1,2,2,m]) <- logit(PHI[1,2,1,m]) + beta.transloc
    logit(PHI[2,2,2,m]) <- logit(PHI[2,2,1,m]) + beta.transloc
  }# months
  
  ## Multi-fragments model
  for(f in 1:n.fragments){
    for(t in 1:n.intervals[f]){
      phi[1,1,t,f] <- prod(PHI[1,1,status[f],start.int[t,f]:end.int[t,f]])
      phi[2,1,t,f] <- prod(PHI[2,1,status[f],start.int[t,f]:end.int[t,f]])
      phi[1,2,t,f] <- prod(PHI[1,2,status[f],start.int[t,f]:end.int[t,f]])
      phi[2,2,t,f] <- prod(PHI[2,2,status[f],start.int[t,f]:end.int[t,f]])
    }# t
    
    for(i in 1:n.individuals[f]){
      z[i,first[i,f],f] ~ dbern(1)  
      for(t in first[i,f]:n.intervals[f]){
        z[i,t+1,f] ~ dbern(z[i,t,f] * phi[sex[i,f],transloc[i,f],t,f])  
      }#t
    }#i
  }#f
  
  
  ## DETECTION PROCESS
  lambda0[1,1] ~ dnorm(0,0.01)    ## Detection Hazard rate female/disturbed
  lambda0[2,1] ~ dnorm(0,0.01)    ## Detection Hazard rate male/disturbed
  lambda0[1,2] ~ dnorm(0,0.01)    ## Detection Hazard rate female/protected
  lambda0[2,2] ~ dnorm(0,0.01)    ## Detection Hazard rate male/protected
  gamma.temp[1,1] ~ dnorm(0,0.01) ## Detection Hazard rate female/disturbed
  gamma.temp[2,1] ~ dnorm(0,0.01) ## Detection Hazard rate male/disturbed
  gamma.temp[1,2] ~ dnorm(0,0.01) ## Detection Hazard rate female/protected
  gamma.temp[2,2] ~ dnorm(0,0.01) ## Detection Hazard rate male/protected
  
  ## Multi-fragments model
  for(f in 1:n.fragments){
    for(t in 1:n.sessions[f]){
      ## Allows for unequal sampling sessions (sessionDuration[t])
      log(lambda[1,t,f]) <- lambda0[1,status[f]] + gamma.temp[1,status[f]] * temp2[t,f]
      log(lambda[2,t,f]) <- lambda0[2,status[f]] + gamma.temp[2,status[f]] * temp2[t,f]
      p[1,t,f] <- 1-exp(-lambda[1,t,f] * duration[t,f]) 
      p[2,t,f] <- 1-exp(-lambda[2,t,f] * duration[t,f]) 
    }# session
    
    for(i in 1:n.individuals[f]){
      for(t in (first[i,f]+1):n.sessions[f]){
        y[i,t,f] ~ dbern(p[sex[i,f],t,f] * z[i,t,f])
      }#t
    }#i
  }#f
})



## ------   2. DATA ------
nimData <- list( y = y)

nimConstants <- list( n.individuals = n.individuals,
                      n.fragments = n.fragments,
                      n.intervals = n.sessions-1,
                      n.sessions = n.sessions,
                      n.months = maxMonth,
                      sex = sex,
                      transloc = transloc,
                      temp = temp,
                      temp2 = temp2,
                      status = status,
                      first = first,
                      start.int = start.int,
                      end.int = end.int,
                      duration = duration)

last <- apply(nimData$y, c(1,3), function(x)max(which(x == 1)))
last[is.infinite(last)] <- NA
z.init <- array(NA, c(maxIDs,maxSessions,n.fragments))
for(f in 1:n.fragments){
  for(i in 1:n.individuals[f]){
    z.init[i,first[i,f]:last[i,f],f] <- 1
    if(last[i,f] < n.sessions[f]){
      z.init[i,(last[i,f]+1):n.sessions[f],f] <- 0
    }
  }#f
}#i

nimInits <- list( z = z.init,
                  logit.phi0 = matrix(2.5,2,2),
                  beta.time = matrix(0,2,2),
                  beta.temp = matrix(0,2,2),
                  beta.transloc = 0,
                  lambda0 = matrix(0,2,2),
                  gamma.temp = matrix(0,2,2))

nimParams <- c("logit.phi0", "beta.temp", "beta.time","beta.transloc",
               "gamma.temp", "lambda0")



## ------   3. SAVE INPUT ------
for(c in 1:4){
  save(nimModel,
       nimData,
       nimConstants,
       nimInits,
       nimParams, 
       file = file.path(analysisDir,
                        modelName,
                        "inFiles",
                        paste0(modelName,c,".RData")))
}#c



## ------   4. FIT ------
for(c in 1:4){
  
  ##-- Create a NIMBLE model object
  Rmodel <- nimbleModel( code = nimModel,
                         constants = nimConstants,
                         data = nimData,
                         inits = nimInits,
                         calculate = F)
  Rmodel$calculate()
  
  ##-- Configure and Build MCMC objects
  conf <- configureMCMC(Rmodel,
                        monitors = nimParams,
                        print = FALSE)
  Rmcmc <- buildMCMC(conf)
  
  ##-- Compile and Run MCMC
  ## Finally, we compile both the model and MCMC objects and execute the 
  ## compiled MCMC for 50 000 iterations and 3 chains.
  # Cmodel <- compileNimble(Rmodel)
  Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
  MCMC_runtime <- system.time(
    nimOutput <- runMCMC( Cmcmc,
                          niter = 3000,
                          nburnin = 0,
                          nchains = 1,
                          thin = 1,
                          samplesAsCodaMCMC = T))
  plot(nimOutput)
  save(nimOutput,
       MCMC_runtime,
       file = file.path(analysisDir,
                        modelName, 
                        paste0("outFiles/output_",c,".RData")))
}#c




## -----------------------------------------------------------------------------
## ------ III. EXPLORE OUTPUTS -----
## ------   1. TRACEPLOTS -----
outputs <- list.files(file.path(analysisDir,modelName,"outFiles"))

nimOutput <- mcmc.list()
for(i in 1:length(outputs)){
  load(file.path(analysisDir,modelName,"outFiles",outputs[i]))
  nimOutput[[i]] <- as.mcmc(myNimbleOutput)
}

pdf(file = file.path(analysisDir,modelName,"Traceplots.pdf"),paper = "A4")
plot(nimOutput)
dev.off()



## ------   2. SURVIVAL ~ TEMPERATURE -----
load(file = file.path(analysisDir,
                      modelName,
                      "inFiles",
                      paste0(modelName,"1.RData")))
n.months <- nimConstants$n.months
temp <- nimConstants$temp
nimMat <- do.call(rbind, nimOutput)
n.iter <- dim(nimMat)[1]
beta.temp <- array(nimMat[ ,1:4], c(n.iter,2,2))
beta.time <- array(nimMat[ ,5:8], c(n.iter,2,2))
logit.phi0 <- array(nimMat[ ,17:20], c(n.iter,2,2))
PHI <- array(NA, c(n.iter,2,2,2,n.months))
for(m in 1:n.months){
  PHI[ ,1,1,1,m] <- ilogit(logit.phi0[ ,1,1,1] + beta.temp[ ,1,1] * temp[m] + beta.time[ ,1,1] * m)
  PHI[ ,2,1,1,m] <- ilogit(logit.phi0[ ,2,1,1] + beta.temp[ ,2,1] * temp[m] + beta.time[ ,2,1] * m)
  PHI[ ,1,2,1,m] <- ilogit(logit.phi0[ ,1,2,1] + beta.temp[ ,1,2] * temp[m] + beta.time[ ,1,2] * m)
  PHI[ ,2,2,1,m] <- ilogit(logit.phi0[ ,2,2,1] + beta.temp[ ,2,2] * temp[m] + beta.time[ ,2,2] * m)
  PHI[ ,1,1,2,m] <- ilogit(logit.phi0[ ,1,1,2] + beta.temp[ ,1,1] * temp[m] + beta.time[ ,1,1] * m)
  PHI[ ,2,1,2,m] <- ilogit(logit.phi0[ ,2,1,2] + beta.temp[ ,2,1] * temp[m] + beta.time[ ,2,1] * m)
  PHI[ ,1,2,2,m] <- ilogit(logit.phi0[ ,1,2,2] + beta.temp[ ,1,2] * temp[m] + beta.time[ ,1,2] * m)
  PHI[ ,2,2,2,m] <- ilogit(logit.phi0[ ,2,2,2] + beta.temp[ ,2,2] * temp[m] + beta.time[ ,2,2] * m)
}# months

mean.PHI <- apply(PHI, c(2,3,4,5), mean)
upper.PHI <- apply(PHI, c(2,3,4,5), function(x)quantile(x,0.975))
lower.PHI <- apply(PHI, c(2,3,4,5), function(x)quantile(x,0.025))


pdf(file.path(analysisDir, modelName,"survivalProbabilities.pdf"),
    width = 12, height = 7)
par(mfrow = c(1,2))
## FEMALES
plot(1, type = "n", xlim = c(0,n.months+1), ylim = c(0, 1), axes = F,
     ylab = "Survival prob.", xlab = "Months", main = "Females")
axis(1, at = seq(0,250,50), labels = seq(0,250,50))
axis(2, at = seq(0,1,0.2), labels = seq(0,1,0.2))
legend( x = 1, y = 0.3,
        title = "Fragment status:",
        legend = c("Degraded", "Protected"),
        bty = "n",
        fill = myCols[c(2,4)])

## Female disturbed
polygon(x = c(1:n.months,n.months:1),
        y = c(upper.PHI[1,1, ],rev(lower.PHI[1,1, ])),
        col = adjustcolor(myCols[2],alpha.f = 0.5), border = F)
points(1:n.months, mean.PHI[1,1,], type = "l", lwd = 2, col = myCols[2])

## Female protected
polygon(x = c(1:n.months,n.months:1),
        y = c(upper.PHI[1,2, ],rev(lower.PHI[1,2, ])),
        col = adjustcolor(myCols[4],alpha.f = 0.5), border = F)
points(1:n.months, mean.PHI[1,2,], type = "l", lwd = 2, col = myCols[4])


## MALES
plot(1, type = "n", xlim = c(0,n.months+1), ylim = c(0, 1), axes = F,
     ylab = "Survival prob.", xlab = "Months", main = "Males")
axis(1, at = seq(0,250,50), labels = seq(0,250,50))
axis(2, at = seq(0,1,0.2), labels = seq(0,1,0.2))

## Male disturbed
polygon(x = c(1:n.months,n.months:1),
        y = c(upper.PHI[2,1, ],rev(lower.PHI[2,1, ])),
        col = adjustcolor(myCols[2],alpha.f = 0.5), border = F)
points(1:n.months, mean.PHI[2,1,], type = "l", lwd = 2, col = myCols[2])

## Male protected
polygon(x = c(1:n.months,n.months:1),
        y = c(upper.PHI[2,2, ],rev(lower.PHI[2,2, ])),
        col = adjustcolor(myCols[4], alpha.f = 0.5), border = F)
points(1:n.months, mean.PHI[2,2,], type = "l", lwd = 2, col = myCols[4])
dev.off()


## For recruitment, take the approach of N.Hostetter
## use monthly betas, with sum(betas[1:n.months]) = 1
## then derive gammas from betas using:
## gamma[t] <- 1-sum(beta[1:t]) ... or something similar
## check in the age model!!!!
## -----------------------------------------------------------------------------