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
modelName <- "microcebe_temp_time_transloc2"
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
##-- Identify protection status for each fragment
STATUS <- c(1,1,2,2,2,2,2,1)

capture_data <- capture_data[capture_data$site %in% fragments, ] 


##-- Manually fix a few problematic individuals
## e.g. individuals captured in multiple fragments
capture_data$date[capture_data$transponder == "/635231D/" & capture_data$date == "2010-08-22"] <- "2010-08-23"
capture_data$site[capture_data$transponder == "/635231D/"] <- "M20"

capture_data$site[capture_data$transponder == "/6395D57/"] <- "M16"

capture_data$date[capture_data$transponder == "/6CBD354/" & capture_data$date == "2014-04-08"] <- "2014-05-08"

capture_data$site[capture_data$transponder == "/6CBE39E/"] <- "M15C"

capture_data$site[capture_data$transponder == "/74C839F/"] <- "M15B"

capture_data$Sexe[capture_data$transponder == "/6CD1185/"] <- "m"


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



## ------   4. Clean data for each individual ------
## Identify first year of capture
minYear <- min(capture_sessions$start.year)
## Identify last year of capture
maxYear <- max(capture_sessions$start.year)
## Identify first month of capture in the first year
minMonth <- min(capture_sessions$start.month[capture_sessions$start.year == minYear])
firstMonth <- minMonth
## Identify last month of capture in the last year
lastMonth <- max(capture_sessions$start.month[capture_sessions$start.year == maxYear])

## Store data for each individual in a list
data_list <- list()
IDs <- unique(capture_data$transponder)
n.individuals <- length(IDs)
# i <- which(IDs == "/6CC3AF9/")

for(i in 1:n.individuals){
  data_list[[i]] <- list()
  
  ##-- Subset data & capture sessions to individual "i"
  print(IDs[i])
  thisData <- capture_data[capture_data$transponder == IDs[i], ]
  thisFragment <- unique(thisData$site)
  data_list[[i]]$fragment <- thisFragment
  print(thisFragment)
  
  ##-- Identify first capture of focal individual
  firstCap <- thisData[which.min(thisData$date), ]
  
  ##-- Identify recapture sessions in this fragment
  thisSessions <- capture_sessions[capture_sessions$site == thisFragment, ]
  
  ##-- Ensure sessions are ordered by start.date
  thisSessions <- thisSessions[order(thisSessions$start.date), ]
  dim(thisSessions)
  
  ##-- Remove recapture sessions before first capture of the focal individual
  thisSessions <- thisSessions[thisSessions$end.date >= firstCap$date, ]
  dim(thisSessions)
  
  ##-- Check if this individual was captured outside "official" sessions
  if(thisSessions$start.date[1] >= firstCap$date){
    thisSessions <- rbind(rep(NA,13),thisSessions)
    thisSessions[1,"start.month"] <- firstCap$month
    thisSessions[1,"start.year"] <- firstCap$year
    thisSessions[1,"site"] <- firstCap$site
    thisSessions[1,"start.date"] <- firstCap$date
    thisSessions[1,"end.date"] <- firstCap$date
    thisSessions[1,"duration"] <- difftime(time1 = thisSessions$end.date[1],
                                           time2 = thisSessions$start.date[1],
                                           units = "days") 
  } 
  
  ##-- Give an index to each recapture session
  thisSessions <- thisSessions[order(thisSessions$start.date), ]
  thisSessions$index <- 1:nrow(thisSessions)
  
  ##-- Identify end months and years of each recapture session
  thisSessions$end.year <- as.numeric(format(thisSessions$end.date,"%Y"))
  thisSessions$end.month <- as.numeric(format(thisSessions$end.date,"%m"))
  
  ##-- Store number of sessions for this fragment
  n.sessions <- nrow(thisSessions)
  data_list[[i]]$n.sessions <- n.sessions
  
  ##-- Calculate months index of each capture and recapture sessions
  for(s in 1:n.sessions){
    thisSessions$start.month.index[s] <- (thisSessions$start.year[s]-minYear)*12 +
      thisSessions$start.month[s] - minMonth + 1
    thisSessions$end.month.index[s] <- (thisSessions$end.year[s]-minYear)*12 +
      thisSessions$end.month[s] - minMonth + 1
  }#s
  data_list[[i]]$start.month.index <- thisSessions$start.month.index
  data_list[[i]]$end.month.index <- thisSessions$end.month.index
  data_list[[i]]$duration <- thisSessions$duration 
  data_list[[i]]$seas <- thisSessions$season 
  data_list[[i]]$temp2 <- thisSessions$T
  data_list[[i]]$precip <- thisSessions$meanPP
  
  ##-- Calculate range of years covered by the study in this fragment
  years <- min(thisSessions$start.year):max(thisSessions$start.year)
  data_list[[i]]$years <- years
  data_list[[i]]$n.years <- length(years)
  
  ##-- Calculate range of months covered by the study in this fragment
  months <- min(thisSessions$start.month.index):max(thisSessions$start.month.index)
  data_list[[i]]$months <- months
  data_list[[i]]$n.months <- length(months)-1
  
  
  
  ## ---------------------------------------------------------------------------
  ## ------ II. CREATE CAPTURE HISTORY ------
  ##-- Identify in which session each individual was recaptured
  for(c in 1:nrow(thisData)){
    temp <- thisSessions$index[thisSessions$start.date <= thisData$date[c] &
                                 thisSessions$end.date >= thisData$date[c]]
    thisData$session[c] <- ifelse(length(temp) == 0,NA,temp)
  }#c
  
  ##-- CHECK IF SOME CAPTURES ARE STILL OUTSIDE THE OFFICIAL CAPTURE SESSIONS
  if(any(is.na(thisData$session))){stop(print(thisData[is.na(thisData$session), ]))}
  
  ##-- Create a dummy dataset
  dummy <- data.frame(transponder  = "dummy", session = thisSessions$index)
  
  ##-- Combine real & dummy datasets
  thisData.dummy <- rbind.fill(thisData, dummy)
  
  ##-- Create the matrix of capture history
  thisCH <- table(thisData.dummy$transponder, thisData.dummy$session)
  
  ##-- Remove dummy individual
  thisCH <- thisCH[-which(dimnames(thisCH)[[1]] == "dummy"), ] 
  
  ##-- Save individual detection history and associated info
  thisCH[thisCH > 0] <- 1
  data_list[[i]]$CH <- thisCH
  data_list[[i]]$ids <- IDs[i]
  
  ##-- Sex
  sex <- unique(thisData$Sexe)
  if(length(sex) != 1)stop(paste("sex issue with id:", IDs[i]))
  data_list[[i]]$sex <- ifelse(sex == "f", 1, 2) 
  
  ##-- Translocation
  transloc <- unique(thisData$residency)
  if(length(transloc) != 1)stop("translocation issue")
  data_list[[i]]$transloc <- ifelse(transloc == "resident", 1, 2) 
  
  ##-- Interval lengths between capture sessions
  data_list[[i]]$start.int <- data_list[[i]]$end.int <- NULL
  for(s in 1:(n.sessions-1)){
    data_list[[i]]$start.int[s] <- thisSessions$start.month.index[s]
    data_list[[i]]$end.int[s] <- thisSessions$start.month.index[s+1]-1
  }#s
}#i
names(data_list) <- IDs[i]



## ------   5. Format data in arrays ------
##-- Extract the number of ids detected 
n.individuals <- length(data_list)

##-- Extract the number of capture sessions per indivdiual
n.sessions <- unlist(lapply(data_list, function(x)x$n.sessions))
maxSessions <- max(n.sessions)

##-- Extract the number of months
months <- unlist(lapply(data_list, function(x)x$months))
maxMonth <- max(months)

##-- Put into ragged arrays
y <- matrix(NA, n.individuals, maxSessions)
sex <- transloc <- status <- rep(NA, n.individuals)
start.int <- end.int <- duration <- seas <- temp2 <- matrix( NA,
                                                             nrow = n.individuals,
                                                             ncol = maxSessions)
for(i in 1:n.individuals){
  y[i,1:n.sessions[i]] <- data_list[[i]]$CH
  sex[i] <- data_list[[i]]$sex
  transloc[i] <- data_list[[i]]$transloc
  start.int[i,1:n.sessions[i]] <- data_list[[i]]$start.month.index
  end.int[i,1:n.sessions[i]] <- data_list[[i]]$end.month.index
  duration[i,1:n.sessions[i]] <- data_list[[i]]$duration
  seas[i,1:n.sessions[i]] <- data_list[[i]]$seas
  temp2[i,1:n.sessions[i]] <- data_list[[i]]$temp2
  status[i] <- STATUS[fragments == data_list[[i]]$fragment]
}#i

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
  
  ## Individual state
  for(i in 1:n.individuals){
    z[i,1] ~ dbern(1)
    for(t in 1:n.intervals[i]){
      phi[i,t] <- prod(PHI[sex[i],status[i],transloc[i],start.int[i,t]:end.int[i,t]])
      z[i,t+1] ~ dbern(z[i,t] * phi[i,t])  
    }#t
  }#i
  
  
  ## DETECTION PROCESS
  lambda0[1,1] ~ dnorm(0,0.01)    ## Detection Hazard rate female/disturbed
  lambda0[2,1] ~ dnorm(0,0.01)    ## Detection Hazard rate male/disturbed
  lambda0[1,2] ~ dnorm(0,0.01)    ## Detection Hazard rate female/protected
  lambda0[2,2] ~ dnorm(0,0.01)    ## Detection Hazard rate male/protected
  
  gamma.temp[1,1] ~ dnorm(0,0.01) ## Detection Hazard rate female/disturbed
  gamma.temp[2,1] ~ dnorm(0,0.01) ## Detection Hazard rate male/disturbed
  gamma.temp[1,2] ~ dnorm(0,0.01) ## Detection Hazard rate female/protected
  gamma.temp[2,2] ~ dnorm(0,0.01) ## Detection Hazard rate male/protected
  
  ## Individual detection
  for(i in 1:n.individuals){
    for(t in 2:n.sessions[i]){
      log(lambda[i,t]) <- lambda0[sex[i],status[i]] + gamma.temp[sex[i],status[i]] * temp2[i,t]
      p[i,t] <- 1-exp(-lambda[i,t] * duration[i,t]) 
      y[i,t] ~ dbern(p[i,t] * z[i,t])
    }#t
  }#i
})



## ------   2. DATA ------
laste <- apply(nimData$y, 1, function(x)max(which(x == 1)))
laste[is.infinite(laste)] <- NA
z.init <- z.data <- matrix(NA, n.individuals, maxSessions)
for(i in 1:n.individuals){
  z.data[i,1:laste[i]] <- 1
  if(laste[i] < n.sessions[i]){
    z.init[i,(laste[i]+1):n.sessions[i]] <- 0
  }
}#i

nimData <- list(y = y,
                z = z.data)

nimConstants <- list( n.individuals = n.individuals,
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


nimInits <- list( z = z.init,
                  logit.phi0 = matrix(2.5,2,2),
                  beta.time = matrix(0,2,2),
                  beta.temp = matrix(0,2,2),
                  beta.transloc = 0,
                  lambda0 = matrix(-1,2,2),
                  gamma.temp = matrix(0,2,2))

nimParams <- c("logit.phi0", "beta.temp", "beta.time","beta.transloc",
               "gamma.temp", "lambda0")



## ------   3. SAVE INPUT ------
# for(c in 1:4){
#   save(nimModel,
#        nimData,
#        nimConstants,
#        nimInits,
#        nimParams, 
#        file = file.path(analysisDir,
#                         modelName,
#                         "inFiles",
#                         paste0(modelName,c,".RData")))
# }#c
# 


## ------   4. FIT ------
#for(c in 1:4){
  # 
  # ##-- Create a NIMBLE model object
  # Rmodel <- nimbleModel( code = nimModel,
  #                        constants = nimConstants,
  #                        data = nimData,
  #                        inits = nimInits,
  #                        calculate = F)
  # Rmodel$calculate()
  # 
  # ##-- Configure and Build MCMC objects
  # conf <- configureMCMC(Rmodel,
  #                       monitors = nimParams,
  #                       print = FALSE)
  # Rmcmc <- buildMCMC(conf)
  # 
  # ##-- Compile and Run MCMC
  # ## Finally, we compile both the model and MCMC objects and execute the 
  # ## compiled MCMC for 50 000 iterations and 3 chains.
  # Cmodel <- compileNimble(Rmodel)
  # Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
  # MCMC_runtime <- system.time(
  #   nimOutput <- runMCMC( Cmcmc,
  #                         niter = 5000,
  #                         nburnin = 0,
  #                         nchains = 1,
  #                         thin = 1,
  #                         samplesAsCodaMCMC = T))
  # plot(nimOutput)
  # save(nimOutput,
  #      MCMC_runtime,
  #      file = file.path(analysisDir,
  #                       modelName, 
  #                       paste0("outFiles/output_",c,".RData")))
#}#c




## -----------------------------------------------------------------------------
## ------ III. PROCESS OUTPUTS -----
## ------   1. TRACEPLOTS -----
outputs <- list.files(file.path(analysisDir, modelName, "outFiles"))

nimOutput2 <- mcmc.list()
for(i in 1:length(outputs)){
  load(file.path(analysisDir,modelName,"outFiles",outputs[i]))
  nimOutput2[[i]] <- as.mcmc(nimOutput)
}
nimOutput <- nimOutput2

pdf(file = file.path(analysisDir, modelName, "Traceplots.pdf"), paper = "A4")
plot(nimOutput)
dev.off()



# ## ------   2. SURVIVAL ~ TEMPERATURE -----
# load(file = file.path(analysisDir,
#                       modelName,
#                       "inFiles",
#                       paste0(modelName,"1.RData")))
# n.months <- nimConstants$n.months
# temp <- nimConstants$temp
# nimMat <- do.call(rbind, nimOutput)
# n.iter <- dim(nimMat)[1]
# beta.temp <- array(nimMat[ ,1:4], c(n.iter,2,2))
# beta.time <- array(nimMat[ ,5:8], c(n.iter,2,2))
# logit.phi0 <- array(nimMat[ ,18:21], c(n.iter,2,2))
# PHI <- array(NA, c(n.iter,2,2,n.months))
# for(m in 1:n.months){
#   PHI[ ,1,1,m] <- ilogit(logit.phi0[ ,1,1] + beta.temp[ ,1,1] * temp[m] + beta.time[ ,1,1] * m)
#   PHI[ ,2,1,m] <- ilogit(logit.phi0[ ,2,1] + beta.temp[ ,2,1] * temp[m] + beta.time[ ,2,1] * m)
#   PHI[ ,1,2,m] <- ilogit(logit.phi0[ ,1,2] + beta.temp[ ,1,2] * temp[m] + beta.time[ ,1,2] * m)
#   PHI[ ,2,2,m] <- ilogit(logit.phi0[ ,2,2] + beta.temp[ ,2,2] * temp[m] + beta.time[ ,2,2] * m)
# }# months



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
beta.transloc <- nimMat[ ,9]
logit.phi0 <- array(nimMat[ ,18:21], c(n.iter,2,2))
PHI <- array(NA, c(n.iter,2,2,2,n.months))
for(m in 1:n.months){
  for(d in 1:2){
    for(s in 1:2){
      for(t in 1:2){
        PHI[ ,s,d,t,m] <- ilogit(logit.phi0[ ,s,d] +
                                   beta.temp[ ,s,d] * temp[m] +
                                   beta.time[ ,s,d] * m +
                                   beta.transloc[ ] * (t-1))
      }#t
    }#s
  }#d
}# months

mean.PHI <- apply(PHI, c(2,3,4,5), mean)
upper.PHI <- apply(PHI, c(2,3,4,5), function(x)quantile(x,0.975))
lower.PHI <- apply(PHI, c(2,3,4,5), function(x)quantile(x,0.025))


pdf(file.path(analysisDir, modelName, "survivalProbabilities.pdf"),
    width = 12, height = 12)
par(mfrow = c(2,2))
sex <- c("females", "males")
degrad <- c("degraded","protected")
for(d in 1:2){
  for(s in 1:2){
    plot(1, type = "n", xlim = c(0,n.months+1), ylim = c(0, 1), axes = F,
         ylab = "Survival prob.",
         xlab = "Months",
         main = paste0(sex[s],"-",degrad[d]))
    axis(1, at = seq(0,250,50), labels = seq(0,250,50))
    axis(2, at = seq(0,1,0.2), labels = seq(0,1,0.2))
    legend( x = 1, y = 0.2,
            title = "Individual status:",
            legend = c("resident", "translocated"),
            bty = "n",
            fill = myCols[c(2,4)])
    
    polygon(x = c(1:n.months,n.months:1),
            y = c(upper.PHI[s,d,1, ],rev(lower.PHI[s,d,1, ])),
            col = adjustcolor(myCols[2],alpha.f = 0.5), border = F)
    points(1:n.months, mean.PHI[s,d,1,], type = "l", lwd = 2, col = myCols[2])
    
    polygon(x = c(1:n.months,n.months:1),
            y = c(upper.PHI[s,d,2, ],rev(lower.PHI[s,d,2, ])),
            col = adjustcolor(myCols[4],alpha.f = 0.5), border = F)
    points(1:n.months, mean.PHI[s,d,2,], type = "l", lwd = 2, col = myCols[4])
  }#s
}#t
dev.off()



## ------   3. DETECTION PROBABILITY ~ TEMPERATURE -----
load(file = file.path( analysisDir,
                       modelName,
                       "inFiles",
                       paste0(modelName,"1.RData")))
temp2 <- seq(min(nimConstants$temp2, na.rm = T),max(nimConstants$temp2, na.rm = T), 0.05)
nimMat <- do.call(rbind, nimOutput)
n.iter <- dim(nimMat)[1]
gamma.temp <- array(nimMat[ ,10:13], c(n.iter,2,2))
lambda0 <- array(nimMat[ ,14:17], c(n.iter,2,2))
P <- array(NA, c(n.iter,2,2,length(temp2)))

## Multi-fragments model
for(t in 1:length(temp2)){
  for(s in 1:2){
    for(f in 1:2){
      P[ ,s,f,t] <- 1-exp(-exp(lambda0[ ,s,f] + gamma.temp[ ,s,f] * temp2[t])*4) 
    }#f
  }#s
}#t

mean.P <- apply(P, c(2,3,4), function(x)mean(x, na.rm = T))
upper.P <- apply(P, c(2,3,4), function(x)quantile(x, 0.975, na.rm = T))
lower.P <- apply(P, c(2,3,4), function(x)quantile(x, 0.025, na.rm = T))

pdf(file.path(analysisDir, modelName, "detectionProbabilities.pdf"),
    width = 20, height = 10)
par(mfrow = c(1,2))
sex <- c("females", "males")
degrad <- c("degraded","protected")
for(d in 1:2){
    plot(1, type = "n",
         xlim = c(min(temp2)-0.4,max(temp2)+0.05),
         ylim = c(0, 1), axes = F,
         ylab = "Detection prob.",
         xlab = "Temperature",
         main = degrad[d])
    axis(1,
         at = seq(min(temp2)-0.4,max(temp2),0.5),
         labels = seq(min(temp2)-0.4,max(temp2),0.5))
    axis(2, at = seq(0,1,0.2), labels = seq(0,1,0.2))
    legend( x = 26, y = 1,
            legend = c("females", "males"),
            bty = "n",
            fill = myCols[c(2,4)])
    
    polygon(x = c(temp2,rev(temp2)),
            y = c(upper.P[1,d, ],rev(lower.P[1,d, ])),
            col = adjustcolor(myCols[2],alpha.f = 0.5), border = F)
    points(temp2, mean.P[1,d, ], type = "l", lwd = 2, col = myCols[2])
    
    polygon(x = c(temp2,rev(temp2)),
            y = c(upper.P[2,d, ],rev(lower.P[2,d, ])),
            col = adjustcolor(myCols[4],alpha.f = 0.5), border = F)
    points(temp2, mean.P[2,d, ], type = "l", lwd = 2, col = myCols[4])
  }#d
dev.off()


## -----------------------------------------------------------------------------