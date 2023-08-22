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
library(xtable)
library(data.table)
library(ggplot2)


## ------ WORKING DIRECTORIES ------
source("workingDirectories.R")
modelName <- "m_phi[status_sex_temp_time_transloc]_p[site_sex_temp]_RJ"
source("C:/My_documents/RovQuant/Temp/PD/FUNCTIONS/FunctionScripts/wildMap.R")
myCols <- wildMap(4)


## -----------------------------------------------------------------------------
## ------ I. LOAD & CLEAN DATA ------
## ------   1. Individual captures ------
##-- Load individual capture data 
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
  
  ##-- Site
  data_list[[i]]$site <- which(fragments == data_list[[i]]$fragment) 
  
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
sex <- transloc <- status <- site <- rep(NA, n.individuals)
start.int <- end.int <- duration <- seas <- temp2 <- matrix( NA,
                                                             nrow = n.individuals,
                                                             ncol = maxSessions)
for(i in 1:n.individuals){
  y[i,1:n.sessions[i]] <- data_list[[i]]$CH
  site[i] <- data_list[[i]]$site
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
scaleTemp <- as.numeric(scale(temp))
scaleMonth <- as.numeric(scale(1:length(scaleTemp)))
precip <- weather$meanPP[weather$index >= 1 & weather$index <= maxMonth]




## -----------------------------------------------------------------------------
## ------ II. NIMBLE MODEL ------
## ------   1. MODEL ------
nimModel <- nimbleCode({
  
  ## DEMOGRAPHIC PROCESS
  psi ~ dunif(0,1)
  
  for(d in 1:2){
    for(s in 1:2){
      z.time[d,s] ~ dbern(psi)
      z.temp[d,s] ~ dbern(psi)
      z.transloc[d,s] ~ dbern(psi)
      logit.phi0[d,s] ~ dnorm(0,0.01)     ## Baseline survival[site,sex]
      beta.time[d,s] ~ dnorm(0,0.01)      ## Temporal effect[site,sex]
      beta.temp[d,s] ~ dnorm(0,0.01)      ## Temperature effect[site,sex]
      beta.transloc[d,s] ~ dnorm(0,0.01)  ## Translocation effect[site,sex]
      for(m in 1:n.months){
        logit(PHI[d,s,1,m]) <- logit.phi0[d,s] +
          beta.temp[d,s] * z.temp[d,s] * temp[m] +
          beta.time[d,s] * z.time[d,s] * month[m]
        logit(PHI[d,s,2,m]) <- logit.phi0[d,s] +
          beta.temp[d,s] * z.temp[d,s] * temp[m] +
          beta.time[d,s] * z.time[d,s] * month[m] +
          beta.transloc[d,s] * z.transloc[d,s]
      }#m
    }#ss
  }#s
  
  
  ## Individual state
  for(i in 1:n.individuals){
    z[i,1] ~ dbern(1)
    for(t in 1:n.intervals[i]){
      phi[i,t] <- prod(PHI[status[i],sex[i],transloc[i],start.int[i,t]:end.int[i,t]])
      z[i,t+1] ~ dbern(z[i,t] * phi[i,t])
    }#t
  }#i
  
  
  ## DETECTION PROCESS
  for(f in 1:n.sites){
    for(s in 1:2){
      lambda0[f,s] ~ dnorm(0,0.01)    ## Detection Hazard rate female/disturbed
      gamma.temp[f,s] ~ dnorm(0,0.01) ## Detection Hazard rate female/disturbed
      z.det[f,s] ~ dbern(psi)
    }
  }
  
  ## Individual detection
  for(i in 1:n.individuals){
    for(t in 2:n.sessions[i]){
      log(lambda[i,t]) <- lambda0[site[i],sex[i]] +
        gamma.temp[site[i],sex[i]] * temp2[i,t] * z.det[site[i],sex[i]]
      p[i,t] <- 1-exp(-lambda[i,t] * duration[i,t])
      y[i,t] ~ dbern(p[i,t] * z[i,t])
    }#t
  }#i
}) 



## ------   2. DATA ------
laste <- apply(y, 1, function(x)max(which(x == 1)))
laste[is.infinite(laste)] <- NA
z.init <- z.data <- matrix(NA, n.individuals, maxSessions)
for(i in 1:n.individuals){
  z.data[i,1:laste[i]] <- 1
  if(laste[i] < n.sessions[i]){
    z.init[i,(laste[i]+1):n.sessions[i]] <- 0
  }
}#i

nimData <- list( y = y,
                 z = z.data)

nimConstants <- list( n.individuals = n.individuals,
                      n.intervals = n.sessions-1,
                      n.sessions = n.sessions,
                      n.sites = length(fragments),
                      n.months = maxMonth,
                      sex = sex,
                      site = site,
                      transloc = transloc,
                      temp = scaleTemp,
                      month = scaleMonth,
                      temp2 = temp2,
                      status = status,
                      first = first,
                      start.int = start.int,
                      end.int = end.int,
                      duration = duration)


nimInits <- list( z = z.init,
                  psi = 0.9,
                  logit.phi0 = matrix(2.5,2,2),
                  beta.time = matrix(0,2,2),
                  beta.temp = matrix(0,2,2),
                  beta.transloc = matrix(0,2,2),
                  z.time = matrix(1,2,2),
                  z.temp = matrix(1,2,2),
                  z.transloc = matrix(1,2,2),
                  
                  lambda0 = matrix(-1,8,2),
                  gamma.temp = matrix(0,8,2),
                  z.det = matrix(1,8,2)
                  )

nimParams <- c( "logit.phi0",
                "beta.temp", "beta.time", "beta.transloc",
                "z.temp", "z.time", "z.transloc", "z.det",
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
##-- Create a NIMBLE model object
Rmodel <- nimbleModel( code = nimModel,
                       constants = nimConstants,
                       data = nimData,
                       inits = nimInits,
                       calculate = F)
Rmodel$calculate()

##-- Configure and Build RJ-MCMC objects
conf <- configureMCMC(Rmodel,
                      monitors = nimParams,
                      print = FALSE)

configureRJ(conf,
            targetNodes = c("beta.temp"),
            indicatorNodes = c('z.temp'),
            control = list(mean = 0, scale = .2))


configureRJ(conf,
            targetNodes = c( "beta.time"),
            indicatorNodes = c('z.time'),
            control = list(mean = 0, scale = .2))

configureRJ(conf,
            targetNodes = c("beta.transloc"),
            indicatorNodes = c('z.transloc'),
            control = list(mean = 0, scale = .2))

configureRJ(conf,
            targetNodes = c("gamma.temp"),
            indicatorNodes = c('z.det'),
            control = list(mean = 0, scale = .2))

## Check the assigned samplers
conf$printSamplers(nimParams)

## Build MCMC
Rmcmc <- buildMCMC(conf)

##-- Compile and Run MCMC
## Finally, we compile both the model and MCMC objects and execute the
## compiled MCMC for 50 000 iterations and 3 chains.
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
MCMC_runtime <- system.time(
  nimOutput <- runMCMC( Cmcmc,
                        niter = 55000,
                        nburnin = 5000,
                        nchains = 4,
                        thin = 1,
                        samplesAsCodaMCMC = T))
#plot(nimOutput)
save(nimOutput,
     MCMC_runtime,
     file = file.path(analysisDir,
                      modelName,
                      paste0("outFiles/output.RData")))




## -----------------------------------------------------------------------------
## ------ III. PROCESS OUTPUTS -----
load(file.path(analysisDir, modelName, "inFiles", paste0(modelName,"1.RData")))
load(file.path(analysisDir, modelName, paste0("outFiles/output.RData")))
load(file.path(analysisDir, modelName, paste0("outFiles/processed_output.RData")))


## ------   1. Process and save MCMC samples -----
# outputs <- list.files(file.path(analysisDir, modelName, "outFiles"))
# parm.index <- c(1,2,9,10,17,18,25,26,33,34,41,42,49:82,89,90)
# nimOutput2 <- mcmc.list()
# for(i in 1:length(outputs)){
#   load(file.path(analysisDir,modelName,"outFiles",outputs[i]))
#   nimOutput2[[i]] <- as.mcmc(nimOutput[, parm.index])
# }
# nimOutput <- nimOutput2
# 
# pdf(file = file.path(analysisDir, modelName, "Traceplots.pdf"))
# plot(nimOutput)
# dev.off()
#
#
# res <- ProcessCodaOutput(nimOutput)
# save(res, nimOutput,
#      file = file.path(analysisDir,
#                       modelName,
#                       paste0("outFiles/processed_output.RData")))



## ------   2. RJ-MCMC PLOTS ------
##---- Get model features
n.chains <- length(nimOutput)
n.iterations <- dim(nimOutput[[1]])[1]

##---- Get target covariates
covNames <- c("lphi0.dist.f","lphi0.dist.m","lphi0.prot.f","lphi0.prot.m",
              "temp.dist.f", "temp.dist.m", "temp.prot.f", "temp.prot.m",
              "time.dist.f","time.dist.m", "time.prot.f", "time.prot.m",
              "transloc.dist.f","transloc.dist.m", "transloc.prot.f", "transloc.prot.m")
# , "temp.det.1.f","temp.det.2.f","temp.det.3.f","temp.det.4.f","temp.det.5.f","temp.det.6.f","temp.det.7.f","temp.det.8.f",
#   "temp.det.1.m","temp.det.2.m","temp.det.3.m","temp.det.4.m","temp.det.5.m","temp.det.6.m","temp.det.7.m","temp.det.8.m")

zRJ.wide <- data.table(cbind(rep(1,200000),rep(1,200000),rep(1,200000),rep(1,200000),
                             res$sims.list$z.temp[ ,1,1:2],res$sims.list$z.temp[ ,2,1:2],
                             res$sims.list$z.time[ ,1,1:2],res$sims.list$z.time[ ,2,1:2],
                             res$sims.list$z.transloc[ ,1,1:2],res$sims.list$z.transloc[ ,2,1:2]))
                       #, res$sims.list$z.det[ , ,1], res$sims.list$z.det[ , ,2]))
dimnames(zRJ.wide) <- list(NULL, covNames)

betas.wide <- data.table(cbind(res$sims.list$logit.phi0[ ,1,1:2],res$sims.list$logit.phi0[ ,2,1:2],
                               res$sims.list$beta.temp[ ,1,1:2],res$sims.list$beta.temp[ ,2,1:2],
                               res$sims.list$beta.time[ ,1,1:2],res$sims.list$beta.time[ ,2,1:2],
                               res$sims.list$beta.transloc[ ,1,1:2],res$sims.list$beta.transloc[ ,2,1:2]))
                         #, res$sims.list$gamma.temp[ , ,1], res$sims.list$gamma.temp[ , ,2]))
dimnames(betas.wide) <- list(NULL, covNames)


##---- List model combinations
mods <- apply(zRJ.wide, 1, function(x){paste(covNames[x == 1], collapse = "+")})

betas.wide$model <- zRJ.wide$model <- gsub("\\(Intercept\\)\\+", "", mods)
betas.wide$chain <- zRJ.wide$chain <- rep(1:n.chains, each = n.iterations)
betas.wide$iteration <- zRJ.wide$iteration <- rep(1:n.iterations, n.chains)

zRJ.df <- melt(zRJ.wide, id.vars = c("iteration", "chain", "model"))
names(zRJ.df) <- c("iteration", "chain", "model", "variable", "value")

betas.df <-  melt(betas.wide, id.vars = c("iteration", "chain", "model"))
names(betas.df) <-  c("iteration", "chain", "model", "variable", "value")

betas.df$value[zRJ.df$value == 0] <- NA

# betas.aggr <- betas.df %>%
#   group_by("var") %>%
#   summarise(p.inclusion = mean(!is.na("value")))
betas.aggr <- data.table(do.call(rbind, lapply(levels(betas.df$variable), function(x){
  tmp <- betas.df[betas.df$variable == x, ]
  out <- cbind("variable" = x,
               "p.inclusion" = mean(!is.na(tmp$value)))
})))

betas.df <- merge(betas.df, betas.aggr)

included <- zRJ.df$value == 1

betas.df <- betas.df[included, ]
zRJ.df <- zRJ.df[included, ]

betas.df <- betas.df[order(betas.df$variable, betas.df$model, betas.df$chain), ]
zRJ.df <- zRJ.df[order(zRJ.df$variable, zRJ.df$model, zRJ.df$chain), ]

myfun1 <- function(x) 1:length(x)

temp <- betas.df %>% group_by(variable, model, chain) %>%
  summarize(iteration.model = myfun1(value))

betas.df$iteration.model <- temp$iteration.model

aggr <- data.frame(table(betas.df$model) / length(betas.df$model))
names(aggr) <- c("model", "weight")
betas.df <- merge(betas.df, aggr)

##---- MODEL TALLY
aggr <- aggr[order(aggr$weight, decreasing = TRUE),]
aggr$model <- factor(aggr$model, levels = aggr$model)
reduced_aggr <- filter(aggr, weight >= 0.01)

pdf(file.path(analysisDir, modelName, "results_RJ-MCMC.pdf"),
    width = 10, height = 12)
ggplot(data = reduced_aggr,
       mapping =  aes(x = model, y = weight, alpha = weight)) +
  geom_col(fill = "magenta") +
  theme(axis.text.x = element_text(
               angle = 90,
               vjust = 1,
               hjust = 1
             )) + ylab("Weight") + xlab("Models")


# ##---- COEFFICIENT TRACE PLOTS (OVERALL)
# ggplot(data = betas.df, aes(
#   x = iteration,
#   y = value,
#   color = factor(chain))) +
#   geom_line() +
#   facet_wrap(~ variable, scales = "free") +
#   xlab("Iteration") +  theme(legend.position = "none")


# ##---- COEFFICIENT TRACE PLOTS (MODEL-SPECIFIC)
# ggplot(data = betas.df, aes(
#   x = iteration.model,
#   y = value,
#   color = factor(chain))) +
#   geom_line() +
#   facet_grid(variable ~ model, margins = FALSE, scales = "free") +
#   xlab("Iteration") +  theme(legend.position = "none")


##---- PLOT COEFFICIENT ESTIMATES (OVERALL)
ggplot(betas.df, aes(value, variable, alpha = p.inclusion)) +
  geom_violin(
    draw_quantiles = c(0.025, 0.5, 0.975),
    fill = "turquoise",
    color = "white") +
  xlim(-10,2) +
  geom_vline(xintercept = 0)


##---- PLOT COEFFICIENT ESTIMATES (MODEL-SPECIFIC)
reduced_betas.df <- betas.df[betas.df$model %in% reduced_aggr$model]
ggplot(reduced_betas.df, aes(value, variable, alpha = weight)) +
  geom_violin(
    draw_quantiles = c(0.025, 0.5, 0.975),
    fill = "magenta",
    color = grey(1)) +
  xlim(-10,2) +
  geom_vline(xintercept = 0) +
  facet_wrap(~ model)

dev.off()



## ------   3. SURVIVAL PLOTS (best model) -----
n.months <- nimConstants$n.months
temp <- nimConstants$temp
month <- nimConstants$month

##-- Subset MCMC samples to the best model
nimMat <- betas.df[betas.df$weight == max(betas.df$weight), ]
table(nimMat$variable)

##-- Females protected ----
logit.phi0 <- nimMat$value[nimMat$variable %in% "lphi0.prot.f"]
beta.time <- nimMat$value[nimMat$variable %in% "time.prot.f"]
n.iterations <- length(logit.phi0)
phi.f.prot <- matrix(NA, n.iterations, n.months)
for(m in 1:n.months){
  phi.f.prot[ ,m] <- ilogit(logit.phi0 + beta.time * month[m])
}#m
mean.phi.f.prot <- mean.phi.f.prot.t <- apply(phi.f.prot, 2, mean)
upper.phi.f.prot <- upper.phi.f.prot.t <- apply(phi.f.prot, 2, function(x)quantile(x,0.975))
lower.phi.f.prot <- lower.phi.f.prot.t <- apply(phi.f.prot, 2, function(x)quantile(x,0.025))


##-- Females disturbed ----
logit.phi0 <- nimMat$value[nimMat$variable %in% "lphi0.dist.f"]
beta.time <- nimMat$value[nimMat$variable %in% "time.dist.f"]
beta.transloc <- nimMat$value[nimMat$variable %in% "transloc.dist.f"]
n.iterations <- length(logit.phi0)
phi.f.dist <- matrix(NA, n.iterations, n.months)
phi.f.dist.t <- matrix(NA, n.iterations, n.months)
for(m in 1:n.months){
  phi.f.dist[ ,m] <- ilogit(logit.phi0 + beta.time * month[m])
  phi.f.dist.t[ ,m] <- ilogit(logit.phi0 + beta.time * month[m] + beta.transloc)
}#m
mean.phi.f.dist <- apply(phi.f.dist, 2, mean)
upper.phi.f.dist <- apply(phi.f.dist, 2, function(x)quantile(x,0.975))
lower.phi.f.dist <- apply(phi.f.dist, 2, function(x)quantile(x,0.025))

mean.phi.f.dist.t <- apply(phi.f.dist.t, 2, mean)
upper.phi.f.dist.t <- apply(phi.f.dist.t, 2, function(x)quantile(x,0.975))
lower.phi.f.dist.t <- apply(phi.f.dist.t, 2, function(x)quantile(x,0.025))



##-- Males protected ----
logit.phi0 <- nimMat$value[nimMat$variable %in% "lphi0.prot.m"]
beta.transloc <- nimMat$value[nimMat$variable %in% "transloc.prot.m"]
n.iterations <- length(logit.phi0)
phi.m.prot <- matrix(NA, n.iterations, n.months)
phi.m.prot.t <- matrix(NA, n.iterations, n.months)
for(m in 1:n.months){
  phi.m.prot[ ,m] <- ilogit(logit.phi0) 
  phi.m.prot.t[ ,m] <- ilogit(logit.phi0 + beta.transloc)
  
}#m
mean.phi.m.prot <- apply(phi.m.prot, 2, mean)
upper.phi.m.prot <- apply(phi.m.prot, 2, function(x)quantile(x,0.975))
lower.phi.m.prot <- apply(phi.m.prot, 2, function(x)quantile(x,0.025))

mean.phi.m.prot.t <- apply(phi.m.prot.t, 2, mean)
upper.phi.m.prot.t <- apply(phi.m.prot.t, 2, function(x)quantile(x,0.975))
lower.phi.m.prot.t <- apply(phi.m.prot.t, 2, function(x)quantile(x,0.025))


##-- Males disturbed ----
logit.phi0 <- nimMat$value[nimMat$variable %in% "lphi0.dist.m"]
beta.transloc <- nimMat$value[nimMat$variable %in% "transloc.dist.m"]
n.iterations <- length(logit.phi0)
phi.m.dist <- matrix(NA, n.iterations, n.months)
phi.m.dist.t <- matrix(NA, n.iterations, n.months)
for(m in 1:n.months){
  phi.m.dist[ ,m] <- ilogit(logit.phi0) 
  phi.m.dist.t[ ,m] <- ilogit(logit.phi0 + beta.transloc)
  
}#m
mean.phi.m.dist <- apply(phi.m.dist, 2, mean)
upper.phi.m.dist <- apply(phi.m.dist, 2, function(x)quantile(x,0.975))
lower.phi.m.dist <- apply(phi.m.dist, 2, function(x)quantile(x,0.025))

mean.phi.m.dist.t <- apply(phi.m.dist.t, 2, mean)
upper.phi.m.dist.t <- apply(phi.m.dist.t, 2, function(x)quantile(x,0.975))
lower.phi.m.dist.t <- apply(phi.m.dist.t, 2, function(x)quantile(x,0.025))



##----
pdf(file.path(analysisDir, modelName, "survivalProbabilities_best model.pdf"),
    width = 10, height = 10)
par(mfrow = c(2,2))
##-- Females protected plot ----
par(mar = c(0,4,6,0))
plot(1, type = "n", xlim = c(0,n.months+1), ylim = c(0, 1), axes = F,
     ylab = "Survival probability",
     xlab = "",
     main = "Females")
axis(1, at = seq(11,260,24), labels = seq(1999,2019,2))
axis(2, at = seq(0,1,0.2), labels = seq(0,1,0.2))
polygon(x = c(1:n.months,n.months:1),
        y = c(upper.phi.f.prot.t,rev(lower.phi.f.prot.t)),
        col = adjustcolor(myCols[2],alpha.f = 0.5), border = F)
points(1:n.months, mean.phi.f.prot.t, type = "l", lwd = 3, col = myCols[2])
polygon(x = c(1:n.months,n.months:1),
        y = c(upper.phi.f.prot,rev(lower.phi.f.prot)),
        col = adjustcolor(myCols[4],alpha.f = 0.5), border = F)
points(1:n.months, mean.phi.f.prot, type = "l", lwd = 3, col = myCols[4])


##-- Males protected plot ----
par(mar = c(0,2,6,2))
plot(1, type = "n", xlim = c(0,n.months+1), ylim = c(0, 1), axes = F,
     ylab = "Survival probability",
     xlab = "",
     main = "Males")
axis(1, at = seq(11,260,24), labels = seq(1999,2019,2))
axis(2, at = seq(0,1,0.2), labels = seq(0,1,0.2))
polygon(x = c(1:n.months,n.months:1),
        y = c(upper.phi.m.prot.t,rev(lower.phi.m.prot.t)),
        col = adjustcolor(myCols[2],alpha.f = 0.5), border = F)
points(1:n.months, mean.phi.m.prot.t, type = "l", lwd = 3, col = myCols[2])
polygon(x = c(1:n.months,n.months:1),
        y = c(upper.phi.m.prot,rev(lower.phi.m.prot)),
        col = adjustcolor(myCols[4],alpha.f = 0.5), border = F)
points(1:n.months, mean.phi.m.prot, type = "l", lwd = 3, col = myCols[4])


##-- Females fragmented plot ----
par(mar = c(4,4,2,0))
plot(1, type = "n", xlim = c(0,n.months+1), ylim = c(0, 1), axes = F,
     ylab = "Survival probability",
     xlab = "",
     main = "")
axis(1, at = seq(11,260,24), labels = seq(1999,2019,2))
axis(2, at = seq(0,1,0.2), labels = seq(0,1,0.2))
polygon(x = c(1:n.months,n.months:1),
        y = c(upper.phi.f.dist.t,rev(lower.phi.f.dist.t)),
        col = adjustcolor(myCols[2],alpha.f = 0.5), border = F)
points(1:n.months, mean.phi.f.dist.t, type = "l", lwd = 3, col = myCols[2])
polygon(x = c(1:n.months,n.months:1),
        y = c(upper.phi.f.dist,rev(lower.phi.f.dist)),
        col = adjustcolor(myCols[4],alpha.f = 0.5), border = F)
points(1:n.months, mean.phi.f.dist, type = "l", lwd = 3, col = myCols[4])


##-- Males fragmented plot ----
par(mar = c(4,2,2,2))
plot(1, type = "n", xlim = c(0,n.months+1), ylim = c(0, 1), axes = F,
     ylab = "Survival probability",
     xlab = "",
     main = "")
axis(1, at = seq(11,260,24), labels = seq(1999,2019,2))
axis(2, at = seq(0,1,0.2), labels = seq(0,1,0.2))
polygon(x = c(1:n.months,n.months:1),
        y = c(upper.phi.m.dist.t,rev(lower.phi.m.dist.t)),
        col = adjustcolor(myCols[2],alpha.f = 0.5), border = F)
points(1:n.months, mean.phi.m.dist.t, type = "l", lwd = 3, col = myCols[2])
polygon(x = c(1:n.months,n.months:1),
        y = c(upper.phi.m.dist,rev(lower.phi.m.dist)),
        col = adjustcolor(myCols[4],alpha.f = 0.5), border = F)
points(1:n.months, mean.phi.m.dist, type = "l", lwd = 3, col = myCols[4])


    
dev.off()



## ------   4. SURVIVAL PLOTS (all iterations) -----
n.months <- nimConstants$n.months
temp <- nimConstants$temp
month <- nimConstants$month
PHI <- array(NA, c(dim(res$sims.list$beta.temp)[1],2,2,2,n.months))
for(d in 1:2){
  for(s in 1:2){
    for(m in 1:n.months){
      PHI[ ,d,s,1,m] <- ilogit(res$sims.list$logit.phi0[ ,d,s] +
        res$sims.list$beta.temp[ ,d,s] * res$sims.list$z.temp[ ,d,s] * temp[m] +
        res$sims.list$beta.time[ ,d,s] * res$sims.list$z.time[ ,d,s] * month[m])
      PHI[ ,d,s,2,m] <- ilogit(res$sims.list$logit.phi0[ ,d,s] +
        res$sims.list$beta.temp[ ,d,s] * res$sims.list$z.temp[ ,d,s] * temp[m] +
        res$sims.list$beta.time[ ,d,s] * res$sims.list$z.time[ ,d,s] * month[m] +
        res$sims.list$beta.transloc[ ,d,s] * res$sims.list$z.transloc[ ,d,s])
    }#m
  }#ss
}#s

mean.PHI <- apply(PHI, c(2,3,4,5), mean)
upper.PHI <- apply(PHI, c(2,3,4,5), function(x)quantile(x,0.975))
lower.PHI <- apply(PHI, c(2,3,4,5), function(x)quantile(x,0.025))
dim(lower.PHI)

##----
pdf(file.path(analysisDir, modelName, "survivalProbabilities_all iterations.pdf"),
    width = 10, height = 10)
par(mfrow = c(2,2))

##-- Females protected plot ----
par(mar = c(0,4,6,0))
plot(1, type = "n", xlim = c(0,n.months+1), ylim = c(0, 1), axes = F,
     ylab = "Survival probability",
     xlab = "",
     main = "Females")
axis(1, at = seq(11,260,24), labels = seq(1999,2019,2))
axis(2, at = seq(0,1,0.2), labels = seq(0,1,0.2))
polygon(x = c(1:n.months,n.months:1),
        y = c(upper.PHI[2,1,2, ],rev(lower.PHI[2,1,2, ])),
        col = adjustcolor(myCols[2],alpha.f = 0.5), border = F)
points(1:n.months, mean.PHI[2,1,2, ], type = "l", lwd = 3, col = myCols[2])
polygon(x = c(1:n.months,n.months:1),
        y = c(upper.PHI[2,1,1, ],rev(lower.PHI[2,1,1, ])),
        col = adjustcolor(myCols[4],alpha.f = 0.5), border = F)
points(1:n.months, mean.PHI[2,1,1, ], type = "l", lwd = 3, col = myCols[4])

##-- Males protected plot ----
par(mar = c(0,2,6,2))
plot(1, type = "n", xlim = c(0,n.months+1), ylim = c(0, 1), axes = F,
     ylab = "",
     xlab = "",
     main = "Males")
axis(1, at = seq(11,260,24), labels = seq(1999,2019,2))
axis(2, at = seq(0,1,0.2), labels = seq(0,1,0.2))
polygon(x = c(1:n.months,n.months:1),
        y = c(upper.PHI[2,2,2, ],rev(lower.PHI[2,2,2, ])),
        col = adjustcolor(myCols[2],alpha.f = 0.5), border = F)
points(1:n.months, mean.PHI[2,2,2, ], type = "l", lwd = 3, col = myCols[2])
polygon(x = c(1:n.months,n.months:1),
        y = c(upper.PHI[2,2,1, ],rev(lower.PHI[2,2,1, ])),
        col = adjustcolor(myCols[4],alpha.f = 0.5), border = F)
points(1:n.months, mean.PHI[2,2,1, ], type = "l", lwd = 3, col = myCols[4])


##-- Females fragmented plot ----
par(mar = c(4,4,2,0))
plot(1, type = "n", xlim = c(0,n.months+1), ylim = c(0, 1), axes = F,
     ylab = "Survival probability",
     xlab = "",
     main = "")
axis(1, at = seq(11,260,24), labels = seq(1999,2019,2))
axis(2, at = seq(0,1,0.2), labels = seq(0,1,0.2))
polygon(x = c(1:n.months,n.months:1),
        y = c(upper.PHI[1,1,2, ],rev(lower.PHI[1,1,2, ])),
        col = adjustcolor(myCols[2],alpha.f = 0.5), border = F)
points(1:n.months, mean.PHI[1,1,2, ], type = "l", lwd = 3, col = myCols[2])
polygon(x = c(1:n.months,n.months:1),
        y = c(upper.PHI[1,1,1, ],rev(lower.PHI[1,1,1, ])),
        col = adjustcolor(myCols[4],alpha.f = 0.5), border = F)
points(1:n.months, mean.PHI[1,1,1, ], type = "l", lwd = 3, col = myCols[4])

##-- Males fragmented plot ----
par(mar = c(4,2,2,2))
plot(1, type = "n", xlim = c(0,n.months+1), ylim = c(0, 1), axes = F,
     ylab = "",
     xlab = "",
     main = "")
axis(1, at = seq(11,260,24), labels = seq(1999,2019,2))
axis(2, at = seq(0,1,0.2), labels = seq(0,1,0.2))
polygon(x = c(1:n.months,n.months:1),
        y = c(upper.PHI[1,2,2, ],rev(lower.PHI[1,2,2, ])),
        col = adjustcolor(myCols[2],alpha.f = 0.5), border = F)
points(1:n.months, mean.PHI[1,2,2, ], type = "l", lwd = 3, col = myCols[2])
polygon(x = c(1:n.months,n.months:1),
        y = c(upper.PHI[1,2,1, ],rev(lower.PHI[1,2,1, ])),
        col = adjustcolor(myCols[4],alpha.f = 0.5), border = F)
points(1:n.months, mean.PHI[1,2,1, ], type = "l", lwd = 3, col = myCols[4])

dev.off()





## ------   5. DETECTION PROBABILITY -----
load(file = file.path( analysisDir,
                       modelName,
                       "inFiles",
                       paste0(modelName,"1.RData")))
## Individual detection
temp2 <- seq(min(nimConstants$temp2, na.rm = T), 
             max(nimConstants$temp2, na.rm = T),
             0.05)
n.sites <- nimConstants$n.sites
n.iter <- dim(nimMat)[1]
P <- array(NA, c(n.iter,n.sites,2,length(temp2)))
for(t in 1:length(temp2)){
  for(s in 1:2){
    for(f in 1:n.sites){
      P[ ,f,s,t] <- 1-exp(-exp( res$sims.list$lambda0[ ,f,s] +
                                  res$sims.list$gamma.temp[ ,f,s] *
                                  res$sims.list$z.det[ ,f,s] *
                                  temp2[t]) * 4) 
    }#f
  }#s
}#t
mean.P <- apply(P, c(2,3,4), function(x)mean(x, na.rm = T))
upper.P <- apply(P, c(2,3,4), function(x)quantile(x, 0.975, na.rm = T))
lower.P <- apply(P, c(2,3,4), function(x)quantile(x, 0.025, na.rm = T))

pdf(file.path(analysisDir, modelName, "detectionProbabilities.pdf"),
    width = 8, height = 12)
par(mfrow = c(4,2))
sex <- c("females", "males")
for(f in 1:n.sites){
  plot(1, type = "n",
       xlim = c(min(temp2)-0.05,max(temp2)+0.05),
       ylim = c(0,1), axes = F,
       ylab = "Detection prob.",
       xlab = "Temperature",
       main = fragments[f])
  axis(1,
       at = seq(min(temp2)-0.05,max(temp2),0.5),
       labels = seq(min(temp2)-0.05,max(temp2),0.5))
  axis(2, at = seq(0,1,0.2), labels = seq(0,1,0.2))
  legend( x = 26, y = 1,
          legend = c("females", "males"),
          bty = "n",
          fill = myCols[c(1,3)])
  
  polygon(x = c(temp2,rev(temp2)),
          y = c(upper.P[f,1, ],rev(lower.P[f,1, ])),
          col = adjustcolor(myCols[1],alpha.f = 0.5), border = F)
  polygon(x = c(temp2,rev(temp2)),
          y = c(upper.P[f,2, ],rev(lower.P[f,2, ])),
          col = adjustcolor(myCols[3],alpha.f = 0.5), border = F)
  
  points(temp2, mean.P[f,1, ], type = "l", lwd = 3, col = myCols[1])
  points(temp2, mean.P[f,2, ], type = "l", lwd = 3, col = myCols[3])
}#d
dev.off()



## ------   2. TABLES -----
##---- Create empty table
survival <-  matrix(NA, nrow = 5, ncol = 4)
rownames(survival) <- c("","phi0", "beta.time", "beta.temp", "beta.transloc")
colnames(survival) <- c("Protected", "Protected", "Degraded", "Degraded")
survival[1, ] <- c("F", "M", "F", "M")
survival[2, ] <- paste0( round(plogis(res$mean$logit.phi0), 2),
                         " (",  round(plogis(res$q2.5$logit.phi0),2),
                         "-", round(plogis(res$q97.5$logit.phi0),2), ")")
survival[3, ] <- paste0( round(res$mean$beta.time, 2),
                         " (",  round(res$q2.5$beta.time,2),
                         "-", round(res$q97.5$beta.time,2), ")")
survival[4, ] <- paste0( round(plogis(res$mean$beta.temp), 2),
                         " (",  round(res$q2.5$beta.temp,2),
                         "-", round(res$q97.5$beta.temp,2), ")")
survival[5, ] <- paste0( round(res$mean$beta.transloc, 2),
                         " (",  round(res$q2.5$beta.transloc,2),
                         "-", round(res$q97.5$beta.transloc,2), ")")

addtorow <- list()
addtorow$pos <- list(0, 0, c(2,4))
addtorow$command <- c(paste0(paste0('& \\multicolumn{2}{c}{',
                                    sort(unique(colnames(survival))),
                                    '}', collapse = ''), '\\\\'),
                      "\\rowcolor[gray]{.85}",
                      "\\rowcolor[gray]{.95}")

print(xtable(survival, type = "latex",
             align = paste(c("l",rep("c",ncol(survival))),collapse = "")),
      floating = FALSE,include.colnames=F, include.rownames = TRUE,
      add.to.row = addtorow,
      file = file.path(analysisDir, modelName, "Params.tex"))

##----- Export as .csv
write.csv( survival,
           file =  file.path(analysisDir, modelName, "Params.csv"))




## -----------------------------------------------------------------------------