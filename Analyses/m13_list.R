################################################
######DATA STRUCTURE MICROCEBUS FRAGMENT M13####
################################################
library(dplyr)
library(tidyr)
library(magrittr)
library(purrr)
library(plyr)


#####Data preparation for MCMC model fragment M13
setwd("C:/Users/anvargas/Dropbox/Mouse lemur CMR data/02_Data")

#--Read capture long
capture_his<-read.csv("cmr_mhc_input_data_microcebus.csv")
capture_his$X<-NULL
capture_his<-capture_his[!duplicated(capture_his[ ]),] %>%droplevels()
nrow(capture_his)

#--Fragment M13
m13<-capture_his[capture_his$site=="M13",]

ids<-unique(m13$transponder)
years<-sort(as.numeric(c(unique(m13$year),"2004"))) #no sampling in 2004
n.years.m<-length(years)

#--Session only of captured individuals
sessionIndex.c<-sort(as.numeric(c(unique(m13$session))))
# length(sessionIndex)
# months<-sort(as.numeric(c(unique(m13$month_start))))
n.months.m <- 12*n.years.m


#--Including time intervals
tinterm13<-read.csv("time_interval_session_frag.csv")
tinm13<-tinterm13[tinterm13$Site=="M13",]
names(tinm13)

tinm13<-tinm13[,c(2,6)]
uniq13<-tinm13[!duplicated(tinm13[1:2]),] %>%droplevels()
sessions.m<-uniq13$new_session_av
nrow(uniq13)

#--Compare the session with and without individuals
uniq13$n0ind<- ifelse(match(sessions.m, sessionIndex.c), "Yes", "No")

#--Exctract sessions without individuals to create dummy individuals
sesionsn0oind<-uniq13[is.na(uniq13$n0ind),]


#--First create dummy individuals for session without captures
dummy <- data.frame (transponder  = c("/a/", "/b/"),
                  month_start = c(11,2),
                  session=sesionsn0oind[,1],
                  season=c("wet","wet"))

datadum<- rbind.fill(list(m13, dummy))
datadum<-datadum[order(datadum$session),]

season.m<-data.frame(session=datadum$session,
                     season=datadum$season, 
                     month_start=datadum$month_start)
season.m<-season.m[order(season.m$session), ]
season.m<-season.m[!duplicated(season.m[ ]),] %>%droplevels()

#--Session index 
sessionIndex.m<-sort(as.numeric(c(unique(datadum$session))))
n.sessions.m<- length(sessionIndex.m)

#testdum[testdum$transponder== "/a/",] #check  if the dummy individuals are present

#--Dataset with dummy individuals
m13.dum<-datadum
class(m13.dum)
class(m13.dum$transponder)


#--Capture history
#--Create the matrix of capture history
m13.det <- table(m13.dum$transponder, m13.dum$session)

ch13.m<- ifelse(m13.det>0,1,0)
ch13.m<-as.matrix(ch13.m)
dim(ch13.m)

#--Set to 0 dummy individuals
ch13.m["/a/", "108"]<-0
ch13.m["/b/", "194"]<-0

#--Indexing each transponder for sex
sex<-unique(m13[c("transponder", "Sexe")])
sex$Sexe<- ifelse(sex$Sexe== "f", 1,0) 
listm13s<-sex[,] %>% purrr::transpose()
listm13s<-as.array(listm13s)
dim(listm13s)

#--Index for each transpoder with each session, season and age
changingt<-unique(m13.dum[c("year","month_start","transponder", "session","season", "age_estimation_2", "effort_days_sess")]) ## all captures
changingt$season<- ifelse(changingt$season== "wet", 1, 0) #season
changingt$age_estimation_2<- ifelse(changingt$age_estimation_2== "Y", 0,
                                    ifelse(changingt$age_estimation_2== "NI", 1, 2))
names(changingt)
changingt<-changingt%>% left_join(uniq13, by=c("session"="new_session_av"))

#--Number of month from the begining of the study 
yearsi<-sort(rep(c(2000:2019),12))
month.m<-rep(c(1:12), 20)
d.months<-as.data.frame(cbind(yearsi, month.m))[-c(1:10),] #here the month when the study start is 11, it's manually excluded the previous month 
d.months$month.in<- 1:nrow(d.months)
d.months$month.m<-as.character(d.months$month.m)

months.ord$month_start<-as.character(months.ord$month_start)
months.ord<-as.data.frame(changingt[,c("year","month_start","session")])
months.ord<-months.ord%>%left_join(d.months, by=c("year"="yearsi","month_start"="month.m"))
months.ord<- months.ord[!duplicated(months.ord[c("year", "month_start", "session", "month.in")]),]

changingt$month_start<-as.character(changingt$month_start)
changingt<-changingt%>% left_join(months.ord, by=c("session", "year", "month_start"))
head(changingt)
nrow(changingt)
#--list index by capture with month when iniciated, time interval, age(state), season, year
M13_chan<-changingt%>% filter(!(transponder== "/a/"| transponder=="/b/" | transponder=="/c/")) %>% purrr::transpose()
M13_chan<-as.array(M13_chan)
dim(M13_chan)

## Keep only detected individuals
detected.m <- apply(ch13.m, 1, function(x)any(x == 1))
nrow(ch13.m)

y.detected.m <- ch13.m[detected.m, ]
nrow(y.detected.m)

## Extract the first detection occasion per individual
f.m <- apply(y.detected.m, 1, function(x)min(which(x == 1)))

## Remove individuals detected for the first time on the last session
y.detected.m <- y.detected.m[which(f.m != dim(y.detected.m)[2]), ]
f.m <- f.m[which(f.m != dim(y.detected.m)[2])]



start.int.m <- end.int.m <- dt1 <- dt2 <- NULL
for(t in 1:(length(sessionIndex.m)-1)){
  start.int.m[t] <- sessionIndex.m[t]
  end.int.m[t] <- sessionIndex.m[t+1]-1
  # dt1[t] <- sum(season[start.int[t]:end.int[t]] == 1)
  # dt2[t] <- sum(season[start.int[t]:end.int[t]] == 2)
}#t


########################################################
##############SIMULATED MODELS PIERRE###################
########################################################



##Pierre script

library(nimble)
library(basicMCMCplots)
library(coda)

#Set simulation characteristics
N0 <- 50                                           ## Initial population size       
n.years <- 10                                       ## Number of years simulated
season <- rep(c(1,1,1,1,1,1,1,1,2,2,2,2), n.years) ## Wet/dry seasons index
phi0 <- c(0.98,0.87)                                ## Time- and season-specific survival
beta <- c(-0.02,-0.01)
r <- 0.08                                          ## Reproduction probability
litterSize <- c(0.5,0.35,0.15)                     ## Vector of litter size probabilities

#parameters for the sampling sessions (% of months sampled and mean duration of each sampling session)
meanDuration <- 4 
lambdaDet <- 0.45
propSession <- 0.46

#Generate CMr dataset:start simulating the CMR dataset and to initialize simulation objects.

POP <- list()                                       ## List of population composition 
POP[[1]] <- rep(1, N0)                              ## Initial population composition
N <- vector()                                       ## List of population size
N[1] <- N0                                          ## Initial population size

#simulate individual- and time-specific survival, recruitment, and litter size. From that we derive monthly population sizes.

n.months <- 12*n.years
phi <- NULL 
start.int <- end.int <- NULL 
for(t in 2:n.months){
  phi[t-1] <- 1/(1+exp(-(logit(phi0[season[t]]) + beta[season[t]] * t)))
  Z <- B <- rep(0,length(POP[[t-1]]))        
  for(i in 1:length(POP[[t-1]])){       
    if(POP[[t-1]][i] == 1){         ## Sample states only if individual is alive
      Z[i] <- rbinom(1, 1, phi[t-1])## Sample individual survival
      repro <- rbinom(1, 1, r*Z[i]) ## Sample individual repro (only if individual survived)
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

#Here, we re-organize individual life-histories in a matrix.

z <- matrix(0, length(POP[[n.months]]), n.months)
for (t in 1:n.months){
  z[1:length(POP[[t]]), t] <- POP[[t]]
}#t

#Then, we sample in which months the trapping sessions occured and how long each one lasts. (First and last sessions must be sampling sessions for the model to make sense)

sessionIndex <- c( 1,
                   sample( x = 2:(n.months-1),
                           size = round(propSession*(n.months-2))),
                   n.months)
sessionIndex <- sessionIndex[order(sessionIndex)]
n.sessions <- length(sessionIndex)

#Before sampling individual observations for each session.

sessionDuration <- rpois(n = n.sessions, lambda = meanDuration)
p <- 1 - exp(-lambdaDet * sessionDuration)

y <- matrix(0, dim(z)[1], n.sessions)
for(i in 1:dim(z)[1]){
  for (t in 1:n.sessions){
    y[i,t] <- rbinom(size = 1, n = 1, prob = p[t] * z[i,t])
  }#t
}#i

#To match real-life CMR datasets, we need to filter out individuals that were never captured, as well as individuals that were detected for the first time on the last sampling session.

## Keep only detected individuals
detected <- apply(y, 1, function(x)any(x == 1))
y.detected <- y[detected, ]


## Extract the first detection occasion per individual
f <- apply(y.detected, 1, function(x)min(which(x == 1)))

## Remove individuals detected for the first time on the last session
y.detected <- y.detected[which(f != dim(y.detected)[2]), ]
f <- f[which(f != dim(y.detected)[2])]

#start.int <- end.int <- dt1 <- dt2 <- NULL
for(t in 1:(length(sessionIndex)-1)){
  start.int[t] <- sessionIndex[t]
  end.int[t] <- sessionIndex[t+1]-1
  # dt1[t] <- sum(season[start.int[t]:end.int[t]] == 1)
  # dt2[t] <- sum(season[start.int[t]:end.int[t]] == 2)
}#t

#Define NIMBLE Model Structure
#Here, we define the NIMBLE model.

modelCode <- nimbleCode({
  ##---------------------
  ## DEMOGRAPHIC PROCESS 
  phi0[1] ~ dunif(0,1)
  phi0[2] ~ dunif(0,1)
  logit.phi0[1] <- logit(phi0[1])
  logit.phi0[2] <- logit(phi0[2])
  beta[1] ~ dnorm(0,0.01)
  beta[2] ~ dnorm(0,0.01)
  for(m in 1:n.months){
    # PHI[m] <- phi0[season[m]] + beta[season[m]] * m
    logit(PHI[m]) <- logit.phi0[season[m]] + beta[season[m]] * m
  }# months
  
  for(t in 1:n.intervals){
    phi[t] <- prod(PHI[start.int[t]:end.int[t]])
    # phi[t] <- getPhi(phiVec = PHI[start.int[t]:end.int[t]])
    # phi[t] <- pow(phi0[1], dt1[t]) * pow(phi0[2], dt2[t]) # Session-specific survival probabilities
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
  
  ##------------------------------------------------------------------------------------------------------------------
})

# Format the data for NIMBLE
# Then, we organize the simulated data in a format useable by NIMBLE; i.e. we create objects constants, data, and inits for later use in the function nimbleModel.

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

# 
# Create a NIMBLE model object
# Now, we can create the nimble model object, using the model structure defined in code, and the constants, data, and initial values.

Rmodel <- nimbleModel( code = modelCode,
                       constants = nimConstants,
                       data = nimData,
                       inits = nimInits,
                       calculate = F)

Rmodel$calculate()

# Configure and Build MCMC objects
# We configure an MCMC algorithm to the Rmodel model object.
# 
# We assign MCMC monitors to ϕ, γ, lambda, and ψ.

conf <- configureMCMC(Rmodel, monitors = c("phi0", "beta", "lambda"), print = FALSE)
Rmcmc <- buildMCMC(conf)

# Compile and Run MCMC
# Finally, we compile both the model and MCMC objects and execute the compiled MCMC for 5 000 iterations and 3 chains.

Cmodel <- compileNimble(Rmodel)


Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

MCMC_runtime <- system.time(
  samples <- runMCMC( Cmcmc,
                      niter = 5000,
                      nburnin = 1000,
                      nchains = 3,
                      samplesAsCodaMCMC = T))

summary(samples)


##############Model with data M13####

# Format the data for NIMBLE
# Then, we organize the simulated data in a format useable by NIMBLE; i.e. we create objects constants, data, and inits for later use in the function nimbleModel.
season <- rep(c(1,1,1,1,1,1,1,1,2,2,2,2), n.years.m) ## Wet/dry seasons index
nimData <- list( y = y.detected.m,
                 sessionDuration = sessionDuration )

nimConstants <- list( n.individuals = dim(y.detected.m)[1],
                      n.months = n.months,
                      n.intervals = dim(y.detected.m)[2]-1,
                      n.sessions = dim(y.detected.m)[2],
                      season = season,
                      f = f.m,
                      start.int.m = start.int.m,
                      end.int.m = end.int.m)


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

# 
# Create a NIMBLE model object
# Now, we can create the nimble model object, using the model structure defined in code, and the constants, data, and initial values.

Rmodel <- nimbleModel( code = modelCode,
                       constants = nimConstants,
                       data = nimData,
                       inits = nimInits,
                       calculate = F)

Rmodel$calculate()

# Configure and Build MCMC objects
# We configure an MCMC algorithm to the Rmodel model object.
# 
# We assign MCMC monitors to ϕ, γ, lambda, and ψ.

conf <- configureMCMC(Rmodel, monitors = c("phi0", "beta", "lambda"), print = FALSE)
Rmcmc <- buildMCMC(conf)

# Compile and Run MCMC
# Finally, we compile both the model and MCMC objects and execute the compiled MCMC for 5 000 iterations and 3 chains.

Cmodel <- compileNimble(Rmodel)


Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

MCMC_runtime <- system.time(
  samples <- runMCMC( Cmcmc,
                      niter = 3000,
                      nburnin = 1000,
                      nchains = 3,
                      samplesAsCodaMCMC = T))





