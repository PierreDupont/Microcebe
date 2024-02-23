##----------------------------------------##
## AV_ Microcebus summary data for manuscript ##
##----------------------------------------##

library(ggplot2)
library(dplyr)
library(wesanderson)
library(dplyr)
library(tidyr)
library(magrittr)
library(purrr)
library(plyr)

source("C:/Users/anvargas/EmptyForests/temp/AFVV/WorkingDirectories.R")

source("C:/Users/anvargas/Microcebe/wildMap.R")
R.utils::sourceDirectory("C:/Users/anvargas/Source",modifiedOnly = FALSE)

modelName <- "m_phi[status_sex_temp_time_transloc]_p[site_sex_temp]_RJ"

#modelName<-"NimbleOutFORm_phi[status_sex_temp_time_transloc]_p[site_sex_temp]_RJ8.RData"

## ------ SOURCE THE REQUIRED FUNCTIONS ------
#source("C:/Users/anvargas/Source")#,modifiedOnly = FALSE)

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

names(capture_data)
unique(capture_data$site)

n.trans.site<-capture_data %>% select(transponder, residency, site) %>%  
  #filter(site %in% c("M15A", "M15B", "M15C", "M15D") )%>%
  filter(site %in% c("M13") )%>%
  #filter(Sexe== "m" ) %>% 
  filter(residency== "translocated" ) %>%  
  distinct() #%>% group_by(site) %>%  
  #dplyr::summarise (n=n())





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

## Column for months
colnames(weather)[3]<-"Temperature"
weather<-weather%>%mutate( month.n= month.abb[MONTH] ) 

## Colors
temp_col<-wes_palette("Darjeeling1")[3]
fill_col<-wes_palette("Royal1")[3]

## Plot
temp_p<-weather %>% ggplot(aes (y= Temperature, x= month.n)) + geom_boxplot(na.rm = TRUE, col= temp_col, fill=fill_col)+
  scale_x_discrete(limits = month.abb) + 
  theme(panel.background = element_blank())+ylab("Temperature (\u00B0C)")+ xlab("Month")


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
  
}
data_list[1]

## Check the total number of captures, and individuals, by sex.
names(capture_data)

# Captures
n.capt <- capture_data %>% select(transponder, Sexe) %>% 
                              #distinct() %>% 
                              count("Sexe")

# Individuals 
n.ind <- capture_data %>% select(transponder, Sexe) %>% 
  distinct() %>% 
  count("Sexe")


## ------ III. PROCESS OUTPUTS -----
load(file.path(analysisDir, modelName, "inFiles", paste0(modelName,"1.RData")))
load(file.path(analysisDir, modelName, paste0("outFiles/output.RData")))
load(file.path(analysisDir, modelName, paste0("outFiles/processed_output.RData")))
load(paste0(analysisDir,"/", modelName, "/outFiles/", "NimbleOutFOR",modelName, "8.RData"))

class(nimOutput)
res <- ProcessCodaOutput(nimOutput)

res$mean$beta.time
res$q2.5$beta.time
res$q97.5$beta.time


res$q97.5$beta.transloc
res$q2.5$beta.transloc
res$mean$beta.transloc


res$mean$beta.temp
res$q2.5$beta.temp
res$q97.5$beta.temp

res$Rhat
res$sims.list$beta.time
res$sims.list$beta.time


dim(survival)
survival
# first dimension : degraded (1), protected (2)
# second dimension: female (1), male (2)
# third dimension: resident (1), translocated (2) 
# fourth dimension:  month 1: length months

## ------   2. SURVIVAL PLOTS (all iterations) -----
n.months <- nimConstants$n.months
temp <- nimConstants$temp
month <- nimConstants$month
PHI <- array(NA, c(dim(res$sims.list$beta.temp)[1],2,2,2,n.months))
for(d in 1:2){ #is this resident vs. translocated
  for(s in 1:2){ # 1 females, 2 males
    for(m in 1:n.months){ # number of months 
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


# first dimension : degraded (1), protected (2)
# second dimension: female (1), male (2)
# third dimension: resident (1), translocated (2) 
# fourth dimension:  month 1: length months


## Extract first and last PHI of females in degraded 
phi.f.d.1m<-  mean.PHI[1,1,1,1]
phi.uci.f.d.1m<- upper.PHI[1,1,1,1]
phi.lci.f.d.1m<-lower.PHI[1,1,1,1]

## Extract first and last PHI of females in protected
phi.f.p.1m<-  mean.PHI[2,1,1,1]
phi.uci.f.p.1m<- upper.PHI[2,1,1,1]
phi.lci.f.p.1m<-lower.PHI[2,1,1,1]


## Last PHI of females in protected
phi.f.p.252<-  mean.PHI[2,1,1,n.months]
phi.uci.f.p.252<- upper.PHI[2,1,1,n.months]
phi.lci.f.p.252<-lower.PHI[2,1,1,n.months]


## First PHI of females in degraded translocated 
phi.f.d.1m<-  mean.PHI[1,1,2,1]
phi.uci.f.d.1m<- upper.PHI[1,1,2,1]
phi.lci.f.d.1m<-lower.PHI[1,1,2,1]


## First PHI of females in protected translocated 
phi.f.d.1m<-  mean.PHI[2,1,2,1]
phi.uci.f.d.1m<- upper.PHI[2,1,2,1]
phi.lci.f.d.1m<-lower.PHI[2,1,2,1]



mean.PHI[1,1,2, ]

## Extract first and last PHI of males in degraded 
phi.m.d.1m<-  mean.PHI[1,2,1,1]
phi.uci.m.d.1m<- upper.PHI[1,2,1,1]
phi.lci.m.d.1m<-lower.PHI[1,2,1,1]

## Extract first and last PHI of males in protected
phi.m.p.1m<-  mean.PHI[2,2,1,1]
phi.uci.m.p.1m<- upper.PHI[2,2,1,1]
phi.lci.m.p.1m<-lower.PHI[2,2,1,1]



## PHI of males in degraded translocated 
phi.m.d.1m<-  mean.PHI[1,2,2,1]
phi.uci.m.d.1m<- upper.PHI[1,2,2,1]
phi.lci.m.d.1m<-lower.PHI[1,2,2,1]


## First PHI of females in protected translocated 
phi.f.m.1m<-  mean.PHI[2,2,2,1]
phi.uci.f.m.1m<- upper.PHI[2,2,2,1]
phi.lci.f.m.1m<-lower.PHI[2,2,2,1]


mean.PHI[1,2,1,1] 
dim(mean.PHI)




