########################################################
##### ------------ MANDENA MICROCEBE ------------- #####
##### SIMULTANEOUS ANALYSIS for MULTIPLE FRAGMENTS #####
########################################################
rm(list = ls())

## ------ LIBRARIES ------
library(dplyr)
#library(tidyr)
library(magrittr)
#library(purrr)
library(plyr)
library(coda)
library(nimble)
library(xtable)


## ------ WORKING DIRECTORIES ------
source("workingDirectories.R")
modelName <- "m_phi[status_sex_temp_time_transloc]_p[site_sex_temp]_RJ"
source("wildMap.R")
myCols <- wildMap(4)
source(file.path(analysisDir, "ProcessCodaOutput_v3.R"))
source(file.path(analysisDir, "PlotViolins.R"))


## -----------------------------------------------------------------------------
## ------ MAKE PLOTS -----
## ------   1. LOAD OUTPUTS -----
outputs <- list.files(file.path(analysisDir, modelName, "outFiles"))
parm.index <- c(1,2,9,10,17,18,25,26,33,34,41,42,49:82,89,90)
nimOutput2 <- mcmc.list()
for(i in 1:length(outputs)){
  load(file.path(analysisDir,modelName,"outFiles",outputs[i]))
  nimOutput2[[i]] <- as.mcmc(nimOutput[, parm.index])
}
nimOutput <- nimOutput2

# pdf(file = file.path(analysisDir, modelName, "Traceplots.pdf"))
# plot(nimOutput)
# dev.off()


## ------   2. PRINT TABLES -----

results <- ProcessCodaOutput_v3(x = nimOutput) 

##---- Create empty table
survival <-  matrix(NA, nrow = 5, ncol = 4)
rownames(survival) <- c("","phi0", "beta.time", "beta.temp", "beta.transloc")
colnames(survival) <- c("Protected", "Protected", "Degraded", "Degraded")
survival[1, ] <- c("F", "M", "F", "M")
survival[2, ] <- paste0( round(plogis(results$mean$logit.phi0), 2),
                         " (",  round(plogis(results$q2.5$logit.phi0),2),
                         "-", round(plogis(results$q97.5$logit.phi0),2), ")")
survival[3, ] <- paste0( round(plogis(results$mean$beta.time), 2),
                         " (",  round(plogis(results$q2.5$beta.time),2),
                         "-", round(plogis(results$q97.5$beta.time),2), ")")
survival[4, ] <- paste0( round(plogis(results$mean$beta.temp), 2),
                         " (",  round(plogis(results$q2.5$beta.temp),2),
                         "-", round(plogis(results$q97.5$beta.temp),2), ")")
survival[5, ] <- paste0( round(plogis(results$mean$beta.transloc), 2),
                         " (",  round(plogis(results$q2.5$beta.transloc),2),
                         "-", round(plogis(results$q97.5$beta.transloc),2), ")")

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



## ------   2. PLOT SURVIVAL -----
load(file = file.path(analysisDir,
                      modelName,
                      "inFiles",
                      paste0(modelName,"1.RData")))
n.months <- nimConstants$n.months
n.sites <- nimConstants$n.sites
temp <- nimConstants$temp
month <- nimConstants$month
nimMat <- do.call(rbind, nimOutput)
n.iter <- dim(nimMat)[1]
beta.temp <- array(nimMat[ ,grep("beta.temp",dimnames(nimMat)[[2]])], c(n.iter,2,2))
beta.time <- array(nimMat[ ,grep("beta.time",dimnames(nimMat)[[2]])], c(n.iter,2,2))
beta.transloc <- array(nimMat[ ,grep("beta.transloc",dimnames(nimMat)[[2]])], c(n.iter,2,2))
logit.phi0 <- array(nimMat[ ,grep("logit.phi0",dimnames(nimMat)[[2]])], c(n.iter,8,2))
PHI <- array(NA, c(n.iter,2,2,2,n.months))
dimnames(PHI) <- list( iteration = 1:n.iter,
                       status = c("protected","degraded"),
                       sex = c("female","male"),
                       translocated = c("no","yes"),
                       month = 1:n.months)
for(f in 1:2){
  for(m in 1:n.months){
    for(s in 1:2){
      for(t in 1:2){
        PHI[ ,f,s,t,m] <- ilogit(logit.phi0[ ,f,s] +
                                   beta.temp[ ,f,s] * temp[m] +
                                   beta.time[ ,f,s] * month[m] +
                                   beta.transloc[ ,f,s] * (t-1))
      }#t
    }#s
  }#m
}# f

##-- Identify first month w/ translocation
test <- cbind.data.frame( 
  transloc = nimConstants$transloc,
  sex = nimConstants$sex,
  status = nimConstants$status,
  startM = nimConstants$start.int[ ,1]) %>%
  filter(., transloc == 2)

start <- matrix(NA,2,2)
start[1,1] <- min(test$startM[test$status == 1 & test$sex == 1])
start[1,2] <- min(test$startM[test$status == 1 & test$sex == 2])
start[2,1] <- min(test$startM[test$status == 2 & test$sex == 1])
start[2,2] <- min(test$startM[test$status == 2 & test$sex == 2])


##-- Plot Survival
mean.PHI <- apply(PHI, c(2,3,4,5), mean)
upper.PHI <- apply(PHI, c(2,3,4,5), function(x)quantile(x,0.975))
lower.PHI <- apply(PHI, c(2,3,4,5), function(x)quantile(x,0.025))
pdf(file.path(analysisDir, modelName, "survivalProbabilities.pdf"),
    width = 10, height = 12)
par(mfrow = c(2,2))
sex <- c("females", "males")
status <- c("degraded", "protected")
for(f in 1:2){
  for(s in 1:2){
    plot(1, type = "n",
         xlim = c(0,n.months+1),
         ylim = c(0, 1),
         axes = F,
         ylab = "Survival prob.",
         xlab = "Months",
         main = paste0(status[f],"-", sex[s]))
    axis(1, at = seq(0,250,50), labels = seq(0,250,50))
    axis(2, at = seq(0,1,0.2), labels = seq(0,1,0.2))
    legend( x = 1, y = 0.2,
            title = "Individual status:",
            legend = c("resident", "translocated"),
            bty = "n",
            fill = myCols[c(2,4)])
    
    polygon(x = c(1:n.months,n.months:1),
            y = c(upper.PHI[f,s,1, ],rev(lower.PHI[f,s,1, ])),
            col = adjustcolor(myCols[2],alpha.f = 0.5), border = F)
    points(1:n.months, mean.PHI[f,s,1,1:n.months], type = "l", lwd = 2, col = myCols[2])
    
    polygon(x = c(start[f,s]:n.months,n.months:start[f,s]),
            y = c( upper.PHI[f,s,2,start[f,s]:n.months],
                  rev(lower.PHI[f,s,2,start[f,s]:n.months])),
            col = adjustcolor(myCols[4],alpha.f = 0.5), border = F)
    points(start[f,s]:n.months, mean.PHI[f,s,2,start[f,s]:n.months], type = "l", lwd = 2, col = myCols[4])
  }#s
}#f
graphics.off()




## ------   3. PLOT DETECTION PROBABILITY -----
load(file = file.path( analysisDir,
                       modelName,
                       "inFiles",
                       paste0(modelName,"1.RData")))
temp2 <- seq(min(nimConstants$temp2, na.rm = T), 
             max(nimConstants$temp2, na.rm = T),
             0.05)
n.sites <- nimConstants$n.sites
nimMat <- do.call(rbind, nimOutput)
n.iter <- dim(nimMat)[1]
gamma.temp <- array(nimMat[ ,grep("gamma.temp",dimnames(nimMat)[[2]])], c(n.iter,8,2))
lambda0 <- array(nimMat[ ,grep("lambda0",dimnames(nimMat)[[2]])], c(n.iter,8,2))
P <- array(NA, c(n.iter,8,2,length(temp2)))
for(t in 1:length(temp2)){
  for(s in 1:2){
    for(f in 1:8){
      P[ ,f,s,t] <- 1-exp(-exp(lambda0[ ,f,s] + gamma.temp[ ,f,s] * temp2[t]) * 4) 
    }#f
  }#s
}#t

mean.P <- apply(P, c(2,3,4), function(x)mean(x, na.rm = T))
upper.P <- apply(P, c(2,3,4), function(x)quantile(x, 0.975, na.rm = T))
lower.P <- apply(P, c(2,3,4), function(x)quantile(x, 0.025, na.rm = T))

pdf(file.path(analysisDir, modelName, "detectionProbabilities.pdf"),
    width = 10, height = 12)
par(mfrow = c(4,2))
sex <- c("females", "males")
for(f in 1:n.sites){
  plot(1, type = "n",
       xlim = c(min(temp2)-0.4,max(temp2)+0.05),
       ylim = c(0,1), axes = F,
       ylab = "Detection prob.",
       xlab = "Temperature",
       main = fragments[f])
  axis(1,
       at = seq(min(temp2)-0.4,max(temp2),0.5),
       labels = seq(min(temp2)-0.4,max(temp2),0.5))
  axis(2, at = seq(0,1,0.2), labels = seq(0,1,0.2))
  legend( x = 26, y = 1,
          legend = c("females", "males"),
          bty = "n",
          fill = myCols[c(2,4)])
  
  polygon(x = c(temp2,rev(temp2)),
          y = c(upper.P[f,1, ],rev(lower.P[f,1, ])),
          col = adjustcolor(myCols[2],alpha.f = 0.5), border = F)
  points(temp2, mean.P[f,1, ], type = "l", lwd = 2, col = myCols[2])
  
  polygon(x = c(temp2,rev(temp2)),
          y = c(upper.P[f,2, ],rev(lower.P[f,2, ])),
          col = adjustcolor(myCols[4],alpha.f = 0.5), border = F)
  points(temp2, mean.P[f,2, ], type = "l", lwd = 2, col = myCols[4])
}#d
dev.off()


## -----------------------------------------------------------------------------