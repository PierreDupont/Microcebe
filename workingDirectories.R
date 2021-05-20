if(Sys.info()['user'] == 'pidu') {
  dataDir <- 'C:\Users\pidu\Dropbox (Personal)\Mouse lemur CMR data\02_Data'  
  figDir <- 'C:\Users\pidu\Dropbox (Personal)\Mouse lemur CMR data\02_Data'            ## Pierre
  dataDir <- 'C:\Users\pidu\Dropbox (Personal)\Mouse lemur CMR data\02_Data'            ## Pierre
  ## Pierre
} else if(Sys.info()['user'] == 'cymi') {
  dataDir <- 'C:/Personal_Cloud/OneDrive/Work/nimbleSCR/'   ## Andrea
} else stop('unknown user')
