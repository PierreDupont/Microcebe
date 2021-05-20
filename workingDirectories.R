if(Sys.info()['user'] == 'pidu') { ## Pierre
  dataDir <- 'C:/Users/pidu/Dropbox (Personal)/Mouse lemur CMR data/02_Data'  
  analysisDir <- 'C:/Users/pidu/Dropbox (Personal)/Mouse lemur CMR data/03_Analysis'
  simDir <- 'C:/Users/pidu/Dropbox (Personal)/Mouse lemur CMR data/04_Simulation'
  figDir <- 'C:/Users/pidu/Dropbox (Personal)/Mouse lemur CMR data/05_Paper/Figures'            
  tableDir <- 'C:/Users/pidu/Dropbox (Personal)/Mouse lemur CMR data/05_Paper/Tables'            
} else if(Sys.info()['user'] == 'anva') {## Andrea
  dataDir <- 'C:/Users/pidu/Dropbox (Personal)/Mouse lemur CMR data/02_Data'  
  analysisDir <- 'C:/Users/pidu/Dropbox (Personal)/Mouse lemur CMR data/03_Analysis'
  simDir <- 'C:/Users/pidu/Dropbox (Personal)/Mouse lemur CMR data/04_Simulation'
  figDir <- 'C:/Users/pidu/Dropbox (Personal)/Mouse lemur CMR data/05_Paper/Figures'            
  tableDir <- 'C:/Users/pidu/Dropbox (Personal)/Mouse lemur CMR data/05_Paper/Tables'    
} else stop('unknown user')
