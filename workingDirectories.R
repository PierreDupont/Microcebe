if(Sys.info()['user'] == 'pidu') { ## Pierre
  dataDir <- 'C:/Users/pidu/Dropbox (Personal)/Mouse lemur CMR data/02_Data'  
  analysisDir <- 'C:/Users/pidu/Dropbox (Personal)/Mouse lemur CMR data/03_Analysis'
  simDir <- 'C:/Users/pidu/Dropbox (Personal)/Mouse lemur CMR data/04_Simulation'
  figDir <- 'C:/Users/pidu/Dropbox (Personal)/Mouse lemur CMR data/05_Paper/Figures'            
  tableDir <- 'C:/Users/pidu/Dropbox (Personal)/Mouse lemur CMR data/05_Paper/Tables'
  resultDir <- 'C:/Users/pidu/Dropbox (Personal)/Mouse lemur CMR data/06_Results'
} else if(Sys.info()['user'] == 'anvargas') {## Andrea
  dataDir <- 'C:/Users/anvargas/Dropbox/Mouse lemur CMR data/02_Data'  
  analysisDir <- 'C:/Users/anvargas/Dropbox/Mouse lemur CMR data/03_Analysis'
  simDir <- 'C:/Users/anvargas/Dropbox/Mouse lemur CMR data/04_Simulation'
  figDir <- 'C:/Users/anvargas/Dropbox/Mouse lemur CMR data/05_Paper/Figures'            
  tableDir <- 'C:/Users/anvargas/Dropbox/Mouse lemur CMR data/05_Paper/Tables'  
  resultDir <- 'C:/Users/anvargas/Dropbox/Mouse lemur CMR data/06_Results'
} else stop('unknown user')
