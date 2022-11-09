measureMatrixHetero<-function(A,RC,DMD){
  #### FOr Testing Only ####
  # library(EAD)
  # data("csd_EAD")
  # A<-CSD_EAD$RES_CONS_PAT
  #### END For testing only ####


  CHECK<-list()
  NUMB_PRO<-NROW(A)
  NUMB_RES<-NCOL(A)
  #### 1. Density, e.g. see Balakrishnan et al. (2011) ####
  CHECK[['DENS']]<-length(which(A!=0))/(NUMB_PRO*NUMB_RES)


  #### 2. Range of Consumption ####
  CHECK[['RANGE_CON']]<-(as.numeric(quantile(A,.95))-as.numeric(quantile(A,.05)))/as.numeric(mean(A))

  DMD_mat <- matrix(DMD,
                    ncol = NCOL(A),
                    byrow = F,
                    nrow = NROW(A))
  TRC<-colSums(A * DMD_mat)
  RCU <- RC/TRC
  PC<-A %*% RCU
  CHECK[['RANGE_PC']]<-(quantile(PC,.95)-quantile(PC,.05))/mean(PC)
  CHECK[['DMD_T10']]<-measure_TOP10(DMD)
  CHECK[['RC_T10']]<-measure_TOP10(RC)
  return(CHECK)
}
