
COR_ANAND<-function(ProductionEnvironment,COR1,COR2){
  ## Original implemented by Anand et al. 2017
  # use the framework of Balakrishnan, Hanse and Labro (2011)
  # in Anand et al. Flowchart 5.1 (b)


  DISP1<-ProductionEnvironment[['DISP1']]
  DISP2<-ProductionEnvironment[['DISP2']]
  NUMB_RES<-ProductionEnvironment[['NUMB_RES']]
  NUMB_PRO<-ProductionEnvironment[['NUMB_PRO']]

  #### ---- STEP BASELINE CONSUMPTOION   ---- ####
  BASE<-rnorm(NUMB_PRO) #Anand Approach
  # BASE<-dqrnorm(NUMB_PRO)

  # GENERATE RANDOM NUMBERS FOR CORRLEATIONS
  # in Anand et al. Flowchart 5.1 (a)

  RES_CONS_PAT_PRE<-matrix(0,nrow = NUMB_PRO,ncol = NUMB_RES)
  RES_CONS_PAT_PRE[,1]<-BASE
  for (i in 2:NUMB_RES) {
    RES_CONS_PAT_PRE[,i]<-rnorm(NUMB_PRO)
  }


  RES_CONS_PAT<-matrix(0,nrow =NUMB_PRO ,ncol = NUMB_RES ) # PREALLOCAITON

  RES_CONS_PAT[,1]<-BASE

  # Generate Correlation of DISP1 rescources
  # in Anand et al. Flowchart 5.1 (b)
  sqrtConstant1<-sqrt(1-COR1^2)

  for(i in 2:DISP1){

    # Anand Version
    RES_CONS_PAT[,i]<-(COR1*BASE)+(sqrtConstant1*RES_CONS_PAT_PRE[,i])
  }

  sqrtConstant2<-sqrt(1-COR2^2)


  for(i in (DISP1+1):NUMB_RES){

    RES_CONS_PAT[,i]<-(COR2*BASE)+(sqrtConstant2*RES_CONS_PAT_PRE[,i])

  }

  return(RES_CONS_PAT)

}
