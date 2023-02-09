processVariety<-function(PrD,
                         PrC,
                         DMD,
                         includedProd=NULL,
                         prop_setupChange=c(0.5,0.5),
                         prop_setupTime=1,
                         prop_demand=0.5){
  require(dplyr)
  #### Input for Testing ####
  # PrD <- matrix(c(10,12,10,0,0,
  #               6,5,2,2,4),ncol = 2)
  # PrC <- c(120,130)
  # prop_setupChange <- c(0.5,0.5)
  # includedProd <- c(4,5)
  # DMD <- c(100,250,300,20,120)
  # prop_setupTime <- 0.1
  # prop_demand <- 0.5
  #### END Input for Testing ####

  #### 1. Calculate Process costs for each process ####
  mean_process_duration<-apply(PrD,2,function(x) mean(x[x!=0]))
  setupCosts <- PrC * runif(length(PrC),min = prop_setupChange[1], max = prop_setupChange[2]) * prop_setupTime * mean_process_duration

  #### 2. Calculate process variants for each process ####
  if(is.null(includedProd)) includedProd<-1:NROW(PrD)
  PrVV<-apply(PrD[includedProd,],2,function(t){
    unique(t[t>0])
  })

  #### 3. Demand for each product variant ####
  PrVD<-lapply(1:length(PrVV),function(x){
    temp<-tibble(PrD_n=PrD[includedProd,x],
           DMD=DMD[includedProd]) %>%
      dplyr::group_by(PrD_n) %>%
      dplyr::summarise(DMD=sum(DMD)) %>%
      dplyr::ungroup() %>%
      dplyr::filter(PrD_n!=0)
    return(temp$DMD[match(PrVV[[x]],temp$PrD_n)])
  })

  #### 4. Costs per Setup changes ####

  PrV_setupCosts<-lapply(1:length(PrVD),function(x){
    setupCosts[x] / PrVD[[x]] * clc_meanSetups(lotSize = ceiling(mean(PrVD[[x]])*prop_demand),
                                               DMD = PrVD[[x]],
                                               runs = 100)
  })

  #### 5. Costs per Changes and Product ####
  PC_setup<-sapply(1:length(PrVD),function(x){
    PC_setups<-PrV_setupCosts[[x]][match(PrD[,x],PrVV[[x]])]
    return(ifelse(is.na(PC_setups),0,PC_setups))
  })
  PC_setup<-rowSums(PC_setup)
  PC_setup[-includedProd] <- 0
  return(PC_setup)
}

clc_processCostRate<-function(EAD,RCU){
  ## calculate costs for each process
  RD<-diag(1,NCOL(EAD[[1]]$P$PrD)) %*% ((EAD[[1]]$DMM$PrD_RD %*% EAD[[1]]$DSM$RD) + EAD[[1]]$DMM$PrD_RD)
  return(RD %*% RCU)
}


clc_meanSetups<-function(lotSize,DMD,runs=50){
  #### INput for Testing ####
  # lotSize<-79
  # DMD<-c(400,250)
  #### END Input for Testing ####

  DMD_restore<-DMD
  probs <- DMD/sum(DMD)
  setups<-t(sapply(1:runs,function(i){
    DMD<-DMD_restore
    setups<-rep(0,length(DMD))
    while (any(DMD>0)) {
      DMD_choose <- which(DMD>0)
      if(length(DMD_choose)>1){
        idx<-sample(DMD_choose,1,prob=probs[DMD_choose])
        setups[idx] <- setups[idx]+1
        DMD[idx]<- DMD[idx] - lotSize
      }else{
        DMD[DMD_choose]<-0
        setups[DMD_choose] <- 1
      }
    }
    return(setups)
  }))
  return(ceiling(colMeans(setups)))
}
