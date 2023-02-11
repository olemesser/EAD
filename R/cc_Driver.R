processVariety<-function(PrD,
                         PrC,
                         DMD,
                         includedProd=NULL,
                         prop_setupChange=c(0.5,0.5),
                         prop_setupTime=1,
                         r_Order_to_Hold=2){
  require(dplyr)
  #### Input for Testing ####
  # PrD <- matrix(c(10,12,10,0,0,
  #               6,5,2,2,4),ncol = 2)
  # PrC <- c(120,130)
  # prop_setupChange <- c(0.5,0.5)
  # includedProd <- c(4,5)
  # DMD <- c(100,250,300,20,120)
  # prop_setupTime <- 0.1
  # r_Order_to_Hold <- 2
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
    ## economic order quantity ##
    if (length(PrVD[[x]])==0) {
      lotSize<-NA
      costs<-0
    }else{
      lotSize<-ceiling(sqrt(2*PrVD[[x]]*r_Order_to_Hold))
      costs<-setupCosts[x] / PrVD[[x]] * clc_meanSetups(lotSize = lotSize,DMD = PrVD[[x]],runs = 50)
    }
    temp<-list(costs = costs,
               mean_LotSize = mean(lotSize))
    return(temp)
  })

  #### 5. Costs per Changes and Product ####
  PC_setup<-sapply(1:length(PrVD),function(x){
    PC_setups<-PrV_setupCosts[[x]]$costs[match(PrD[,x],PrVV[[x]])]
    return(ifelse(is.na(PC_setups),0,PC_setups))
  })
  mean_LotSize<-mean(na.omit(sapply(PrV_setupCosts,function(x) x$mean_LotSize)))
  PC_setup<-rowSums(PC_setup)
  PC_setup[-includedProd] <- 0

  out<-list(PC_setup=PC_setup,
            lotsize=mean_LotSize)
  return(out)
}

clc_processCostRate<-function(EAD,RCU){
  ## calculate costs for each process
  RD<-diag(1,NCOL(EAD[[1]]$P$PrD)) %*% ((EAD[[1]]$DMM$PrD_RD %*% EAD[[1]]$DSM$RD) + EAD[[1]]$DMM$PrD_RD)
  return(RD %*% RCU)
}


clc_meanSetups<-function(lotSize,DMD,runs=50){
  #### INput for Testing ####
  # lotSize<-c(40,1)
  # DMD<-c(100,50)
  # runs=50
  #### END Input for Testing ####

  if(length(lotSize)==1) lotSize<-rep(lotSize,length(DMD))
  DMD_restore<-DMD
  probs <- DMD/sum(DMD)

  setups<-t(sapply(1:runs,function(i){
    DMD<-DMD_restore
    setups<-rep(0,length(DMD))
    idx_old<-0
    while (any(DMD>0)) {
      DMD_choose <- which(DMD>0)
        if(length(DMD_choose)==1){
          idx<-DMD_choose
        }else{
          idx<-sample(DMD_choose,1,prob = probs[DMD_choose])
        }
        if(idx_old!=idx & DMD[idx]>0){
          setups[idx] <- setups[idx]+1
        }
        DMD[idx]<- DMD[idx] - lotSize[idx]
        idx_old<-idx
    }
    return(setups)
  }))
  return(ceiling(colMeans(setups)))
}
