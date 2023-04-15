#' @title Substitutes two components by an overdesigned one
#' @description This function randomly chooses two component and replace those by one overdesigned one.
#' The first component is chosen randomly. The second component is selected by finding the nearest neighbor of component one.
#' Since components are substituted this function changes the following matrices \eqn{DMM_{FD,PD}}, \eqn{DSM_{PD}} and \eqn{DMM_{PD,PrD}}.
#' @param EAD An EAD object created via \link[EAD]{crt_EAD_MC}.
#' @return Returns the manipulated EAD.
#' @details For further details see the corresponding vignette by running \code{utils::vignette(package = ‘EAD’)} and \insertCite{Meerschmidt.2024;textual}{EAD}
#' @references
#' \insertAllCited{}
#' @examples
#'
#' EAD <- smallEAD
#' overdesign_EAD(EAD)
overdesign_EAD<-function(EAD){
  require(EAD)
  #### Testing ####
  # EAD <- smallEAD
  # EAD$measures$SYSTEM$SDC_n$FD_PD
  # EAD$measures$SYSTEM$SDC$FD_PD
  # el_substitute<-4
  # el_replaced_by<-2
  #### End Testing ####

    available_components <- which(colSums(EAD$DMM$FD_PD)>0)
    if(length(available_components)<2) stop("Not enough elements for overdesign available.")
    el_substitute <- sample(available_components,1)
    replacement_candidates <- setdiff(available_components,el_substitute)
    neighbors <- as.matrix(dist(t(EAD$DMM$FD_PD)))
    neighbors <- neighbors[el_substitute,replacement_candidates,drop=F]
    el_replaced_by <- replacement_candidates[which.min(neighbors)]


    #### Manipulate DMM ####
    EAD$DMM$FD_PD[,el_replaced_by] <- pmax(EAD$DMM$FD_PD[,el_replaced_by],EAD$DMM$FD_PD[,el_substitute])
    EAD$DMM$FD_PD[,el_substitute] <- 0

    #### Adjust DSM ####
      ## since el_substitute is replaced the dependencies are moved to the el_replaced_by element
      EAD$DSM$PD[,el_replaced_by] <- pmax(EAD$DSM$PD[,el_substitute], EAD$DSM$PD[,el_replaced_by])
      EAD$DSM$PD[,el_substitute] <- 0

      ## since the el_substitute can may require additional elements in DSM_PD we have to remap their dependencies as well
      EAD$DSM$PD[el_replaced_by,] <- pmax(EAD$DSM$PD[el_substitute,],EAD$DSM$PD[el_replaced_by,])
      EAD$DSM$PD[el_substitute,] <- 0

      ## an exception is made when el_substitute is required by el_replaced_by
      ## in such cases the the resulting DSM_PD would contain diagonal entries which is not allowed
      diag(EAD$DSM$PD)<-0

      if(EAD$DSM$PD[el_replaced_by,el_substitute]>0){
        EAD$DSM$PD[el_replaced_by,el_substitute]<-0
      }else{
        EAD$DSM$PD[el_replaced_by,] <- pmax(EAD$DSM$PD[el_substitute,],EAD$DSM$PD[el_replaced_by,])
        EAD$DSM$PD[el_substitute,] <- 0
      }

    #### Adjust Process domain ####
    EAD$DMM$PD_PrD[el_replaced_by,] <- EAD$DMM$PD_PrD[el_replaced_by,] + EAD$DMM$PD_PrD[el_substitute,]
    EAD$DMM$PD_PrD[el_substitute,] <- 0

  EAD <- suppressWarnings(update_EAD(EAD))
  return(EAD)
}

#' @title Calculates Components' Development and Part Administration Costs
#' @description This function calculates components' development and part administration costs.
#' @param P_PD The product mix in the physical domain containing products in rows and components in columns.
#' @param DEMAND A numeric demand vector. If products are excluded from the product mix, the entries are zero.
#' @param PDC_fix A numeric vector holding the fixed costs for a component. These costs are calculated by tracing resource costs back into the physical domain.
#' @param R_dvl Proportion of development costs on \code{PDC_fix}.
#' @param R_PA Proportion of part administration costs on \code{PDC_fix}.
#' @return Returns the summed development costs and part administration costs
#' @details For further details see \insertCite{Meerschmidt.2024;textual}{EAD}
#' @references
#' \insertAllCited{}
#' @examples
#'
#' require(EAD)
#' EAD <- smallEAD
#' P_PD <- EAD$P$PD
#' DEMAND<-c(10,50,100,30)
#' PDC_fix<-c(700,1819,4589,1367)
#' R_dvl<-0.3
#' R_PA<-0.2
#'
#' development_costs(P_PD,
#'                   DEMAND,
#'                   PDC_fix,
#'                   R_dvl,
#'                   R_PA)
development_costs<-function(P_PD,DEMAND,PDC_fix,R_dvl,R_PA){
  require(EAD)
  #### Testing ####
  # EAD <- smallEAD
  # P_PD <- EAD$P$PD
  # DEMAND<-c(0,50,100,30)
  # PDC_fix<-c(700,1819,4589,1367)
  # R_dvl<-0.3
  # R_PA<-0.2
  #### End Testing ####

  total_consumption <- as.numeric(DEMAND %*% P_PD)

  non_zero <- as.numeric(total_consumption > 0)
  TC_dvl <- PDC_fix * R_dvl * non_zero
  TC_PA <- PDC_fix * R_PA * non_zero

  return(list(TC_dvl = sum(TC_dvl),
              TC_PA = sum(TC_PA)))
}



clc_domainCosts<-function(EAD){
  #### Testing ####
  # EAD <- smallEAD
  #### End Testing ####

  TRC <- colSums(EAD$P$RD * EAD$DEMAND)
  RCU_fix <- EAD$RC$fix / TRC
  RCU_var <- (EAD$RC$direct + EAD$RC$var) / TRC

  #### Trace costs ####
    ### in PrD ###
      RCU_PrD_fix <- as.numeric(EAD$DMM$PrD_RD %*% RCU_fix + (EAD$DMM$PrD_RD %*% EAD$DSM$RD) %*% RCU_fix)
      PrDC_fix <- RCU_PrD_fix * colSums(EAD$P$PrD * EAD$DEMAND)
      RCU_PrD_var <- as.numeric(EAD$DMM$PrD_RD %*% RCU_var + (EAD$DMM$PrD_RD %*% EAD$DSM$RD) %*% RCU_var)

    ### in PD ###
      RCU_PD_fix <- as.numeric(EAD$DMM$PD_PrD %*% RCU_PrD_fix + (EAD$DMM$PD_PrD %*% EAD$DSM$PrD) %*% RCU_PrD_fix)
      PDC_fix <- RCU_PD_fix * colSums(EAD$P$PD * EAD$DEMAND)
      RCU_PD_var <- as.numeric(EAD$DMM$PD_PrD %*% RCU_PrD_var + (EAD$DMM$PD_PrD %*% EAD$DSM$PrD) %*% RCU_PrD_var)

  return(list(fix = list(FD = NA,
                               PD = PDC_fix,
                               PrD = PrDC_fix,
                               RD = EAD$RC$fix),
              var = list(FD = NA,
                               PD = RCU_PD_var,
                               PrD = RCU_PrD_var,
                               RD = RCU_var)))
}


#' @title Calculates the setup costs
#' @description This function calculates the setup costs for a specific EAD design. The exact procedure is described in \insertCite{Meerschmidt.2024;textual}{EAD}.
#' See details for more information.
#' @param EAD An EAD object created by \link[EAD]{crt_EAD}.
#' @param R_hold the ratio of component's holding costs on total component costs.
#' @param R_order the ratio of component's order costs on total costs.
#' @param R_setup the proportion of setup costs on process manufacturing costs for a given lot.
#' @param PrCU_var The variable process unit costs.
#' @return Returns the summed setup costs as well as the lot sizes for the individual components
#' @details For further details see \insertCite{Meerschmidt.2024;textual}{EAD} as well as the sub models of \insertCite{Thonemann.2000;textual}{EAD} and \insertCite{Zhang.2020;textual}{EAD}.
#' A vignette is available by running \code{utils::vignette('setupCosts',package='EAD')}.
#' @references
#' \insertAllCited{}
#' @examples
#'
#' require(EAD)
#' EAD <- smallEAD
#' EAD$DEMAND <- c(10,50,100,30)
#' R_hold <- rep(0.05,4)
#' R_order <- rep(0.15,4)
#' R_setup <- rep(0.2,4)
#' PrCU_var <- c(100,26.9,2.19,35.7)
#'
#' clc_setupCosts(EAD,
#'                R_hold,
#'                R_order,
#'                R_setup,
#'                PrCU_var)
setupCosts<-function(EAD,C_hold,C_order,R_setup,PrCU_var){
  require(EAD)
  #### Input for testing ####
  # EAD <- smallEAD
  # EAD$DEMAND <- c(10,50,100,30)
  # C_hold <- rep(5,4)
  # C_order <- rep(15,4)
  # C_setup <- rep(0.2,4)
  # PrCU_var <- c(100,26.9,2.19,35.7)
  #### END Input for Testing ####

  #### 1. Calculate Component Demand ####
  DMD_component <- as.numeric(EAD$DEMAND %*% EAD$P$PD)

  #### 2. Calculate the consumption of processes in total ####
  DMD_PD <- EAD$DMM$PD_PrD + EAD$DMM$PD_PrD %*% EAD$DSM$PrD
  DMD_PD[DMD_PD>0] <- 1
  DMD_PD <- DMD_PD * DMD_component

  #### 3. Calculate Lot Sizes ####
  C_order <- sapply(1:NROW(DMD_PD),function(x) C_order)
  C_hold <- sapply(1:NROW(DMD_PD),function(x) C_hold)
  LZM <- sqrt(2*DMD_PD * C_order / C_hold)

  #### 4. Calculate the Tasks ####
  TM <- ceiling(DMD_PD / LZM)
  TM[is.nan(TM)] <- 0

  #### 5. Sample execution order ####
  execution <- apply(TM,2,function(x){
    idx <- which(x!=0)
    execution <- unlist(sapply(idx,function(y) rep(y,x[y])))
    execution <- execution[sample(1:length(execution))]
    return(execution)
  })

  ## plus one is added since there is a initial setup
  n_setups <-sapply(execution,function(x) sum(diff(x) != 0)+1)


  #### 6. Calculate Lot Process Consumption ####
  LPC <- LZM * (EAD$DMM$PD_PrD + EAD$DMM$PD_PrD %*% EAD$DSM$PrD)

    ## Average Process Consumption per Lot (APCONS)
        APCONS <- apply(LPC,2,function(x) mean(x[x!=0]))

    ## Average lot process costs (ALPC)
        ALPC <- APCONS * PrCU_var

    return(list(TC_setup = sum(ALPC * R_setup * n_setups),
                lotSize = ceiling(apply(LZM, 1, max)),
                DMD_component = DMD_component))

}


orderCosts<-function(DMD_PD,
                     N_lot,
                     C_order){
  TC_order <- ceiling(DMD_PD/N_lot) * C_order
  return(sum(TC_order))
}

supplyCosts <- function(P_PD,
                          C_supply){
  return(sum(apply(P_PD,2,function(x) sum(x)>0)) * C_supply)
}








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
  # includedProd <- c(1,2,3,4,5)
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
      mean_setup<-NA
    }else{
      lotSize<-ceiling(sqrt(2*PrVD[[x]]*r_Order_to_Hold))
      mean_setup<-clc_meanSetups(lotSize = lotSize,DMD = PrVD[[x]],runs = 50)
      costs<-setupCosts[x] / PrVD[[x]] * mean_setup
    }
    temp<-list(costs = costs,
               mean_LotSize = mean(lotSize),
               mean_setup = mean_setup)
    return(temp)
  })

  #### 5. Costs per Changes and Product ####
  PC_setup<-sapply(1:length(PrVD),function(x){
    PC_setups<-PrV_setupCosts[[x]]$costs[match(PrD[,x],PrVV[[x]])]
    return(ifelse(is.na(PC_setups),0,PC_setups))
  })
  mean_LotSize<-sapply(PrV_setupCosts,function(x) x$mean_LotSize)
  PC_setup<-rowSums(PC_setup)
  PC_setup[-includedProd] <- 0

  out<-list(PC_setup=PC_setup,
            lotsize=mean_LotSize,
            setups=sapply(PrV_setupCosts,function(x) x$mean_setup))
  return(out)
}

clc_setupCosts_old<-function(EAD,
                         RCU_var,
                         includedProd = NULL,
                         prop_setupChange=c(0.5,0.5),
                         prop_setupTime=1,
                         r_Order_to_Hold=2){

  DMM_PD_PrD <- EAD[[1]]$DMM$PD_PrD
  DMM_PrD_RD <- EAD[[1]]$DMM$PrD_RD
  DSM_PrD <- EAD[[1]]$DSM$PrD
  DSM_RD <- EAD[[1]]$DSM$RD
  DMD <- EAD[[1]]$DEMAND
  P_PD <- EAD[[1]]$P$PD
  if(is.null(includedProd)) includedProd<-1:NROW(PrD)


  #### Calculate Resource Usage for producing each component ####
  P_PD_PD<-diag(1,nrow=ncol(P_PD))
  P_PrD<-P_PD_PD %*% DMM_PD_PrD + (P_PD_PD %*% DMM_PD_PrD) %*% DSM_PrD
  P_RD<-P_PrD %*% DMM_PrD_RD + (P_PrD %*% DMM_PrD_RD) %*% DSM_RD

  #### calculate lot sizes ####
  ### economic order quantity model
  DMD_PD<-colSums(P_PD[includedProd,]*DMD[includedProd])
  lotSize<-ceiling(sqrt(2*DMD_PD*r_Order_to_Hold))

  #### Calculate Setup Changes ####
  mean_SetupMatrix<-apply(P_RD,2,function(x){
    DMD_temp <- DMD_PD
    DMD_temp[x==0]<-0
    clc_meanSetups(lotSize = lotSize,DMD =DMD_temp,runs = 20)
  })

  c_Setup <- RCU_var * runif(length(RCU_var),min=prop_setupChange[1],max = prop_setupChange[2]) * apply(P_RD,2,function(x) mean(x[x>0]))

  TC_setup <- mean_SetupMatrix %*% c_Setup

  PC_setup <- P_PD %*% ifelse(is.nan(TC_setup / DMD_PD),0,TC_setup / DMD_PD)
  PC_setup[-includedProd] <- 0

  out<-list(PC_setup=PC_setup,
            lotsize=mean(lotSize[lotSize>0]),
            setups=apply(mean_SetupMatrix,1,function(x) mean(x[x>0])))
  return(out)
}

clc_meanSetups<-function(lotSize,DMD,runs=50){
  #### INput for Testing ####
  # lotSize<-c(40,1)
  # DMD<-c(100,50)
  # runs=50
  #### END Input for Testing ####

  if(length(lotSize)==1) lotSize<-rep(lotSize,length(DMD))
  DMD_restore<-DMD
  # probs <- DMD/sum(DMD)
  probs_vec<-ifelse(is.nan( (DMD /lotSize)),0, (DMD /lotSize))
  probs <- probs_vec / sum(probs_vec)
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




clc_processCostRate<-function(EAD,RCU){
  ## calculate costs for each process
  RD<-diag(1,NCOL(EAD[[1]]$P$PrD)) %*% ((EAD[[1]]$DMM$PrD_RD %*% EAD[[1]]$DSM$RD) + EAD[[1]]$DMM$PrD_RD)
  return(RD %*% RCU)
}


