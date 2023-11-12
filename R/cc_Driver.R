#' @title Substitutes two components by an overdesigned one
#' @description This function randomly chooses two component and replace those by one overdesigned one.
#' There are two methods available for manipulating the degree of overdesign (see \code{method}).
#' Since components are substituted this function changes the following matrices \eqn{DMM_{FD,PD}}, \eqn{DSM_{PD}} and \eqn{DMM_{PD,PrD}}.
#' @param EAD An EAD object created via \link[EAD]{crt_EAD_MC}.
#' @param bounds A numerical vector specifying the lower an upper bound for an increased in material costs due to overdesign. Default \code{bounds=c(1,1)}
#' @param method The method which is used to overdesign the \eqn{DMM_{FD,PD}} matrix.
#' \describe{
#'   \item{random}{The first component is chosen randomly.
#'   The second component is selected by finding the nearest neighbor of component one by using the Euclidean distance.}
#'    \item{optimized}{As a first component, the one with the smallest demand rate is chosen.
#'    The second component, is selected based on two criteria.
#'    First, the procedures searches if there is a component available which overfills the requirements.
#'    If this is the case, there are no additional costs since the component already exists.
#'    If there is no component which overfill the requirements of component one, the second component is selected by finding the nearest neighbor of component one by using the Euclidean distance.}
#'    }
#' @return Returns the manipulated EAD.
#' @details For further details see the corresponding vignette by running \code{utils::vignette('cc-experiment', package = ‘EAD’)} and \insertCite{Meerschmidt.2024;textual}{EAD}
#' @references
#' \insertAllCited{}
#' @examples
#'
#' EAD <- smallEAD
#' overdesign_EAD(EAD,bounds=c(1,1),method='random')
overdesign_EAD<-function(EAD,bounds=c(1,1),method=c('random','optimized')){
  require(EAD)
  #### Testing ####
  # EAD <- smallEAD
  # EAD$measures$SYSTEM$SDC_n$FD_PD
  # EAD$measures$SYSTEM$SDC$FD_PD
  # el_substitute<-2
  # el_replaced_by<-3
  # bounds=c(0.5,0.5)
  # bounds=c(1,1)
  # method <- "optimized"
  #### End Testing ####

  if(length(method)>1) method <- method[1]

    costs <- clc_domainCosts(EAD)

    available_components <- which(colSums(EAD$DMM$FD_PD)>0)
    if(length(available_components)<2) stop("Not enough elements for overdesign available.")

    if(method == "random"){
       el_substitute <- sample(available_components,1)
       replacement_candidates <- setdiff(available_components,el_substitute)
       neighbors <- as.matrix(dist(t(EAD$DMM$FD_PD)))
       neighbors <- neighbors[el_substitute,replacement_candidates,drop=F]
       el_replaced_by <- replacement_candidates[which.min(neighbors)]
    }else if(method == "optimized"){
      DMD_component <- as.numeric(EAD$DEMAND %*% EAD$P$PD)
      el_substitute <- which.min(ifelse(DMD_component==0,NA,DMD_component))
      replacement_candidates <- setdiff(available_components,el_substitute)
      ## check if there is a component which overfills the requirements
      overfill <- EAD$DMM$FD_PD[,el_substitute] %*% EAD$DMM$FD_PD[,replacement_candidates,drop=F]
      overfill <- as.numeric(overfill[1,])
      el_replaced_by <- replacement_candidates[which(overfill>=sum(EAD$DMM$FD_PD[,el_substitute]))]
      if(length(el_replaced_by)>0){
        el_replaced_by <-  el_replaced_by[which.min(costs$var$PD[el_replaced_by])]
      }else{
        neighbors <- as.matrix(dist(t(EAD$DMM$FD_PD)))
        neighbors <- neighbors[el_substitute,replacement_candidates,drop=F]
        el_replaced_by <- replacement_candidates[which.min(neighbors)]
      }
    }else{
      stop("You selected an incorrect method for the overdesign! Options are: 'random', 'optmized'. ")
    }



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
      max_costsEL <-  apply(EAD$DMM$PD_PrD[c(el_substitute,el_replaced_by),],2,function(x) x[which.max(x)])
      min_costsEL <-apply(EAD$DMM$PD_PrD[c(el_substitute,el_replaced_by),],2,function(x) x[which.min(x)])

      min_bound <- bounds[1]
      repeat{
        EAD_temp <- EAD
        EAD_temp$DMM$PD_PrD[el_replaced_by,] <-max_costsEL + min_costsEL * min_bound
        EAD_temp$DMM$PD_PrD[el_substitute,] <- 0
        #### Update EAD (Design & Costs) ####
        EAD_temp <- suppressWarnings(update_EAD(EAD_temp))
        costs_new <- clc_domainCosts(EAD_temp)
        EAD_temp$RC$direct <-costs$var$RD * colSums(EAD_temp$P$RD * EAD_temp$DEMAND)

        if(sum(EAD_temp$RC$direct)<sum(EAD$RC$direct)){
          min_bound <- min_bound+0.05
        }else{
          break
        }
      }

    EAD$DMM$PD_PrD[el_replaced_by,] <-max_costsEL + min_costsEL * runif(length(min_costsEL),
                                                                         min = min_bound,
                                                                         max = ifelse(min_bound>bounds[2],min_bound,bounds[2]))
    EAD$DMM$PD_PrD[el_substitute,] <- 0

    #### Update EAD (Design & Costs) ####
    EAD <- suppressWarnings(update_EAD(EAD))
    costs_new <- clc_domainCosts(EAD)
    EAD$RC$direct <-costs$var$RD * colSums(EAD$P$RD * EAD$DEMAND)

    return(list(EAD = EAD,
                overdesign = list(substitute = el_substitute,
                                  replaced_by = el_replaced_by)))
}

materialCosts <- function(TC_var_S0,
                          P_RD,
                          DMD,
                          RCU){
  TC_var_S1 <- sum(as.numeric(P_RD %*% RCU) * DMD)
  return(TC_var_S1 - TC_var_S0)
}

overdesign_costs <- function(PDUC_var,
                             design,
                             bounds = c(0,1)){

  #### Input for Testing ####
  # PDUC_var <- c(1,30,4,5,6)
  # design <- c(2,1)
  # bounds <- c(1,1)
  #### End Input for testing ####

  PDUC_var_new <- PDUC_var
  if(!is.null(design)){
    PDUC_var_new[design[2]] <- max(PDUC_var[design]) + min(PDUC_var[design]) * runif(1,
                                                                                     min = bounds[1],
                                                                                     max = bounds[2])
    PDUC_var_new[design[1]] <- 0
  }

  return(list(PDUC = list(var = PDUC_var_new,
                          fix = NA)))
}


#' @title Calculates Components' Development  Costs
#' @description This function calculates components' development costs as defined by \insertCite{Meerschmidt.2024;textual}{EAD}.
#' @param P_PD The product mix in the physical domain containing products in rows and components in columns.
#' @param DEMAND A numeric demand vector. If products are excluded from the product mix, the entries are zero.
#' @param C_dvl The development costs for each component.
#' @param overdesign_Change A vector of length two, specifying the overdesign. The first entry indicates the index of the component which is substituted by the second indexed component.
#' @param bounds A vector specifying the bounds for a uniform distribution.
#' @return Returns the summed development costs.
#' @details For further details see \insertCite{Meerschmidt.2024;textual}{EAD}.
#' @references
#' \insertAllCited{}
#' @examples
#'
#' require(EAD)
#' EAD <- smallEAD
#' P_PD <- EAD$P$PD
#' DEMAND<-c(10,50,100,30)
#' C_dvl<-c(210,117,303,210)
#'
#' development_costs(P_PD,
#'                   DEMAND,
#'                   PDC_fix,
#'                   R_dvl)
development_costs<-function(P_PD,
                            DEMAND,
                            C_dvl,
                            overdesign_Change=NULL,
                            bounds = c(1,1)){
  require(EAD)
  #### Testing ####
  # EAD <- smallEAD
  # P_PD <- EAD$P$PD
  # DEMAND<-c(0,50,100,30)
  # C_dvl <- c(10,7,12,14)
  # overdesign_Change <- c(3,2)
  # bounds = c(1,1)
  #### End Testing ####

  TCC <- as.numeric(DEMAND %*% P_PD)
  non_zero <- as.numeric(TCC > 0)
  if(!is.null(overdesign_Change)){
    C_dvl[overdesign_Change[2]] <-  max(C_dvl[overdesign_Change]) +  min(C_dvl[overdesign_Change]) * runif(1,
                                                                                                       min = bounds[1],
                                                                                                       max = bounds[2])
    C_dvl[overdesign_Change[1]] <- 0
  }
  TC_dvl <- C_dvl * non_zero

  return(list(TC_dvl = sum(TC_dvl),
              C_dvl = C_dvl))
}

#' @title Calculates Components' Part Administration Costs
#' @description This function calculates components' part administration costs as defined by \insertCite{Meerschmidt.2024;textual}{EAD}.
#' @param P_PD The product mix in the physical domain containing products in rows and components in columns.
#' @param DEMAND A numeric demand vector. If products are excluded from the product mix, the entries are zero.
#' @param PDC_fix A numeric vector holding the fixed costs for a component. These costs are calculated by tracing resource costs back into the physical domain.
#' @param R_PA Proportion of part administration costs on \code{PDC_fix}.
#' @return Returns the summed development costs.
#' @details For further details see \insertCite{Meerschmidt.2024;textual}{EAD}.
#' @references
#' \insertAllCited{}
#' @examples
#'
#' require(EAD)
#' EAD <- smallEAD
#' P_PD <- EAD$P$PD
#' DEMAND<-c(10,50,100,30)
#' C_PA<-c(140,78,202,140)
#'
#' development_costs(P_PD,
#'                   DEMAND,
#'                   PDC_fix,
#'                   C_PA)
partmgmtCosts <- function(P_PD,DEMAND,C_PA){
  require(EAD)
  TCC <- as.numeric(DEMAND %*% P_PD)
  non_zero <- as.numeric(TCC > 0)
  TC_PA <- C_PA * non_zero
  return(sum(TC_PA))
}


clc_domainCosts<-function(EAD){
  #### Testing ####
  # EAD <- smallEAD
  #### End Testing ####

  TRC <- colSums(EAD$P$RD * EAD$DEMAND)
  RCU_fix <- EAD$RC$fix_d + EAD$RC$fix_i / TRC
  RCU_var <- ( EAD$RC$var_d + EAD$RC$var_i) / TRC

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
#' @param C_hold Component's holding costs.
#' @param C_order Component's order costs.
#' @param C_setup The average setup costs.
#' @return Returns the summed setup costs as well as the lot sizes for the individual components
#' @details For further details see \insertCite{Meerschmidt.2024;textual}{EAD} as well as the sub models of \insertCite{Thonemann.2000;textual}{EAD} and \insertCite{Zhang.2020;textual}{EAD}.
#' A vignette is available by running \code{utils::vignette('setupCosts',package='EAD')}.
#' @references
#' \insertAllCited{}
#' @examples
#'
#' require(EAD)
#' set.seed(1234)
#' EAD <- smallEAD
#' EAD$DEMAND <- c(10,50,100,30)
#' C_hold <- rep(0.05,4)
#' C_order <- rep(0.15,4)
#' C_setup <- 444.4
#'
#' clc_setupCosts(EAD,
#'                C_hold,
#'                C_order,
#'                C_setup)
setupCosts<-function(P_PD,
                     DMM_PD_PrD,
                     DSM_PrD,
                     DEMAND,
                     C_hold,
                     C_order,
                     C_setup){
  require(EAD)
  #### Input for testing ####
  # EAD <- smallEAD
  # P_PD <- EAD$P$PD
  # DMM_PD_PrD <- EAD$DMM$PD_PrD
  # DSM_PrD <- EAD$DSM$PrD
  # DEMAND <- c(10,50,100,30)
  # C_hold <- rep(5,4)
  # C_order <- rep(15,4)
  # C_setup <- 444.4
  #### END Input for Testing ####

  #### 1. Calculate Component Demand ####
  DMD_component <- as.numeric(DEMAND %*% P_PD)

  #### 2. Calculate the consumption of processes in total ####
  DMD_PD <- DMM_PD_PrD + DMM_PD_PrD %*% DSM_PrD
  DMD_PD[DMD_PD>0] <- 1
  DMD_PD <- DMD_PD * DMD_component

  #### 3. Calculate Lot Sizes ####
  C_order <- sapply(1:NROW(DMD_PD),function(x) C_order)
  C_hold <- sapply(1:NROW(DMD_PD),function(x) C_hold)
  LZM <- sqrt(2*DMD_PD * C_order / C_hold)

  #### 4. Calculate the Tasks ####
  TM <- ceiling(DMD_PD / LZM)
  TM[is.nan(TM)] <- 0


  #### 5. Sample execution order  and calculate the number of setups ####
  runs<-500
  n_setups <- matrix(NA,nrow = runs,ncol = NCOL(TM))
  for(i in 1:runs){
    for (j in 1:NCOL(TM)) {
      x <- TM[,j]
      execution <- rep(which(x != 0), x[x != 0])
      execution <- sample(execution)
      n_setups[i,j] <- sum(diff(execution) != 0) + 1
    }
  }
  n_setups <- ceiling(colMeans(n_setups))

  TM_bin <- TM
  TM_bin[TM_bin>0] <- 1

     return(list(TC_setup = sum(n_setups * C_setup),
                lotSize = ceiling(apply(LZM, 1, max)),
                DMD_component = DMD_component,
                n_setups = n_setups,
                n_tasks = sum(TM),
                n_taks_distinct = sum(TM_bin),
                DNS_TM = measure_DENS(TM[rowSums(TM)>0,]),
                TM = TM,
                rm_TM = mean(rowMeans(TM_bin)[rowMeans(TM_bin)>0])))

}


clc_numbSetups<-function(TM,rep = 500){
  n_setups <- t(sapply(1:rep,function(y){
    n_setups <- apply(TM,2,function(x){
      execution <- rep(which(x != 0), x[x != 0])
      execution <- execution[sample(1:length(execution))]
      ## plus one is added since there is a initial setup
      return(sum(diff(execution) != 0)+1)
    })
    return(n_setups)
  }))
 return(n_setups)
}


orderCosts<-function(DMD_PD,
                     N_lot,
                     C_order){
  N_order <- ifelse(!is.nan(ceiling(DMD_PD/N_lot)),ceiling(DMD_PD/N_lot),0)
  TC_order <- ifelse(is.nan(N_order),0,N_order) * C_order
  return(list(TC_order = sum(TC_order),
              N_order = N_order))
}

supplyCosts <- function(P_PD,
                        DEMAND,
                        C_supply){
  P_PD[DEMAND==0,] <- 0
  return(sum(apply(P_PD,2,function(x) sum(x)>0)) * C_supply)
}


stockCosts <-function(C_hold,
                      lotSize){
  Units_stock <- sum(ceiling(3/2 *lotSize))
  return( list(TC_stock = C_hold * Units_stock,
               Units_stock = Units_stock))
}


toolingCosts <- function(DMM_PD_PrD,
                         DSM_PrD,
                         DMD_component,
                         C_tooling){
  PVM <- DMM_PD_PrD + DMM_PD_PrD %*% DSM_PrD
  PVM[DMD_component==0,] <- 0
  C_tooling[PVM==0] <- 0
  N_processVariety <- sum(PVM>0)
  return(list(C_tooling = sum(C_tooling),
              N_processVariety = N_processVariety))
}


clc_initialCosts <- function(EAD,
                             R_dvl,
                             R_PA,
                             R_tooling,
                             R_hold,
                             R_order,
                             R_setup,
                             R_supply){
  ## checks ###
  if(any(c(R_dvl,R_PA,R_tooling,R_hold,R_order,R_setup,R_supply)<0)) stop("At least on of the cost ratios is negative. Please pass positive values.")
  if(sum(R_dvl,R_PA,R_tooling,R_supply)>1) stop("R_dvl + R_PA + R_tooling + R_supply > 1. Please reduce at least one of these values.")
  if(sum(R_hold,R_order,R_setup)>1) stop("R_hold + R_order + R_setup > 1. Please reduce at least one of these values.")

  costs <- clc_domainCosts(EAD)
  TC_var <- sum(EAD$RC$direct+EAD$RC$var)
  TC_fix <- sum(EAD$RC$fix)
  out<-list()
  #### Fixed Cost Drivers ####
    ## calculate the initial development costs ##
    out[['C_dvl']] <- costs$fix$PD * R_dvl

    ## calculate the initial part administration costs ##
    out[['C_pa']] <- costs$fix$PD * R_PA

    ## calculate the initial tooling costs ##
    PVM <- EAD$DMM$PD_PrD + EAD$DMM$PD_PrD %*% EAD$DSM$PrD
    N_proVar <- sum(PVM>0)
    C_tooling <- matrix(0,nrow = NROW(PVM),ncol = NCOL(PVM))
    C_tooling[PVM>0] <- runif(N_proVar,
                              min = 0,
                              max = 1)
    C_tooling <- C_tooling /sum(C_tooling)
    out[['C_tooling']] <- C_tooling * (TC_fix * R_tooling)

    ## calculate initial supply costs ##
    N_supplier <- supplyCosts(P_PD = EAD$P$PD,
                              DEMAND = EAD$DEMAND,
                              C_supply = 1)
    out[['C_supply']] <- TC_fix * R_supply / N_supplier


  #### Variable Cost Drivers ###
  ## calculate costs for a single setup ##
  setups <- setupCosts(P_PD = EAD$P$PD,
                           DMM_PD_PrD = EAD$DMM$PD_PrD,
                           DSM_PrD = EAD$DSM$PrD,
                           DEMAND = EAD$DEMAND,
                           C_hold = R_hold,
                           C_order = R_order,
                           C_setup = 0)
  out[['C_setup']] <- TC_var * R_setup / sum(setups$n_setups)

  ## calculate initial costs per order ##
  N_orders <- sum(ceiling(setups$DMD_component/setups$lotSize))
  out[['C_order']] <- TC_var * R_order / N_orders

  ## calculate initial holding costs per component unit ##
  out[['C_hold']] <- TC_var * R_hold / sum(3/2 * setups$lotSize)

  out[['TC_fix_CC']] <- TC_fix * (R_supply + R_tooling + R_PA + R_dvl)
  out[['TC_var_CC']] <- TC_var * (R_hold + R_order +  R_setup)


  return(out)
}


clc_variableComponentCosts <- function(EAD,RCU){
  require(EAD)
  #### Input For Testing ####
  # EAD <- smallEAD
  # costs <- clc_PCB(EAD$P$RD,
  #                  DMD = EAD$DEMAND,
  #                  RC_direct = EAD$RC$direct + EAD$RC$var,
  #                  RC_var = rep(0,length(EAD$RC$var)),
  #                  RC_fix = EAD$RC$fix)
  # RCU <- costs$RCU_direct
  #### End Input for Testing ####

  P_PrD <- diag(1,nrow = NROW(EAD$DMM$PD_PrD)) %*% ((EAD$DMM$PD_PrD %*% EAD$DSM$PrD) + EAD$DMM$PD_PrD)
  P_RD <- P_PrD %*% ((EAD$DMM$PrD_RD %*% EAD$DSM$RD) + EAD$DMM$PrD_RD)
  costs <- clc_PCB(RES_CONS_PAT = P_RD,
                  RCU_direct = RCU,
                  RCU = rep(0,length(RCU)),
                  RC_fix = rep(0,length(RCU)))

  PDUC_var <- costs$PC_direct
  return(list(var = PDUC_var,
              fix = NA))
}
