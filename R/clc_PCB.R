#' @title Calculate (true) Product Costs
#' @description This function calculates the benchmark (true) product costs (PC_B).
#' There are two options to use this function:
#' \describe{
#'  \item{\code{!is.null(RC_var) & is.null(RCU)}}{This returns the product costs for a given EAD. \code{sum(PC_B*DMD)==TC}
#'  See the the first example for details.}
#'  \item{\code{is.null(RC_var) & !is.null(RCU)}}{This returns the product costs for an EAD with an reduced product mix. \code{sum(PC_B*DMD)<TC}
#'  See the second example for details.}
#' }
#' @param RES_CONS_PAT The resource consumption pattern matrix.
#' @param DMD A demand vector. If no demand vector is supplied, each product is assumed to have a demand of one.
#' @param RC_var A resource cost vector holding the variable costs.
#' @param RC_fix A resource cost vector holding the fixed costs.
#' @param RCU A vector holding the resource costs per unit for the variable indirect costs.
#' @param RCU_direct A vector holding the resource costs per unit for the variable direct costs.
#' @return Calculates the true product costs (PC_B) and returns the resource costs per unit (RCU)
#' @examples
#' ### Calculate product costs for the designed EAD ###
#' data("exampleEAD")
#' TC =  sum(exampleEAD[[1]][[1]]$RC$fix +
#'           exampleEAD[[1]][[1]]$RC$var +
#'           exampleEAD[[1]][[1]]$RC$direct)
#'
#' costs<- clc_PCB(RES_CONS_PAT = exampleEAD[[1]][[1]]$P$RD,
#'                 DMD = exampleEAD[[1]][[1]]$DEMAND,
#'                 RC_direct = exampleEAD[[1]][[1]]$RC$direct,
#'                 RC_var = exampleEAD[[1]][[1]]$RC$var,
#'                 RC_fix = exampleEAD[[1]][[1]]$RC$fix)
#'
#' costs$PC_B
#'
#' ### Calculate product costs for an EAD with an reduced product mix
#' ## lets assume we drop the first ten products
#'
#' DMD<-exampleEAD[[1]][[1]]$DEMAND
#' DMD[1:10]<-0
#' costs<-clc_PCB(RES_CONS_PAT = exampleEAD[[1]][[1]]$P$RD,
#'                DMD = DMD,
#'                RC_direct = exampleEAD[[1]][[1]]$RC$direct,
#'                RC_fix = exampleEAD[[1]][[1]]$RC$fix,
#'                RCU=costs$RCU)
clc_PCB<-function(RES_CONS_PAT,
                  DMD = NULL,
                  RC_direct = NULL,
                  RC_var = NULL,
                  RC_fix = NULL,
                  RCU = NULL,
                  RCU_direct = NULL){
  #### 0. Initial Checks ####
  if(is.null(DMD)) DMD<-rep(1,NROW(RES_CONS_PAT))
  if(is.null(RC_fix)) RC_fix<-rep(0,length(RC_var))
  if(!is.null(RC_var) & !is.null(RCU)) stop("RC_var and RCU supplied. \nPlease insert either RC_var or RCU.")
  if(is.null(RC_var) & is.null(RCU)) stop("RC_var and RCU are not supplied. \nPlease insert either RC_var or RCU.")

  #### 1. Check if Demand is zero for some products ####
  RES_CONS_PAT[DMD==0,]<-0

  #### 2. Calculate TRC and RCU ####
  TRC<-colSums(RES_CONS_PAT * DMD)
  if(is.null(RCU) & !is.null(RC_var)) RCU <- RC_var/TRC
  RCU<-ifelse(is.nan(RCU) | is.infinite(RCU),0,RCU)

  if(is.null(RC_direct)) RC_direct <- rep(0,NCOL(RES_CONS_PAT))
  if(is.null(RCU_direct)) RCU_direct <- (RC_direct/TRC)
  RCU_direct<-ifelse(is.nan(RCU_direct) | is.infinite(RCU),0,RCU_direct)
  PC_direct <- as.numeric(RES_CONS_PAT %*% RCU_direct)

  #### 3. Calculate (true) Benchmark Costs ####
  PC_B_var <- as.numeric(RES_CONS_PAT %*% RCU)
  PC_B_fixed <- RES_CONS_PAT %*% ifelse(is.infinite(RC_fix/TRC),0,(RC_fix/TRC))

  #### 4. Calculate Total Costs ####
  TC_var <- sum((PC_B_var + PC_direct) * DMD)
  TC_fix <- sum(PC_B_fixed * DMD)


  return(list(PC_B=as.numeric(PC_B_var + PC_B_fixed + PC_direct),
              RCU=RCU,
              RCU_direct=RCU_direct,
              TRC = TRC,
              PC_direct = as.numeric(PC_direct),
              PC_B_indirect = as.numeric(PC_B_var + PC_B_fixed),
              PC_B_fixed=as.numeric(PC_B_fixed),
              PC_B_var=as.numeric(PC_B_var),
              TC_var = TC_var,
              TC_fix = TC_fix,
              TC = TC_var + TC_fix))
}



clc_costingERROR<-function(PC_B,PC_H,DMD=NULL){
  if(is.null(DMD)) DMD<-rep(1,length(PC_B))
  PC_B<-rep(PC_B,DMD)
  PC_H<-rep(PC_H,DMD)
  PE <- (PC_B - PC_H)/PC_B
  output<-list(EUCD = sqrt(sum((PC_B - PC_H)^2)),
               MPE = as.numeric(mean(abs(PE))),
               PE = as.numeric(PE),
               APE = as.numeric(abs(PE)))

  return(output)
}
