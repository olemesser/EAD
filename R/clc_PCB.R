#' @title Calculate (true) Product Costs
#' @description This function calculates the benchmark (true) product costs (PC_B).
#' There are two options to use this function:
#' \describe{
#'  \item{\code{!is.null(RC_var) & is.null(RCU_indirect)}}{This returns the product costs for a given EAD. \code{sum(PC_B*DMD)==TC}
#'  See the the first example for details.}
#'  \item{\code{is.null(RC_var) & !is.null(RCU_indirect)}}{This returns the product costs for an EAD with an reduced product mix. \code{sum(PC_B*DMD)<TC}
#'  See the second example for details.}
#' }
#' @param RES_CONS_PAT The resource consumption pattern matrix.
#' @param DMD A demand vector. If no demand vector is supplied, each product is assumed to have a demand of one.
#' @param RC_var_d A resource cost vector holding the direct variable costs.
#' @param RC_var_i A resource cost vector holding the indirect variable costs.
#' @param RC_fix_d A resource cost vector holding the direct fixed costs.
#' @param RC_fix_i A resource cost vector holding the indirect fixed costs.
#' @param RCU_indirect A vector holding the resource costs per unit for the variable indirect costs.
#' @param RCU_direct A vector holding the resource costs per unit for the variable direct costs.
#' @return Calculates the true product costs (PC_B) and returns the resource costs per unit (RCU_indirect)
#' @examples
#' ### Calculate product costs for the designed EAD ###
#' data("smallEAD")
#' TC =  sum(smallEAD[[1]]$RC$var_d +
#'           smallEAD[[1]]$RC$var_i +
#'           smallEAD[[1]]$RC$fix_d +
#'           smallEAD[[1]]$RC$fix_i)
#'
#' costs<- clc_PCB(RES_CONS_PAT = smallEAD[[1]]$P$RD,
#'                 DMD = smallEAD[[1]]$DEMAND,
#'                 RC_var_d = smallEAD[[1]]$RC$var_d,
#'                 RC_var_i = smallEAD[[1]]$RC$var_i ,
#'                 RC_fix_i = smallEAD[[1]]$RC$fix_i,
#'                 RC_fix_d = smallEAD[[1]]$RC$fix_d)
#'
#' costs$PC_B
#'
#' ### Calculate product costs for an EAD with an reduced product mix
#' ## lets assume we drop the first ten products
#'
#' DMD<-smallEAD[[1]]$DEMAND
#' DMD[1:10]<-0
#' costs<-clc_PCB(RES_CONS_PAT = smallEAD[[1]]$P$RD,
#'                DMD = DMD,
#'                RC_fix_d = smallEAD[[1]]$RC$fix_d
#'                RC_fix_i = smallEAD[[1]]$RC$fix,
#'                RCU_indirect = costs$RCU_indirect,
#'                RCU_direct = costs$RC_direct)
clc_PCB<-function(RES_CONS_PAT,
                  DMD = NULL,
                  RC_var_d = NULL,
                  RC_var_i = NULL,
                  RC_fix_i = NULL,
                  RC_fix_d = NULL,
                  RCU_indirect = NULL,
                  RCU_direct = NULL){

  #### Input for Testing ####
  # RES_CONS_PAT = EAD[[1]]$P$RD
  # DMD = EAD[[1]]$DEMAND
  # RC_var_d = EAD[[1]]$RC$var_d
  # RC_fix_i =  EAD[[1]]$RC$fix_i
  # RC_var_i = EAD[[1]]$RC$var_i
  # RC_fix_d =  EAD[[1]]$RC$fix_d
  # RCU_indirect = NULL
  # RCU_direct = NULL
  #### End input for Testing ####

  #### 0. Initial Checks ####
  if(is.null(DMD)) DMD<-rep(1,NROW(RES_CONS_PAT))
  if(is.null(RC_fix_i)) RC_fix_i<-rep(0,length(RC_var_d))
  if(is.null(RC_fix_d)) RC_fix_d<-rep(0,length(RC_var_d))
  if(!is.null(RC_var_i) & !is.null(RCU_indirect)) stop("RC_var and RCU_indirect supplied. \nPlease insert either RC_var or RCU_indirect.")
  if(!is.null(RC_var_d) & !is.null(RCU_direct)) stop("RC_var and RCU_indirect supplied. \nPlease insert either RC_var or RCU_indirect.")


  #### 1. Check if Demand is zero for some products ####
  RES_CONS_PAT[DMD==0,]<-0

  #### 2. Calculate TRC and RCU for variable costs ####
  TRC<-colSums(RES_CONS_PAT * DMD)
  if(is.null(RCU_indirect) & !is.null(RC_var_i)) RCU_indirect <- RC_var_i/TRC
  RCU_indirect<-ifelse(is.nan(RCU_indirect) | is.infinite(RCU_indirect),0,RCU_indirect)

  if(is.null(RC_var_d)) RC_var_d <- rep(0,NCOL(RES_CONS_PAT))
  if(is.null(RCU_direct)) RCU_direct <- (RC_var_d/TRC)
  RCU_direct<-ifelse(is.nan(RCU_direct) | is.infinite(RCU_direct),0,RCU_direct)
  PC_B_var_d <- as.numeric(RES_CONS_PAT %*% RCU_direct)
  PC_B_var_i <- as.numeric(RES_CONS_PAT %*% RCU_indirect)

  #### 3. Calculate Fixed Costs ####
  PC_B_fix_i <- RES_CONS_PAT %*% ifelse(is.infinite(RC_fix_i/TRC),0,(RC_fix_i/TRC))
  PC_B_fix_d <- RES_CONS_PAT %*% ifelse(is.infinite(RC_fix_d/TRC),0,(RC_fix_d/TRC))

  #### 4. Calculate Total Costs ####
  TC_var <- sum((PC_B_var_d + PC_B_var_i) * DMD)
  TC_fix <- sum((PC_B_fix_i + PC_B_fix_d)  * DMD)


  return(list(PC_B=as.numeric(PC_B_var_i + PC_B_var_d + PC_B_fix_i + PC_B_fix_d),
              RCU_indirect = RCU_indirect,
              RCU_direct = RCU_direct,
              TRC = TRC,
              PC_direct = as.numeric(PC_B_var_d + PC_B_fix_d),
              PC_B_indirect = as.numeric(PC_B_var_i + PC_B_fix_i),
              PC_B_var_i = as.numeric(PC_B_var_i),
              PC_B_var_d = as.numeric(PC_B_var_d),
              PC_B_fix_i = as.numeric(PC_B_fix_i),
              PC_B_fix_d = as.numeric(PC_B_fix_d),
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
