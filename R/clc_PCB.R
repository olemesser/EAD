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
#' @param RCU A vector holding the resource costs per unit.
#' @return Calculates the true product costs (PC_B) and returns the resource costs per unit (RCU)
#' @examples
#' ### Calculate product costs for the designed EAD ###
#' data("exampleEAD")
#' costs<- clc_PCB(RES_CONS_PAT = exampleEAD[[1]][[1]]$P$RD,
#'                 DMD = exampleEAD[[1]][[1]]$DEMAND,
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
#'                RC_fix = exampleEAD[[1]][[1]]$RC$fix,
#'                RCU=costs$RCU)
clc_PCB<-function(RES_CONS_PAT,DMD=NULL,RC_var=NULL,RC_fix=NULL,RCU=NULL){
  #### 0. Initial Checks ####
  if(is.null(DMD)) DMD<-rep(1,NROW(RES_CONS_PAT))
  if (is.null(RC_fix)) RC_fix<-rep(0,length(RC_var))
  if(!is.null(RC_var) & !is.null(RCU)) stop("RC_var and RCU supplied. \nPlease insert either RC_var or RCU.")
  if(is.null(RC_var) & is.null(RCU)) stop("RC_var and RCU are not supplied. \nPlease insert either RC_var or RCU.")

  #### 1. Check if Demand is zero for some products ####
  RES_CONS_PAT[DMD==0,]<-0


  #### 2. Calculate TRC and RCU ####
  TRC<-colSums(RES_CONS_PAT * DMD)
  if(is.null(RCU) & !is.null(RC_var)) RCU <- RC_var/TRC


  #### 3. Calculate (true) Benchmark Costs ####
  PC_B_var <- RES_CONS_PAT %*% RCU
  PC_B_fixed <-RES_CONS_PAT %*% (RC_fix/TRC)
  PC_B <- PC_B_var + PC_B_fixed

  return(list(PC_B=PC_B,RCU=RCU))
}
