#' @title Create a resource cost object
#' @description This function creates three resource cost vectors according to traditional costing systems.
#' In a first step, a direct cost vector is created with the total costs of  \code{TC*(1-r_in)} and a coefficient of variation of \code{cv}.
#' The remaining \code{TC*r_in} costs are indirect costs.
#' Both, the variable and fixed indirect cost vector are created via sampling of values. The final vectors have a correlation of \code{cor_fix} and \code{cor_var} regarding the direct cost vector.
#' The total fixed an indirect costs are defined as \code{(TC*r_in*r_fix)}.
#' @param N_RD Number of resources.
#' @param TC Total costs.
#' @param r_in Proportion of indirect costs on total costs.
#' @param r_fix Proportion of fixed costs on total indirect costs.
#' @param sdlog The log-normal standard deviation for the total resource cost vector.
#' @return A list containing the three cost vectors. For each vector the coefficient of variation as well as the top 10% largest resource costs for each vector are calculated.
#' Since the procedure is not able to match the exact input values, the final correlation value are measured. \code{cor_fix} measures the correlation between indirect fixed costs and direct cost vector.
#' \code{cor_var} between indirect variable costs and direct cost vector and \code{cor_fix_var} between both indirect cost vectors.
#' @examples
#' set.seed(1234)
#'
#' RC<-crt_RC(N_RD = 50,
#'            TC = 10^6,
#'            r_in = 0.5,
#'            r_fix = 0.4,
#'            sdlog = 0.2)
#'
crt_RC<-function(N_RD,
                 TC=10^6,
                 r_in,
                 r_fix,
                 sdlog){

  #### Input for Testing ####
  # N_RD = 50
  # TC = 10^6
  # r_in = 0.5
  # r_fix = 0.5
  # sdlog = 1
  #### End Input for Testing ####

  #### 0. Initial Checks ####
  if(r_in>0.9) stop("You selected more than 90% of the costs as indirect. The maximum is proportion is 90%. Please reduce r_in!")
  if(r_in<0) stop("The proportion of indirect costs cannot be negative. Please set r_in>=0!")
  if(r_fix<0 | r_fix>1) stop("r_fix must be within the bound of 0<=r_fix<=1")
  if(sdlog<=0) stop("The log-normal standard deviation cannot be zero or negative. Please increase the input value.")

  #### 1. Split Costs according the input ratios ####
  TC_direct <- TC * (1-r_in)
  TC_indirect <- TC * (r_in)
  TC_fix <- TC * r_fix
  TC_var <- TC - TC_fix

  #### 2. Create Total variable and fixed cost vector ####
  ### 2.1 fixed costs ###
  RC_fix <- crt_RCid(l = N_RD,
                  TC = TC_fix,
                  sdlog = sdlog)

  ### 2.2 variable costs ###
  RC_var <- crt_RCid(l = N_RD,
                  TC = TC_var,
                  sdlog = sdlog)

  #### 3. Split RC into direct and indirect costs ####
  ## calculate possible ratios for indirect fixed and variable costs
  ## the total sum of indirect costs must be TC_indirect
  r_in_fix_min <- ifelse(sum(RC_var$RC) < TC_indirect,
                         abs(TC_indirect - sum(RC_var$RC)) / sum(RC_fix$RC),
                         0)
  r_in_max <- ifelse(sum(RC_fix$RC) > TC_indirect,
                     TC_indirect/sum(RC_fix$RC) * 0.9,
                     0.9)
  ## 3.1 Split fixed cost vector into direct and indirect costs ##
  RC_fix <- crt_RC_indirect(RC = RC_fix$RC, r_in = runif(1,min = r_in_fix_min,max = r_in_max))

  ## 3.2 Split variable cost vector into direct and indirect costs ##
  r_in_var <- (TC_indirect - sum(RC_fix$RC_i$RC)) / sum(RC_var$RC)
  RC_var <- crt_RC_indirect(RC = RC_var$RC, r_in = r_in_var)


  #### 4. Create Output Object ####
  output<-list(RC_var = RC_var,
               RC_fix = RC_fix,
               r_fix =   sum(RC_fix$RC_i$RC +   RC_fix$RC_d$RC) / TC,
               r_in = sum(RC_fix$RC_i$RC +   RC_var$RC_i$RC) / TC,
               TC = sum(RC_fix$RC_i$RC + RC_fix$RC_d$RC + RC_var$RC_i$RC +   RC_var$RC_d$RC))

  return(output)

}


crt_RCid<-function(l,TC,sdlog){
  RC <- rlnorm(l,meanlog = 0, sdlog = sdlog)
  RC <- (RC/sum(RC)*TC)
  cv_is<-sd(RC)/mean(RC)
  return(list(RC=RC,
              cv=cv_is))
}

crt_RC_indirect <- function(RC,r_in,tol=0.005){
  N_RD <- length(RC)
  repeat{
    prop_in <- runif(N_RD,min = r_in -  min(c(r_in, 1 - r_in)),
                     max = r_in +  min(c(r_in, 1 - r_in)))
    RC_i <- RC * prop_in
    RC_d <- RC - RC_i
    r_in_is <- sum(RC_i)/sum(c(RC_d,RC_i))
    if(abs(r_in - r_in_is) < tol){
      break
    }else if(r_in_is > r_in){
      idx <- sample(1:N_RD,1)
      prop_in[idx] <- prop_in[idx] * runif(1,min=0.9,max=1)
    }else if(r_in_is < r_in){
      idx <- sample(1:N_RD,1)
      prop_in[idx] <- prop_in[idx] * runif(1,min=1,max=1.1)
    }
  }
  if(any(RC_i < 0) |  any(RC_d < 0)) stop("Invalid cost vectors created")

  out <- list(RC_i = list(RC = RC_i,
                          top10 = measure_TOP10(RC_i)),
              RC_d = list(RC = RC_d,
                          top10 = measure_TOP10(RC_d)))

  return(out)
}
