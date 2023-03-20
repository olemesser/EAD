#' @name crt_DEMAND
#' @title Creates the demand of products
#' @description Creates a demand vector using the \link[stats]{rlnorm} function.
#' @param n_PROD The number of products.
#' @param TOTAL_DEMAND the total demand the products in the portfolio.
#' @param Q_VAR A numeric vector of length one or two which specifies the logarithmic standard deviation.
#' This values is passed as \code{sdlog} to \link[stats]{rlnorm}.
#' @return Returns the demand vector and statistics
#' @examples
#' # example
#' TOTAL_DEMAND<-1000
#' crt_DEMAND(n_PROD=120,TOTAL_DEMAND,Q_VAR=c(0.1,1))
#' @rdname crt_DEMAND
crt_DEMAND<-function(n_PROD,TOTAL_DEMAND,Q_VAR){
  #### For Testing ####
  # n_PROD<-100
  # TOTAL_DEMAND<-10000
  # Q_VAR<-0
  #### End Input Testing ####

  if(length(Q_VAR)==1) Q_VAR <- rep(Q_VAR,2)
  preDemand <- rlnorm(n_PROD,meanlog = 0, sdlog = runif(1,Q_VAR[1],Q_VAR[2]))
  DEMAND <- ceiling((preDemand/sum(preDemand))*TOTAL_DEMAND)

  CHECK<-list(CV=sd(DEMAND)/mean(DEMAND),
              DMD_T10=measure_TOP10(DEMAND),
              Q_VAR = sd(log(DEMAND)))
  out<-list(DEMAND=DEMAND,CHECK=CHECK)

  return(out)

}
