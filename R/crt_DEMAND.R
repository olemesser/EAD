#' @name crt_DEMAND
#' @title Creates the demand of products
#' @description Creates a demand vector with a given level of heterogeneity which is specified by cv.
#' @param n_PROD The number of products.
#' @param TOTAL_DEMAND the total demand the products in the portfolio.
#' @param cv Desired coefficient of variation as a vector representing the lower or upper bounds or as a double.
#' @return Returns the demand vector and statistics
#' @examples
#' # example
#' TOTAL_DEMAND<-1000
#' crt_DEMAND(n_PROD=120,TOTAL_DEMAND,cv=c(0.1,1))
#' @rdname crt_DEMAND
crt_DEMAND<-function(n_PROD,TOTAL_DEMAND,cv){
  #### For Testing ####
  # n_PROD<-30
  # TOTAL_DEMAND<-10000
  # cv=c(3,3)
  #### End Input Testing ####

  if(length(cv)==1) cv<-rep(cv,2)
  cv<-runif(1,min=cv[1],max=cv[2])
  cv_temp<-cv
  tries<-0
  repeat{
    tries<-tries+1
    DEMAND<-with_timeout(rand_vect(n_PROD,TOTAL_DEMAND,sd=cv_temp*TOTAL_DEMAND/n_PROD),
                 timeout =  5)
    if(DEMAND$message=="error"){
      cv_temp<-1
      DEMAND<-rand_vect(n_PROD,TOTAL_DEMAND,sd=cv_temp*TOTAL_DEMAND/n_PROD)
    }else{
      DEMAND<-DEMAND$res
    }
    cv_is<-(sd(DEMAND)/mean(DEMAND))
    err<-abs(cv_is-cv)/cv
    err<-ifelse(is.nan(err),0,err)
    err<-ifelse(is.infinite(err),0,err)
    if(tries>15) break
    if(err<=0.1 & all(DEMAND>0) & sum(DEMAND)==TOTAL_DEMAND){
      break
    }else if(cv_is<cv){
      cv_temp<-cv_temp*(1+err/2)
    }else if(cv_is>cv){
      cv_temp<-cv_temp*(1+err/2)
    }
  }

  CHECK<-list(CV=cv_is,
              DMD_T10=measure_TOP10(DEMAND))
  out<-list(DEMAND=DEMAND,CHECK=CHECK)

  return(out)

}
