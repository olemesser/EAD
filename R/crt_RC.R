#' @title Create a resource cost object
#' @description This function creates two resource cost vectors.
#' A fixed resource cost vector representing the total fixed costs as proportion of total costs defined by \code{TC*ratio_fixedC}.
#' Variable costs are given as \code{TC*(1-ratio_fixedC)}.
#' @param N_RD Number of resources.
#' @param TC Total costs.
#' @param ratio_fixedC Proportion of fixed costs on total costs
#' @param RC_cor Correlation between variable cost vector and fixed cost vector
#' @param RC_cv Coefficient of variation for resource cost distribution
#' @return A list containing the fixed and variable cost vector and the measured coefficient of variation for each cost vector as well as
#' the top 10% largest resource costs for each vector.
#' @examples
#' set.seed(1234)
#'
#' RC<-crt_RC(N_RD=50,TC=10^6,ratio_fixedC=0.4,RC_cor=0.5,RC_cv=0.2)
#'
crt_RC<-function(N_RD,TC=10^6,
                 ratio_fixedC,
                 RC_cor,RC_cv,
                 max_tries=15){
  suppressWarnings(require(faux))
  suppressWarnings(require(GA))
  #### Input for Testing ####
  # N_RD=50
  # TC=10^6
  # ratio_fixedC=0
  # RC_cor=-0.5
  # RC_cv=0.5
  # max_tries=100
  #### End Input for Testing ####


  TC_fix <- TC * ratio_fixedC
  TC_var <- TC - TC_fix


  x<-ga(type = "real-valued",
        fitness = obj_DMD,
        N_RD = N_RD,
        TC_var = TC_var,
        RC_cv = RC_cv,
        lower = c(0,
                  0,
                  -1),
        upper = c(TC_var/N_RD*10000,10^7,4),
        popSize = 50,
        suggestions = c(TC_var/N_RD,10^7,1),
        run=75,
        maxiter = 200,
        monitor = F,
        keepBest = T,
        pcrossover = 0.4,
        pmutation = 0.2,
        maxFitness = -0.05)

    RC_var<-lapply(1:100,function(i){
      RC_var<-rsnorm(N_RD ,
                     mean =x@solution[1,1],
                     sd = x@solution[1,2],
                     xi=x@solution[1,3])
      RC_var<-(RC_var+abs(min(RC_var))+TC_var/N_RD*.2)
      RC_var<-RC_var/sum(RC_var)*TC_var
      cv_var<-sd(RC_var)/mean(RC_var)
      out<-list(RC_var=RC_var,
                cv=cv_var,
                err=abs(cv_var-RC_cv))
      return(out)
    })
    RC_var<-RC_var[[which.min(sapply(RC_var, function(x) x$err))]]$RC_var
    cv_var<-sd(RC_var)/mean(RC_var)


#### Generate Fixed Cost Vector ####
  ## use RC_var to generate RC_fixed vector with a given correlation ##
    tries<-0
    max_tries<-50
    add<-0
    if(TC_fix==0){
      RC_fix<-rep(0,N_RD)
      cv_fix<-0
    }else{
      repeat{
        tries<-tries+1
        RC_fix<-faux::rnorm_pre(RC_var,r=RC_cor+add,mu=TC_fix/N_RD,sd=cv_var*TC_fix/N_RD)
        RC_fix<-abs(RC_fix)
        RC_fix<-RC_fix/sum(RC_fix)*TC_fix
        cv_fix<-sd(RC_fix)/mean(RC_fix)
        cor_is<-cor(RC_var,RC_fix)
        if((cor_is-RC_cor<=0.1 & cor_is-RC_cor>=-0.1) | max_tries==tries){
          break
        }else{
          add<-add+0.05
          if(RC_cor+add>=1) add<-1-RC_cor
        }
      }
    }


#### Create Output Object ####
  output<-list(RC_var = list(RC = RC_var,
                           cv = cv_var,
                           top10 = measure_TOP10(RC_var)),
               RC_fix = list(RC = RC_fix,
                           cv = cv_fix,
                           top10 = measure_TOP10(RC_fix)),
               cor_check=cor(RC_fix,RC_var),
               ratio_fixedC=ratio_fixedC)

  return(output)

}



obj_DMD<-function(x,N_RD,TC_var,RC_cv){
  suppressWarnings(library(fGarch))
  cv_var<-sapply(1:50,function(i){
    RC_var<-rsnorm(N_RD ,mean = x[1], sd = x[2], xi=x[3])
    RC_var<-abs(RC_var)
    RC_var<-RC_var/sum(RC_var)*TC_var
    cv_var<-sd(RC_var)/mean(RC_var)
    return(cv_var)
  })
  cv_var<-mean(cv_var)
  return(abs(RC_cv-cv_var)*-1)
}
