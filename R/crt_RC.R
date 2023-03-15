#' @title Create a resource cost object
#' @description This function creates three resource cost vectors according to traditional costing systems.
#' In a first step, a direct cost vector is created by generate with the total costs of  \code{TC*(1-r_in)} and a coefficient of variation of \code{cv}.
#' The remaining \code{TC*r_in} costs are indirect costs.
#' Both, the variable and fixed indirect cost vector are created via sampling of values. The final vectors have a correlation of \code{cor_fix} and \code{cor_var} regarding the direct cost vector.
#' The total fixed an indirect costs are defined as \code{(TC*r_in*r_fix)}.
#' @param N_RD Number of resources.
#' @param TC Total costs.
#' @param r_in Proportion of indirect costs on total costs.
#' @param r_fix Proportion of fixed costs on total indirect costs.
#' @param cor_var Correlation between direct cost vector and variable indirect cost vector.
#' @param cor_fix Correlation between direct cost vector and fixed indirect cost vector.
#' @param cv Coefficient of variation for resource cost distribution.
#' @return A list containing the three cost vectors. For each vector the coefficient of variation as well as the top 10% largest resource costs for each vector are calculated.
#' Since the procedure is not able to match the exact input values, the final correlation value are measured. \code{cor_fix} measures the correlation between indirect fixed costs and direct cost vector.
#' \code{cor_var} between indirect variable costs and direct cost vector and \code{cor_fix_var} between both indirect cost vectors.
#' @examples
#' set.seed(1234)
#'
#' RC<-crt_RC(N_RD=50,
#'            TC=10^6,
#'            r_in = 0.5,
#'            r_fix=0.4,
#'            cor_var=0.5,
#'            cor_fix=0.5,
#'            cv=0.2)
#'
crt_RC<-function(N_RD,
                 TC=10^6,
                 r_in,
                 r_fix,
                 cor_var,
                 cor_fix,
                 cv,
                 max_tries=15){
  suppressMessages(suppressWarnings(require(faux)))
  suppressMessages(suppressWarnings(require(fGarch)))
  suppressMessages(suppressWarnings(require(GA)))
  #### Input for Testing ####
  # N_RD=50
  # TC=10^6
  # r_in=0.3
  # r_fix=0.5
  # cor_var=0
  # cor_fix=0
  # cv=2
  # max_tries=2000
  #### End Input for Testing ####

  #### 0. Initial Checks ####
  if(r_in>0.9) stop("You selected more than 90% of the costs as indirect. The maximum is proportion is 90%. Please reduce r_in!")
  if(r_in<0) stop("The proportion of indirect costs cannot be negative. Please set r_in>=0!")
  if(r_fix<0 | r_fix>1) stop("r_fix must be within the bound of 0<=r_fix<=1")
  if(cv==0) cv<-0.1

  #### 1. Split Costs according the input ratios ####
  TC_direct <- TC * (1-r_in)
  TC_indirect <- TC * (r_in)
  TC_fix <- TC_indirect * r_fix
  TC_var <- TC_indirect - TC_fix

  #### 2. Create Direct Cost Vector ####
  RC_direct<- sample_RC(N_RD = N_RD,
                        TC = TC_direct,
                        cv = cv,
                        max_run = max_tries)

  cv_direct_is <- sd(RC_direct)/mean(RC_direct)

  #### 3. Create Indirect cost vector ####
    ### 3.1 fixed costs ###
    fix <- crt_corVEC(TC = TC_fix,
                     RC_base = RC_direct,
                     cv = cv,
                     cor = cor_fix,
                     max_tries = 50)

    ### 3.2 variable costs ###
    var <- crt_corVEC(TC = TC_var,
                      RC_base = RC_direct,
                      cv = cv,
                      cor = cor_var,
                      max_tries = 50)

  #### 4. Create Output Object ####
    output<-list(RC_var = list(RC = var$RC,
                             cv = var$cv,
                             top10 = measure_TOP10( var$RC)),
                 RC_fix = list(RC = fix$RC,
                             cv = fix$cv,
                             top10 = measure_TOP10(fix$RC)),
                 RC_direct = list(RC = RC_direct,
                                  cv = cv_direct_is,
                                  top10 = measure_TOP10(RC_direct)),
                 cor_fix = suppressWarnings(cor(fix$RC,RC_direct)),
                 cor_var = suppressWarnings(cor(RC_direct,var$RC)),
                 cor_fix_var = suppressWarnings(cor(fix$RC,var$RC)),
                 cor_indirect = suppressWarnings(cor((fix$RC+var$RC),RC_direct)),
                 r_fix = r_fix,
                 r_in = r_in,
                 TC = sum(var$RC+fix$RC+RC_direct))

  return(output)

}

sample_RC<-function(N_RD,TC,cv,max_run=1000){
  x<-rep(0,N_RD)
  x[-1]<-TC/length(x)*.001
  x[1]<-TC-sum(x)
  x<-x
  run<-0
  repeat{
    rm_costs<-x[1]*runif(1,min=0.001,max=0.01)
    if(x[1]-rm_costs<0) break
    idx<-sample(2:length(x),1)
    x[idx]<-x[idx]+rm_costs
    x[1]<-x[1]-rm_costs
    cv_is <- sd(x)/mean(x)
    check_cv<-abs(cv_is-cv)/cv
    run<-run+1
    if((check_cv<=0.05 | run>max_run) & all(x>0)) break
  }
  x<-x[sample(1:length(x))]
  return(x)
}

obj_DMD<-function(x,N_RD,TC,cv){
  cv_is<-sapply(1:1,function(i){
    set.seed(1)
    RC<-fGarch::rsnorm(N_RD ,mean = x[1], sd = x[2], xi=x[3])
    RC<-abs(RC)
    RC<-RC/sum(RC)*TC
    cv_is<-sd(RC)/mean(RC)
    set.seed(NULL)
    return(cv_is)
  })
  cv_is<-mean(cv_is)
  return(abs(cv-cv_is)*-1)
}


crt_corVEC<-function(TC,
                     RC_base,
                     cv,
                     cor,
                     max_tries=50){

  tries<-0
  add<-0
  N_RD<-length(RC_base)
  if(TC==0){
    RC<-rep(0,N_RD)
    cv_is<-0
    cor_is<-0
  }else{
    repeat{
      tries<-tries+1
      RC<-faux::rnorm_pre(RC_base,r=cor+add,mu=TC/N_RD,sd=cv*TC/N_RD)
      RC<-abs(RC)
      RC<-RC/sum(RC)*TC
      cv_is<-sd(RC)/mean(RC)
      cor_is<-cor(RC_base,RC)
      if((cor_is-cor<=0.1 & cor_is-cor>=-0.1) | max_tries==tries){
        break
      }else{
        add<-add+0.05
        if(cor+add>=1) add<-1-cor
      }
    }
  }
  return(list(RC=RC,
              cv=cv_is,
              cor=cor_is))
}
