#' @title Create product mix (P_FD)
#' @description Creates a product mix based on a given number of functional requirements \code{N_FR} and optional input values which depend on the \code{method} argument.
#' @param N_FR Number of functional requirements. It is recommend to keep the value between 3-15.
#' Product mixes with much more than 15 functional requirements require much RAM as all combinations are generated. Therefore, the function is stopped if \code{N_FR>15}
#' @param DNS A numeric value between \code{0<DNS<=1} specifying the desired density. This argument is only used if either \code{method='DNS'} or \code{method='combination'} are selected.
#' @param N_PROD The number of desired products. Note that \code{N_FR<=N_PROD<=2^(N_FR)-1}. This argument is only used if \code{method='DNS'} is selected.
#' @param prop_PROD The proportion of products to sample from the free product mix This argument is only used if \code{method='random'} is selected.
#' @param method Character which specifies the method, used to generate the product mix. There are two options:
#' \describe{
#'   \item{combination}{The product mix is created by using the DSM_FD matrix which contains the configuration constraints in conjuntive normal form.}
#'   \item{random}{Randomly choose \code{prop_PROD} percentage of products from the product mix. Is is ensured, that each functional requirement is used at least ones.}
#'   \item{DNS}{Creates a product matrix with a given density.}
#' }
#' @return A list object containing the product mix. This list consists of:
#' \describe{
#'   \item{P_FD}{The free product mix matrix}
#'   \item{P_FD_const}{The constraint product mix matrix}
#'   \item{DSM_FD}{The DSM_FD matrix containing the entries to constraint the variety.}
#'   \item{measures}{A named list containing the following measures: number of product variants \code{N_P}, uniform Diversification index \code{D_u},
#'     local outlier factor \code{LOF}, the density \code{DENS_FD} and the product line commonality index \code{PCI_FD}.
#'     For further information on the measures see \link[EAD]{measure_diversificationINDEX}, \link[EAD]{measure_LOF}, \link[EAD]{measure_DENS}, \link[EAD]{measure_PCI}.}
#' }
#' @examples
#' set.seed(1234)
#' prodMIX<-create_ProductMix(N_FR=7,prop_PROD=0.2,method="random")
#' prodMIX$measures
#' prodMIX$P_FD_const
#'
#' prodMIX<-create_ProductMix(N_FR=15,DNS=0.1,N_PROD=100,method="DNS")
#' prodMIX$measures
#' prodMIX$P_FD_const
create_ProductMix<-function(N_FR=13,
                            DNS=0.1,
                            prop_PROD=1,
                            N_PROD=100,
                            method=c("combination","random","DNS")){

  suppressMessages(suppressWarnings(require(dplyr)))
  suppressMessages(suppressWarnings(require(tidyr)))


  #### 1.1 Create Free Product Mix ####
  if(method %in% c("random","combination")){
    if(N_FR>20) stop("The maximum N_FR values is hardcoded to 20. Please reduce N_FR!")
    P_FD<-expand.grid(lapply(1:N_FR,function(x) c(0,1)))
    P_FD<-P_FD[rowSums(P_FD)!=0,]
  }


  if(length(DNS)==2) DNS<-runif(1,min=DNS[1],max=DNS[2])
  if(length(prop_PROD)==2) prop_PROD<-runif(1,min=prop_PROD[1],max=prop_PROD[2])

  if(method=="combination"){
    require(rpicosat)
    #### 1.2 Create Constraint Matrix ####
    ## ensure that there are not unused features
    repeat{
      ## create matrix ##
      DSM_FD <- matrix(
        sample(
          c(0,1),
          N_FR^2
          ,replace = T,
          prob = c(1-DNS,DNS)
        ),
        nrow = N_FR
      )
      diag(DSM_FD)<-0
      ## create conjunctive normal form
      formula<-data.frame(merge(c(1:N_FR),c(1:N_FR)),entry=as.numeric(DSM_FD)) %>%
        mutate(x_arg=ifelse(entry==1,sample(c(-1,1),n(),replace = T),NA),
               y_arg=ifelse(entry==1,sample(c(-1,1),n(),replace = T),NA)) %>%
        mutate(x_arg=x_arg*x,
               y_arg=y_arg*y) %>%
        dplyr::select(x_arg,y_arg) %>%
        tidyr::drop_na() %>%
        asplit(1)

      formula<-lapply(formula,function(x){
        names(x)<-NULL
        x<-c(x[1],x[2])
        return(x)
      })

      ## add missing literals
      formula<-c(formula,lapply(1:N_FR,function(x) c(-x,x)))
      formula<-unique(formula)

      ## check if satisfiable
      check<-apply(P_FD,1,function(x){
        c(which(x==1),-which(x==0))
      })
      check<-as.list(as.data.frame(check))

      res <- sapply(check,function(x){
        temp<-picosat_sat(formula,x)
        temp<-picosat_solution_status(temp)
        if(temp=="PICOSAT_UNSATISFIABLE"){
          return(FALSE)
        }else if(temp=="PICOSAT_SATISFIABLE"){
          return(TRUE)
        }
      })
      P_FD_const <- P_FD[as.logical(res),]
      if(all(colSums(P_FD_const)>0)) break
    }
  }else if(method=="random"){
    tries<-0
    if(prop_PROD==0) prop_PROD<-0.01
    repeat{
      selection<-sample(1:NROW(P_FD),ceiling(NROW(P_FD)*prop_PROD))
      P_FD_const<-P_FD[selection,]
      tries<-tries+1
      if(all(colSums(P_FD_const)>0) & NROW(P_FD_const)>2){
        DSM_FD<-diag(0,nrow = NCOL(P_FD))
        break
      }else if(tries>50){
        prop_PROD<-prop_PROD*1.01
      }
    }
  }else if(method=="DNS"){
    if(N_PROD<N_FR) N_PROD <- N_FR
    DSM_FD<-diag(0,nrow = N_PROD)
    product_mix <- list()
    prod <- 1
    while(prod<(N_PROD+1)){
      repeat{
        candidate <- sample(c(0,1),N_FR,replace = T,prob = c(1-DNS,DNS))
        existing_prod <- list(candidate) %in% product_mix
        if(!sum(candidate)==0 & existing_prod == FALSE){
          product_mix[[prod]] <- candidate
          prod <- prod + 1
          break
        }
      }
    }
    P_FD_const <- do.call(rbind,product_mix)
  }

  colnames(P_FD_const)<-paste0("FD_",1:N_FR)
  measures<-list(N_P=NROW(P_FD_const), # number of product variants
                 D_u=measure_diversificationINDEX(P_FD_const), # Diversification Index
                 LOF_10=measure_LOF(P_FD_const), # Local Outlier Factor
                 DENS_FD=measure_DENS(P_FD_const),# density of DSM_FD
                 PCI_FD = measure_PCI(P_FD_const)) # product mix commonality


  productMIX<-list(
                    P_FD = NA,
                    P_FD_const = P_FD_const,
                    DSM_FD = DSM_FD,
                    measures = measures)
  gc()
  return(productMIX)

}


