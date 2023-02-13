#' @title Create product mix (P_FD)
#' @description Creates a product mix based on a given number of functional requirements \code{N_FR} and a density \code{PARAM}.
#' @param N_FR Number of functional requirements. It is recommend to keep the value between 3-15.
#' Product mixes with much more than 15 functional requirements require much RAM as all combinations are generated. Therefore, the function is stopped if \code{N_FR>15}
#' @param PARAM A numeric value between \code{0<PARAM<=1} where its meaning depend on the selected \code{'method'}.
#' @param method Character which specifies the method, used to generate the product mix. There are two options:
#' \describe{
#'   \item{combination}{The product mix is created by using the DSM_FD matrix which contains the configuration constraints in conjuntive normal form.}
#'   \item{random}{Randomly choose \code{PARAM} percentage of products from the product mix. Is is ensured, that each functional requirement is used at least ones.}
#'   \item{PCI}{Creates a product mix based on a pre defined product line commonality index (PCI). The density must be within the bounds of \code{0<PARAM<=0.5}.}
#'   \item{DNS}{Creates a product matrix with a given density.}
#' }
#' @return A list object containing the product mix. This list consists of:
#' \describe{
#'   \item{P_FD}{The free product mix matrix}
#'   \item{P_FD_const}{The constraint product mix matrix}
#'   \item{DSM_FD}{The DSM_FD matrix containing the entries to constraint the variety.}
#'   \item{measures}{A named list containing the following measures: number of product variants (N_P), uniform diversification index (D_u),
#'     local outlier factor (LOF) and Product Variety (PV) }
#' }
#' @examples
#' set.seed(1234)
#' prodMIX<-create_ProductMix(N_FR=7,PARAM=0.2,method="random")
#' prodMIX$P_FD_const
create_ProductMix<-function(N_FR=13,
                            PARAM=0.05,
                            method=c("combination","random","PCI","DNS")){

  suppressWarnings(require(dplyr))
  suppressWarnings(require(tidyr))

  if(N_FR>15) stop("The maximum N_FR values is hardcoded to 15. Please reduce N_FR!")
  #### 1.1 Create Free Product Mix ####
  P_FD<-expand.grid(lapply(1:N_FR,function(x) c(0,1)))
  P_FD<-P_FD[rowSums(P_FD)!=0,]
  colnames(P_FD)<-paste0("FD_",1:N_FR)
  if(length(PARAM)==2) PARAM<-runif(1,min=PARAM[1],max=PARAM[2])

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
          prob = c(1-PARAM,PARAM)
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
    if(PARAM==0) PARAM<-0.01
    repeat{
      selection<-sample(1:NROW(P_FD),ceiling(NROW(P_FD)*PARAM))
      P_FD_const<-P_FD[selection,]
      tries<-tries+1
      if(all(colSums(P_FD_const)>0) & NROW(P_FD_const)>2){
        DSM_FD<-diag(0,nrow = NCOL(P_FD))
        break
      }else if(tries>50){
        PARAM<-PARAM*1.01
      }
    }
  }else if(method=="PCI"){
    suppressWarnings(require(GA))
    DSM_FD<-diag(0,nrow = NCOL(P_FD))
    x_init<-sample(c(0,1),NROW(P_FD),replace = T)
    x<-ga(type = "binary",
          fitness =  objFct_productMix,
          P = P_FD,
          PARAM = PARAM,
          nBits = length(x_init),
          suggestions=x_init,
          popSize = 150,
          run=100,
          maxiter = ifelse(PARAM<0.25,4000,2000),
          monitor = F,
          keepBest = T,
          pcrossover = 0.02,
          pmutation = 0.05,
          maxFitness = -0.01)
    P_FD_const<-P_FD[as.logical(x@solution[1,]),]
  }else if(method=="DNS"){
    if(PARAM<=0.07) PARAM<-0.07
    suppressMessages(suppressWarnings(require(GA)))
    DSM_FD<-diag(0,nrow = NCOL(P_FD))
    dns_temp<-as.numeric(rowSums(P_FD)/NCOL(P_FD))
    range<-function(x){
      y<-(x^2-(1)*x+0.39)
      return(y*8)
    }
    P_FD<-P_FD[dns_temp<=PARAM*range(PARAM) & dns_temp>=PARAM*1/range(PARAM),]
    max_p <- 200
    prob<-ifelse(max_p>NROW(P_FD),NROW(P_FD),max_p)/NROW(P_FD)
    x_init<-sample(c(0,1),NROW(P_FD),replace = T,prob = c(1-prob,prob))
    prod<-sum(x_init)
    if(prod<=15) x_init[sample(which(x_init==0),15-prod)]<-1
    if(prod>=200) x_init[sample(which(x_init==1),prod-max_p)]<-0
    x<-ga(type = "binary",
          fitness =  function(x,P_FD,PARAM,min_p,max_p){
            prod<-sum(x)
            temp<-P_FD[as.logical(x),,drop = FALSE]
            temp<-measure_DENS(temp)*ifelse(all(colSums(temp)>0),1,2)
            obj<-abs(mean(temp)-PARAM)*-1
            pen<-max(c(prod-max_p,min_p-prod))
            pen<-ifelse(pen>0,pen*1.1,1)
            obj<-ifelse(prod<min_p | prod>max_p,obj*pen,obj)
            return(obj)
          },
          P = P_FD,
          PARAM = PARAM,
          min_p = N_FR,
          max_p = max_p,
          nBits = length(x_init),
          suggestions=x_init,
          popSize = 40,
          run=300,
          maxiter = 2000,
          monitor = F,
          keepBest = T,
          pcrossover = 0.02,
          pmutation = 0.05,
          maxFitness = -0.02)
    P_FD_const<-P_FD[as.logical(x@solution[1,]),]
  }

  measures<-list(N_P=NROW(P_FD_const), # number of product variants
                 D_u=measure_diversificationINDEX(P_FD_const), # Diversification Index
                 LOF_10=measure_LOF(P_FD_const), # Local Outlier Factor
                 DENS_FD=sum(P_FD_const>0)/prod(dim(P_FD_const)),# density of DSM_FD
                 PCI_FD = measure_PCI(P_FD_const)) # product mix commonality

  productMIX<-list(
                    P_FD = P_FD,
                    P_FD_const = P_FD_const,
                    DSM_FD = DSM_FD,
                    measures = measures)
  gc()
  return(productMIX)

}



objFct_productMix<-function(x,P,PARAM){
  prod<-sum(x)
  if(prod<2) x[sample(which(x==0),3-prod)]<-1
  temp<-P[as.logical(x),,drop = FALSE]
  temp<-ifelse(all(colSums(temp)>0),measure_PCI(temp),measure_PCI(temp)*10)
  return(abs(temp-PARAM)*-1)
}
