#' @title Create a domain dependency matrix (DSM)
#' @description Creates a \code{N_el} x \code{N_el}  domain dependency matrix (DSM) based on different methods.
#' @param N_el Number of rows and columns for the DSM.
#' @param method A character specifying the method used for DSM generation. The following methods are available:
#' \describe{
#'   \item{DNS}{Generates a DSM by sampling random matrices to achieve a given density. Only  \code{PARAM$DNS} is needed}
#'   \item{modular}{Generates a DSM using two conditions. The procedure tries to achieve a certain density specified via should be achieved \code{PARAM$DNS} and, second,
#'   tries to create a certain distribution for the off-diagonal elements which is specified via \code{PARAM$cv}. Low values of cv indicate a random distribution while high values indicate that values are more clustered around the DSM main diagonal.
#'   The latter represents modular designs.}
#'   \item{SC}{Creates a matrix in order to achieved a specified normalized structural complexity measure. For more details see: \link[EAD]{measure_structuralcomplexity}.}
#' }
#' @param PARAM A named list object with the following entries:
#' \describe{
#'   \item{DNS}{Density for the matrix}
#'   \item{cv}{Coefficient of variation for off-diagonal entries.}
#'   \item{SC}{Specifies the system design complexity which should be generated. For further details see \link[EAD]{measure_structuralcomplexity}.}
#' }
#' @param upper_Bound Integer, default \code{upper_Bound=20}, specifies the maximum number in the DMM if \code{binary=FALSE}
#' @param DENS A numeric value specifying the density within the DSM. This is an alternative method to generate the DSM.
#' @param forced A numeric vector containing indices of elements which cannot be zero since otherwise a invalid design is created later on.
#' @return A list object containing the product mix. This list consists of:
#' \describe{
#'   \item{DSM}{The created DSM}
#'   \item{measures}{A list of measures for the given DSM. The measures are: \link[EAD]{measure_DENS}, \link[EAD]{measure_structuralcomplexity},
#'   \link[EAD]{measure_modularity}, \link[EAD]{measure_neumannEntropy}, \link[EAD]{measure_MCC}, \link[EAD]{measure_HVM}, \link[EAD]{measure_HIC}}
#' }
#' @examples
#' set.seed(1234)
#' ## using the DSM method ##
#' crt_DSM(N_el=7,
#'         PARAM=list(DNS=0.01),
#'         method='DNS',
#'         upper_Bound=5)
crt_DSM<-function(N_el,
                  method=c('DNS','modular','SC'),
                  PARAM=list(DNS=0.1,
                             SC=0.1,
                             cv=0.1),
                  upper_Bound=5,
                  forced=NULL){

  #### Input Testing ###
  # N_el<-10
  # upper_Bound<-20
  # PARAM=list(DNS=0.,
  #            SC=0.1,
  #            cv=0.1)
  # forced=c(1,4)
  # method<-'modular'
  # P<-expand.grid(list(c(0,1),
  #                  c(0,1),
  #                  c(0,1),
  #                  c(0,1)))
  # P<-as.matrix(P)
  # DMM<-cbind(diag(1,4),matrix(0,nrow=4,ncol = 6))
  # forced<-which(colSums(DMM)==0)
  # P %*% (( DMM %*% DSM) + DMM)
  #### END Input for Testing ####
  if(length(forced)==0) forced<-NULL

  if(method=="SC"){
    require(GA)
    x_init<-rep(0,N_el^2)
    x<-ga(type = "binary",
          fitness = dsm_objFCT,
          N_el = N_el,
          PARAM = PARAM$SC,
          nBits = length(x_init),
          suggestions = x_init,
          popSize = 200,
          maxiter = 50,
          run=5,
          monitor = F,
          keepBest = T,
          pcrossover = 0.05,
          pmutation = 0.02,
          maxFitness = 0.9999)
    DSM<-matrix(as.numeric(x@solution[1,]),nrow = N_el,ncol = N_el)
  }else if(method=="DNS"){
    DSM<-matrix(sample(c(0,1),N_el^2,prob = c(1-PARAM$DNS,PARAM$DNS),replace = T),nrow = N_el,ncol = N_el)
  }else if(method=="modular"){
      DSM<-crt_DSMmod(N_el = N_el,
                      DNS = PARAM$DNS,
                      cv = PARAM$cv)
  }
  diag(DSM)<-0
  DSM<-force_DSMentries(DSM,forced = forced)
  DSM_bin<-DSM
  DSM[DSM>0]<-sample(1:upper_Bound,size=sum(DSM_bin>0),replace = T)

  measures<-list(DENS_DSM=measure_DENS(DSM),
                 SC_adj_n=measure_structuralcomplexity(DSM,norm=T),
                 Q=measure_modularity(DSM),
                 NE=measure_neumannEntropy(DSM),
                 MCC=measure_MCC(DSM),
                 HVM=measure_HVM(DSM),
                 HIC=measure_HIC(DSM))

  output<-list(DSM = DSM,
               measures = measures)
  return(output)
}


dsm_objFCT<-function(x,N_el,PARAM){
  A<-matrix(x,nrow = N_el,ncol = N_el)
  err<-1-abs(measure_structuralcomplexity(A,norm=T)-PARAM)
  return(err)
}



crt_DSMmod<-function(N_el,DNS,cv){
  require(Matrix)
  require(GA)
  #### INput for testing ####
  # N_el<-10
  # cv<-0.5
  # DNS<-0.2
  # n<-c(0,1)
  #### END input testing ####

  x<-ga(type = "real-valued",
        fitness = dsm_objFCT_mod,
        N_el = N_el,
        DNS = DNS,
        cv = cv,
        lower=c(0,0),
        upper=c(1,10),
        popSize = 40,
        maxiter = 100,
        run=10,
        monitor = F,
        keepBest = T,
        pcrossover = 0.05,
        pmutation = 0.02,
        maxFitness = 0.01)

  DSM<-dsm_objFCT_mod(x@solution[1,],
                      N_el,
                      DNS=DNS,
                      cv = cv,
                      opt = F)
  return(DSM)
}

dsm_objFCT_mod<-function(n,
                         N_el,
                         DNS,
                         cv,
                         opt=T){
  diag <- lapply(1:N_el,function(x){
    prob<-1/((x+n[1])^n[2])
    sample(c(0,1),size=N_el,replace = T,prob = c(1-prob,prob))
  })
  diag[[1]]<-rep(0,N_el)
  d<-sapply(2:length(diag),function(x) sum(diag[[x]]))
  cv_is<-ifelse(mean(d)==0,0,sd(d)/mean(d))

  k<-c(0:(N_el-1))
  k<-k*sample(c(-1,1),size=length(k),replace=T)
  DSM<-bandSparse(N_el,
                  k=k,
                  diagonals=diag,
                  symmetric = F)
  DSM<-as.matrix(DSM)
  if(opt){
    obj<-abs(measure_DENS(DSM)-DNS)*+abs(cv_is-cv)/5
    obj<-obj*(-1)
  }else{
    obj<-DSM
  }
  return(obj)
}


force_DSMentries<-function(DSM,forced=NULL){
  if(!is.null(forced)){
    repeat{
      z<-forced[1]
      cond<-T
      repeat{
        row<-sample(setdiff(1:NROW(DSM),forced),1)
        if(row!=z) break
      }
      DSM[row,z]<-1
      forced<-forced[-1]
      if(length(forced)==0) break
    }
  }
  return(DSM)
}
