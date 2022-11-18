#' @title Create a domain dependency matrix (DSM)
#' @description Creates a \code{N_el} x \code{N_el}  domain dependency matrix (DSM) by trying to reach a given SC_adj_n or a DENS value.
#' @param N_el Number of rows and columns for the DSM.
#' @param SC_adj_n_in Double with \code{0<=SC_adj_n_in<1}. Specifies the system design complexity which should be generated. For further details see \link[EAD]{measure_structuralcomplexity}.
#' @param upper_Bound Integer, default \code{upper_Bound=20}, specifies the maximum number in the DMM if \code{binary=FALSE}
#' @param DENS A numeric value specifying the density within the DSM. This is an alternative method to generate the DSM.
#' @param forced A numeric vector containing indices of elements which cannot be zero since otherwise a invalid design is created later on.
#' @return A list object containing the product mix. This list consists of:
#' \describe{
#'   \item{DSM}{The created DSM}
#'   \item{DENS_DSM}{Density of the DSM}
#'   \item{SC_adj_n}{The achieved adjusted and normalized structural complexity measure.}
#' }
#' @examples
#' set.seed(1234)
#'
#' crt_DSM(N_el=7,SC_adj_n_in=0.01,upper_Bound=5)
crt_DSM<-function(N_el,
                  SC_adj_n_in=NULL,
                  DENS=NULL,
                  upper_Bound=5,
                  forced=NULL){

  #### Input Testing ###
  # N_el<-10
  # SC_adj_n_in=NULL
  # upper_Bound<-20
  # DENS<-0.1
  # P<-expand.grid(list(c(0,1),
  #                  c(0,1),
  #                  c(0,1),
  #                  c(0,1)))
  # P<-as.matrix(P)
  # DMM<-cbind(diag(1,4),matrix(0,nrow=4,ncol = 6))
  # forced<-which(colSums(DMM)==0)
  # P %*% (( DMM %*% DSM) + DMM)
  #### END Input for Testing ####

  if(is.null(DENS) & is.null(SC_adj_n_in)) stop("Both, DENS and SC_adj_n_in are null. Please specify at least one input!")
  if(is.null(DENS) & !is.null(SC_adj_n_in)){
    require(GA)
    x_init<-rep(0,N_el^2)
    x<-ga(type = "binary",
          fitness = dsm_objFCT,
          N_el = N_el,
          SC_adj_n_in = SC_adj_n_in,
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
  }else if(is.null(SC_adj_n_in) & !is.null(DENS)){
    DSM<-matrix(sample(c(0,1),N_el^2,prob = c(1-DENS,DENS),replace = T),nrow = N_el,ncol = N_el)
    if(!is.null(forced)){
       repeat{
          z<-forced[1]
          DSM[sample(setdiff(1:NROW(DSM),forced),1),z]<-1
          forced<-forced[-1]
          if(length(forced)==0) break
        }
    }
  }else{
    stop("Both arguments (DENS,SC_adj_n_in) are supplied. Please use either DENS or SC_adj_n_in!")
  }

  diag(DSM)<-0
  DSM_bin<-DSM
  DSM[DSM>0]<-sample(1:upper_Bound,size=sum(DSM_bin>0),replace = T)

  output<-list(DSM=DSM,
               DENS_DSM=sum(DSM_bin)/N_el^2,
               SC_adj_n=measure_structuralcomplexity(DSM,norm=T),
               Q=measure_modularity(DSM),
               NE=measure_neumannEntropy(DSM),
               GCI=measure_GCI(DSM),
               MCC=measure_MCC(DSM),
               HVM=measure_HVM(DSM),
               HIC=measure_HIC(DSM))
  return(output)
}


dsm_objFCT<-function(x,N_el,SC_adj_n_in){
  A<-matrix(x,nrow = N_el,ncol = N_el)
  err<-1-abs(measure_structuralcomplexity(A,norm=T)-SC_adj_n_in)
  return(err)
}
