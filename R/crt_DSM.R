#' @title Create a domain dependency matrix (DSM)
#' @description Creates a \code{N_el} x \code{N_el}  domain dependency matrix (DSM) based on different methods.
#' @param N_el Number of rows and columns for the DSM.
#' @param method A character specifying the method used for DSM generation. The following methods are available:
#' \describe{
#'   \item{DNS}{Generates a DSM by sampling random matrices to achieve a given density. Only  \code{PARAM$DNS} is needed}
#'   \item{modular}{Generates a DSM using two conditions. The procedure tries to achieve a certain density specified via should be achieved \code{PARAM$DNS} and, second,
#'   tries to create a certain distribution for the off-diagonal elements which is specified via \code{PARAM$cv}. Low values of cv indicate a random distribution while high values indicate that values are more clustered around the DSM main diagonal.
#'   The latter represents modular designs.}
#' }
#' @param PARAM A named list object with the following entries:
#' \describe{
#'   \item{DNS}{Density for the matrix}
#'   \item{cv}{Coefficient of variation for off-diagonal entries.}
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
                  method=c('DNS','modular'),
                  PARAM=list(DNS=0.1,
                             SC=0.1,
                             modular=0),
                  upper_Bound=5,
                  forced=NULL){

  #### Input Testing ###
  # set.seed(12)
  # N_el<-6
  # upper_Bound<-20
  # PARAM=list(DNS=0.2,
  #            SC=0.1,
  #            cv=0.1,
  #            modular=0.5)
  # forced=c(1,4)
  # method<-'modular'
  # P<-expand.grid(list(c(0,1),
  #                  c(0,1),
  #                  c(0,1),
  #                  c(0,1)))
  # P<-as.matrix(P)
  # DMM<-crt_DMM(4,N_el,method='SDC',
  #              DMM_PAR=0.1,
  #              upper_Bound=1,
  #              allowZero = T)$DMM
  # forced<-which(colSums(DMM)==0)
  # P %*% ((DMM %*% DSM) + DMM)
  #### END Input for Testing ####
  if(length(forced)==0) forced<-NULL

  if(method=="DNS"){
    DSM<-matrix(sample(c(0,1),N_el^2,prob = c(1-PARAM$DNS,PARAM$DNS),replace = T),nrow = N_el,ncol = N_el)
  }else if(method=="modular"){
      DSM<-crt_DSMmod(N_el = N_el,
                      DNS = PARAM$DNS,
                      modular = PARAM$modular)
  }
  diag(DSM)<-0
  DSM<-force_DSMentries(DSM,forced = forced)
  DSM_bin<-DSM
  DSM[DSM>0]<-sample(1:upper_Bound,size=sum(DSM_bin>0),replace = T)

  measures<-list(DENS_DSM=measure_DENS(DSM),
                 SC=measure_structuralcomplexity(DSM),
                 Q=measure_modularity(DSM),
                 NE=measure_neumannEntropy(DSM),
                 MCC_n=measure_MCC(DSM,norm=T),
                 HVM_n=measure_HVM(DSM,norm=T),
                 HIC_n=measure_HIC(DSM,norm=T),
                 MCC=measure_MCC(DSM,norm = F),
                 HVM=measure_HVM(DSM,norm = F),
                 HIC=measure_HIC(DSM,norm = F))

  output<-list(DSM = DSM,
               measures = measures)
  return(output)
}



crt_DSMmod<-function(N_el,DNS,modular=0){
  require(Matrix)
  #### INput for testing ####
  # N_el<-10
  # modular<-0
  # DNS<-0.2
  #### END input testing ####

  size_bM<-rand_vect(3,N_el)
  bm<-lapply(size_bM,function(x){
    matrix(sample(0:1,x^2,replace = T,prob = c(1-modular,modular)),ncol=x)
  })
  DSM<-bdiag(bm)
  DSM<-as.matrix(DSM)
  DNS_is<-measure_DENS(DSM)
  if(DNS_is>DNS){
    el_rmv<-sum(DSM)-ceiling(NROW(DSM)^2*DNS)
    DSM[DSM>0][sample(1:sum(DSM>0),el_rmv)]<-0
  }else if(DNS_is<DNS){
    el_add<-ceiling(NROW(DSM)^2*DNS)-sum(DSM)
    DSM[DSM==0][sample(1:sum(DSM==0),el_add)]<-1
  }
  return(DSM)
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
