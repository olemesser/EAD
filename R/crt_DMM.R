#' @title Create a domain mapping matrix (DMM)
#' @description Creates a \code{N_src} x \code{N_tgt}  domain mapping matrix (DMM) by reaching a certain given SDC.
#' @param N_src Number of rows for the DMM.
#' @param N_tgt Number of columns for the DMM.
#' @param DMM_PAR Double within \code{0<=DMM_PAR<=1}. Specifies the desired standardized system design complexity which should be generated. For further details on measure calculation see: \link[EAD]{measure_designComplexity}.
#' @param binary Boolean, default=F. If the output DMM should be an binary or non binary matrix. If \code{binary=TRUE} then the input \code{'upper_Bound'} will be ignored.
#' @param upper_Bound Integer, default \code{upper_Bound=20}, specifies the maximum number in the DMM if \code{binary=FALSE}
#' @param allowZero Boolean, default \code{allowZero=T}, specifies whether zero columns are allowed in the DMM matrices.
#' @return A list object containing the product mix. This list consists of:
#' \describe{
#'   \item{DMM}{The created DMM}
#'   \item{DMM_bin}{The DMM in a binary version}
#'   \item{SDC}{The created SDC.}
#'   \item{SDC_err}{Percentage error between input SDC and created SDC.}
#' }
#' @examples
#' set.seed(1234)
#'
#' crt_DMM(N_src=7,
#'         N_tgt=16,
#'         DMM_PAR=0.4,
#'         upper_Bound=1)
#'
crt_DMM<-function(N_src,
                  N_tgt,
                  DMM_PAR,
                  binary=F,
                  upper_Bound=20,
                  allowZero=T){
  require(EAD)
  #### Testing ####
  # N_src<-7
  # N_tgt<-10
  # binary=F
  # upper_Bound=20
  # DMM_PAR<-.1
  # allowZero=F
  #### END Testing ###


  repeat{
    DMM<-matrix(0,N_src,N_tgt)
    diag(DMM)<-1
    if(allowZero==F){
      empty_col <- which(colSums(DMM)==0)
      if(length(empty_col)>0){
        DMM[,empty_col] <- apply(DMM[,empty_col,drop=F],2,function(x){
          temp<-rep(0,NROW(DMM))
          temp[sample(1:NROW(DMM),1)]<-1
          return(temp)
        })
      }
    }
    sdc_candidate<-measure_designComplexity(DMM,norm=T)
    while (sdc_candidate<DMM_PAR){
      DMM <- DMM
      DMM[sample(which(DMM==0),1)] <- 1
      sdc_candidate <- measure_designComplexity(DMM,norm=T)
      if(sdc_candidate>=DMM_PAR) break
    }
      cond_1<-all(rowSums(DMM)>0)
      if(allowZero){
        cond_2<-T
        cond_3<-T
      }else{
        cond_2<-all(colSums(DMM)>0)
        cond_3<-T
      }

      ## break only if each row and column has at least one entry
      if(cond_1 & cond_2 & cond_3) break
  }

  DMM_bin<-DMM
  if(binary) upper_Bound<-1
  DMM[DMM>0]<-sample(c(1:upper_Bound),size=sum(DMM>0),replace = T)
  SDC_n<-measure_designComplexity(DMM_bin,norm = T)
  SDC_err <- (SDC_n-DMM_PAR)
  SDC_err<-ifelse(is.nan(SDC_err),0,SDC_err)


  #### Write Output object ####
  output<-list(
    DMM = DMM,
    DMM_bin = DMM_bin,
    SDC = measure_designComplexity(DMM_bin,norm = F),
    SDC_n = SDC_n,
    JSDC = measure_JSDC(DMM_bin),
    R = measure_reangularity(DMM_bin),
    S = measure_semiangularity(DMM_bin),
    SDC_err= SDC_err
  )

  gc()
  return(output)
}
