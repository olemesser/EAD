#' @title Create a domain mapping matrix (DMM)
#' @description Creates a \code{N_src} x \code{N_tgt}  domain mapping matrix (DMM) by reaching a certain given SDC.
#' @param N_src Number of rows for the DMM.
#' @param N_tgt Number of columns for the DMM.
#' @param DMM_PAR Double within \code{0<=DMM_PAR<=1}. Specifies either the system design complexity which should be generated \code{method='SDC'}.
#' If \code{method='DENS'} is selected, DMM_PAR percentage of entries are non-zero. If \code{DMM_PAR=0} a uncoupled design is created under both methods.
#' @param method Character which specifies the method used for DMM generation. The options are:
#' \describe{
#' \item{SDC}{Uses a genetic algorithm to generate a DMM where \code{measure_designComplexity(DMM)=DMM_PAR}. For further details see \link[EAD]{measure_designComplexity}}
#' \item{DENS}{Randomly samples DENS percentage of entries within the DMM.}
#' }
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
#'         method='SDC',
#'         DMM_PAR=0.4,
#'         upper_Bound=1)
#'
crt_DMM<-function(N_src,
                  N_tgt,
                  DMM_PAR,
                  method=c('SDC','DNS'),
                  binary=F,
                  upper_Bound=20,
                  allowZero=T){
  require(EAD)
  require(dplyr)
  require(data.table)
  #### Testing ####
  # N_src<-7
  # N_tgt<-16
  # binary=F
  # upper_Bound=20
  # DMM_PAR<-0.2
  # allowZero=T
  # method="DENS"
  #### END Testing ###


  if(method=='DNS'){
    tries<-0
    DMM_start<-DMM_PAR
    repeat{
      tries<-tries+1
      if(DMM_PAR==0){
        DMM_PAR<-(DMM_PAR+0.02)
      }
      # each 100 tries increase DMM_PAR by 5%
      if(tries %% 500 ==0) DMM_PAR<-DMM_PAR*1.05
      if(DMM_PAR>1) DMM_PAR<-1
      if(DMM_start==0 & N_src<=N_tgt){
          DMM_left<-diag(1,N_src)
          if(allowZero){
            DMM_right<-DMM_right<-matrix(0,nrow = N_src,ncol = N_tgt-N_src)
          }else{
            DMM_right<-matrix(sample(c(0,1),
                               N_src*N_tgt,
                               replace = T,
                               prob = c(1-DMM_PAR,DMM_PAR)),
                        nrow=N_src,ncol=N_tgt)
          }
          DMM<-cbind(DMM_left,DMM_right)
        }else{
          DMM<-matrix(sample(c(0,1),
                             N_src*N_tgt,
                             replace = T,
                             prob = c(1-DMM_PAR,DMM_PAR)),
                      nrow=N_src,ncol=N_tgt)
        }
      cond1<-all(rowSums(DMM)>0)
      cond2<-ifelse(allowZero,TRUE,all(colSums(DMM)>0))
      if(cond1 & cond2) break
    }
  }else if(method=='SDC'){
    require(GA)
    require(Matrix)
    ## estimate initial suggestion ##
    if(allowZero){
      DMM_left<-diag(0,N_src)
      DMM_left<-t(apply(DMM_left,1,function(x){
        x[sample(c(1:NCOL(DMM_left)),1)]<-1
        return(x)
      }))
    }else{
      DMM_left<-diag(1,N_src)
    }

    ## run optimization ##
    if(N_src>=N_tgt){
      x_init<-DMM_left[,1:N_tgt]
    }else if(N_src<N_tgt){
      DMM_right<-matrix(0,nrow = N_src,ncol = N_tgt-N_src)
      x_fill<-sample(1:NROW(DMM_right),replace = T)
      if(NCOL(DMM_right)-NROW(DMM_right)>0 & allowZero==F){
        x_missing<-sample(1:NROW(DMM_right),size=NCOL(DMM_right)-NROW(DMM_right),replace = T)
      }else{
        x_missing<-0
      }
      x_fill<-c(x_fill,x_missing)
      remove(x_missing)
      for(j in 1:NCOL(DMM_right)) DMM_right[x_fill[j],j]<-1
      x_init<-cbind(DMM_left,DMM_right)
    }
    x_init<-as.numeric(x_init)
    dens<-(exp(DMM_PAR)-1)/(exp(1)-1)
    x_init[x_init==0]<-sample(0:1,size=sum(x_init==0),replace = T,prob = c(1-dens,dens))
    repeat{
      x<-ga(type = "binary",
            fitness = dmm_objFCT,
            N_src = N_src,
            N_tgt = N_tgt,
            DMM_PAR = DMM_PAR,
            allowZero = allowZero,
            nBits = length(x_init),
            suggestions=x_init,
            popSize = 200,
            run=75,
            maxiter = 200,
            monitor = F,
            keepBest = T,
            pcrossover = 0.5,
            pmutation = 0.1,
            maxFitness = -0.01)

      DMM<-matrix(as.numeric(x@solution[1,]),nrow = N_src,ncol = N_tgt)

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
  }

  DMM_bin<-DMM
  if(binary) upper_Bound<-1
  DMM[DMM>0]<-sample(c(1:upper_Bound),size=sum(DMM>0),replace = T)
  SDC_n<-measure_designComplexity(DMM_bin,norm = T)
  SDC_err <- (SDC_n-DMM_PAR)/DMM_PAR
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

dmm_objFCT<-function(x,N_src,N_tgt,DMM_PAR,allowZero){
  require(Matrix)
  A<-matrix(x,nrow = N_src,ncol = N_tgt)
  err<-abs(measure_designComplexity(A)-DMM_PAR)
  cond_1<-all(rowSums(A)>0)
  if(allowZero){
    cond_2<-T
  }else{
    cond_2<-all(colSums(A)>0)
  }
  # cond_3<-rankMatrix(A)==min(dim(A))
  cond_3<-T
  if(cond_1==F | cond_2==F | cond_3==F) err<-(err+1)*10^3
  return(err*(-1))
}

