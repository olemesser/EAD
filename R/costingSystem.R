#' @title Calculate the Reported Product Costs (PC_H)
#' @description This function calculates the reported product costs using ABC product costing systems.
#' Methods which are used in \insertCite{BLH2011;textual}{EAD} are implimented.
#' @param RES_CONS_PAT The resource consumption pattern matrix.
#' @param DMD A demand vector. If no demand vector is supplied, each product is assumed to have a demand of one.
#' @param RC_direct The direct resource cost vector.
#' @param ACP The number of activity cost pools.
#' @param RD The resource driver method to allocate resource costs to activity cost pools. The following methods are available:
#' \describe{
#'  \item{random}{Resources are randomly assigned to activity cost pools.}
#'  \item{size-random}{Focus on the resources with the largest monetary value. First, the ACP largest resources are assigned to individual cost pools. The remaining resources are random-ly distributed across ACP.}
#'  \item{size-misc}{Focus on the resources with the largest monetary value. First, the ACP-1 largest resources are assigned to individual cost pools. The remaining resources are assigned to a miscellaneous cost pool.}
#'  \item{correl-random}{First ACP resources are randomly assigned to cost pools. The remaining resources are assigned to the ACP with the highest corre-lation.}
#'  \item{correl-size}{Such as heuristic four, except the initial assignment is done by resource’s size rather than random choice.}
#' }
#' @param AD The activity driver used in the second stage allocation.
#' \describe{
#'  \item{big-pool}{The big-pool heuristic uses the largest resource as a driver for each activity cost pool}
#'  \item{indexed-N}{More information demanding is the indexed method which builds a composite driver for the N-largest resources. See examples how to use this allocation method.}
#' }
#' @return Returns the reported product costs vector (PC_H)
#' @references
#' \insertAllCited{}
#' @examples
#'
#' ### Calculate the reported product costs ###
#' data("exampleEAD")
#'   EAD <- exampleEAD
#'   RES_CONS_PAT <- EAD[[4]][[1]]$P$RD
#'   RC_indirect <- EAD[[4]][[1]]$RC$var + EAD[[4]][[1]]$RC$fix
#'   DMD <-  EAD[[4]][[1]]$DEMAND
#'   ACP<-1
#'   RD<-"random"
#'   AD<-"indexed-2"
#'
#'   costingSystem_ABC(RES_CONS_PAT,
#'                    RC_indirect,
#'                    DMD,
#'                    RD=RD,
#'                    AD=AD,
#'                    ACP=10)
costingSystem_ABC<-function(RES_CONS_PAT,
                        RC_indirect,
                        DMD,
                        RD=c("random","size-random","size-misc","correl-random","correl-size"),
                        AD=c("big-pool","indexed-2","indexed-3"),
                        ACP=10){

  #### Input for Testing ####
  # EAD <- exampleEAD
  # RES_CONS_PAT <- EAD[[4]][[1]]$P$RD
  # RC_indirect <- EAD[[4]][[1]]$RC$var + EAD[[4]][[1]]$RC$fix
  # DMD <-  EAD[[4]][[1]]$DEMAND
  # ACP<-1
  # RD<-"random"
  # AD<-"big-pool"
  # AD<-"indexed-3"
  #### END Testing ####

  #### Stage (1/2) - Allocation to Cost Center (CC) ####
    if(RD=="random"){
      CC <- RCP_ACP_random(RC = RC_indirect,ACP = ACP)
    }else if(RD=="size-random"){
      CC <- RCP_ACP_sizerandom(RC = RC_indirect,ACP = ACP)
    }else if(RD=="size-misc"){
      CC <- RCP_ACP_sizemisc(RC = RC_indirect,ACP = ACP)
    }else if(RD=="correl-random"){
      CC <- RCP_ACP_correl(RC = RC_indirect,
                           ACP = ACP,
                           RES_CONS_PAT = RES_CONS_PAT,
                           method="random")
    }else if(RD=="correl-size"){
      CC <- RCP_ACP_correl(RC = RC_indirect,
                           ACP = ACP,
                           RES_CONS_PAT = RES_CONS_PAT,
                           method="size")
    }

  #### Stage (2/2) - Allocation across products ####
  ACT_CONS_PAT<-stage_twoAllocation(CC,RES_CONS_PAT,AD)
  ACP <- ACT_CONS_PAT$ACP
  ACT_CONS_PAT <- ACT_CONS_PAT$ACT_CONS_PAT
  TAC <- colSums(ACT_CONS_PAT * DMD)
  PC_H<-as.numeric(ACT_CONS_PAT %*% (ACP/TAC))

  # sum(PC_H * DMD) == sum(RC_indirect)
  return(PC_H)
}

RCP_ACP_random<-function(RC,ACP,exclude=NULL){
  NUMB_RC<-length(RC)
  repeat{
    RC_to_ACP<-split(sample(c(1:NUMB_RC),NUMB_RC),
                     c(1:ACP,sample(c(1:ACP),NUMB_RC-ACP,replace = T)))
    CC<-lapply(RC_to_ACP,function(x){
      return(list(mapping=x,
                  size=RC[x]))
    })
    names(CC)<-NULL

    # break if every ACP is greater then 0
    if(all(sapply(CC,function(x) sum(x$size)>0 & length(x$mapping)>0)) & length(RC_to_ACP)==ACP){
      break
    }
  }
  if(!is.null(exclude)){
    CC <- lapply(CC,function(x){
            x$size <- x$size[!(x$mapping %in% exclude)]
            x$mapping <- x$mapping[!(x$mapping %in% exclude)]
            return(x)
          })
  }
  return(CC)
}

RCP_ACP_sizerandom<-function(RC,ACP){
  ## largest ACP resources assigned to number of activity pools chosen by the system designer (i.e., ACP).
  ## The costs of remaining resources are assigned randomly to the activity pools.

  #### Pre assignment of first ACP resources ####
  RCs<-sort(RC,decreasing = TRUE,index.return=TRUE)   # sorted Resource cost vector
  output <- lapply(1:ACP,function(x) list(mapping=RCs$ix[x],
                                          size=RCs$x[x]))
  exclude <- RCs$ix[1:ACP]


  #### Assign remaining resources randomly ####
  CC <- RCP_ACP_random(RC,ACP,exclude = exclude)
  output <- lapply(1:ACP,function(x){
              list(mapping = c(output[[x]]$mapping,CC[[x]]$mapping),
                   size = c(output[[x]]$size,CC[[x]]$size))
            })

  return(output)
}

RCP_ACP_sizemisc<-function(RC,ACP){

  #### Pre assignment of first ACP resources ####
  RCs<-sort(RC,decreasing = TRUE,index.return=TRUE)   # sorted Resource cost vector
  if(ACP>1){
    output <- lapply(1:(ACP-1),function(x) list(mapping=RCs$ix[x],
                                            size=RCs$x[x]))
    output[[ACP]]<-list(mapping = RCs$ix[ACP:length(RCs$ix)],
                        size = RCs$x[ACP:length(RCs$x)])
  }else{
    output<-list(list(mapping=RCs$ix,
                 size=RCs$x))
  }

  return(output)

}

RCP_ACP_correl<-function(RC,ACP = ACP,RES_CONS_PAT,method=c("random","size")){
  # Pick ACP resources randomly and allocate one each to the number of activity pools chosen by the system designer (ACP).
  # For the first activity pool, select those resources with the highest correlation with the resource in the pool.
  # Assign a total of INT(RCP/ACP) resources to this pool. Repeat for the second activity pool and so on.
  #### Randomly Assigned ####
  if(method=="random"){
    assigned <- sample(1:length(RC),ACP)
    CC <- lapply(1:length(assigned),function(x){
      list(mapping=assigned[x],
           size=RC[assigned[x]])
    })
  }else if(method=="size"){
    #### Pre assignment of first ACP resources ####
    RCs<-sort(RC,decreasing = TRUE,index.return=TRUE)   # sorted Resource cost vector
    CC <- lapply(1:ACP,function(x) list(mapping=RCs$ix[x],
                                            size=RCs$x[x]))
    assigned <- RCs$ix[1:ACP]
  }

  #### Correlation Based Assigned ####
  n_assigned<-setdiff(1:length(RC),assigned)
  if(length(n_assigned)>0){
    for (x in n_assigned) {
      cor_vec <- cor(RES_CONS_PAT[,x],RES_CONS_PAT)
      cor_vec <- sort(cor_vec,decreasing = TRUE,index.return=TRUE)
      cor_vec$x<-cor_vec$x[cor_vec$ix %in% assigned]
      cor_vec$ix<-cor_vec$ix[cor_vec$ix %in% assigned]
      idx <- which(assigned==cor_vec$ix[1])
      CC[[idx]] <- list(mapping=c(CC[[idx]]$mapping,x),
                        size=c(CC[[idx]]$size,RC[x]))
    }
  }
  return(CC)
}

stage_twoAllocation<-function(CC,RES_CONS_PAT,AD=c("big-pool","indexed-2","indexed-3")){

  ## reorder cost pools ##
  CC<-lapply(CC,function(x){
    x$mapping<-x$mapping[order(x$size,decreasing = T)]
    x$size<-x$size[order(x$size,decreasing = T)]
    return(x)
  })
  ACP<-sapply(CC, function(x) sum(x$size))


  if(AD=="big-pool"){
    ACT_CONS_PAT<-sapply(CC,function(x){
        RES_CONS_PAT[,x$mapping[1]] # big pool
    })
  }else if (startsWith(AD,"indexed")){
    i<-as.numeric(strsplit(AD,"indexed-")[[1]][2])
    ACT_CONS_PAT<-sapply(CC,function(x){
      i_max<-ifelse(length(x$mapping)<i,length(x$mapping),i)
      rowMeans(RES_CONS_PAT[,x$mapping[1:i_max],drop = F])
    })
  }

  return(list(ACT_CONS_PAT = ACT_CONS_PAT,
              ACP = ACP))
}


#' @title Calculate the Reported Product Costs (PC_H)
#' @description This function calculates the reported product costs using a simple volume-based product costing system.
#' @param RES_CONS_PAT The resource consumption pattern matrix.
#' @param DMD A demand vector. If no demand vector is supplied, each product is assumed to have a demand of one.
#' @param RC_direct The direct resource cost vector.
#' @param RC_indirect The indirect resource cost vector.
#' @param method The allocation method. The following methods are available:
#' \describe{
#'  \item{PU}{Indirect costs are allocated equally across products.}
#'  \item{DC-0.2}{Indirect costs are allocated based on the direct costs. The index 0.2 indicates that the 20 percentage largest direct costs are used as an estimation basis}
#' }
#' @return Returns the reported product costs vector (PC_H)
#' @examples
#'
#' ### Calculate the reported product costs ###
#' data("exampleEAD")
#'   EAD <- exampleEAD
#'   RES_CONS_PAT <- EAD[[4]][[1]]$P$RD
#'   RC_indirect <- EAD[[4]][[1]]$RC$var + EAD[[4]][[1]]$RC$fix
#'   RC_direct <- EAD[[4]][[1]]$RC$direct
#'   DMD <-  EAD[[4]][[1]]$DEMAND
#'
#'   costingSystem_VD(RES_CONS_PAT,
#'                    RC_direct,
#'                    RC_indirect,
#'                    DMD,
#'                    method="PU")
costingSystem_VD<-function(RES_CONS_PAT,
                           RC_direct,
                           RC_indirect,
                           DMD,
                           method=c("DIV","DLH-0.8")){
  #### Input for Testing ####
  # EAD <- exampleEAD
  # RES_CONS_PAT <- EAD[[4]][[1]]$P$RD
  # RC_indirect <- EAD[[4]][[1]]$RC$var + EAD[[4]][[1]]$RC$fix
  # RC_direct <- EAD[[4]][[1]]$RC$direct
  # DMD <-  EAD[[4]][[1]]$DEMAND
  # method <- "DC-0.1"
  # PC_B<-clc_PCB(RES_CONS_PAT,DMD,RC_indirect)$PC_B
  #### END Testing ####

  #### Method Selection ####
  if(method=="DIV"){
    ### PU - indirect costs are allocated equally across products
    PC_H <- rep(sum(RC_indirect)/sum(DMD),length(DMD))
  }else if(startsWith(as.character(method),"DLH")){
    ### DC - indirect costs are allocated based on the direct costs
    p_res<-as.numeric(strsplit(as.character(method),"DLH-")[[1]][2])
    p_res<-ceiling(NCOL(RES_CONS_PAT)*p_res)
    TRC<-colSums(RES_CONS_PAT * DMD)
    RCU <- RC_direct/TRC
    RCU<-ifelse(is.nan(RCU) | is.infinite(RCU),0,RCU)
    RCU[-order(RC_direct,decreasing = T)[1:p_res]]<-0
    PC_direct <- (RES_CONS_PAT %*% RCU) * DMD
    PC_direct <- PC_direct/sum(PC_direct)
    PC_H <- PC_direct * sum(RC_indirect) / DMD
  }

  # sum(PC_H * DMD) == sum(RC_indirect)
  return(as.numeric(PC_H))

}
