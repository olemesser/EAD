fixedCost_allocation<-function(RES_CONS_PAT,
                               RC_fix,
                               DMD,
                               method=c("VBA","GRP","IDX","BP"),
                               NUMB_CC=10){

  #### Stage (1/2) - Allocation to Cost Center (CC) ####
    ## random assignment of RCP into CC ##
    CC<-RCP_ACP_random(RC_fix = RC_fix,
                       NUMB_CC = NUMB_CC)
  if(method=="VBA"){
    ACT_CONS_PAT<-matrix(1,ncol = NUMB_CC,nrow = nrow(RES_CONS_PAT))
  }else if(method=="GRP"){
    ACT_CONS_PAT<-sapply(CC,function(x){
      if(length(x$mapping)>1){
        y<-rowMeans(RES_CONS_PAT[,x$mapping])
      }else{
        y<-RES_CONS_PAT[,x$mapping]
      }
      return(as.numeric(cut(y,3)))
    })
  }else if(method=="BP"){
    ACT_CONS_PAT<-sapply(CC,function(x){
      if(length(x$mapping)>1){
        RES_CONS_PAT[,x$mapping[which.max(RC_fix[x$mapping])]] # big pool
      }else{
        RES_CONS_PAT[,x$mapping]
      }
    })
  }else if(method=="IDX"){
    ACT_CONS_PAT<-sapply(CC,function(x){
      if(length(x$mapping)>1){
        rowMeans(RES_CONS_PAT[,x$mapping]) # indexed driver full average
      }else{
        RES_CONS_PAT[,x$mapping]
      }
    })
  }else{
    stop("You choose a non specified method.")
  }

  #### Stage (2/2) - Allocation across products ####
    ACP<-sapply(CC,function(x) x$size)
    DMD_mat <- matrix(DMD,
                      ncol = NCOL(ACT_CONS_PAT),
                      byrow = F,
                      nrow = NROW(ACT_CONS_PAT))
    TRC<-colSums(ACT_CONS_PAT * DMD_mat)
    RCU_fix <- ACP/TRC
    PC_H<-ACT_CONS_PAT %*% RCU_fix
  return(PC_H)
}


RCP_ACP_random<-function(RC_fix,NUMB_CC){
  NUMB_RC<-length(RC_fix)
  repeat{
    RC_to_ACP<-split(sample(c(1:NUMB_RC),NUMB_RC),
                     c(1:NUMB_CC,sample(c(1:NUMB_CC),NUMB_RC-NUMB_CC,replace = T)))
    CC<-lapply(RC_to_ACP,function(x){
      return(list(mapping=x,
                  size=sum(RC_fix[x])))
    })
    names(CC)<-NULL
    # break if every ACP is greater then 0
    if(all(sapply(CC,function(x) x$size)>0) & length(RC_to_ACP)==NUMB_CC){
      break
    }
  }
  return(CC)
}
