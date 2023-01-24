processVariety<-function(EAD,
                         includedProd=NULL,
                         prop_setupChange=c(0.05,0.1)){

  if(is.null(includedProd)) includedProd<-1:NROW(EAD[[1]]$P$PrD)

  variety<-apply(EAD[[1]]$P$PrD[includedProd,],2,function(t){
    length(unique(t[t>0]))
  })

  RD<-diag(1,ncol = NCOL(EAD[[1]]$P$PrD),nrow = NCOL(EAD[[1]]$P$PrD)) %*% (( EAD[[1]]$DMM$PrD_RD %*% EAD[[1]]$DSM$RD) + EAD[[1]]$DMM$PrD_RD)

  costs<-clc_PCB(RD,DMD = 1,RC_var = EAD[[1]]$RC$var)
  setupCosts <- costs$PC_B %*% runif(10,min = prop_setupChange[1], max = prop_setupChange[2])

  i<-1
  lapply(1:NCOL(setupCosts),function(i){
    ## costs for a single setup change
    costs_perChange<-setupCosts[,i] / variety

    ## distribute costs based on the product variety
    pr<-1
    sapply(1:NCOL(EAD[[1]]$P$PrD),function(pr){
      vec<-EAD[[1]]$P$PrD[includedProd,pr]
      temp<-dplyr::as_data_frame(table(vec)) %>%
        mutate(costProd=costs_perChange[pr]/n)
        return(temp$costProd[match(vec,temp$vec)]/EAD[[1]]$DEMAND)
      })

  })


  tapply(seq_along(vec), vec, identity)[unique(vec)]




}

