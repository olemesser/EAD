processVariety<-function(EAD,
                         includedProd=NULL,
                         prop_setupChange=c(0.05,0.1),
                         samples=10){

  if(is.null(includedProd)) includedProd<-1:NROW(EAD[[1]]$P$PrD)

  ## calculate process variants
  variety<-apply(EAD[[1]]$P$PrD[includedProd,],2,function(t){
    length(unique(t[t>0]))
  })

  ## calculate costs for each process
  RD<-diag(1,ncol = NCOL(EAD[[1]]$P$PrD),nrow = NCOL(EAD[[1]]$P$PrD)) %*% (( EAD[[1]]$DMM$PrD_RD %*% EAD[[1]]$DSM$RD) + EAD[[1]]$DMM$PrD_RD)
  costs<-clc_PCB(RD,DMD = 1,RC_var = EAD[[1]]$RC$var)
  ## derive setup costs
  setupCosts <- costs$PC_B %*% runif(samples,min = prop_setupChange[1], max = prop_setupChange[2])

  i<-1
  PC_setup<-sapply(1:NCOL(setupCosts),function(i){
    ## costs for a single setup change
    ## -1 is added since if there is only one process variant no changes are needed
    costs_perChange<-setupCosts[,i] / (variety-1)
    costs_perChange<-ifelse(is.infinite(costs_perChange),0,costs_perChange)
    ## distribute costs based on the product variety
    pr<-1
    df<-sapply(1:NCOL(EAD[[1]]$P$PrD),function(pr){
      vec<-EAD[[1]]$P$PrD[includedProd,pr]
      temp<-dplyr::as_tibble(table(vec)) %>%
        mutate(costspCh=costs_perChange[pr],
               costspCh=ifelse(vec==0,0,costspCh))
      temp$costspCh[match(vec,as.numeric(temp$vec))]
      # %>%
      #   left_join(tibble(vec=as.character(EAD[[1]]$P$PrD[includedProd,pr]),
      #                    DMD=EAD[[1]]$DEMAND[includedProd]) %>%
      #               group_by(vec) %>%
      #               summarise(DMD=sum(DMD)),
      #             by="vec") %>%
      #   mutate(costPROD=costspCh/DMD)
        return(temp$costPROD[match(vec,as.numeric(temp$vec))])
      })
    return(rowSums(df))
  })
  return(PC_setup)
}

