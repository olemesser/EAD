indirectCost_allocation<-function(RES_CONS_PAT,
                                  DMD,
                                  RC_var,
                                  RC_fix,
                                  CSD_method=c("VBA","GRP","BP","IDX"),
                                  NUMB_CC){

  #### INput for Testing ####
  # data("csd_EAD")
  # RES_CONS_PAT <- CSD_EAD$RES_CONS_PAT
  # DMD <- CSD_EAD$DMD
  # RC_var <- CSD_EAD$RC_var_nc
  # RC_fix <- CSD_EAD$RC_fix_nc
  # remove(CSD_EAD)
  # NUMB_CC=1
  # CSD_method<-"VBA"
  #### END Input Testing ####


  #### 1. Calculate TRC and RCU
    DMD_mat <- matrix(DMD,
                      ncol = NCOL(RES_CONS_PAT),
                      byrow = F,
                      nrow = NROW(RES_CONS_PAT))
    TRC<-colSums(RES_CONS_PAT * DMD_mat)
    RCU_var <- RC_var/TRC
    RCU_fix <- RC_fix/TRC

  #### 2. Calculate Benchmark Costs ####
    PC_B_var <- RES_CONS_PAT %*% RCU_var
    PC_B_fixed <-RES_CONS_PAT %*% RCU_fix
    PC_B = PC_B_var + PC_B_fixed

  #### 3. Calculate Heuristic Costs ####
    PC_H_var <- RES_CONS_PAT %*% RCU_var # assumed to be error free for this experiment
    PC_H_fixed <-fixedCost_allocation(RES_CONS_PAT,
                                      RC_fix,
                                      DMD,
                                      method=CSD_method,
                                      NUMB_CC=NUMB_CC)
    PC_H=PC_H_var + PC_H_fixed


  #### 4. Calculate Measures ####
    errors<-allocation_Errors(PC_B = PC_B,
                      PC_H = PC_H,
                      DMD = DMD)
    sum(errors$ME_we)

  #### 5. Create output object
  measures<-data.frame(EUCD_uw=errors$EUCD_uw,
                 EUCD_we = errors$EUCD_we,
                 MPE_uw=mean(abs(errors$MPE_uw)),
                 MPE_we = mean(abs(errors$MPE_we)),
                 ME_uw= mean(abs(errors$ME_uw)),
                 ME_we = mean(abs(errors$ME_we)),
                 UC_p_uw = errors$UC_p_uw,
                 OC_p_uw = errors$OC_p_uw,
                 UC_p_we = errors$UC_p_we,
                 OC_p_we = errors$OC_p_we,
                 UC_top10_uw = errors$UC_top10_uw,
                 OC_top10_uw = errors$OC_top10_uw,
                 UC_top10_we = errors$UC_top10_we,
                 OC_top10_we = errors$OC_top10_we,
                 UC_top10_uw_abs = errors$UC_top10_uw_abs,
                 OC_top10_uw_abs = errors$OC_top10_uw_abs,
                 UC_top10_we_abs = errors$UC_top10_we_abs,
                 OC_top10_we_abs = errors$OC_top10_we_abs,
                 MPE_range_uw = errors$MPE_range_uw,
                 MPE_range_we = errors$MPE_range_we)
    out<-list(PC_B=PC_B,
              PC_H=PC_H,
              measures = measures)
  return(out)

}


allocation_Errors<-function(PC_B,PC_H,DMD=NULL){
  #### Testing ####
  # PC_B<-sample(1:100,50)
  # PC_H<-sample(1:100,50)
  # DMD<-NULL
  #### End Testing ###
  if(is.null(DMD)) DMD<-rep(1,length(PC_B))
  out<-list()
  PC_H_tilde<-PC_H[rep(1:length(PC_H),DMD)]
  PC_B_tilde<-PC_B[rep(1:length(PC_B),DMD)]


  #### Total Error ####
  out[['EUCD_uw']]<-sqrt(sum((PC_B-PC_H)^2)) # total euclidean error
  out[['MPE_uw']] <- as.numeric((PC_B-PC_H)/PC_B) # mean percentage error
  out[['ME_uw']] <- as.numeric((PC_B-PC_H)) # mean error
  out[['EUCD_we']]<-sqrt(sum((PC_B_tilde-PC_H_tilde)^2)) # total euclidean error
  out[['MPE_we']] <- (PC_B_tilde-PC_H_tilde)/PC_B_tilde # mean percentage error
  out[['ME_we']] <- (PC_B_tilde-PC_H_tilde) # mean error


  #### Distribution of Errors ####
  out[['UC_p_uw']]<-1/length(PC_B)*sum(PC_H>PC_B*1.05) # proportion of products being under-costed
  out[['OC_p_uw']]<-1/length(PC_B)*sum(PC_H<PC_B*0.95) # proportion of products being over-costed
  out[['UC_p_we']]<-1/length(PC_B_tilde)*sum(PC_H_tilde>PC_B_tilde*1.05) # proportion of products being under-costed
  out[['OC_p_we']]<-1/length(PC_B_tilde)*sum(PC_H_tilde<PC_B_tilde*0.95) # proportion of products being over-costed

  ### MPE ###
  MPE_UC_uw<-out[['MPE_uw']][out[['MPE_uw']]>0]
  MPE_OC_uw<-out[['MPE_uw']][out[['MPE_uw']]<0]
  MPE_UC_we<-out[['MPE_we']][out[['MPE_we']]>0]
  MPE_OC_we<-out[['MPE_we']][out[['MPE_we']]<0]

  ### ME ###
  ME_UC_uw<-out[['ME_uw']][out[['ME_uw']]>0]
  ME_OC_uw<-out[['ME_uw']][out[['ME_uw']]<0]
  ME_UC_we<-out[['ME_we']][out[['ME_we']]>0]
  ME_OC_we<-out[['ME_we']][out[['ME_we']]<0]


  ## relative values ##
  out[['UC_top10_uw']]<-mean(sort(MPE_UC_uw,decreasing = T)[1:ceiling(length(MPE_UC_uw)*0.1)]) # mean of top 10%  largest under-costed products
  out[['OC_top10_uw']]<-mean(sort(MPE_OC_uw,decreasing = F)[1:ceiling(length(MPE_OC_uw)*0.1)]) # mean of top 10%  largest over-costed products
  out[['UC_top10_we']]<-mean(sort(MPE_UC_we,decreasing = T)[1:ceiling(length(MPE_UC_we)*0.1)]) # mean of top 10%  largest under-costed products
  out[['OC_top10_we']]<-mean(sort(MPE_OC_we,decreasing = F)[1:ceiling(length(MPE_OC_we)*0.1)]) # mean of top 10%  largest over-costed products

  ## absolute values ##
  out[['UC_top10_uw_abs']]<-mean(sort(ME_UC_uw,decreasing = T)[1:ceiling(length(ME_UC_uw)*0.1)]) # mean of top 10%  largest under-costed products
  out[['OC_top10_uw_abs']]<-mean(sort(ME_OC_uw,decreasing = F)[1:ceiling(length(ME_OC_uw)*0.1)]) # mean of top 10%  largest over-costed products
  out[['UC_top10_we_abs']]<-mean(sort(ME_UC_we,decreasing = T)[1:ceiling(length(ME_UC_we)*0.1)]) # mean of top 10%  largest under-costed products
  out[['OC_top10_we_abs']]<-mean(sort(ME_OC_we,decreasing = F)[1:ceiling(length(ME_OC_we)*0.1)]) # mean of top 10%  largest over-costed products


  ## range of errors ##
  out[['MPE_range_uw']]<-abs(out[['UC_top10_uw']]-out[['OC_top10_uw']])
  out[['MPE_range_we']]<-abs(out[['UC_top10_we']]-out[['OC_top10_we']])

  return(out)

}


EXP_indirectAllocation<-function(x,NUMB_CC,CSD_method=c("VBA","GRP","BP","IDX")){
  require(dplyr)
  ### define costing system
  CSD<-tidyr::expand_grid(CSD_method=CSD_method,
                   NUMB_CC=NUMB_CC) %>%
    filter(!(CSD_method=="VBA" & NUMB_CC>first(NUMB_CC)))

  out<-lapply(1:NROW(CSD),function(cs){
    ### Apply Volume Based Costing ###
    costs<-indirectCost_allocation(RES_CONS_PAT = x[[1]]$P$RD,
                                   DMD = x[[1]]$DEMAND,
                                   RC_var = x[[1]]$RC$var,
                                   RC_fix = x[[1]]$RC$fix,
                                   CSD_method=CSD$CSD_method[cs],
                                   NUMB_CC = CSD$NUMB_CC[cs])
    colnames(costs$measures) <- paste0("cm.",colnames(costs$measures))

    data<-tibble(ID=x[[1]]$ID,
                 as.data.frame(x[[1]]$measures),
                 CSD_method=CSD$CSD_method[cs],
                 NUMB_CC = CSD$NUMB_CC[cs],
                 costs$measures,
                 DENS_RCP = sum(x[[1]]$P$RD>0)/prod(dim(x[[1]]$P$RD)))
    data$DISS_norm<-measure_DISS(x[[1]]$P$FD)


    MPE = (costs$PC_B-costs$PC_H)/costs$PC_B

    pc_raw<-tibble(ID=x[[1]]$ID,
                   PC_B = as.numeric(costs$PC_B),
                   PC_H = as.numeric(costs$PC_H),
                   MPE = as.numeric(MPE),
                   MPE_abs = as.numeric(abs(MPE)),
                   DMD = x[[1]]$DEMAND,
                   DMD_perc = x[[1]]$DEMAND/sum(x[[1]]$DEMAND),
                   E = measure_EXOTICNESS(x[[1]]$P$RD),
                   MEAN_DIST_ww =NA,
                   NN_dist = measure_NN(x[[1]]$P$RD,prop=0.2),
                   INTRA = measure_INTRA(x[[1]]$P$RD),
                   INTRA_norm = measure_INTRA(x[[1]]$P$RD,norm = T),
                   INTER = as.numeric(measure_INTER(x[[1]]$P$RD)),
                   CSD_method=CSD$CSD_method[cs],
                   NUMB_CC = CSD$NUMB_CC[cs],
                   RC_fix = x[[1]]$measures$RC$ratio_fixedC)
    out<-list(data = data,
              pc_raw = pc_raw)
  })

  out<-list(data = data.table::rbindlist(lapply(out,function(x) x$data)),
            pc_raw = data.table::rbindlist(lapply(out,function(x) x$pc_raw)))
  return(out)

}
