update_EAD<-function(EAD){

  ### Calculate Physical Domain ###
  P_DD_1 <-EAD$P$FD %*% EAD$DMM$FD_PD
  P_DD_2 <- P_DD_1 %*% EAD$DSM$PD
  ### second multiply with DSM
  EAD$P$PD <- pmax(P_DD_1,P_DD_2)
  remove(P_DD_1,P_DD_2)

  ### Calculate Process Domain ###
  if(length(EAD$DMM$FD_PrD)==0) EAD$DMM$FD_PrD <- matrix(0,nrow=NCOL(EAD$P$FD),ncol=NCOL(EAD$P$PrD))
  EAD$P$PrD <- EAD$P$PD %*% (( EAD$DMM$PD_PrD %*% EAD$DSM$PrD) + EAD$DMM$PD_PrD) + EAD$P$FD  %*% EAD$DMM$FD_PrD

  ### Calculate Resource Domain ###
  if(length(EAD$DMM$FD_RD)==0) EAD$DMM$FD_RD <- matrix(0,nrow=NCOL(EAD$P$FD),ncol=NCOL(EAD$P$RD))
  if(length( EAD$DMM$PD_RD)==0)  EAD$DMM$PD_RD <- matrix(0,nrow=NCOL(EAD$P$PD),ncol=NCOL(EAD$P$RD))
  EAD$P$RD <- EAD$P$PrD %*% ((  EAD$DMM$PrD_RD %*%  EAD$DSM$RD) +  EAD$DMM$PrD_RD) + EAD$P$FD %*%  EAD$DMM$FD_RD + EAD$P$PD  %*%  EAD$DMM$PD_RD

  #### Update Measures ####
  DSM_PD_keep <- !(colSums(EAD$DSM$PD)==0 & rowSums(EAD$DSM$PD)==0 & colSums(EAD$DMM$FD_PD)==0)
  DSM_PrD_keep <- !(colSums(EAD$DSM$PrD)==0 & rowSums(EAD$DSM$PrD)==0 & colSums(EAD$DMM$PD_PrD)==0)
  DSM_RD_keep <- !(colSums(EAD$DSM$RD)==0 & rowSums(EAD$DSM$RD)==0 & colSums(EAD$DMM$PrD_RD)==0)


  EAD$measures <- list(SYSTEM = list(DMD_cv = sd(EAD$DEMAND[EAD$DEMAND>0])/mean(EAD$DEMAND[EAD$DEMAND>0]),
                                 DMD_T10 = measure_TOP10(EAD$DEMAND[EAD$DEMAND>0]),
                                 Q_VAR = sd(log(EAD$DEMAND[EAD$DEMAND>0])),
                                 SDC = list(FD_PD=measure_designComplexity(EAD$DMM$FD_PD,norm = F),
                                            FD_PrD=NA,
                                            FD_RD=NA,
                                            PD_PrD=measure_designComplexity(EAD$DMM$PD_PrD,norm = F),
                                            PD_RD=NA,
                                            PrD_RD=measure_designComplexity(EAD$DMM$PrD_RD,norm = F)),
                                 SDC_n = list(FD_PD=measure_designComplexity(EAD$DMM$FD_PD[rowSums(EAD$DMM$FD_PD)>0,colSums(EAD$DMM$FD_PD)>0,drop=F]),
                                              # FD_PrD=NA,
                                              # FD_RD=NA,
                                              PD_PrD=measure_designComplexity(EAD$DMM$PD_PrD[rowSums(EAD$DMM$PD_PrD)>0,colSums(EAD$DMM$PD_PrD)>0,drop=F]),
                                              # PD_RD=NA,
                                              PrD_RD=measure_designComplexity(EAD$DMM$PrD_RD[rowSums(EAD$DMM$PrD_RD)>0,colSums(EAD$DMM$PrD_RD)>0,drop=F])),
                                 R = list(FD_PD=NA, #measure_reangularity(EAD$DMM$FD_PD),
                                          PD_PrD=NA, #measure_reangularity(EAD$DMM$PD_PrD),
                                          PrD_RD=NA), #measure_reangularity(EAD$DMM$PrD_RD)),
                                 S = list(FD_PD=NA, #measure_semiangularity(EAD$DMM$FD_PD),
                                          FD_PrD=NA,
                                          FD_RD=NA,
                                          PD_PrD=NA, #measure_semiangularity(EAD$DMM$PD_PrD),
                                          PD_RD=NA,
                                          PrD_RD=NA), #measure_semiangularity(EAD$DMM$PrD_RD)),
                                 JSDC = list(FD_PD=measure_JSDC(EAD$DMM$FD_PD[rowSums(EAD$DMM$FD_PD)>0,colSums(EAD$DMM$FD_PD)>0,drop=F]),
                                             # FD_PrD=measure_JSDC(EAD$DMM$FD_PrD[rowSums(EAD$DMM$FD_PrD)>0,colSums(EAD$DMM$FD_PrD)>0,drop=F]),
                                             # FD_RD=measure_JSDC(EAD$DMM$FD_RD[rowSums(EAD$DMM$FD_RD)>0,colSums(EAD$DMM$FD_RD)>0,drop=F]),
                                             PD_PrD=measure_JSDC(EAD$DMM$PD_PrD[rowSums(EAD$DMM$PD_PrD)>0,colSums(EAD$DMM$PD_PrD)>0,drop=F]),
                                             # PD_RD=measure_JSDC(EAD$DMM$PD_RD[rowSums(EAD$DMM$PD_RD)>0,colSums(EAD$DMM$PD_RD)>0,drop=F]),
                                             PrD_RD=measure_JSDC(EAD$DMM$PrD_RD[rowSums(EAD$DMM$PrD_RD)>0,colSums(EAD$DMM$PrD_RD)>0,drop=F])),
                                 SC = list(PD=measure_designComplexity(EAD$DSM$PD[DSM_PD_keep, DSM_PD_keep, drop=F],norm = F),
                                           PrD=measure_designComplexity(EAD$DSM$PrD[DSM_PrD_keep, DSM_PrD_keep, drop=F],norm = F),
                                           RD=measure_designComplexity(EAD$DSM$RD[DSM_RD_keep, DSM_RD_keep, drop=F],norm = F)),
                                 Q = list(PD=measure_modularity(EAD$DSM$PD[DSM_PD_keep, DSM_PD_keep, drop=F]),
                                          PrD=measure_modularity(EAD$DSM$PrD[DSM_PrD_keep, DSM_PrD_keep, drop=F]),
                                          RD=measure_modularity(EAD$DSM$RD[DSM_RD_keep, DSM_RD_keep, drop=F])),
                                 NE = list(PD=measure_neumannEntropy(EAD$DSM$PD[DSM_PD_keep, DSM_PD_keep, drop=F]),
                                           PrD=measure_neumannEntropy(EAD$DSM$PrD[DSM_PrD_keep, DSM_PrD_keep, drop=F]),
                                           RD=measure_neumannEntropy(EAD$DSM$RD[DSM_RD_keep, DSM_RD_keep, drop=F])),
                                 DENS_DSM = list(PD=measure_DENS(EAD$DSM$PD[DSM_PD_keep, DSM_PD_keep, drop=F]),
                                                 PrD=measure_DENS(EAD$DSM$PrD[DSM_PrD_keep, DSM_PrD_keep, drop=F]),
                                                 RD=measure_DENS(EAD$DSM$RD[DSM_RD_keep, DSM_RD_keep, drop=F])),
                                 MCC = list(PD=measure_MCC(EAD$DSM$PD[DSM_PD_keep, DSM_PD_keep, drop=F],norm=F),
                                            PrD=measure_MCC(EAD$DSM$PrD[DSM_PrD_keep, DSM_PrD_keep, drop=F],norm=F),
                                            RD=measure_MCC(EAD$DSM$RD[DSM_RD_keep, DSM_RD_keep, drop=F],norm=F)),
                                 MCC_n = list(PD=measure_MCC(EAD$DSM$PD[DSM_PD_keep, DSM_PD_keep, drop=F]),
                                              PrD=measure_MCC(EAD$DSM$PrD[DSM_PrD_keep, DSM_PrD_keep, drop=F]),
                                              RD=measure_MCC(EAD$DSM$RD[DSM_RD_keep, DSM_RD_keep, drop=F])),
                                 HVM = list(PD=measure_HVM(EAD$DSM$PD[DSM_PD_keep, DSM_PD_keep, drop=F],norm=F),
                                            PrD=measure_HVM(EAD$DSM$PrD[DSM_PrD_keep, DSM_PrD_keep, drop=F],norm=F),
                                            RD=measure_HVM(EAD$DSM$RD[DSM_RD_keep, DSM_RD_keep, drop=F],norm=F)),
                                 HVM_n = list(PD=measure_HVM(EAD$DSM$PD[DSM_PD_keep, DSM_PD_keep, drop=F]),
                                              PrD=measure_HVM(EAD$DSM$PrD[DSM_PrD_keep, DSM_PrD_keep, drop=F]),
                                              RD=measure_HVM(EAD$DSM$RD[DSM_RD_keep, DSM_RD_keep, drop=F])),
                                 HIC = list(PD=measure_HIC(EAD$DSM$PD[DSM_PD_keep, DSM_PD_keep, drop=F],norm=F),
                                            PrD=measure_HIC(EAD$DSM$PrD[DSM_PrD_keep, DSM_PrD_keep, drop=F],norm=F),
                                            RD=measure_HIC(EAD$DSM$RD[DSM_RD_keep, DSM_RD_keep, drop=F],norm=F)),
                                 HIC_n = list(PD=measure_HIC(EAD$DSM$PD[DSM_PD_keep, DSM_PD_keep, drop=F]),
                                              PrD=measure_HIC(EAD$DSM$PrD[DSM_PrD_keep, DSM_PrD_keep, drop=F]),
                                              RD=measure_HIC(EAD$DSM$RD[DSM_RD_keep, DSM_RD_keep, drop=F])),
                                 DNS = list(FD = measure_DENS(EAD$P$FD[EAD$DEMAND>0,,drop=F]),
                                            PD = measure_DENS(EAD$P$PD[EAD$DEMAND>0,,drop=F]),
                                            PrD = measure_DENS(EAD$P$PrD[EAD$DEMAND>0,,drop=F]),
                                            RD = measure_DENS(EAD$P$RD[EAD$DEMAND>0,,drop=F])),
                                 PCI = list(FD = measure_PCI(EAD$P$FD[EAD$DEMAND>0,,drop=F]),
                                            PD = measure_PCI(EAD$P$PD[EAD$DEMAND>0,,drop=F]),
                                            PrD = measure_PCI(EAD$P$PrD[EAD$DEMAND>0,,drop=F]),
                                            RD = measure_PCI(EAD$P$RD[EAD$DEMAND>0,,drop=F])),
                                 D_u = list(FD = NA, #measure_diversificationINDEX(EAD$P$FD[EAD$DEMAND>0,,drop=F]),
                                            PD = NA, #measure_diversificationINDEX(EAD$P$FD[EAD$DEMAND>0,,drop=F]),
                                            PrD = NA, #measure_diversificationINDEX(EAD$P$PrD[EAD$DEMAND>0,,drop=F]),
                                            RD = measure_diversificationINDEX(EAD$P$RD[EAD$DEMAND>0,,drop=F])), #),
                                 NPV = list(FD = measure_NPV(EAD$P$FD[EAD$DEMAND>0,,drop=F])),
                                 CI = list(FD = measure_CI(EAD$P$FD[EAD$DEMAND>0,,drop=F]),
                                           PD = measure_CI(EAD$P$PD[EAD$DEMAND>0,,drop=F]),
                                           PrD = measure_CI(EAD$P$PrD[EAD$DEMAND>0,,drop=F]),
                                           RD = measure_CI(EAD$P$RD[EAD$DEMAND>0,,drop=F])),
                                 OV = list(FD = NA, #measure_OV(EAD$P$FD[EAD$DEMAND>0,,drop=F]),
                                           PD = NA, #measure_OV(EAD$P$PD[EAD$DEMAND>0,,drop=F]),
                                           PrD = NA, #measure_OV(EAD$P$PrD[EAD$DEMAND>0,,drop=F]),
                                           RD = NA), #measure_OV(EAD$P$RD[EAD$DEMAND>0,,drop=F])),
                                 # RES_COR = NA,
                                 D = list(FD = NA, #measure_diversificationINDEX(EAD$P$FD[EAD$DEMAND>0,,drop=F],DMD = EAD$DEMAND[EAD$DEMAND>0]),
                                          PD = NA, #measure_diversificationINDEX(EAD$P$PD[EAD$DEMAND>0,,drop=F],DMD = EAD$DEMAND[EAD$DEMAND>0]),
                                          PrD = NA, #measure_diversificationINDEX(EAD$P$PrD[EAD$DEMAND>0,,drop=F],DMD = EAD$DEMAND[EAD$DEMAND>0]),
                                          RD = measure_diversificationINDEX(EAD$P$RD[EAD$DEMAND>0,,drop=F],DMD = EAD$DEMAND[EAD$DEMAND>0])), #),
                                 N_FD = sum(colSums(EAD$P$FD[EAD$DEMAND>0,,drop=F])>0),
                                 N_PD = sum(colSums(EAD$P$PD[EAD$DEMAND>0,,drop=F])>0),
                                 N_PrD = sum(colSums(EAD$P$PrD[EAD$DEMAND>0,,drop=F])>0),
                                 N_RD = sum(colSums(EAD$P$RD[EAD$DEMAND>0,,drop=F])>0),
                                 N_PROD = sum(EAD$DEMAND[EAD$DEMAND>0]>0),
                                 TSS = sum(colSums(EAD$P$FD[EAD$DEMAND>0,,drop=F])>0) + sum(colSums(EAD$P$PD[EAD$DEMAND>0,,drop=F])>0) + sum(colSums(EAD$P$PrD[EAD$DEMAND>0,,drop=F])>0) + sum(colSums(EAD$P$RD[EAD$DEMAND>0,,drop=F])>0),
                                 RC = list(RC_vari_top10 = measure_TOP10(EAD$RC$var_i),
                                           RC_vard_top10 = measure_TOP10(EAD$RC$var_d),
                                           RC_fixi_top10 = measure_TOP10(EAD$RC$fix_i),
                                           RC_fixd_top10 = measure_TOP10(EAD$RC$fix_d),
                                           r_fix = sum(EAD$RC$fix_d,EAD$RC$fix_i) / sum(EAD$RC$var_d,EAD$RC$var_i,EAD$RC$fix_d,EAD$RC$fix_i),
                                           r_in = sum(EAD$RC$var_i,EAD$RC$fix_i) / sum(EAD$RC$var_d,EAD$RC$var_i,EAD$RC$fix_d,EAD$RC$fix_i),
                                           TC = sum(EAD$RC$var_d,EAD$RC$var_i,EAD$RC$fix_d,EAD$RC$fix_i))),
                   PRODUCT = list(INTER = list(FD = measure_INTER(EAD$P$FD[EAD$DEMAND>0,,drop=F]),
                                               PD = measure_INTER(EAD$P$PD[EAD$DEMAND>0,,drop=F]),
                                               PrD = measure_INTER(EAD$P$PrD[EAD$DEMAND>0,,drop=F]),
                                               RD = measure_INTER(EAD$P$RD[EAD$DEMAND>0,,drop=F])),
                                  INTRA = list(FD = measure_INTRA(EAD$P$FD[EAD$DEMAND>0,,drop=F]),
                                               PD = measure_INTRA(EAD$P$PD[EAD$DEMAND>0,,drop=F]),
                                               PrD = measure_INTRA(EAD$P$PrD[EAD$DEMAND>0,,drop=F]),
                                               RD = measure_INTRA(EAD$P$RD[EAD$DEMAND>0,,drop=F])),
                                  LOF = list(FD = NA, #measure_LOF(EAD$P$FD[EAD$DEMAND>0,,drop=F]),
                                             PD = NA, #measure_LOF(EAD$P$PD[EAD$DEMAND>0,,drop=F]),
                                             PrD = NA, #measure_LOF(EAD$P$PrD[EAD$DEMAND>0,,drop=F]),
                                             RD = NA))) #measure_LOF(EAD$P$RD[EAD$DEMAND>0,,drop=F]))))
  return(EAD)
}
