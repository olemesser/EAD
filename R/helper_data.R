#' @title Creates sample data for CSD
#' @description Creates a sample data object containing RES_CONS_PAT, RC_var_nc, RC_var_cc, RC_fix_nc, RC_fix_cc and DMD.
#' @return The calculated modularity measure
#' @noRd
create_RES_CONS_PAT<-function(n_RES=50,
                              n_PROD=70){

  TC<-10^6
  repeat{
    RES_CONS_PAT<-matrix(sample(0:50,n_RES*n_PROD,replace = T),nrow=n_PROD)
    if(all(colSums(RES_CONS_PAT)>0) & all(rowSums(RES_CONS_PAT)>0)) break
  }

  DMD<-crt_DEMAND(n_PROD,
             TOTAL_DEMAND = 1000,
             Q_VAR = runif(1,0,3))$DEMAND

  ratio <- runif(4)
  ratio <- ratio/sum(ratio)
  RC_var_nc<-genRC(n_RES = n_RES,
        DISP2 = 0.5,
        TC=ratio[1]*TC)$RC
  RC_var_cc<-genRC(n_RES = n_RES,
                   DISP2 = 0.5,
                   TC=ratio[2]*TC)$RC
  RC_fix_nc<-genRC(n_RES = n_RES,
                   DISP2 = 0.5,
                   TC=ratio[3]*TC)$RC
  RC_fix_cc<-genRC(n_RES = n_RES,
                   DISP2 = 0.5,
                   TC=ratio[4]*TC)$RC



  CSD_EAD<-list(RES_CONS_PAT=RES_CONS_PAT,
                DMD=DMD,
                RC_var_nc=RC_var_nc,
                RC_var_cc=RC_var_cc,
                RC_fix_nc=RC_fix_nc,
                RC_fix_cc=RC_fix_cc)


  path<-getwd()
  setwd(paste0(path,"/data"))
  save(CSD_EAD,file="csd_EAD.RData")
  setwd(path)


}
