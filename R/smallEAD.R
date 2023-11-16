crt_smallEAD <- function(){
  require(tidyr)
  set.seed(1234)
  DOE<-expand_grid(N_FR = 4, # number of functional requirements
                   N_DD = 4, # number of physical domain elements
                   N_PrD = 8, # number of process domain elements
                   N_RD = 10, # number of resource domain elements
                   prop_PROD = 1, # density of the functional product mix
                   method_FD = "random", # method to create the product mix
                   N_PROD = list(c(20,25)), # number of products
                   TOTAL_DEMAND = list(c(100,10000)), # total demand
                   Q_VAR = list(c(0,1.7)), # demand heterogeneity
                   DMM_PAR = expand_grid(FD_PD=list(c(0,0.1)),
                                         PD_PrD=list(c(0,0.05)),
                                         PrD_RD=list(c(0,0.05))), # desired system design complexity
                   uB_DMM = 10,
                   allowZero = T,
                   ut_DMM = F, # if the upper triangle DMMs should be generated too (DMM_FD_PrD,DMM_FD_RD,DMM_PD_RD)
                   DSM_param=expand_grid(PD=list(c(0,0.1,0,1)),
                                         PrD=list(c(0,0,0,1)),
                                         RD=list(c(0,0,0,1))), # first two entries refer to the density of the DSMs and the second pair to the cv if the DSM_method='modular' is used.
                   DSM_method='modular',
                   ub_DSM = 1,
                   TC = 10^6, # total costs
                   r_in = list(c(0,0.9)),
                   r_fix = list(c(0,1)), # proportion of fixed costs on total costs
                   RC_sdlog = list(c(0.1,2)), # skewness of resource costs
                   N_RUN = 1:1 # number of runs
  )
  smallEAD <- crt_EAD(DOE)[[1]]
}


crt_sampleEAD <- function(){
  require(digest)
  P <- list(FD = matrix(c(1,0,0,
                            0,1,0,
                            0,0,1,
                            0,1,1),byrow=T,nrow=4))
  DSM <- list(PD = matrix(c(0,0,0,1,
                            0,0,0,0,
                            0,0,0,0,
                            0,0,0,0),byrow=T,nrow=4),
              PrD = diag(0,4),
              RD = diag(0,4))

  DMM <- list(FD_PD = matrix(c(1,0,0,0,
                               0,1,1,0,
                               0,0,1,1),byrow=T,nrow=3),
              PD_PrD = matrix(c(1,0,0,0,
                                0,1,3,0,
                                0,2,0,0,
                                0,0,0,1),byrow=T,nrow=4),
              PrD_RD = matrix(c(1,0,0,0,
                                0,3,5,0,
                                0,0,1,0,
                                0,0,0,2),byrow=T,nrow=4),
              FD_PrD = matrix(0,nrow=3,ncol=4),
              FD_RD = matrix(0,nrow=3,ncol=4),
              PD_RD = matrix(0,nrow=4,ncol=4))


  P_DD_1 <-P$FD %*% DMM$FD_PD
  P_DD_2 <- P_DD_1 %*% DSM$PD
  ### second multiply with DSM
  P[["PD"]] <- pmax(P_DD_1,P_DD_2)
  remove(P_DD_1,P_DD_2)

  ### 2.8.2 Calculate Process Domain ###
  P[["PrD"]] <- P[["PD"]] %*% ((DMM$PD_PrD %*% DSM$PrD) + DMM$PD_PrD)

  ### 2.8.3 Calculate Resource Domain ###
  P[["RD"]] <- P[["PrD"]] %*% (( DMM$PrD_RD %*% DSM$RD) + DMM$PrD_RD)

  EAD<-list(P=P,
            DSM=DSM,
            DMM=DMM,
            DEMAND=c(10,5,2,4),
            RC=list(var_d = c(100,800,600,500),
                    fix_d = rep(700,4),
                    var_i = rep(0,4),
                    fix_i = rep(0,4)),
            DOE=NA,
            message = "sucess")
  EAD$measures <- update_EAD(EAD)$measures

  EAD$ID<-digest(EAD, algo="crc32", serialize =TRUE)
  return(EAD)
}
