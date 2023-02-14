sim_SetupCosts<-function(DOE){
  require(EAD)
  suppressWarnings(require(digest))
  suppressWarnings(require(tidyr))

  #### Input for Testing ####
  # DOE<-expand_grid(N_FR = list(c(13)), # number of functional requirements
  #                  N_DD = list(c(26)), # number of physical domain elements
  #                  N_PrD = list(c(52)), # number of process domain elements
  #                  N_RD = list(c(104)), # number of resource domain elements
  #                  PARAM_FD = seq(0.07,0.5,0.02), # generate free combination
  #                  method_FD = "DNS",
  #                  TOTAL_DEMAND = 10000, # total demand
  #                  DMD_cv = list(c(0,3)), # coefficient of variation for demand distribution
  #                  DMM_PAR = list(c(0,0.115)), # desired design complexity
  #                  DMM_method="SDC", # method for generating the DMM
  #                  ut_DMM = F, # if the upper triangle DMMs should be generated too (DMM_FD_PrD,DMM_FD_RD,DMM_PD_RD)
  #                  DSM_param=list(c(0,0.14,0,1)), # first two entries refer to the density of the DSMs and the second pair to the cv if the DSM_method='modular' is used.
  #                  DSM_method='modular',
  #                  TC = 10^6, # total costs
  #                  ratio_fixedC = list(c(0,0)), # proportion of fixed costs on total costs
  #                  RC_cor = list(c(0,1)), # correlation between variable cost vector and fixed cost vector
  #                  RC_cv = list(c(0,1.5)), # coefficient of variation for resource cost distribution
  #                  N_RUN = 1:2 # number of runs
  # )
  # x<-1
  #### END Input for Testing ####

  #### For each DOE ####
    DOE<-DOE %>%
      group_by_all() %>%
      group_split()

  res<-lapply(1:length(DOE),function(x){
      #### Create an EAD ###
        EAD<-crt_EAD(DOE[[x]])
        P_PD<-EAD[[1]]$P$PD
        P_PD[P_PD>0] <- 1

        #### Calculate Costs ###
        costs<- clc_PCB(RES_CONS_PAT = EAD[[1]]$P$RD,
                        DMD = EAD[[1]]$DEMAND,
                        RC_var = EAD[[1]]$RC$var,
                        RC_fix = EAD[[1]]$RC$fix)
        RCU_initial<-costs$RCU

      #### Manipulate Product Design ####
        res<-list()
        i<-2
        while (sum(colSums(P_PD)>0)>2) {
          if(i>1){
            PD_left <- which(colSums(P_PD)>0)
            PD_remove <- sample(PD_left,1)
            PD_replace <- sample(setdiff(PD_left,PD_remove),1)

            ## remove element from DMM_FD_PD
            DMM_repl <- EAD[[1]]$DMM$FD_PD[,PD_remove] + EAD[[1]]$DMM$FD_PD[,PD_replace]
            DMM_repl[DMM_repl>0] <- 1
            EAD[[1]]$DMM$FD_PD[,PD_replace] <- DMM_repl
            EAD[[1]]$DMM$FD_PD[,PD_remove] <- 0

            ## remove element from DSM_PD in rows
            DSM_row <- EAD[[1]]$DSM$PD[PD_remove,] + EAD[[1]]$DSM$PD[PD_replace,]
            DSM_row[DSM_row>0] <- 1
            EAD[[1]]$DSM$PD[PD_replace,] <- DSM_row
            EAD[[1]]$DSM$PD[PD_remove,] <- 0

            ## remove elemtns from column in DSM_PD
            DSM_col <- EAD[[1]]$DSM$PD[,PD_remove] + EAD[[1]]$DSM$PD[,PD_replace]
            DSM_col[DSM_col>0] <- 1
            EAD[[1]]$DSM$PD[,PD_replace] <- DSM_col
            EAD[[1]]$DSM$PD[,PD_remove] <- 0

            P_DD_1 <- EAD[[1]]$P$FD %*% EAD[[1]]$DMM$FD_PD
            P_DD_2 <- P_DD_1 %*% EAD[[1]]$DSM$PD
            ### second multiply with DSM
            P_PD <- pmax(P_DD_1,P_DD_2)
            P_PD[P_PD>0] <- 1
          }
          remove(P_DD_1)
          remove(P_DD_2)
          remove(DSM_col)
          remove(DSM_row)
          remove(DMM_repl)

          EAD_temp <- EAD
          EAD_temp[[1]]$P$PD <- P_PD
          PC_setupChange<- clc_setupCosts(EAD_temp,
                                          RCU_var = RCU_initial,
                                          includedProd = 1:NROW(P_PD),
                                          prop_setupChange=c(0.5,0.5),
                                          prop_setupTime=1,
                                          r_Order_to_Hold=2)

          res[[i]]<-list(id=EAD[[1]]$ID,
                         NPV=NROW(P_PD),
                          DMD_total = sum(EAD_temp[[1]]$DEMAND),
                          TC = list(TC=sum(costs$PC_B*EAD_temp[[1]]$DEMAND),
                                    TC_setup=sum(PC_setupChange$PC_setup * EAD_temp[[1]]$DEMAND)),
                          PCI_PD = measure_PCI(P_PD),
                          PCI_PD_w = measure_PCI(P_PD,DMD = EAD_temp[[1]]$DEMAND),
                          SDC_n_PD_PrD = EAD[[1]]$measures$SYSTEM$SDC_n$PD_PrD,
                          SDC_n_PrD_RD = EAD[[1]]$measures$SYSTEM$SDC_n$PrD_RD,
                          HIC_n_PrD = EAD[[1]]$measures$SYSTEM$HIC_n$PrD,
                          HIC_n_RD = EAD[[1]]$measures$SYSTEM$HIC_n$RD,
                          mean_lotsize = PC_setupChange$lotsize,
                          mean_setups = sum(na.omit(PC_setupChange$setups)))
          res[[i]]<-as.data.frame(res[[i]])
          i<-i+1
        }
        res<-data.table::rbindlist(res)
    return(res)
    })
  res<-data.table::rbindlist(res)
  return(res)
}




sim_SetupCosts_MC<-function(DOE=NULL,
                            NUMB_CORES=4,
                            cluster=F,
                            logfile="",
                            extMC_lib=F,
                            ehNodes="remove"){
  suppressWarnings(require(parallel))
  suppressWarnings(require(doSNOW))
  suppressWarnings(require(foreach))
  suppressWarnings(require(dplyr))

    DOE_list<-DOE %>%
      group_split(split=1:n())%>%
      as.list()

  if(extMC_lib){
    library(odegoparallel)
    cl <- odegoparallel::initMC(NUMB_CORES = NUMB_CORES,
                                cluster = cluster,
                                logfile = logfile)
    parallel::clusterExport(cl, envir = globalenv(),c("time_limit"))
    print(cl)
    res <- odegoparallel::run_MC(cl, X = DOE_list,
                                 FUN = function(ip, ...) {
                                   res<-list()
                                     res<-with_timeout(sim_SetupCosts(DOE=ip),timeout = time_limit)
                                   gc()
                                   return(res)
                                 }, packages = c("EAD","odegoparallel",
                                                 "dplyr", "rpicosat", "tidyr", "GA",
                                                 "Matrix", "digest", "faux","data.table",
                                                 "plyr", "DescTools", "igraph", "R.utils","fGarch"),
                                 errorhandlingNodes = ehNodes)
  }else{
    cl<-EAD::setupMC(NUMB_CORES=NUMB_CORES,logfile=logfile)
    registerDoSNOW(cl)
    on.exit(stopCluster(cl))
    print(cl)
    print(paste0("Time  Limit set to: ",time_limit, " seconds per EAD realization."))
    snow::clusterExport(cl,"time_limit")
    res<-par_apply(cl,X=DOE_list,FUN=sim_SetupCosts)
  }
  res<-res[sapply(res,length)>0]
  return(res)
}
