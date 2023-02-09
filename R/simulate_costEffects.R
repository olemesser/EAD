simulate_costEffects<-function(DOE){
  require(EAD)
  suppressWarnings(require(digest))
  suppressWarnings(require(tidyr))

  #### Input for Testing ####
  # DOE<-expand_grid(N_FR = list(c(7)), # number of functional requirements
  #                  N_DD = list(c(20)), # number of physical domain elements
  #                  N_PrD = list(c(40)), # number of process domain elements
  #                  N_RD = list(c(80)), # number of resource domain elements
  #                  PARAM_FD =1, # generate free combination
  #                  method_FD = "random",
  #                  TOTAL_DEMAND = 10000, # total demand
  #                  DMD_cv = list(c(0,3)), # coefficient of variation for demand distribution
  #                  DMM_PAR = list(c(0,0.115)), # desired design complexity
  #                  DMM_method="SDC", # method for generating the DMM
  #                  ut_DMM = F, # if the upper triangle DMMs should be generated too (DMM_FD_PrD,DMM_FD_RD,DMM_PD_RD)
  #                  DSM_param=list(c(0,0.14,0,1)), # first two entries refer to the density of the DSMs and the second pair to the cv if the DSM_method='modular' is used.
  #                  DSM_method='modular',
  #                  TC = 10^6, # total costs
  #                  ratio_fixedC = list(c(0.5,0.8)), # proportion of fixed costs on total costs
  #                  RC_cor = list(c(0,1)), # correlation between variable cost vector and fixed cost vector
  #                  RC_cv = list(c(0,1.5)), # coefficient of variation for resource cost distribution
  #                  N_RUN = 1:2 # number of runs
  # )
  # x<-1
  #### END Input for Testing ####

  #### For each DOE ####
  res<-lapply(1:NROW(DOE),function(x){
    EAD<-crt_EAD(DOE[x,])

    ## order product quantity in decreasing order
    order_prod_quantity<-order(EAD[[1]]$DEMAND,decreasing = T)

    costs<- clc_PCB(RES_CONS_PAT = EAD[[1]]$P$RD,
                    DMD =EAD[[1]]$DEMAND,
                    RC_var = EAD[[1]]$RC$var,
                    RC_fix = EAD[[1]]$RC$fix)
    RCU_initial<-costs$RCU
    #### 2. For each new product ###
    qty<-10
    res<-lapply(10:length(order_prod_quantity),function(qty){
      ### 2.1 Select the qty largest products ###
      idx<-order_prod_quantity[1:qty]

      ### 2.2 calculate Product costs and total costs ###
      DMD_temp<-EAD[[1]]$DEMAND
      DMD_temp[setdiff(1:length(DMD_temp),idx)]<-0
      costs<-clc_PCB(RES_CONS_PAT = EAD[[1]]$P$RD,
              DMD = DMD_temp,
              RCU = RCU_initial,
              RC_fix = EAD[[1]]$RC$fix)
      TC<-sum(costs$PC_B*DMD_temp)

      #### complexity costs due to additional process variety
      PrC<-clc_processCostRate(EAD,
                          RCU = RCU_initial)
      PC_setupChange<-processVariety(PrD = EAD[[1]]$P$PrD,
                                     PrC = PrC,
                                     DMD = EAD[[1]]$DEMAND,
                                     includedProd=idx,
                                     prop_setupChange=c(0.5,0.5),
                                     prop_setupTime=1)
      TC_setup<-sum(PC_setupChange * DMD_temp)

      out<-list(NPV=qty,
                DMD_total = sum(DMD_temp),
                TC = list(TC=TC,
                          TC_setup=TC_setup),
                # PC_B = data.frame(PC_B= costs$PC_B,
                #             PC_B_setup = PC_setupChange),
                PCI_PD = measure_PCI(EAD[[1]]$P$PD[idx,]),
                DMD_T10 = measure_TOP10(EAD[[1]]$DEMAND[idx]),
                SDC_n_FD_PD = EAD[[1]]$measures$SYSTEM$SDC_n$FD_PD,
                SDC_n_PD_PrD = EAD[[1]]$measures$SYSTEM$SDC_n$PD_PrD)
      return(out)
    })
    df<-lapply(res,function(t){
      data.frame(NPV=t$NPV,
                 DMD_total = t$DMD_total,
                 PCI_PD = t$PCI_PD,
                 DMD_T10 = t$DMD_T10,
                 SDC_n_FD_PD = t$SDC_n_FD_PD,
                 SDC_n_PD_PrD = t$SDC_n_PD_PrD,
                 TC=t$TC$TC,
                 TC_setup=t$TC$TC_setup)
    })
    df<-data.table::rbindlist(df)
    return(df)
  })
return(res)
}



simulate_costEffects_MC<-function(DOE,
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
                                 FUN = function(DOE, ...) {
                                   res<-list()
                                   res<-with_timeout(simulate_costEffects(DOE[1,]),timeout = time_limit)
                                   if(res$message=="error"){
                                     res<-list()
                                   }else if(res$message=="success"){
                                     res<-res$res
                                   }
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
    res<-par_apply(cl,X=DOE_list,FUN=simulate_costEffects)
  }
  res<-res[sapply(res,length)>0]
  return(res)
}
