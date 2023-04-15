simulate_costEffects<-function(DOE){
  require(EAD)
  suppressWarnings(require(digest))
  suppressWarnings(require(tidyr))

  #### Input for Testing ####
  set.seed(1243)
  DOE<-expand_grid(N_FR = list(c(7)), # number of functional requirements
                   N_DD = list(c(15,30)), # number of physical domain elements
                   N_PrD = list(c(30,60)), # number of process domain elements
                   N_RD = list(c(60)), # number of resource domain elements
                   prop_PROD  = 1,
                   method_FD = "random",
                   TOTAL_DEMAND = 10000, # total demand
                   Q_VAR = list(c(0,2)), # coefficient of variation for demand distribution
                   DMM_PAR = expand_grid(FD_PD=list(c(0,0.1)),
                                         PD_PrD=list(c(0,0.0)),
                                         PrD_RD=list(c(0,0.0))), # desired system design complexity
                   uB_DMM = 5,
                   allowZero = F,
                   ut_DMM = F, # if the upper triangle DMMs should be generated too (DMM_FD_PrD,DMM_FD_RD,DMM_PD_RD)
                   DSM_param=expand_grid(PD=list(c(0,0.1,0,1)),
                                         PrD=list(c(0,0,0,1)),
                                         RD=list(c(0,0,0,1))), # first two entries refer to the density of the DSMs and the second pair to the cv if the DSM_method='modular' is used.
                   DSM_method='modular',
                   ub_DSM = 5,
                   TC = 10^6, # total costs
                   r_in = list(c(0.1,0.3)),
                   r_fix = list(c(1,1)), # proportion of fixed costs on total indirect costs
                   RC_cv = list(c(0,3)), # coefficient of variation for resource cost distribution
                   R_dvl = list(c(0.1,0.6)),
                   R_PA = list(c(0.1,0.6)),
                   C_order = list(c(50,100)), # order costs
                   C_hold = list(c(500,1000)), # holding costs
                   C_setup = list(c(10^6*0.05,10^6*0.5)), # 5%-40%
                   C_supply = list(c(100,1000)),
                   N_RUN = 1:1 # number of runs
  )
  x<-1
  #### END Input for Testing ####

  #### For each DOE ####
  res<-lapply(1:NROW(DOE),function(x){
    #### 1. Initialize ####
      ### 1.1 Generate the EAD ###
      EAD<-crt_EAD(DOE[x,])[[1]]



    #### 2. Create Scenarios for Overdesign ####
      ### for each overdesign scenario the product variety is increased starting with one product until all products are introduced
      DMD_sort <- sort(EAD$DEMAND,decreasing = T,index.return=TRUE)
      while(sum(colSums(EAD$DMM$FD_PD)>0)>1){
        ### 2.2 Calculate Initial Costs ###
          costs <- clc_domainCosts(EAD)

        #### 3. Product Variety ####
        p<-10
        R_dvl <- runif(1,min = DOE$R_dvl[x][[1]][1],max = DOE$R_dvl[x][[1]][2])
        R_PA <- runif(1,min = DOE$R_PA[x][[1]][1],max = DOE$R_PA[x][[1]][2])
        C_hold <- runif(1,DOE$C_hold[x][[1]][1],DOE$C_hold[x][[1]][2])
        C_order <- runif(1,DOE$C_order[x][[1]][1],DOE$C_order[x][[1]][2])
        C_setup <- runif(1,DOE$C_setup[x][[1]][1],DOE$C_setup[x][[1]][2])
        C_supply <- runif(1,DOE$C_supply[x][[1]][1],DOE$C_supply[x][[1]][2])


          lapply(1:length(DMD_sort$x),function(p){
                    out<-list()
                    ### 3.1 Exclude Products ###
                    products_in <- DMD_sort$ix[1:p]
                    products_out <- setdiff(DMD_sort$ix,products_in)
                    DEMAND_temp <- rep(0,length(EAD$DEMAND))
                    DEMAND_temp[products_in] <- EAD$DEMAND[products_in]

                    ### 3.2 Calculate Cost by going through the individual cost driver ###
                    ccDriver_develPA <- development_costs(P_PD = EAD$P$PD,
                                                          DEMAND = DEMAND_temp,
                                                          PDC_fix = costs$fix$PD,
                                                          R_dvl = R_dvl,
                                                          R_PA = R_PA)
                    out['dvl_developmentCosts'] <- ccDriver_develPA$TC_dvl
                    out['dvl_partadminCosts'] <-  ccDriver_develPA$TC_PA

                    ### 3.3 Calculate the Setup Costs ###
                    setup<- setupCosts(EAD,
                                   C_hold = C_hold,
                                   C_order = C_order,
                                   C_setup = C_setup,
                                   PrCU_var = costs$var$PrD)
                    out['man_setupCosts'] <- setup$TC_setup

                    ### 3.4. Tooling Costs ###


                    ### 3.5 Purchasing Orders ###
                    out['pur_orderCosts'] <- orderCosts(DMD_PD = setup$DMD_component,
                                                    N_lot = setup$lotSize,
                                                    C_order = C_order)

                    ### 3.6 Supplier Management Costs ###
                    out['pur_supplyCosts'] <- supplyCosts(P_PD = EAD$P$PD,
                                                C_supply = C_supply)



                    return(out)
          })



        EAD<-overdesign_EAD(EAD)

      }







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

      #### complexity costs due to smaller lot sizes
      PC_setupChange<- clc_setupCosts(EAD,
                           RCU_var = RCU_initial,
                           includedProd = idx,
                           prop_setupChange=c(0.5,0.5),
                           prop_setupTime=1,
                           r_Order_to_Hold=2)


      TC_setup<-sum(PC_setupChange$PC_setup * DMD_temp)

      out<-list(NPV=qty,
                DMD_total = sum(DMD_temp),
                TC = list(TC=TC,
                          TC_setup=TC_setup),
                # PC_B = data.frame(PC_B= costs$PC_B,
                #             PC_B_setup = PC_setupChange),
                PCI_PD = measure_PCI(EAD[[1]]$P$PD[idx,]),
                PCI_PD_w = measure_PCI(EAD[[1]]$P$PD[idx,],DMD=DMD_temp[idx]),
                DMD_T10 = measure_TOP10(EAD[[1]]$DEMAND[idx]),
                SDC_n_FD_PD = EAD[[1]]$measures$SYSTEM$SDC_n$FD_PD,
                SDC_n_PD_PrD = EAD[[1]]$measures$SYSTEM$SDC_n$PD_PrD,
                mean_lotsize = PC_setupChange$lotsize,
                mean_setups = sum(PC_setupChange$setups))
      out<-as.data.frame(out)
      return(out)
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
