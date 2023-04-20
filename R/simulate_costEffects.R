simulate_costEffects<-function(DOE){
  suppressWarnings(require(EAD))
  suppressWarnings(require(digest))
  suppressWarnings(require(tidyr))
  suppressWarnings(require(data.table))
  suppressWarnings(require(dplyr))

  #### Input for Testing ####
  # set.seed(1243)
  # DOE<-expand_grid(N_FR = list(c(7)), # number of functional requirements
  #                  N_DD = list(c(10)), # number of physical domain elements
  #                  N_PrD = list(c(30)), # number of process domain elements
  #                  N_RD = list(c(60)), # number of resource domain elements
  #                  prop_PROD  = 1,
  #                  method_FD = "random",
  #                  TOTAL_DEMAND = 10000, # total demand
  #                  Q_VAR = list(c(0,2)), # coefficient of variation for demand distribution
  #                  DMM_PAR = expand_grid(FD_PD=list(c(0,0.1)),
  #                                        PD_PrD=list(c(0,0.1)),
  #                                        PrD_RD=list(c(0,0.1))), # desired system design complexity
  #                  uB_DMM = 5,
  #                  allowZero = F,
  #                  ut_DMM = F, # if the upper triangle DMMs should be generated too (DMM_FD_PrD,DMM_FD_RD,DMM_PD_RD)
  #                  DSM_param=expand_grid(PD=list(c(0,0.1,0,1)),
  #                                        PrD=list(c(0,0,0,1)),
  #                                        RD=list(c(0,0,0,1))), # first two entries refer to the density of the DSMs and the second pair to the cv if the DSM_method='modular' is used.
  #                  DSM_method='modular',
  #                  ub_DSM = 5,
  #                  TC = 10^6, # total costs
  #                  r_in = list(c(0.1,0.3)),
  #                  r_fix = list(c(1,1)), # proportion of fixed costs on total indirect costs
  #                  RC_cv = list(c(0,3)), # coefficient of variation for resource cost distribution
  #                  R_dvl = list(c(0.1,0.4)),
  #                  R_PA = list(c(0.1,0.3)),
  #                  R_order = list(c(0.01,0.1)), # total order costs
  #                  R_hold = list(c(0.1,0.2)), # total holding costs
  #                  R_setup = list(c(0.05,0.4)), # total setup costs 5%-40%
  #                  R_tooling = list(c(0.05,0.2)),
  #                  R_supply = list(c(0.05,0.2)),
  #                  N_RUN = 1:1 # number of runs
  # )
  # x<-1
  #### END Input for Testing ####

  #### For each DOE ####
  res<-lapply(1:NROW(DOE),function(x){
    #### 1. Initialize ####
      ### 1.1 Generate the EAD ###
      EAD<-crt_EAD(DOE[x,])[[1]]

      ### 1.2 Sample proportions based on input bounds ###
      R_dvl <- runif(1,min = DOE$R_dvl[x][[1]][1],max = DOE$R_dvl[x][[1]][2])
      R_PA <- runif(1,min = DOE$R_PA[x][[1]][1],max = DOE$R_PA[x][[1]][2])
      R_tooling <- runif(1,DOE$R_tooling[x][[1]][1],DOE$R_tooling[x][[1]][2])
      R_hold <- runif(1,DOE$R_hold[x][[1]][1],DOE$R_hold[x][[1]][2])
      R_order <- runif(1,DOE$R_order[x][[1]][1],DOE$R_order[x][[1]][2])
      R_setup <- runif(1,DOE$R_setup[x][[1]][1],DOE$R_setup[x][[1]][2])
      R_supply <- runif(1,DOE$R_supply[x][[1]][1],DOE$R_supply[x][[1]][2])

      ### 1.3 Calculate initial values for complexity cost drivers (CCD) ###
      costDriver_inital <- clc_initialCosts(EAD,
                                            R_dvl = R_dvl,
                                            R_PA = R_PA,
                                            R_tooling = R_tooling,
                                            R_hold = R_hold,
                                            R_order = R_order,
                                            R_setup = R_setup,
                                            R_supply = R_supply)

      ### 1.4 Calculate costs for non complexity driven costs ###
      nonCC <- clc_PCB(EAD$P$RD,
                       DMD=EAD$DEMAND,
                       RC_direct = EAD$RC$direct,
                       RC_var = EAD$RC$var,
                       RC_fix = EAD$RC$fix)




    #### 2. Create Scenarios for Overdesign ####
      ### for each overdesign scenario the product variety is increased starting with one product until all products are introduced
        ## method, starting with the high-volume and than going into niches
        # order_introduction <- sort(EAD$DEMAND,decreasing = T,index.return=TRUE)$ix
        ## method: high volume products higher chance but not neceassry introduced first
        ix <- 1:length(EAD$DEMAND)
        products_out <- ix
        order_introduction <- vector(mode="numeric")
        repeat{
          order_introduction <- c(order_introduction,sample(ix[products_out],
                                                            size = 1,
                                                            prob = EAD$DEMAND[products_out]))
          products_out <- setdiff(ix,order_introduction)
          if(length(order_introduction)==(length(ix)-1)){
            order_introduction <- c(order_introduction,ix[products_out])
            break
          }
        }

      #### 2.1 Overdesign components step by step until only one component is left ####
        i <- 1
        out <- list()
      repeat{
        #### 3. Product Variety ####
          # p<-10
          conceptCosts <- lapply(1:length(order_introduction),function(p){
                                out<-list()
                                ### 3.0 Exclude Products ###
                                products_in <-order_introduction[1:p]
                                products_out <- setdiff(order_introduction,products_in)
                                DEMAND_temp <- rep(0,length(EAD$DEMAND))
                                DEMAND_temp[products_in] <- EAD$DEMAND[products_in]

                                ### 3.1 Increased Material Costs ###
                                TC_var_S0 <- sum((nonCC$PC_direct + nonCC$PC_B_var) * DEMAND_temp)
                                out['dvl_materialCosts'] <- materialCosts(TC_var_S0 = TC_var_S0,
                                                                          P_RD = EAD$P$RD,
                                                                          DMD = DEMAND_temp,
                                                                          RCU = nonCC$RCU + nonCC$RCU_direct)

                                ### 3.2 Calculate Cost by going through the individual cost driver ###
                                ### starting with development costs ###
                                out['dvl_developmentCosts'] <- development_costs(P_PD = EAD$P$PD,
                                                                      DEMAND = DEMAND_temp,
                                                                      C_dvl = costDriver_inital$C_dvl)

                                ### 3.3 Part Management Costs ###
                                out['dvl_partadminCosts'] <- partmgmtCosts(P_PD = EAD$P$PD,
                                              DEMAND = DEMAND_temp,
                                              C_PA = costDriver_inital$C_pa)

                                ### 3.3 Setup Costs ###
                                setup<- setupCosts(P_PD = EAD$P$PD,
                                                   DMM_PD_PrD = EAD$DMM$PD_PrD,
                                                   DSM_PrD = EAD$DSM$PrD,
                                                   DEMAND = DEMAND_temp,
                                                   C_hold = R_hold,
                                                   C_order = R_order,
                                                   C_setup = costDriver_inital$C_setup)
                                out['man_setupCosts'] <- setup$TC_setup

                                ### 3.4. Tooling Costs ###
                                tooling <- toolingCosts(DMM_PD_PrD = EAD$DMM$PD_PrD,
                                                                        DSM_PrD = EAD$DSM$PrD,
                                                                        DMD_component = setup$DMD_component,
                                                                        C_tooling = costDriver_inital$C_tooling)
                                out['man_toolingCosts'] <- tooling$C_tooling

                                ### 3.5 Purchasing Order Costs ###
                                order <- orderCosts(DMD_PD = setup$DMD_component,
                                                                N_lot = setup$lotSize,
                                                                C_order = costDriver_inital$C_order)
                                out['pur_orderCosts'] <- order$TC_order

                                ### 3.6 Supplier Management Costs ###
                                out['pur_supplyCosts'] <- supplyCosts(P_PD = EAD$P$PD,
                                                                      DEMAND = DEMAND_temp,
                                                                      C_supply = costDriver_inital$C_supply)

                                ### 3.7. Stock costs ###
                                out['pur_stockCosts'] <- stockCosts(C_hold = costDriver_inital$C_hold,
                                                                     lotSize = setup$lotSize)

                                ### 3.8 Calculate total complexity costs ###
                                out['TC_CC'] <- sum(unlist(out))
                                out['TC_CC_var'] <- sum(out$man_setupCosts,out$pur_orderCosts,out$pur_stockCosts,out$dvl_materialCosts)
                                out['TC_CC_fix'] <- sum(out$dvl_developmentCosts,out$dvl_partadminCosts,out$man_toolingCosts,out$pur_supplyCosts)

                                ### 3.9 Calculate non complexity costs ###
                                nonCC_scenario <- clc_PCB(EAD$P$RD,
                                                          DMD = DEMAND_temp,
                                                          RCU = nonCC$RCU,
                                                          RCU_direct = nonCC$RCU_direct,
                                                          RC_fix = EAD$RC$fix)
                                out['TC_NC'] <- nonCC_scenario$TC
                                out['TC_NC_var'] <- nonCC_scenario$TC_var
                                out['TC_NC_fix'] <- nonCC_scenario$TC_fix
                                out['N_PROD_step'] <- sum(DEMAND_temp>0)
                                out['DMD_perc'] <- sum(DEMAND_temp) / sum(EAD$DEMAND)
                                out['LZM'] <- mean(setup$lotSize)
                                out['N_setups'] <- sum(setup$n_setups)
                                out['N_proVar'] <- tooling$N_processVariety
                                out['N_order'] <- order$N_order

                                return(as.data.frame(out))
                      }) # end product variety loop
          #### 4. Measure Output ####
            conceptCosts <- data.table::rbindlist(conceptCosts)
            systemMeasures <- t(data.frame(unlist(EAD$measures$SYSTEM)))
            rownames(systemMeasures) <- NULL
            out[[i]]<-data.frame(step=i-1,
                                 systemMeasures,
                                 conceptCosts)

          if(sum(colSums(EAD$DMM$FD_PD)>0)==1){
            break
          }else{
            EAD<-overdesign_EAD(EAD)
            i <- i + 1
            # print(i)
          } # end break condition
      } # end overdesign loop
      res<-data.table::rbindlist(out,use.names = T)
    return(res)
  }) # end DOE
  res<-data.table::rbindlist(res) %>%
    as_tibble() %>%
    select_if(~sum(!is.na(.)) > 0)
  gc()
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
                                                 "dplyr", "rpicosat", "tidyr",
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
