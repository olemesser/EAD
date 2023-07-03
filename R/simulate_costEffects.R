simulate_costEffects<-function(DOE){
  suppressMessages(suppressWarnings(require(EAD)))
  suppressMessages(suppressWarnings(require(digest)))
  suppressMessages(suppressWarnings(require(tidyr)))
  suppressMessages(suppressWarnings(require(data.table)))
  suppressMessages(suppressWarnings(require(dplyr)))


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
      bounds <- DOE$bounds[x][[1]]

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
        ## product costs and resource unit costs ##
        nonCC <- clc_PCB(EAD$P$RD,
                         DMD=EAD$DEMAND,
                         RC_direct = EAD$RC$direct,
                         RC_var = EAD$RC$var,
                         RC_fix = EAD$RC$fix)
        ## component costs ##
        PDUC_noOD <- clc_variableComponentCosts(EAD = EAD,
                                           RCU = nonCC$RCU_direct)
        PDUC <- PDUC_noOD
        EAD_PD_noOD <- EAD$P$PD


    #### 2. Create Scenarios for Overdesign ####
      ### for each overdesign scenario the product variety is increased starting with one product until all products are introduced
        ## method, starting with the high-volume and than going into niches
        # order_introduction <- sort(EAD$DEMAND,decreasing = T,index.return=TRUE)$ix
        ## method: high volume products higher chance but not necessary introduced first
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

        order_introduction_index <- 1:length(order_introduction)
        variety_step <- split(order_introduction_index, ceiling(seq_along(order_introduction_index) / DOE$prod_step_width[x]))
        order_introduction <- lapply(variety_step,function(x){
          order_introduction[x]
        })
      #### 2.1 Overdesign components step by step until only one component is left ####
        i <- 1
        out <- list()
        overdesign_Change <- NULL
        C_dvl <- costDriver_inital$C_dvl
      repeat{
        #### 3. Product Variety ####
          ### increase product variety within each step ###
          # p<-1
          conceptCosts <- lapply(1:length(order_introduction),function(p){
                                out<-list()
                                ### 3.0 Exclude Products ###
                                products_in <- as.numeric(unlist(order_introduction[1:p]))
                                products_out <- setdiff(as.numeric(unlist(order_introduction)),
                                                        products_in)
                                DEMAND_temp <- rep(0,length(EAD$DEMAND))
                                DEMAND_temp[products_in] <- EAD$DEMAND[products_in]

                                  ## calculate component demand rate for the scenario of no overdesign
                                    DMD_component_noOD <- as.numeric(DEMAND_temp %*% EAD_PD_noOD)
                                  ## calculate component demand rate for the current scenario
                                    DMD_component <- as.numeric(DEMAND_temp %*% EAD$P$PD)

                                ### 3.1 Increased Material Costs ###
                                    variableComponenCosts_noOD <- sum(PDUC_noOD$var * DMD_component_noOD)
                                  ## calculate costs of overdesign
                                  out['dvl_materialCosts'] <- sum(PDUC$var * DMD_component) - variableComponenCosts_noOD


                                ### 3.2 Calculate Cost by going through the individual cost driver ###
                                ### starting with development costs ###
                                  out['dvl_developmentCosts'] <- development_costs(P_PD = EAD$P$PD,
                                                                      DEMAND = DEMAND_temp,
                                                                      C_dvl = C_dvl)$TC_dvl

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
                                stock <- stockCosts(C_hold = costDriver_inital$C_hold,
                                                                     lotSize = setup$lotSize)
                                out['pur_stockCosts'] <- stock$TC_stock

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

                                out['TC_NC'] <- nonCC_scenario$TC - out$dvl_materialCosts
                                out['TC_NC_var'] <- nonCC_scenario$TC_var - out$dvl_materialCosts
                                out['TC_NC_fix'] <- nonCC_scenario$TC_fix
                                out['TC'] <- out$TC_NC +  out$TC_CC
                                out['N_PROD_step'] <- sum(DEMAND_temp>0)
                                out['DMD_perc'] <- sum(DEMAND_temp) / sum(EAD$DEMAND)
                                out['LZM'] <- mean(setup$lotSize[setup$lotSize>0])
                                out['DMD_PD_mean'] <- mean(DMD_component)
                                out['DMD_PD_sum'] <- sum(DMD_component)
                                out['Units_stock'] <- stock$Units_stock
                                out['N_setups'] <- sum(setup$n_setups)
                                out['N_proVar'] <- tooling$N_processVariety
                                out['N_order'] <- sum(order$N_order)
                                out['driver_dvl'] <- mean(costDriver_inital$C_dvl[!is.nan(costDriver_inital$C_dvl)])
                                out['driver_pa'] <- mean(costDriver_inital$C_pa[!is.nan(costDriver_inital$C_pa)])
                                out['driver_setups'] <- costDriver_inital$C_setup
                                out['driver_tooling'] <- mean(costDriver_inital$C_tooling[costDriver_inital$C_tooling>0])
                                out['driver_order'] <- costDriver_inital$C_order
                                out['driver_supply'] <- costDriver_inital$C_supply
                                out['driver_stock'] <- costDriver_inital$C_hold
                                out['ADC'] <-   as.numeric(out['dvl_developmentCosts']) / sum(DEMAND_temp)
                                out['order_to_hold'] <- R_order / R_hold

                                return(as.data.frame(out))
                      }) # end product variety loop
          #### 4. Measure Output ####
            conceptCosts <- data.table::rbindlist(conceptCosts,use.names = T)
            systemMeasures <- t(data.frame(unlist(EAD$measures$SYSTEM)))
            rownames(systemMeasures) <- NULL
            out[[i]]<-data.frame(ID=EAD$ID,
                                 step=i-1,
                                 systemMeasures,
                                 conceptCosts)

          if(sum(colSums(EAD$DMM$FD_PD)>0)<=ceiling(NROW(EAD$DMM$FD_PD)/2)){
            break
          }else{
            EAD<-overdesign_EAD(EAD,bounds=DOE$bounds[x][[1]],method = DOE$method_overdesign[x])
            overdesign_Change <- c(EAD$overdesign$substitute,EAD$overdesign$replaced_by)
            EAD <- EAD$EAD
            PDUC <- overdesign_costs(PDUC_var = PDUC$var,
                                     design = overdesign_Change,
                                     bounds = bounds)$PDUC

            C_dvl <- development_costs(P_PD = EAD$P$PD,
                                       DEMAND = EAD$DEMAND,
                                       C_dvl = C_dvl,
                                       overdesign_Change=overdesign_Change,
                                       bounds = bounds)$C_dvl

            i <- i + 1
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
                                                 "dplyr", "tidyr",
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
  # res<-res[sapply(res,length)>0]
  return(res)
}
