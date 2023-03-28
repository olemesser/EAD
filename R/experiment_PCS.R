experiment_PCS<-function(DOE){
  suppressWarnings(suppressMessages(require(EAD)))
  suppressWarnings(suppressMessages(require(data.table)))
  df<-lapply(1:NROW(DOE), function(i){
    #### Create EAD ####
      output<-list()
      EAD <- crt_EAD(DOE[i,])
      measures_system <- data.frame(t(unlist(EAD[[1]]$measures$SYSTEM)))
      output[['measures_system']] <- t(na.omit(measures_system))

      costs <- clc_PCB(
                        RES_CONS_PAT = EAD[[1]]$P$RD,
                        DMD = EAD[[1]]$DEMAND,
                        RC_fix =  EAD[[1]]$RC$fix,
                        RC_direct = EAD[[1]]$RC$direct,
                        RC_var = EAD[[1]]$RC$var)
      RCU<-costs$RCU
      TRC<-costs$TRC
      RCU_direct <- costs$RCU_direct

    #### 2. Create Variety Scenario ####
      ### start with 30% of variety and increase 2% in each step
      DMD_sort <- sort(EAD[[1]]$DEMAND,decreasing = T,index.return=TRUE)
      step_width<-ceiling(length(DMD_sort$ix)*0.02)
      products_in <- DMD_sort$ix[1:ceiling(length(DMD_sort$ix)*0.3)]
      products_out <- setdiff(DMD_sort$ix,products_in)
      ### for each scenario calculate product costs
      ### update product level measures (PCI, DNS)
      system_level_data<-list()
      l<-1
      repeat{
        if(l>1){
          step_width_max <- ifelse(length(products_out)>step_width,step_width,length(products_out))
          products_in <- c(products_in,products_out[1:step_width_max])
          products_out <- products_out[-c(1:step_width_max)]
        }

        P_FD <- EAD[[1]]$P$FD[products_in,,drop=F]
        P_PD <- EAD[[1]]$P$PD[products_in,,drop=F]
        P_RD <- EAD[[1]]$P$RD[products_in,,drop=F]
        resources_to_keep <- which(colSums(P_RD)>0)
        DMD <- EAD[[1]]$DEMAND[products_in]

        DMD_zero <- rep(0,length(EAD[[1]]$DEMAND))
        DMD_zero[products_in] <- EAD[[1]]$DEMAND[products_in]
        #### Calculate True Costs ####
        benchmark <- clc_PCB(
                             RES_CONS_PAT = EAD[[1]]$P$RD[,resources_to_keep],
                             DMD = DMD_zero,
                             RCU = RCU[resources_to_keep],
                             RC_fix =  ifelse(colSums(EAD[[1]]$P$RD)==0,0,EAD[[1]]$RC$fix)[resources_to_keep],
                             RCU_direct = costs$RCU_direct[resources_to_keep])
        PC_B <- benchmark$PC_B
        PC_direct <- benchmark$PC_direct
        # # check fixed_costs
        # sum(ifelse(colSums(EAD[[1]]$P$RD)==0,0,EAD[[1]]$RC$fix)[resources_to_keep]) == sum(benchmark$PC_B_fixed *DMD_zero)
        # ##check variable costs
        # sum((EAD[[1]]$P$RD[,resources_to_keep] %*% RCU[resources_to_keep])*DMD_zero) == sum(benchmark$PC_B_var *DMD_zero)
        # ## check direct costs
        # sum((benchmark$PC_B-benchmark$PC_B_indirect)*DMD_zero) == sum((EAD[[1]]$P$RD[,resources_to_keep] %*% costs$RCU_direct[resources_to_keep])*DMD_zero)



        #### Create Measure Object ####
        measures_system_temp<-data.frame(PCI.FD = measure_PCI(P_FD),
                                         PCI.PD = measure_PCI(P_PD),
                                         DNS.RD = measure_DENS(P_RD[,resources_to_keep]),
                                         DMD_T10 = measure_TOP10(DMD),
                                         Q_VAR = sd(log(DMD)),
                                         NPV.FD = measure_NPV(P_FD),
                                         N_RD = NCOL(P_RD[,resources_to_keep]),
                                         N_PROD = NROW(P_RD[,resources_to_keep]),
                                         SDC_n.FD_PD = measures_system$SDC_n.FD_PD,
                                         SDC_n.PD_PrD = measures_system$SDC_n.PD_PrD,
                                         SDC_n.PrD_RD = measures_system$SDC_n.PrD_RD,
                                         HIC_n.PD = measures_system$HIC_n.PD,
                                         HIC_n.PrD = measures_system$HIC_n.PrD,
                                         HIC_n.RD = measures_system$HIC_n.RD,
                                         RC_cor.indirect = measures_system$RC.cor_indirect,
                                         r_indirect = sum(benchmark$PC_B_indirect*DMD_zero) /sum(PC_B*DMD_zero),
                                         TC = sum(benchmark$PC_B * DMD_zero),
                                         RC_indirect = sum(benchmark$PC_B_indirect*DMD_zero),
                                         TC_fix = sum(benchmark$PC_B_fixed * DMD_zero),
                                         TC_var = sum(benchmark$PC_B_var * DMD_zero))

        #### Calculate Reported Costs #####
        cs_DOE<-expand.grid(ACP=c(1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,25,30),
                            method = c("random","correl-random","PU","DC-0.2"))
        reported<-lapply(1:NROW(cs_DOE),function(c){
          if(cs_DOE$method[c] %in% c("PU","DC") | startsWith(as.character(cs_DOE$method[c]),"DC")){
            PC_H <- costingSystem_VD(RES_CONS_PAT = P_RD[,resources_to_keep],
                                     RC_indirect = RCU[resources_to_keep] * benchmark$TRC +
                                       ifelse(colSums(P_RD)==0,0,EAD[[1]]$RC$fix)[resources_to_keep],
                                     RC_direct = RCU_direct[resources_to_keep] * benchmark$TRC,
                                     DMD = DMD,
                                     method = cs_DOE$method[c])
          }else if (cs_DOE$method[c] %in% c("random","size-misc","correl-random")){
            PC_H <- costingSystem_ABC(RES_CONS_PAT = P_RD[,resources_to_keep],
                                      DMD = DMD,
                                      RC_indirect = (costs$RCU[resources_to_keep] * benchmark$TRC) +
                                                        ifelse(colSums(P_RD)==0,0,EAD[[1]]$RC$fix)[resources_to_keep],
                                      RD = cs_DOE$method[c],
                                      AD = "big-pool",
                                      ACP = cs_DOE$ACP[c])

            # sum(RC_indirect) == sum(benchmark$PC_B_indirect[products_in]*DMD)
            # sum(PC_H*DMD) == sum(benchmark$PC_B_indirect[products_in]*DMD)
            # PC_H <- PC_H + PC_direct[products_in]
            # clc_costingERROR(PC_B[products_in],PC_H)$EUCD

          }
          PC_H <- PC_H + PC_direct[products_in]
          costingError<-clc_costingERROR(PC_B[products_in],PC_H)
          costingError_w<-clc_costingERROR(PC_B[products_in],PC_H,DMD=DMD)
          return(list(EUCD = costingError$EUCD,
                      EUCD_w = costingError_w$EUCD,
                      MPE = costingError$MPE,
                      MPE_w = costingError_w$MPE,
                      PC_H = PC_H,
                      PC_B = PC_B[products_in],
                      TC = sum(PC_B[products_in] * DMD)))
        })

        #### Create Output Object ####
          ### On System level for each variety step ###
            error <- lapply(reported,function(x){
                      data.frame(EUCD = x$EUCD,
                                 MPE = x$MPE,
                                 EUCD_w = x$EUCD_w,
                                 MPE_w = x$MPE_w,
                                 TC = x$TC)
                      })
            error <- data.table::rbindlist(error)
            system_level_data[[l]]<-data.frame(id = EAD[[1]]$ID,
                                               merge(measures_system_temp,cbind(cs_DOE,error)))
            l <- l+1

        if(length(products_out)==0){
          ### On product level only under full variety ###
          ### select some ACPs
          product_level_data <- lapply(which(cs_DOE$ACP %in% c(4)),function(s){
            data.frame(id= EAD[[1]]$ID,
                       LOF.RD =  EAD[[1]]$measures$PRODUCT$LOF$RD[products_in],
                       INTER.RD = EAD[[1]]$measures$PRODUCT$INTER$RD[products_in],
                       INTRA.RD = EAD[[1]]$measures$PRODUCT$INTRA$RD[products_in],
                       DMD_perc = EAD[[1]]$DEMAND[products_in]/sum(EAD[[1]]$DEMAND[products_in]),
                       ACP = cs_DOE$ACP[s],
                       method = cs_DOE$method[s],
                       PC_B = reported[[s]]$PC_B,
                       PC_H = reported[[s]]$PC_H,
                       PC_direct = PC_direct[products_in],
                       PC_B_var = benchmark$PC_B_var[products_in],
                       PC_B_fix = benchmark$PC_B_fixed[products_in])
          })
          product_level_data <- data.table::rbindlist(product_level_data)
          system_level_data <- data.table::rbindlist(system_level_data)
          break
        }
      }
    return(list(system_level = system_level_data,
                product_level= product_level_data))
  })

  output<-list(system_level=data.table::rbindlist(lapply(df,function(x) x$system_level)),
               product_level=data.table::rbindlist(lapply(df,function(x) x$product_level)))
  return(output)
}



experiment_PCS_MC<-function(DOE=NULL,
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
                                   res<-with_timeout(experiment_PCS(DOE = ip),
                                                     timeout = time_limit)
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
    res<-par_apply(cl,X=DOE_list,FUN=experiment_PCS,errorhandling=ehNodes)
  }
  res<-res[sapply(res,length)>0]
  return(res)

}
