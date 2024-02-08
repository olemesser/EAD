experiment_PCS<-function(DOE,
                         cs_DOE = expand.grid(ACP=c(1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,25,30),
                                              method = c("random","correl-random","DIV","DLH-0.2")),
                         ACP_productLevel = c(5)
                         ){
  #### Input for Testing ####
  # require(EAD)
  # require(tidyr)
  # i<-1
  # DOE<-expand_grid(N_FR = list(c(30)), # number of functional requirements
  #                  N_DD = list(c(50)), # number of physical domain elements
  #                  N_PrD = list(c(50)), # number of process domain elements
  #                  N_RD = list(c(50)), # number of resource domain elements
  #                  DNS = list(c(0.07,0.2)), # density of the P_FD matrix.
  #                  method_FD = "DNS",
  #                  N_PROD = list(c(200,200)), # number of products
  #                  prod_step_width = 20,
  #                  TOTAL_DEMAND = list(c(200,25000)), # total demand
  #                  Q_VAR = list(c(0,1.7)), # coefficient of variation for demand distribution
  #                  DMM_PAR = expand_grid(FD_PD = list(c(0,0.08)),
  #                                        PD_PrD = list(c(0,0.08)),
  #                                        PrD_RD = list(c(0,0.08))), # desired design complexity
  #                  allowZero = F,
  #                  uB_DMM = 2,
  #                  ut_DMM = F, # if the upper triangle DMMs should be generated too (DMM_FD_PrD,DMM_FD_RD,DMM_PD_RD)
  #                  DSM_param=expand_grid(PD = list(c(0,0.14,0,1)),
  #                                        PrD = list(c(0,0.0,0,1)),
  #                                        RD = list(c(0,0.0,0,1))), # first two entries refer to the density of the DSMs and the second pair to the cv if the DSM_method='modular' is used.
  #                  DSM_method='modular',
  #                  ub_DSM = 2,
  #                  TC = 10^6, # total costs
  #                  r_in = list(c(0.2,0.8)),
  #                  r_fix = list(c(0.2,0.6)), # proportion of fixed costs on total indirect costs
  #                  RC_sdlog = list(c(0,3)), # coefficient of variation for resource cost distribution # coefficient of variation for resource cost distribution
  #                  N_RUN = 1:2 # number of runs
  # )
  # cs_DOE = expand.grid(ACP=c(1,2,3,4,5,10,50),
  #                      method = c("random","correl-random","DIV"))
  # ACP_productLevel = c(1,5)
  #### End INput for Testing ####


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
                        RC_var_i = EAD[[1]]$RC$var_i,
                        RC_var_d = EAD[[1]]$RC$var_d,
                        RC_fix_i =  EAD[[1]]$RC$fix_i,
                        RC_fix_d =  EAD[[1]]$RC$fix_d)
      TRC<-costs$TRC
      # sum((costs$PC_B_fix_d + costs$PC_B_fix_i + costs$PC_B_var_d + costs$PC_B_var_i)* EAD[[1]]$DEMAND)
      # sum(RCU_indirect * TRC)
      # sum((costs$PC_B_var_i) * EAD[[1]]$DEMAND)

    #### 2. Create Variety Scenario ####
      ### start with 30% of variety and increase 2% in each step
      ix <- 1:length(EAD[[1]]$DEMAND)
      products_out <- ix
      order_introduction <- vector(mode="numeric")
      repeat{
        order_introduction <- c(order_introduction,
                                sample(ix[products_out],
                                      size = 1,
                                      prob = EAD[[1]]$DEMAND[products_out]
                                      )
                                )
        products_out <- setdiff(ix,order_introduction)
        if(length(order_introduction)==(length(ix)-1)){
          order_introduction <- c(order_introduction,ix[products_out])
          break
        }
      }

      order_introduction_index <- 1:length(order_introduction)
      variety_step <- split(order_introduction_index, ceiling(seq_along(order_introduction_index) / DOE$prod_step_width[i]))
      order_introduction <- lapply(variety_step,function(x){
        order_introduction[x]
      })


      system_level_data<-list()
      l<-1
      p <- 0
      repeat{
        p <- p +1
        products_in <- as.numeric(unlist(order_introduction[1:p]))
        products_out <- setdiff(as.numeric(unlist(order_introduction)),
                                products_in)


        P_FD <- EAD[[1]]$P$FD[products_in,,drop=F]
        P_PD <- EAD[[1]]$P$PD[products_in,,drop=F]
        P_PrD <- EAD[[1]]$P$PrD[products_in,,drop=F]
        P_RD <- EAD[[1]]$P$RD[products_in,,drop=F]
        resources_to_keep <- which(colSums(P_RD)>0)
        P_FD <- P_FD[,which(colSums(P_FD)>0),drop=F]
        P_PD <- P_PD[,which(colSums(P_PD)>0),drop=F]
        P_PrD <- P_PrD[,which(colSums(P_PrD)>0),drop=F]
        P_RD <- P_RD[,resources_to_keep,drop=F]
        DMD <- EAD[[1]]$DEMAND[products_in]

        #### Calculate True Costs ####
        benchmark <- clc_PCB(
                             RES_CONS_PAT = P_RD,
                             DMD = DMD,
                             RCU_indirect = costs$RCU_indirect[resources_to_keep],
                             RCU_direct = costs$RCU_direct[resources_to_keep],
                             RC_fix_i =  EAD[[1]]$RC$fix_i[resources_to_keep],
                             RC_fix_d =  EAD[[1]]$RC$fix_d[resources_to_keep])
        PC_direct <- benchmark$PC_B_var_d + benchmark$PC_B_fix_d
        PC_B <-  benchmark$PC_B_indirect + PC_direct
        #benchmark$TC
        #sum((benchmark$PC_B_var_i + benchmark$PC_B_var_d + benchmark$PC_B_fix_i + benchmark$PC_B_fix_d) * DMD)

        # # check fixed_costs
        #sum(EAD[[1]]$RC$fix[resources_to_keep]) == sum(benchmark$PC_B_fixed *DMD)
        # ##check variable costs
        #sum((P_RD %*% costs$RCU_indirect[resources_to_keep])*DMD) == sum(benchmark$PC_B_var_i *DMD)
        # ## check direct costs
        #sum((benchmark$PC_B_var_d)*DMD) == sum((P_RD %*% costs$RCU_direct[resources_to_keep])*DMD)


        #### Create Measure Object ####
        RC_fix_i <-  EAD[[1]]$RC$fix_i[resources_to_keep]
        RC_var_i <- costs$RCU_indirect[resources_to_keep] * colSums(P_RD * DMD)
        RC_i <- RC_fix_i + RC_var_i
        cor_mat <- cor(EAD[[1]]$P$RD[products_in,resources_to_keep])
        measures_system_temp<-data.frame(PCI.FD = measure_PCI(P_FD),
                                         PCI.PD = measure_PCI(P_PD),
                                         DNS.RD = measure_DENS(P_RD),
                                         DMD_T10 = measure_TOP10(DMD),
                                         TOTAL_DMD = sum(DMD),
                                         RC_i_T10 = measure_TOP10(RC_i),
                                         COR = mean(cor_mat),
                                         Q_VAR = sd(log(DMD)),
                                         NPV.FD = measure_NPV(P_FD),
                                         N_RD = NCOL(P_RD),
                                         N_PROD = NROW(P_RD),
                                         SDC_n.FD_PD = measures_system$SDC_n.FD_PD,
                                         SDC_n.PD_PrD = measures_system$SDC_n.PD_PrD,
                                         SDC_n.PrD_RD = measures_system$SDC_n.PrD_RD,
                                         HIC_n.PD = measures_system$HIC_n.PD,
                                         HIC_n.PrD = measures_system$HIC_n.PrD,
                                         HIC_n.RD = measures_system$HIC_n.RD,
                                         D_RD = measure_diversificationINDEX(P_RD, DMD = DMD),
                                         D_RD_noDemand = measure_diversificationINDEX(P_RD),
                                         INTER.FD = mean(measure_INTER(P_FD)),
                                         INTRA.FD = mean(measure_INTRA(P_FD)),
                                         INTER.PD = mean(measure_INTER(P_PD)),
                                         INTRA.PD = mean(measure_INTRA(P_PD)),
                                         INTER.PrD = mean(measure_INTER(P_PrD)),
                                         INTRA.PrD = mean(measure_INTRA(P_PrD)),
                                         INTER.RD = mean(measure_INTER(P_RD)),
                                         INTRA.RD = mean(measure_INTRA(P_RD)),
                                         R_id = sum(benchmark$PC_B_indirect * DMD) /sum(PC_B * DMD),
                                         TC = sum(benchmark$PC_B * DMD),
                                         TC_indirect = sum(benchmark$PC_B_indirect * DMD),
                                         TC_fix_i = sum((benchmark$PC_B_fix_i) * DMD),
                                         TC_var_i = sum((benchmark$PC_B_var_i) * DMD),
                                         TC_fix = sum((benchmark$PC_B_fix_i + benchmark$PC_B_fix_d) * DMD),
                                         TC_var = sum((benchmark$PC_B_var_i + benchmark$PC_B_var_d) * DMD)) %>%
          mutate(UNIT_SHARE_i =  TC_var_i / (TC_var_i + TC_fix_i),
                 UNIT_SHARE = TC_var / (TC_var + TC_fix))

        #### Calculate Reported Costs #####
        reported<-lapply(1:NROW(cs_DOE),function(c){
          if(cs_DOE$method[c] %in% c("DIV","DLH") | startsWith(as.character(cs_DOE$method[c]),"DLH")){
            PC_H_indirect <- costingSystem_VD(RES_CONS_PAT = P_RD,
                                     RC_indirect = costs$RCU_indirect[resources_to_keep] * benchmark$TRC + EAD[[1]]$RC$fix_i[resources_to_keep],
                                     RC_direct = rep(0,length(costs$RCU_indirect[resources_to_keep])),  # costs$RCU_direct[resources_to_keep] * benchmark$TRC + EAD[[1]]$RC$fix_d[resources_to_keep],
                                     DMD = DMD,
                                     method = cs_DOE$method[c])
          }else if (cs_DOE$method[c] %in% c("random","size-misc","correl-random")){
            PC_H_indirect <- costingSystem_ABC(RES_CONS_PAT = P_RD,
                                      DMD = DMD,
                                      RC_var_i = costs$RCU_indirect[resources_to_keep] * benchmark$TRC,
                                      RC_fix_i = EAD[[1]]$RC$fix_i[resources_to_keep],
                                      RD = cs_DOE$method[c],
                                      AD = "big-pool",
                                      ACP = cs_DOE$ACP[c])

            #sum(PC_H_indirect*DMD) == sum(benchmark$PC_B_indirect*DMD)
          }
          #RC_indirect = costs$RCU_indirect[resources_to_keep] * benchmark$TRC + ifelse(colSums(P_RD)==0,0,EAD[[1]]$RC$fix_i)[resources_to_keep]
          #cond_1 = (sum(RC_indirect) == sum(benchmark$PC_B_indirect*DMD))
          costingError<-clc_costingERROR(PC_B = benchmark$PC_B_indirect,
                                         PC_H = PC_H_indirect,
                                         DMD = DMD)
          return(list(EUCD = costingError$EUCD,
                      MAPE = costingError$MAPE,
                      EUCD_w = costingError$EUCD_w,
                      MAPE_w = costingError$MAPE_w,
                      PC_H = PC_H_indirect + PC_direct,
                      PC_B =  PC_B,
                      PC_H_i = PC_H_indirect,
                      PC_B_i = benchmark$PC_B_indirect,
                      TC = sum(PC_B * DMD)))
        })

        #### Create Output Object ####
          ### On System level for each variety step ###
            error <- lapply(reported,function(x){
                      data.frame(EUCD = x$EUCD,
                                 MAPE = x$MAPE,
                                 EUCD_w = x$EUCD_w,
                                 MAPE_w = x$MAPE_w)
                      })
            error <- data.table::rbindlist(error)
            system_level_data[[l]]<-tibble(id = EAD[[1]]$ID,
                                           measures_system_temp,
                                           cbind(cs_DOE,error))

            l <- l+1

        if(length(products_out)==0){
          ### On product level only under full variety ###
          ### select some ACPs
          product_level_data <- lapply(which(cs_DOE$ACP %in% ACP_productLevel),function(s){
            data.frame(id = EAD[[1]]$ID,
                       INTER.PD = EAD[[1]]$measures$PRODUCT$INTER$PD[products_in],
                       INTRA.PD = EAD[[1]]$measures$PRODUCT$INTRA$PD[products_in],
                       LOF.PD = EAD[[1]]$measures$PRODUCT$LOF$PD[products_in],
                       LOF.RD =  EAD[[1]]$measures$PRODUCT$LOF$RD[products_in],
                       INTER.RD = EAD[[1]]$measures$PRODUCT$INTER$RD[products_in],
                       INTRA.RD = EAD[[1]]$measures$PRODUCT$INTRA$RD[products_in],
                       DMD_perc = EAD[[1]]$DEMAND[products_in]/sum(EAD[[1]]$DEMAND[products_in]),
                       ACP = cs_DOE$ACP[s],
                       method = cs_DOE$method[s],
                       PC_B = reported[[s]]$PC_B,
                       PC_H = reported[[s]]$PC_H,
                       PC_direct = PC_direct[products_in],
                       PC_B_var = benchmark$PC_B_var_i[products_in],
                       PC_B_fix = benchmark$PC_B_fix_i[products_in])
          })
          product_level_data <- data.table::rbindlist(product_level_data)
          system_level_data <- data.table::rbindlist(system_level_data)
          break
        }
      }
    return(list(system_level = system_level_data,
                product_level= product_level_data))
  })

  output<-list(system_level = data.table::rbindlist(lapply(df,function(x) x$system_level)),
               product_level = data.table::rbindlist(lapply(df,function(x) x$product_level)))
  return(output)
}



experiment_PCS_MC<-function(DOE,
                            cs_DOE,
                            ACP_productLevel = 0,
                            NUMB_CORES = 4,
                            logfile = "",
                            time_limit = Inf,
                            cl = NULL,
                            ehNodes = "remove"){
  suppressMessages(suppressWarnings(require(parallel)))
  suppressMessages(suppressWarnings(require(doSNOW)))
  suppressMessages(suppressWarnings(require(foreach)))
  suppressMessages(suppressWarnings(require(dplyr)))
  suppressMessages(suppressWarnings(require(R.utils)))

  DOE_list<-DOE %>%
    group_split(split=1:n())%>%
    as.list()

  if(is.null(cl)) cl<-EAD::setupMC(NUMB_CORES=NUMB_CORES,logfile=logfile)
  on.exit(stopCluster(cl))

  parallel::clusterExport(cl, envir = environment(),c("time_limit","cs_DOE"))
  print(cl)
  res <- run_MC(cl, X = DOE_list,
                               FUN = function(ip, ...) {
                                 res<-list()
                                 res<-with_timeout(experiment_PCS(DOE = ip,cs_DOE = cs_DOE, ACP_productLevel = ACP_productLevel),
                                                   timeout = time_limit)
                                 gc()
                                 return(res)
                               }, packages = c("EAD",
                                               "dplyr", "tidyr",
                                               "Matrix", "digest", "data.table",
                                               "plyr", "DescTools", "igraph", "R.utils"),
                               errorhandlingNodes = ehNodes)


  res<-res[sapply(res,length)>0]
  return(res)

}
