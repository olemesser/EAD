crt_EAD<-function(DOE,
                  file_output=F){
  suppressMessages(suppressWarnings(require(digest)))
  suppressMessages(suppressWarnings(require(tidyr)))

  #### Input Testing ####
  # x<-1
  # DOE<-expand_grid(N_FR = list(c(15)), # number of functional requirements
  #                  N_DD = list(c(18,26)), # number of physical domain elements
  #                  N_PrD = list(c(36,52)), # number of process domain elements
  #                  N_RD = list(c(72,104)), # number of resource domain elements
  #                  DNS = seq(0.07,0.5,0.05), #density within the DSM_FD matrix. Creates product mixes where not all products are included
  #                  N_PROD = 50,
  #                  method_FD = "DNS",
  #                  TOTAL_DEMAND = 10000, # total demand
  #                  Q_VAR = list(c(0,3)), # demand heterogeneity
  #                  DMM_PAR = expand_grid(FD_PD=list(c(0,0.08)),
  #                                        PD_PrD=list(c(0,0.1)),
  #                                        PrD_RD=list(c(0,0.1))), # desired system design complexity
  #                  uB_DMM = 1,
  #                  allowZero = F,
  #                  ut_DMM = F, # if the upper triangle DMMs should be generated too (DMM_FD_PrD,DMM_FD_RD,DMM_PD_RD)
  #                  DSM_param=expand_grid(PD=list(c(0,0.05,0,1)),
  #                                        PrD=list(c(0,0.05,0,1)),
  #                                        RD=list(c(0,0.05,0,1))), # first two entries refer to the density of the DSMs and the second pair to the cv if the DSM_method='modular' is used.
  #                  DSM_method='modular',
  #                  ub_DSM = 1,
  #                  TC = 10^6, # total costs
  #                  r_in = list(c(0,0.9)),
  #                  r_fix = list(c(0,1)), # proportion of fixed costs on total costs
  #                  cor_var = list(c(-1,1)), # correlation between indirect variable cost vector and direct cost vector
  #                  cor_fix = list(c(-1,1)), # correlation between indirect fixed cost vector and direct cost vector
  #                  RC_cv = list(c(0,1.5)), # coefficient of variation for resource cost distribution
  #                  N_RUN = 1:1 # number of runs
  # )
  # DOE<-DOE[1:4,]
  #### END Input Testing ####


  #### For each DOE Setting ####
  EAD<-lapply(1:NROW(DOE),function(x){
    uB_DMM <- DOE$uB_DMM[x]
    ub_DSM <- DOE$ub_DSM[x]
    allowZero <- DOE$allowZero[x]
    measures<-list()
    DMM<-list()
    DSM<-list()
    message<-list()
    P<-list()
    DOE$N_FR[x]<-ifelse(length(DOE$N_FR[x][[1]])>1,sample(DOE$N_FR[[x]][1]:DOE$N_FR[[x]][2],1),DOE$N_FR[[x]][1])
    DOE$N_DD[x]<-ifelse(length(DOE$N_DD[x][[1]])>1,sample(DOE$N_DD[[x]][1]:DOE$N_DD[[x]][2],1),DOE$N_DD[[x]][1])
    DOE$N_PrD[x]<-ifelse(length(DOE$N_PrD[x][[1]])>1,sample(DOE$N_PrD[[x]][1]:DOE$N_PrD[[x]][2],1),DOE$N_PrD[[x]][1])
    DOE$N_RD[x]<-ifelse(length(DOE$N_RD[x][[1]])>1,sample(DOE$N_RD[[x]][1]:DOE$N_RD[[x]][2],1),DOE$N_RD[[x]][1])
    #### 1. Create Product Mix & Demand ####
      ### 1.1 Product Mix Generation ###
      if('DNS' %in% colnames(DOE)){
        DNS_set <- runif(1,DOE$DNS[x][[1]][1],DOE$DNS[x][[1]][2])
      }else{
        DNS_set <- 1
      }
      prodMIX<-create_ProductMix(N_FR = DOE$N_FR[x][[1]],
                        DNS = DNS_set,
                        prop_PROD =ifelse('prop_PROD' %in% colnames(DOE),DOE$prop_PROD[x],NA),
                        N_PROD = ifelse('N_PROD' %in% colnames(DOE),DOE$N_PROD[x],1),
                        method = DOE$method_FD[x])
      P[["FD"]]<-as.matrix(prodMIX$P_FD_const)


      ### 1.2 Demand Generation ###
      DMD<-crt_DEMAND(n_PROD = NROW(P[["FD"]]),
                 TOTAL_DEMAND = DOE$TOTAL_DEMAND[x],
                 Q_VAR = DOE$Q_VAR[x][[1]])
      measures[['SYSTEM']]<-c(measures[['SYSTEM']],
                              DMD_cv = DMD$CHECK$CV,
                              DMD_T10 = DMD$CHECK$DMD_T10,
                              Q_VAR = DMD$CHECK$Q_VAR)
      DMD<-DMD$DEMAND

    #### 2. Create Domains ####
      DOM<-crt_DOMAINS(P_FD = P[["FD"]],
                       N_DD = DOE$N_DD[x][[1]],
                       N_PrD = DOE$N_PrD[x][[1]],
                       N_RD = DOE$N_RD[x][[1]],
                       DMM_PAR = DOE$DMM_PAR[x,],
                       DSM_method = DOE$DSM_method[x],
                       DSM_param = DOE$DSM_param[x,],
                       ut_DMM=F,
                       uB_DMM=uB_DMM,
                       ub_DSM=ub_DSM,
                       allowZero=allowZero)

      DOM$measures$SYSTEM[['D']]<-list(FD=measure_diversificationINDEX(DOM$P$FD,DMD=DMD),
                                       PD=measure_diversificationINDEX(DOM$P$PD,DMD=DMD),
                                       PrD=measure_diversificationINDEX(DOM$P$PrD,DMD=DMD),
                                       RD=measure_diversificationINDEX(DOM$P$RD,DMD=DMD))

      measures[['SYSTEM']]<-c(measures[['SYSTEM']],DOM$measures$SYSTEM,
                              N_FD = NCOL(DOM$P$FD),
                              N_PD = NCOL(DOM$P$PD),
                              N_PrD = NCOL(DOM$P$PrD),
                              N_RD = NCOL(DOM$P$RD),
                              N_PROD = NROW(DOM$P$FD),
                              TSS=DOE$N_FR[x][[1]]+DOE$N_DD[x][[1]]+DOE$N_PrD[x][[1]]+DOE$N_RD[x][[1]])
      measures[['PRODUCT']]<-c(measures[['PRODUCT']],DOM$measures$PRODUCT)
      DSM<-DOM$DSM
      DMM<-DOM$DMM
      P<-DOM$P
      message[['DOMAIN']]<-DOM$message


    #### 3. Create Costs ####
          RC<-crt_RC(N_RD = DOE$N_RD[x][[1]],
                     TC = DOE$TC[x],
                     r_fix = runif(1,min=DOE$r_fix[x][[1]][1],max=DOE$r_fix[x][[1]][2]),
                     r_in = runif(1,min=DOE$r_in[x][[1]][1],max=DOE$r_in[x][[1]][2]),
                     # cor_var = runif(1,min=DOE$cor_var[x][[1]][1],max=DOE$cor_var[x][[1]][2]),
                     # cor_fix = runif(1,min=DOE$cor_fix[x][[1]][1],max=DOE$cor_fix[x][[1]][2]),
                     cv = runif(1,min=DOE$RC_cv[x][[1]][1],max=DOE$RC_cv[x][[1]][2]),
                     max_tries = 2000)


          measures[['SYSTEM']][['RC']]<-list(RC_var_cv =  RC$RC_var$cv,
                                             RC_var_top10 = RC$RC_var$top10,
                                             RC_fix_cv =  RC$RC_fix$cv,
                                             RC_fix_top10 = RC$RC_fix$top10,
                                             cor_fix = RC$cor_fix,
                                             cor_var = RC$cor_var,
                                             cor_fix_var = RC$cor_fix_var,
                                             cor_indirect = RC$cor_indirect,
                                             r_fix = RC$r_fix,
                                             r_in = RC$r_in)
          RC<-list(direct = RC$RC_direct$RC,
                   var=RC$RC_var$RC,
                   fix=RC$RC_fix$RC,
                   P_RD_fix = DOM$P_RD_fix)

    #### 4. Write Output object ####
    EAD<-list(P=P,
              DSM=DSM,
              DMM=DMM,
              DEMAND=DMD,
              RC=RC,
              measures=measures,
              DOE=DOE[x,],
              message = message)

    EAD$ID<-digest(EAD, algo="crc32", serialize =TRUE)

    #### 5. Clean up Workspace ####
    remove(prodMIX)
    remove(DOM)
    gc()
    return(EAD)
  })

  if(file_output){
    save(EAD,file=paste0(Sys.Date(),"_",ids::random_id(n=1,bytes = 4),".RData"))
  }
  return(EAD)
}

#' @title Create EADs using a design of experiments (DOE)
#' @description This function creates EAD realizations using a design of experiments (DOE). The design is specified by \code{DOE}
#' @param DOE A data.frame, where each row contain a parameter set for an EAD realization. The data frame contain the following columns.
#' \describe{
#'   \item{N_FR}{The number of elements in the functional domain (FD)}
#'   \item{N_DD}{The number of elements in the physcial domain (PD)}
#'   \item{N_PrD}{The number of elements in the process domain (PrD)}
#'   \item{N_RD}{The number of elements in the resource domain (RD)}
#'   \item{DENS_FD}{Densitiy of the DSM_FD matrix.}
#'   \item{TOTAL_DEMAND}{Total demand}
#'   \item{DMD_cv}{Coefficient of variation for the demand vector.}
#'   \item{SDC_in}{System design complexity for the DMMs. In the current version this values is applied for all domains.}
#'   \item{ut_DMM}{Boolean. If \code{ut_DMM=TRUE}, the upper triangle DMMs are created representing indirect inter domain dependencies.}
#'   \item{SC_adj_n_in}{Adjusted and normalized structural complexity measure which should be created. In the current version this values is applied for all domains.}
#'   \item{TC}{NTotal costs}
#'   \item{ratio_fixedC}{Proportion of fixed costs on total costs}
#'   \item{RC_cor}{Correlation between variable cost vector and fixed cost vector}
#'   \item{RC_cv}{Coefficient of variation for resource cost distribution}
#'   \item{N_RUN}{Number of runs per unique setting}
#' }
#' @param NUMB_CORES Number of cores to run the simulation. Default: \code{NUMB_CORES=4}
#' @param logfile Optional argument, specifying the name of the logfile which is created in the working directory.
#' @param time_limit Time in seconds after which each EAD creation is interrupted. The interruption produces an error.
#' Therefore runs longer than \code{time_limit} are excluded from the output. The default value is infinity which means that no interruption is done.
#' @return A list of EAD objects
#' @examples
#' set.seed(1234)
#'
#' N_FR=7
#' N_DD=15
#' N_PrD=20
#' N_RD=50
#' DENS_FD=c(0,0.05,0.1,0.2)
#' TOTAL_DEMAND=1000
#' DMD_cv=c(0,0.5,1,3)
#' DMM_PAR=c(0,0.05,0.1,0.2,0.3,0.4,0.5)
#' ut_DMM=c(F)
#' SC_adj_n_in=c(0,0.05,0.01,0.015,0.02)
#' N_RUN=1:1
#' x<-36
#' DOE<-expand.grid(list(N_FR = N_FR,
#'                       N_DD = N_DD,
#'                       N_PrD = N_PrD,
#'                       N_RD = N_RD,
#'                       DENS_FD = DENS_FD,
#'                       TOTAL_DEMAND = TOTAL_DEMAND,
#'                       DMD_cv = DMD_cv,
#'                       DMM_PAR = DMM_PAR,
#'                       ut_DMM = ut_DMM,
#'                       SC_adj_n_in = SC_adj_n_in,
#'                       N_RUN = 1:N_RUN))
#'
#' crt_EAD_MC(DOE,NUMB_CORES=4,logfile="log.txt")
crt_EAD_MC<-function(DOE,
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
    EAD <- odegoparallel::run_MC(cl, X = DOE_list,
                                    FUN = function(DOE, ...) {
                                      res<-list()
                                      res<-with_timeout(crt_EAD(DOE[1,]),timeout = time_limit)
                                      if(res$message=="error"){
                                        res<-list()
                                      }else if(res$message=="success"){
                                        res<-res$res
                                      }
                                      gc()
                                      return(res)
                                    }, packages = c("EAD","odegoparallel",
                                                    "dplyr", "tidyr", "GA",
                                                    "Matrix", "digest", "faux",
                                                    "plyr", "DescTools", "igraph", "R.utils","fGarch"),
                                 errorhandlingNodes = ehNodes)
  }else{
    cl<-EAD::setupMC(NUMB_CORES=NUMB_CORES,logfile=logfile)
    registerDoSNOW(cl)
    on.exit(stopCluster(cl))
    print(cl)
    print(paste0("Time  Limit set to: ",time_limit, " seconds per EAD realization."))
    snow::clusterExport(cl,"time_limit")
    EAD<-par_apply(cl,X=DOE_list,FUN=crt_EAD,errorhandling=ehNodes)
  }
  EAD<-EAD[sapply(EAD,length)>0]
  return(EAD)
}

