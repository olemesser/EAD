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
  #                  DNS = list(c(0.07,0.5)), #density within the DSM_FD matrix. Creates product mixes where not all products are included
  #                  N_PROD = list(c(50,100)), # number of products
  #                  method_FD = "DNS", # method for generating the product mix
  #                  TOTAL_DEMAND = list(c(100,12300)), # total demand
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
  #                  RC_sdlog = list(c(0,1.5)), # coefficient of variation for resource cost distribution
  #                  N_RUN = 1:4 # number of runs
  # )
  # DOE<-DOE[1:4,]
  #### END Input Testing ####

  if(min(unlist(DOE$TOTAL_DEMAND)) < min(unlist(DOE$N_PROD))) stop("TOTAL_DEMAND < N_PROD. Please increase TOTAL_DEMAND or reduce N_PRDO!")

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
    dmd_candidate <- ifelse(length(DOE$TOTAL_DEMAND[x][[1]])>1,sample(DOE$TOTAL_DEMAND[[x]][1]:DOE$TOTAL_DEMAND[[x]][2],1),DOE$TOTAL_DEMAND[[x]][1])
    DOE$N_PROD[x] <- ifelse(length(DOE$N_PROD[x][[1]])>1,sample(DOE$N_PROD[[x]][1]:DOE$N_PROD[[x]][2],1),DOE$N_PROD[[x]][1])

    while(dmd_candidate[[1]] < DOE$N_PROD[x][[1]]){
      dmd_candidate <- ifelse(length(DOE$TOTAL_DEMAND[x][[1]])>1,sample(DOE$TOTAL_DEMAND[[x]][1]:DOE$TOTAL_DEMAND[[x]][2],1),DOE$TOTAL_DEMAND[[x]][1])
    }
    DOE$TOTAL_DEMAND[x] <- dmd_candidate
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
                        N_PROD = DOE$N_PROD[x][[1]],
                        method = DOE$method_FD[x])
      P[["FD"]]<-as.matrix(prodMIX$P_FD_const)


      ### 1.2 Demand Generation ###
      DMD<-crt_DEMAND(n_PROD = NROW(P[["FD"]]),
                 TOTAL_DEMAND = DOE$TOTAL_DEMAND[x][[1]],
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
                     sdlog = runif(1,min=DOE$RC_sdlog[x][[1]][1],max=DOE$RC_sdlog[x][[1]][2]))


          measures[['SYSTEM']][['RC']]<-list(RC_vari_top10 = RC$RC_var$RC_i$top10,
                                             RC_vard_top10 = RC$RC_var$RC_d$top10,
                                             RC_fixi_top10 = RC$RC_fix$RC_i$top10,
                                             RC_fixd_top10 = RC$RC_fix$RC_d$top10,
                                             r_fix = RC$r_fix,
                                             r_in = RC$r_in,
                                             TC = sum(RC$RC_var$RC_d$RC, RC$RC_fix$RC_d$RC,RC$RC_var$RC_i$RC,RC$RC_fix$RC_i$RC))
          RC<-list(var_d = RC$RC_var$RC_d$RC,
                   fix_d = RC$RC_fix$RC_d$RC,
                   var_i=RC$RC_var$RC_i$RC,
                   fix_i=RC$RC_fix$RC_i$RC)

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
#'   \item{DNS}{Density of the product mix \code{P_FD}.}
#'   \item{N_PROD}{Number of products (rows) in \code{P_FD}.}
#'   \item{method_FD}{Method for generating the product mix. For available methods see: \link[EAD]{create_ProductMix}.}
#'   \item{TOTAL_DEMAND}{Total demand}
#'   \item{Q_VAR}{Demand heterogeneity.}
#'   \item{DMM_PAR}{Parameter for creating the DMM matrices. Named list for each of the three DMMs.
#'   The input is list object, where each list element is a 2d-vector specifying the lower and upper bound of the desired standardized system design complextiy.
#'   For further details see: \link[EAD]{crt_DMM}.}
#'   \item{uB_DMM}{The upper bound for entries in the DMM matrices. Passed to: \link[EAD]{crt_DMM}.}
#'   \item{allowZero}{Specifies whether zero columns are allowed in the DMM matrices. Passed to: \link[EAD]{crt_DMM}.}
#'   \item{ut_DMM}{Boolean. If \code{ut_DMM=TRUE}, the upper triangle DMMs (DMM_FD_PrD,DMM_FD_RD,DMM_PD_RD) are created representing indirect inter domain dependencies.}
#'   \item{DSM_param}{The first two entries refer to the density of the DSMs and the second pair to the cv if the DSM_method='modular' is used. For details see: \link[EAD]{crt_DSM}.}
#'   \item{DSM_method}{The method for DSM creation. For details see: \link[EAD]{crt_DSM}.}
#'   \item{ub_DSM}{Upper bounds for entries in the DSMs.}
#'   \item{TC}{Total costs}
#'   \item{r_in}{Proportion of indirect costs on total costs}
#'   \item{r_fix}{Proportion of fixed costs on total indirect costs}
#'   \item{RC_cor}{Correlation between variable cost vector and fixed cost vector}
#'   \item{RC_cv}{Coefficient of variation for resource cost distribution}
#'   \item{N_RUN}{Number of runs per unique setting}
#' }
#' @param NUMB_CORES Number of cores to run the simulation. Default: \code{NUMB_CORES=4}
#' @param logfile Optional argument, specifying the name of the logfile which is created in the working directory.
#' @param time_limit Time in seconds after which each EAD creation is interrupted. This is necessary since there are certain input combinations which result in long computational time.
#' Therefore runs longer than \code{time_limit} are excluded from the output. The default value is infinity which means that no interruption is done.
#' @param cl An optional cluster object created by \link[parallel]{makeCluster}. This allows to distribute the computation across different nodes within a network.
#' For further details see the documentation of the 'parallel' package. By default the function \link[EAD]{setupMC} is called.
#' @param ehNodes A character specifying the parallel node behavior on errors. This argument is passed to the '.errorhandling' of \link[foreach]{foreach}. The default value is \code{'remove'}.
#' @return A list of EAD objects
#' @examples
#' set.seed(1234)
#'
#' library(tidyr)
#' DOE<-expand_grid(N_FR = list(c(15)),
#'                  N_DD = list(c(18,26)),
#'                  N_PrD = list(c(36,52)),
#'                  N_RD = list(c(72,104)),
#'                  DNS = list(c(0.07,0.5)),
#'                  N_PROD = list(c(50)),
#'                  method_FD = "DNS",
#'                  TOTAL_DEMAND = list(c(1000,1000)),
#'                  Q_VAR = list(c(0,3)),
#'                  DMM_PAR = expand_grid(FD_PD = list(c(0,0.08)),
#'                                        PD_PrD = list(c(0,0.1)),
#'                                        PrD_RD = list(c(0,0.1))),
#'                  uB_DMM = 1,
#'                  allowZero = F,
#'                  ut_DMM = F,
#'                  DSM_param=expand_grid(PD = list(c(0,0.05,0,1)),
#'                                        PrD = list(c(0,0.05,0,1)),
#'                                        RD = list(c(0,0.05,0,1))), # first two entries refer to the density of the DSMs and the second pair to the cv if the DSM_method='modular' is used.
#'                  DSM_method='modular',
#'                  ub_DSM = 1,
#'                  TC = 10^6,
#'                  r_in = list(c(0,0.9)),
#'                  r_fix = list(c(0,1)),
#'                  RC_sdlog = list(c(0,1.5)),
#'                  N_RUN = 1:4 # number of runs
#' )
#'
#' EAD <- crt_EAD_MC(DOE,NUMB_CORES=4,logfile="log.txt")
#' EAD[[1]]
#'
crt_EAD_MC<-function(DOE,
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


  parallel::clusterExport(cl, envir = environment(),c("time_limit"))
  print(cl)
  print(paste0("Time  Limit set to: ",time_limit, " seconds per EAD realization."))


  EAD <- run_MC(cl, X = DOE_list,
                FUN = function(DOE, ...) {
                  res<-list()
                  crt_EAD(DOE[1,])
                  res<-with_timeout(crt_EAD(DOE[1,]),timeout = time_limit)
                  if(res$message=="error"){
                    res<-list()
                  }else if(res$message=="success"){
                    res<-res$res
                  }
                  gc()
                  return(res)
                }, packages = c("EAD",
                                "dplyr", "tidyr",
                                "Matrix", "digest", "faux",
                                "plyr", "DescTools", "igraph", "R.utils","fGarch"),
             errorhandlingNodes = ehNodes)

  EAD<-EAD[sapply(EAD,length)>0]
  return(EAD)
}

