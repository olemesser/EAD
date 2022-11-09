#' @name genFirm
#' @title Generates a firm
#' @description Generates a firm by Anand et al. (2017) procedure.
#' @param ProductionEnvironment A structured list. See 'examples' for further details.  \cr
#' @param DISP1 The number of top rescources which should allocate DISP2 percent of total costs.
#' Therefore 1<=DISP1<NUMB_RES.
#' @param DISP2 list containing the elements \code{DISP2$lb} for the lower and \code{DISP2$ub} for the upper bound of DISP2.
#' The function chooses DISP2 from a uniform distribution. If lb=ub the exact value is choosen.
#' DISP2 specifies how mutch percentage of total costs are assigned to the top DISP1 resources? DISP2 0<DISP2<1 \cr
#' @param DENS The density (sparsity) of the RES_CONS_PAT. A list containing the elements \code{DENS$lb} for the lower and \code{DENS$ub} for the upper bound of DENS.
#' @param COR1  A list containing the elements \code{COR1$lb} for the lower and \code{COR1$ub} for the upper bound of COR1.
#' @param exogenious Default=TRUE. If set to true the product mix decision is exognious otherwise endogenious.
#' @examples
#' ProductionEnvironment<-list()
#' ProductionEnvironment[['NUMB_RES']]<-50 #Amount of resources
#' ProductionEnvironment[['NUMB_PRO']]<-50 #Amount of products
#' ProductionEnvironment[['TC']]<-10^6 # Total
#'
#'
#' DISP1<-10
#' DISP2<-list(lb=0.3,ub=0.3)
#' DENS<-list(lb=0.2,ub=0.2)
#' COR1<-list(lb=0.2,ub=0.2)
#' COR2<-list(lb=0.2,ub=0.2)
#'
#'genFirm(ProductionEnvironment,DISP1,DISP2,DENS,COR1,COR2,exogenious=TRUE)
#' @export
NULL
genFirm<-function(ProductionEnvironment,DISP1,DISP2,DENS,COR1,COR2,exogenious=TRUE){
#### Pre Processing ####


  ## Some Checks before ##

  # RCP>=CO
  # stopifnot(ProductionEnvironment[['NUMB_RES']]>=ProductionEnvironment[['NUMB_PRO']])

  # Production Environment
  stopifnot(!is.null(ProductionEnvironment[['NUMB_RES']]),!is.null(ProductionEnvironment[['NUMB_PRO']]),!is.null(ProductionEnvironment[['TC']]))

  # DISP2
  stopifnot(DISP2$lb>0 & DISP2$ub<1 & DISP2$lb<=DISP2$ub)

  # DISP1
  stopifnot(DISP1>=1,DISP1<ProductionEnvironment[['NUMB_RES']])

  ## End of Checks ##






  ## Initialize variables ##
  NUMB_RES<-ProductionEnvironment[['NUMB_RES']]
  NUMB_PRO<-ProductionEnvironment[['NUMB_PRO']]
  TC<-ProductionEnvironment[['TC']]
  CHECK<-list()




  #### A2.1.2-1. Chose random Values ####

    DISP2<-runif(1,min=DISP2$lb,max=DISP2$ub) # random value for G/DISP2
    DENS<-runif(1,min=DENS$lb,max=DENS$ub) # random value for D/DNS
    COR1<-runif(1,min=COR1$lb,max=COR1$ub)
    COR2<-runif(1,min=COR2$lb,max=COR2$ub)



  #### A2.1.2-2. Generate true product margins (MAR) and a decision vector (DECT0) ####
    # ratio of selling price/cost
    ## TODO QT for DECT0
    # MAR<-runif(NUMB_PRO,min = MARLB,max=MARUB)
    # DECT0<-ifelse(MAR<1,0,1)
    DECT0<-rep(1,NUMB_PRO)

  #### A2.1.2-3. Maximum Production Quantities (MXQ) ####

    MXQ<-sample(10:40,NUMB_PRO,replace = T)
    DEMAND_BASE<-DECT0*MXQ
    ProductionEnvironment[['DEMAND']]<-DEMAND_BASE # equivalent to QT

  #### A2.1.2-4. Generate Rescource Consumption pattern matrix (RES_CONS_PAT) ####
    ProductionEnvironment[['DISP1']]<-DISP1
    ProductionEnvironment[['DISP2']]<-DISP2

    pattern<-genRES_CONS_PAT(ProductionEnvironment,DENS,COR1,COR2,RES_METHOD="ANAND") # generate res_cons_pat
    RES_CONS_PAT<-pattern$RES_CONS_PAT

  #### A2.1.2-5. I Generate Vector of resource consumption (TRU, MAXRU) ####

    # Total Resource Units -> Amount Needed to produce mix QT = ProductionEnvironment[['DEMAND']]
    MAXRU<-t(RES_CONS_PAT) %*% MXQ
    TRU<-t(RES_CONS_PAT) %*% ProductionEnvironment[['DEMAND']]


  #### A2.1.2-5. II Generate Vector of Total Rescource Costs (RCC) ####
    RCC<-genRC(ProductionEnvironment$NUMB_RES,DISP1 = DISP1,DISP2 = DISP2,TC=TC)


  #### A2.1.2-5. III Compute Vector of unit resource costs (RCU) ####
    # Unit Resource Costs (RCU)
    RCU<-RCC$RC/MAXRU



  #### A2.1.2-5. IV Reclaculate RCC ####
    if(exogenious==FALSE){
      stop("Not implemented yet")
      # Recompute RCC
      # only needed when exogenious
    }

  #### A2.1.2-5. V Compute Benchmark Costs (PC_B) ####
    CostSystem<-list()
    ## A-2.1.2.5
    CostSystem[['PC_B']] <- RES_CONS_PAT%*%RCU #BENCHMARK PRODUCT COSTS (TOTAL)


  #### A2.1.2-5. VI Compute vector of selling prices (SP) and profit ####
    ## if exogenious this step can be omitted ##

  #### A2.1.2-5. VII Compute Rank vector (RANK) and matrix of percentage resource consumption ####
    ## rank is computed in the cost system function ##

    ## Multiply RES_CONS_PAT with quantities QT
    # RES_CONS_PATp<-RES_CONS_PAT*DEMAND_BASE #RES_CONS_PATp is the vector v (see A-2.1.2.7) or Flowchart 5.9 (b)
    RES_CONS_PATp<-matrix(0,nrow = NROW(RES_CONS_PAT),ncol = NCOL(RES_CONS_PAT))

    for(i in 1:NCOL(RES_CONS_PATp)){
      RES_CONS_PATp[,i]<- RES_CONS_PAT[,i]*DEMAND_BASE/TRU[i]
    }

  #### A2.1.2-5. VIII Compute correlation in resource consumption ####
    ProductionEnvironment[['PEARSONCORR']] = cor(RES_CONS_PATp)


  #### Data Logging ####
    ## Input Parameter Check ##
    CHECK[['INPUT']]<-list(NUMB_RES=NUMB_RES,NUMB_PRO=NUMB_PRO,TC=TC,
                           DISP1=DISP1,DISP2=DISP2,COR1=COR1,COR2=COR2,DENS=DENS)

    ## COR Parameter Check ##


    ## CHECK COR
    # Averange range between lowest and highest consumption.
    CORAP_pre=max(RES_CONS_PATp)-min(RES_CONS_PATp)

    CHECK[['RES_CONS_PAT']]<-c(pattern$CHECK,CORAP=mean(CORAP_pre)*100)


    ## RCP Check ##
    RC_sort = sort(RCC$RC,decreasing = TRUE)
    endRes<-floor(NUMB_RES*.2)
    RC_20p = sum(RC_sort[1:endRes])

    CHECK[['RCP']]<-list(perc_largestPool=RCC$CHECK$cost_largestRCP,
                         perc_top10Pool=RCC$CHECK$cost_topTEN,
                         RC_20p=RC_20p/TC*100)




    # Computing heterogeneity
    # CHECK<- measure_heterogeneity(RES_CONS_PATp,ProductionEnvironment,CHECK)
    ProductionEnvironment[['RES_CONS_PATp']] = RES_CONS_PATp
    ProductionEnvironment[['RES_CONS_PAT']] = RES_CONS_PAT
    ProductionEnvironment[['CostSystem']] = list(RCC=RCC$RC,PC_B=CostSystem[['PC_B']],RCU=RCU,TRU=TRU)



  out<-list(ProductionEnvironment=ProductionEnvironment,
            CHECK=CHECK)

  return(out)
}
