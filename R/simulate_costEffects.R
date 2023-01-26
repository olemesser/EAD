simulate_costEffects<-function(){
  require(EAD)
  suppressWarnings(require(digest))
  suppressWarnings(require(tidyr))

  #### Input for Testing ####
  DOE<-expand_grid(N_FR = list(c(7)), # number of functional requirements
                   N_DD = list(c(20)), # number of physical domain elements
                   N_PrD = list(c(40)), # number of process domain elements
                   N_RD = list(c(80)), # number of resource domain elements
                   PARAM_FD =1, # generate free combination
                   method_FD = "random",
                   TOTAL_DEMAND = 10000, # total demand
                   DMD_cv = list(c(0,3)), # coefficient of variation for demand distribution
                   DMM_PAR = list(c(0,0.115)), # desired design complexity
                   DMM_method="SDC", # method for generating the DMM
                   ut_DMM = F, # if the upper triangle DMMs should be generated too (DMM_FD_PrD,DMM_FD_RD,DMM_PD_RD)
                   DSM_param=list(c(0,0.14,0,1)), # first two entries refer to the density of the DSMs and the second pair to the cv if the DSM_method='modular' is used.
                   DSM_method='modular',
                   TC = 10^6, # total costs
                   ratio_fixedC = list(c(0.5,0.8)), # proportion of fixed costs on total costs
                   RC_cor = list(c(0,1)), # correlation between variable cost vector and fixed cost vector
                   RC_cv = list(c(0,1.5)), # coefficient of variation for resource cost distribution
                   N_RUN = 1:1 # number of runs
  )
  x<-1
  #### END Input for Testing ####

  #### For each DOE ####
  lapply(1:NROW(DOE),function(x){
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

      DMD_temp<-EAD[[1]]$DEMAND
      DMD_temp[setdiff(1:length(DMD_temp),idx)]<-0
      costs<-clc_PCB(RES_CONS_PAT = EAD[[1]]$P$RD,
              DMD = DMD_temp,
              RCU = RCU_initial,
              RC_fix = EAD[[1]]$RC$fix)

      TC<-sum(costs$PC_B*DMD_temp)

      #### complexity costs due to additional process variety
      PC_setupChange<-processVariety(EAD,
                     includedProd=idx,
                     prop_setupChange=c(0,0.05),
                     samples=1)
      TC_setup<-colSums(PC_setupChange * DMD_temp[DMD_temp>0])
      out<-list(NPV=qty,
                DMD_total = sum(DMD_temp),
                TC=TC,
                TC_setup=TC_setup,
                costs_perProd = TC/sum(DMD_temp))
      return(out)
    })

    df<-lapply(res,function(t){
      data.frame(NPV=t$NPV,
                 DMD_total = t$DMD_total,
                 TC=t$TC,
                 TC_setup=t$TC_setup,
                 costs_perProd = t$costs_perProd)
    })
    df<-data.table::rbindlist(df) %>%
      mutate(TC_nc=TC,
             TC=TC+TC_setup,
             costs_perProd=TC/sum(DMD_total))


    plot(df$NPV,df$TC)
    plot(df$NPV,df$costs_perProd)



  })





}
