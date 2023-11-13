crt_smallEAD <- function(){
  require(tidyr)
  set.seed(1234)
  DOE<-expand_grid(N_FR = 4, # number of functional requirements
                   N_DD = 4, # number of physical domain elements
                   N_PrD = 8, # number of process domain elements
                   N_RD = 10, # number of resource domain elements
                   prop_PROD = 1, # density of the functional product mix
                   method_FD = "random", # method to create the product mix
                   N_PROD = list(c(20,25)), # number of products
                   TOTAL_DEMAND = list(c(100,10000)), # total demand
                   Q_VAR = list(c(0,1.7)), # demand heterogeneity
                   DMM_PAR = expand_grid(FD_PD=list(c(0,0.1)),
                                         PD_PrD=list(c(0,0.05)),
                                         PrD_RD=list(c(0,0.05))), # desired system design complexity
                   uB_DMM = 10,
                   allowZero = T,
                   ut_DMM = F, # if the upper triangle DMMs should be generated too (DMM_FD_PrD,DMM_FD_RD,DMM_PD_RD)
                   DSM_param=expand_grid(PD=list(c(0,0.1,0,1)),
                                         PrD=list(c(0,0,0,1)),
                                         RD=list(c(0,0,0,1))), # first two entries refer to the density of the DSMs and the second pair to the cv if the DSM_method='modular' is used.
                   DSM_method='modular',
                   ub_DSM = 1,
                   TC = 10^6, # total costs
                   r_in = list(c(0,0.9)),
                   r_fix = list(c(0,1)), # proportion of fixed costs on total costs
                   RC_sdlog = list(c(0.1,2)), # skewness of resource costs
                   N_RUN = 1:1 # number of runs
  )
  smallEAD <- crt_EAD(DOE)[[1]]
}
