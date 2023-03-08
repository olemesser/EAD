
<!-- README.md is generated from README.Rmd. Please edit that file -->

# The Extended Axiomatic Design (EAD)

<!-- badges: start -->
<!-- badges: end -->

The numerical EAD framework generates firm environments according to the
Extended Axiomatic Design (EAD) introduced by Mertens (2020). The
simulation model is explained in Meßerschmidt (forthcoming). The
framework adapts some function from the cost system design framework by
Anand et al. (2019). The original code in C# is available under
<https://github.com/vanand74/CostSystemSim>.

## Installation

You can install the development version of EAD via any github install
function. For example:

``` r
library(devtools)
install_github("olemesser/EAD",build_vignettes = TRUE)
```

## Create the EAD

First, load all relevant libraries:

``` r
#### Load Libraries ####
library(EAD)
library(dplyr)
library(tidyr)
library(parallel)
library(doSNOW)
library(foreach)
library(R.utils)
library(ids)
```

Second, define global parameter such as the output path, seed and a
timeout. Important note: the `time_limit` variable must be defined as a
global variable via double arrows as shown below. In this example, the
EAD generation process is interrupted after 20 seconds. The individual
time limit may vary from hardware configuration and input values defined
by the design of experiments

``` r
set.seed(1234) # set random seed to reproduce the results
time_limit <<- 50 # global time limit in seconds
# setwd()
```

Third, define a design of experiment (DoE). For each

``` r
 DOE<-expand_grid(N_FR = 8, # number of functional requirements
                   N_DD = 14, # number of physical domain elements
                   N_PrD = 30, # number of process domain elements
                   N_RD = 60, # number of resource domain elements
                   PARAM_FD = seq(0.1,1,0.1), # percentage of products included
                   method_FD = "random",
                   TOTAL_DEMAND = 10000, # total demand
                   DMD_cv = list(c(0,3)), # demand skewness
                   DMM_PAR = expand_grid(FD_PD=list(c(0,0.1)),
                                PD_PrD=list(c(0,0.05)),
                                PrD_RD=list(c(0,0.05))), # desired system design complexity
                   DMM_method="SDC", # method for generating the DMM
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
                   cor_var = list(c(-1,1)), # correlation between indirect variable cost vector and direct cost vector
                   cor_fix = list(c(-1,1)), # correlation between indirect fixed cost vector and direct cost vector
                   RC_cv = list(c(0,1)), # skewness of resource costs
                   N_RUN = 1:100 # number of runs
                    )
```

Now you are ready to run the model. In this example, the EAD is created
only for the first four rows of `DOE`. The code is optimized to
multi-core processing. Where `NUMB_CORES` defines the number of cores
used.

``` r
## run EAD generation procedure
  EAD<-crt_EAD_MC(DOE[1:4,],
                  NUMB_CORES=c(4),
                  logfile="")
```

Finally, write the created EAD environments to a file to conduct further
experiments

``` r
  ## write output to file
  save(EAD,file=paste0(Sys.Date(),"_EAD.RData"))
```

## Calculate (true) Benchmark Costs

The following chunk calculates the product benchmark costs for the given
design.

``` r
# variable indirect costs and the direct costs which are always variable are summed up
RC_var <- EAD[[1]][[1]]$RC$var + EAD[[1]][[1]]$RC$direct
costs<- clc_PCB(RES_CONS_PAT = EAD[[1]][[1]]$P$RD,
               DMD = EAD[[1]][[1]]$DEMAND,
               RC_var = RC_var,
               RC_fix = EAD[[1]][[1]]$RC$fix)

## the product costs multiplied by the demand equals the total costs
sum(costs$PC_B*EAD[[1]][[1]]$DEMAND)==DOE$TC[1]
#> [1] FALSE

PC_B_full<-costs$PC_B
PC_B_full
#>           [,1]
#> 7    52.095851
#> 255 178.655879
#> 104 152.624248
#> 59  111.921369
#> 196 151.288262
#> 151  94.970152
#> 234 172.499959
#> 201  93.291248
#> 153  79.037897
#> 97   56.864679
#> 197  69.403810
#> 66   64.983404
#> 249 133.233277
#> 55   89.368580
#> 58  117.100042
#> 251 171.718321
#> 226 141.674962
#> 21   13.362545
#> 118 109.193691
#> 5     8.965509
#> 247 151.834831
#> 131  86.304858
#> 192 203.669210
#> 113  58.592415
#> 17    4.397035
#> 241  99.738980
```

If the product mix (available products) or the demand varies, the costs
are calculated as:

``` r
## lets assume we drop the first ten products 
DMD<-EAD[[1]][[1]]$DEMAND
DMD[1:10]<-0
costs<-clc_PCB(RES_CONS_PAT = EAD[[1]][[1]]$P$RD,
                DMD = DMD,
                RC_fix = EAD[[1]][[1]]$RC$fix,
                RCU=costs$RCU)

PC_B_reduced<-costs$PC_B

## compare the differences
cbind(PC_B_full,PC_B_reduced)
```

# Full Documentation

For the full documentation use the included vignettes as:

``` r
## list available vignettes
utils::vignette(package = "EAD")

## open main documentation
utils::vignette("documentation",package ="EAD")
```

# References

Anand, V., Balakrishnan, R., & Labro, E. (2019). A Framework for
Conducting Numerical Experiments on Cost System Design. *Journal of
Management Accounting Research*, *31*(1), 41–61.
<https://doi.org/10.2308/jmar-52057>

Mertens, K. G. (2020). *Measure and manage your product costs right –
development and use of an extended axiomatic design for cost modeling.*
TUHH Universitätsbibliothek. <https://doi.org/10.15480/882.2888>

Meyer, M., Meßerschmidt, O., & Mertens, K. G. (2019). *How much does
variety-induced complexity actually cost? Linking axiomatic design with
cost modelling*. In M. Schröder & K. Wegner (Eds.), *Logistik im Wandel
der Zeit – Von der Produktionssteuerung zu vernetzten Supply Chains*
(pp. 813–827). Wiesbaden: Springer Fachmedien Wiesbaden.
<https://doi.org/10.1007/978-3-658-25412-4_39>

Meßerschmidt, O. (forthcoming). *Variety and Costs - The Effects of
Product Variety on Costs and Costing Systems.*
