
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
                   DMM_PAR = list(c(0,0.05)), # desired design complexity
                   DMM_method="SDC", # method for generating the DMM
                   ut_DMM = F, # if the upper triangle DMMs should be generated
                   DSM_param=list(c(0,0.14,0,1)), # the first to entries refer to the density range of the dsm, the second on the weight of the modular structure
                   DSM_method="modular", # method for generating DSMs
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
#>          [,1]
#> 166  93.97619
#> 30   85.16543
#> 153  85.71759
#> 247 138.66876
#> 33   30.33795
#> 45   79.96588
#> 78  100.07538
#> 211 101.53840
#> 62  115.50337
#> 204 152.49754
#> 138 102.17045
#> 120 125.70433
#> 140 128.86535
#> 44  126.68362
#> 213  69.54327
#> 48  133.47603
#> 3    38.78754
#> 2    26.81526
#> 119 117.58072
#> 125 112.68616
#> 122 126.11002
#> 202 125.80265
#> 74   93.28297
#> 150  75.40278
#> 57   92.25810
#> 215 108.33081
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
