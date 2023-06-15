
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
                   DNS = list(c(0.05,0.5)), # density of the functional product mix
                   method_FD = "DNS", # method to create the product mix
                   N_PROD=20, # number of products
                   TOTAL_DEMAND = 10000, # total demand
                   Q_VAR = list(c(0,3)), # demand heterogeneity
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
TC <- sum(EAD[[1]][[1]]$RC$var + EAD[[1]][[1]]$RC$direct + EAD[[1]][[1]]$RC$fix)
costs<- clc_PCB(RES_CONS_PAT = EAD[[1]][[1]]$P$RD,
                DMD = EAD[[1]][[1]]$DEMAND,
                RC_direct = EAD[[1]][[1]]$RC$direct,
                RC_var = EAD[[1]][[1]]$RC$var,
                RC_fix = EAD[[1]][[1]]$RC$fix)

## the product costs multiplied by the demand equals the total costs
sum(costs$PC_B*EAD[[1]][[1]]$DEMAND)==TC
#> [1] FALSE
```

To access the total product costs use the following arguments:

``` r
## total product costs
costs$PC_B
#>  [1] 170.95082  44.52441 262.80348 128.15053  65.33652  46.41945 196.48416
#>  [8] 237.29883 121.73387 200.09699  25.60735 131.75275 104.63150 110.94065
#> [15] 132.83825 176.38468 123.79977 160.89892  92.75856  78.19226
```

For the indirect benchmark costs use:

``` r
## indirect benchmark product costs
costs$PC_B_indirect
#>  [1] 35.531399  8.864546 54.460040 27.476248 13.223336  9.470844 41.329289
#>  [8] 49.428238 25.993841 41.612147  5.112054 26.836826 22.400648 22.478036
#> [15] 26.824494 37.106391 25.181734 32.917573 18.494689 17.036711
```

If the product mix (available products) or the demand varies, the costs
are calculated as:

``` r
## lets assume we drop the first ten products 
DMD<-EAD[[1]][[1]]$DEMAND
DMD[1:10]<-0
costs_reduced<-clc_PCB(RES_CONS_PAT = EAD[[1]][[1]]$P$RD,
               DMD = DMD,
               RC_direct = EAD[[1]][[1]]$RC$direct,
               RC_fix = EAD[[1]][[1]]$RC$fix,
               RCU=costs$RCU)

PC_B_reduced<-costs_reduced$PC_B

## compare the differences
costs<- data.frame(DMD_f = EAD[[1]][[1]]$DEMAND,
                   DMD_r = DMD,
                   PC_B_f = costs$PC_B,
                   PC_B_r = costs_reduced$PC_B,
                   PC_var_f = costs$PC_B_var,
                   PC_var_r = costs_reduced$PC_B_var,
                   PC_fix_f = costs$PC_B_fixed,
                   PC_fix_r = costs_reduced$PC_B_fixed)

head(costs)
```

# Full Documentation

For the full documentation use the included vignettes as:

``` r
## list available vignettes
utils::vignette(package = "EAD")

## open main documentation
utils::vignette("documentation",package ="EAD")
```

# Acknowledgement

Special thanks to [Prof. Matthias
Meyer](https://www.researchgate.net/profile/Matthias-Meyer-12 "Visit on Research Gate"),
[Kai G.
Mertens](https://www.researchgate.net/profile/Kai-G-Mertens "Visit on Research Gate")
and [Mark
Schmidt](https://www.researchgate.net/profile/Mark-Schmidt-19 "Visit on Research Gate")
for detailed discussions during the development of this framework.

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
