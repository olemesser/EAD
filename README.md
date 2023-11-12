
<!-- README.md is generated from README.Rmd. Please edit that file -->

# The Extended Axiomatic Design (EAD)

<!-- badges: start -->
<!-- badges: end -->

The numerical EAD framework generates firm environments according to the
Extended Axiomatic Design (EAD) introduced by Mertens (2020),
Meßerschmidt et al. (2020) and Schmidt et al. (2021). The simulation
model is explained in Meßerschmidt (forthcoming). The framework adapts
some function from the cost system design framework by Anand et
al. (2019). The original code in C# is available under
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
```

Third, define a design of experiment (DoE). For each

``` r
 DOE<-expand_grid(N_FR = 8, # number of functional requirements
                   N_DD = 14, # number of physical domain elements
                   N_PrD = 30, # number of process domain elements
                   N_RD = 60, # number of resource domain elements
                   DNS = list(c(0.05,0.5)), # density of the functional product mix
                   method_FD = "DNS", # method to create the product mix
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
                  NUMB_CORES=4,
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
costs<- clc_PCB(RES_CONS_PAT = EAD[[1]][[1]]$P$RD,
                DMD = EAD[[1]][[1]]$DEMAND,
                RC_var = EAD[[1]][[1]]$RC$var_i + EAD[[1]][[1]]$RC$var_d,
                RC_fix = EAD[[1]][[1]]$RC$fix_i + EAD[[1]][[1]]$RC$fix_d)

## the product costs multiplied by the demand equals the total costs
sum(costs$PC_B*EAD[[1]][[1]]$DEMAND)
#> [1] 1e+06
```

To access the total product costs use the following arguments:

``` r
## total product costs
costs$PC_B
#>  [1] 446.70642 681.04569 225.31326 616.12150 531.15134 318.79079 359.19526
#>  [8] 140.18669 253.81275 436.46515 489.82912 445.58026 226.04489 124.85598
#> [15] 531.87617 542.11744 466.95197 171.22446 147.42956  52.01854 584.50852
#> [22] 614.68510 392.08401
```

For the indirect benchmark costs use:

``` r
## indirect benchmark product costs
costs$PC_B_indirect
#>  [1] 446.70642 681.04569 225.31326 616.12150 531.15134 318.79079 359.19526
#>  [8] 140.18669 253.81275 436.46515 489.82912 445.58026 226.04489 124.85598
#> [15] 531.87617 542.11744 466.95197 171.22446 147.42956  52.01854 584.50852
#> [22] 614.68510 392.08401
```

If the product mix (available products) or the demand varies, the costs
are calculated as:

``` r
## lets assume we drop the first ten products 
DMD<-EAD[[1]][[1]]$DEMAND
DMD[1:10]<-0
costs_reduced<-clc_PCB(RES_CONS_PAT = EAD[[1]][[1]]$P$RD,
               DMD = DMD,
               RC_fix = EAD[[1]][[1]]$RC$fix_i +  EAD[[1]][[1]]$RC$fix_d,
               RCU = costs$RCU)

PC_B_reduced<-costs_reduced$PC_B

## compare the differences
summary<- data.frame(DMD_f = EAD[[1]][[1]]$DEMAND,
                   DMD_r = DMD,
                   PC_B_f = costs$PC_B,
                   PC_B_r = costs_reduced$PC_B,
                   PC_var_f = costs$PC_B_var,
                   PC_var_r = costs_reduced$PC_B_var,
                   PC_fix_f = costs$PC_B_fixed,
                   PC_fix_r = costs_reduced$PC_B_fixed)

summary
#>    DMD_f DMD_r    PC_B_f    PC_B_r  PC_var_f  PC_var_r  PC_fix_f  PC_fix_r
#> 1     20     0 446.70642   0.00000 408.40110   0.00000 38.305316   0.00000
#> 2    103     0 681.04569   0.00000 626.18423   0.00000 54.861458   0.00000
#> 3     47     0 225.31326   0.00000 207.29953   0.00000 18.013733   0.00000
#> 4    191     0 616.12150   0.00000 568.28894   0.00000 47.832557   0.00000
#> 5    422     0 531.15134   0.00000 490.21444   0.00000 40.936908   0.00000
#> 6    293     0 318.79079   0.00000 293.17806   0.00000 25.612731   0.00000
#> 7      7     0 359.19526   0.00000 329.56122   0.00000 29.634032   0.00000
#> 8    271     0 140.18669   0.00000 130.08119   0.00000 10.105501   0.00000
#> 9     12     0 253.81275   0.00000 231.38385   0.00000 22.428899   0.00000
#> 10   222     0 436.46515   0.00000 400.64155   0.00000 35.823600   0.00000
#> 11    13    13 489.82912 576.14027 450.49423 450.49423 39.334887 125.64604
#> 12   138   138 445.58026 516.32090 406.68764 406.68764 38.892623 109.63326
#> 13    17    17 226.04489 260.74770 208.54303 208.54303 17.501854  52.20467
#> 14    31    31 124.85598 146.27844 111.77809 111.77809 13.077891  34.50035
#> 15    27    27 531.87617 614.22797 488.25157 488.25157 43.624600 125.97640
#> 16     4     4 542.11744 628.49091 496.01112 496.01112 46.106316 132.47978
#> 17     8     8 466.95197 536.84994 430.35627 430.35627 36.595699 106.49367
#> 18   460   460 171.22446 193.76936 159.40971 159.40971 11.814754  34.35965
#> 19    16    16 147.42956 170.46955 135.93471 135.93471 11.494847  34.53484
#> 20   117   117  52.01854  59.26326  48.32469  48.32469  3.693847  10.93857
#> 21    58    58 584.50852 685.43503 536.86075 536.86075 47.647765 148.57428
#> 22    99    99 614.68510 722.41871 562.27233 562.27233 52.412778 160.14639
#> 23   145   145 392.08401 453.97818 358.64034 358.64034 33.443670  95.33785
```

# Full Documentation

For the full documentation use the included vignettes as:

``` r
## list available vignettes
utils::vignette(package = "EAD")

## opens main documentation
utils::vignette("createEAD",package ="EAD")

## opens the setup for product costing system experiments
utils::vignette("pcs-experiment",package ="EAD")

## opens the setup forthe complexity cost experiment
utils::vignette("cc-experiment",package ="EAD")
```

This frameworks further comes with some example data sets reported in
literature. To get on overview of available data sets call:

``` r
data(package = "EAD")
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

Meßerschmidt, O., Gumpinger, T., Meyer, M. & Mertens, K. G. Reviewing
Complexity Costs – What Practice Needs and What Research Contributes
Proceedings of the Design Society: DESIGN Conference, 2020. Cambridge
University Press, 647-656. <https://doi.org/10.1017/dsd.2020.152>

Meßerschmidt, O. (forthcoming). Variety and Costs - The Effects of
Product Variety on Costs and Costing Systems.

Meyer, M., Meßerschmidt, O., & Mertens, K. G. (2019). *How much does
variety-induced complexity actually cost? Linking axiomatic design with
cost modelling*. In M. Schröder & K. Wegner (Eds.), *Logistik im Wandel
der Zeit – Von der Produktionssteuerung zu vernetzten Supply Chains*
(pp. 813–827). Wiesbaden: Springer Fachmedien Wiesbaden.
<https://doi.org/10.1007/978-3-658-25412-4_39>

Schmidt M, Mertens KG, Meyer M (2023) Cost hierarchies and the pattern
of product cost cross-subsidization: Extending a computational model of
costing system design. PLOS ONE 18(9): e0290370.
<https://doi.org/10.1371/journal.pone.0290370>
