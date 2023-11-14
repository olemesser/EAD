
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
                   N_RUN = 1:4 # number of runs
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
                RC_var_i = EAD[[1]][[1]]$RC$var_i,
                RC_var_d = EAD[[1]][[1]]$RC$var_d,
                RC_fix_i = EAD[[1]][[1]]$RC$fix_i,
                RC_fix_d = EAD[[1]][[1]]$RC$fix_d)

## the product costs multiplied by the demand equals the total costs
sum(costs$PC_B*EAD[[1]][[1]]$DEMAND)
#> [1] 1e+06
```

To access the total product costs use the following arguments:

``` r
## total product costs
costs$PC_B
#>  [1]  61.36247 148.82923  83.54024 132.32693 168.25209  93.38689  79.70161
#>  [8] 331.47739 228.53084 127.53286 176.92713 117.93086 111.68369  23.13322
#> [15] 164.35135  70.96446 266.76010 162.88119 182.64816 116.52010 191.38530
#> [22] 276.03504 219.44958  56.56839 210.19170
```

For the indirect benchmark costs use:

``` r
## indirect benchmark product costs
costs$PC_B_indirect
#>  [1]  23.16847  74.43491  38.79758  62.93660  90.78427  51.74357  38.49441
#>  [8] 176.59798 112.92931  68.15747  90.54115  51.55780  62.39494  10.10507
#> [15]  91.51170  39.76813 125.99271  81.63696 102.16307  61.84864 100.88934
#> [22] 153.90664 110.02629  28.38933  97.60338
```

If the product mix (available products) or the demand varies, the costs
are calculated as:

``` r
## lets assume we drop the first ten products 
DMD<-EAD[[1]][[1]]$DEMAND
DMD[1:10]<-0
costs_reduced<-clc_PCB(RES_CONS_PAT = EAD[[1]][[1]]$P$RD,
               DMD = DMD,
               RC_fix_i = EAD[[1]][[1]]$RC$fix_i,
               RC_fix_d = EAD[[1]][[1]]$RC$fix_d,
               RCU_indirect = costs$RCU_indirect,
               RCU_direct = costs$RCU_direct)

PC_B_reduced<-costs_reduced$PC_B

## compare the differences
summary<- data.frame(DMD_f = EAD[[1]][[1]]$DEMAND,
                   DMD_r = DMD,
                   PC_B_f = costs$PC_B,
                   PC_B_r = costs_reduced$PC_B,
                   PC_var_f = costs$PC_B_var_d + costs$PC_B_var_i,
                   PC_var_r = costs_reduced$PC_B_var_d + costs_reduced$PC_B_var_i,
                   PC_fix_f = costs$PC_B_fix_d + costs$PC_B_fix_i,
                   PC_fix_r = costs_reduced$PC_B_fix_d + costs_reduced$PC_B_fix_i)

summary
#>    DMD_f DMD_r    PC_B_f    PC_B_r  PC_var_f  PC_var_r PC_fix_f   PC_fix_r
#> 1     55     0  61.36247   0.00000  59.63782   0.00000 1.724649  0.0000000
#> 2     23     0 148.82923   0.00000 144.73247   0.00000 4.096766  0.0000000
#> 3     47     0  83.54024   0.00000  81.30817   0.00000 2.232076  0.0000000
#> 4   1010     0 132.32693   0.00000 128.22674   0.00000 4.100195  0.0000000
#> 5     79     0 168.25209   0.00000 163.66661   0.00000 4.585474  0.0000000
#> 6    181     0  93.38689   0.00000  90.96573   0.00000 2.421158  0.0000000
#> 7    376     0  79.70161   0.00000  77.42543   0.00000 2.276179  0.0000000
#> 8    161     0 331.47739   0.00000 322.04540   0.00000 9.431991  0.0000000
#> 9     47     0 228.53084   0.00000 222.15790   0.00000 6.372945  0.0000000
#> 10    21     0 127.53286   0.00000 123.53152   0.00000 4.001340  0.0000000
#> 11   220   220 176.92713 179.15512 172.27390 172.27390 4.653234  6.8812256
#> 12    38    38 117.93086 119.61639 114.58042 114.58042 3.350443  5.0359720
#> 13    17    17 111.68369 112.77223 108.72401 108.72401 2.959680  4.0482182
#> 14   110   110  23.13322  23.36096  22.48283  22.48283 0.650385  0.8781312
#> 15   180   180 164.35135 168.04524 159.55465 159.55465 4.796704  8.4905938
#> 16    31    31  70.96446  73.75886  68.58892  68.58892 2.375546  5.1699419
#> 17    16    16 266.76010 270.73119 259.31289 259.31289 7.447209 11.4182988
#> 18   201   201 162.88119 164.95874 158.55483 158.55483 4.326359  6.4039145
#> 19    23    23 182.64816 186.53109 177.31293 177.31293 5.335226  9.2181601
#> 20    71    71 116.52010 117.64734 113.44856 113.44856 3.071543  4.1987831
#> 21  3047  3047 191.38530 193.27775 186.14944 186.14944 5.235859  7.1283125
#> 22    23    23 276.03504 280.81747 268.27866 268.27866 7.756384 12.5388120
#> 23    27    27 219.44958 222.10331 213.49743 213.49743 5.952154  8.6058776
#> 24    62    62  56.56839  57.14456  54.94260  54.94260 1.625795  2.2019632
#> 25     3     3 210.19170 213.58663 204.37029 204.37029 5.821414  9.2163356
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
