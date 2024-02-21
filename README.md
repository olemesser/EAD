
<!-- README.md is generated from README.Rmd. Please edit that file -->

# The Extended Axiomatic Design (EAD)

<!-- badges: start -->
<!-- badges: end -->

The numerical EAD framework generates firm environments according to the
Extended Axiomatic Design (EAD) introduced by Mertens (2020),
Meßerschmidt et al. (2020) and Schmidt et al. (2021). The framework
adapts some function from the cost system design framework by Anand et
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
                   N_DD = 8, # number of physical domain elements
                   N_PrD = 8, # number of process domain elements
                   N_RD = 10, # number of resource domain elements
                   DNS = list(c(0.05,0.5)), # density of the functional product mix
                   method_FD = "DNS", # method to create the product mix
                   N_PROD = list(c(20,25)), # number of products
                   TOTAL_DEMAND = list(c(100,10000)), # total demand
                   Q_VAR = list(c(0,1.7)), # demand heterogeneity
                   DMM_PAR = expand_grid(FD_PD=list(c(0,0.1)),
                                PD_PrD=list(c(0,0.05)),
                                PrD_RD=list(c(0,0.05))), # desired system design complexity
                   uB_DMM = 10,
                   allowZero = F,
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
#>  [1]  188.63854  418.85422  554.71796  190.17190  373.57268  457.06087
#>  [7] 1122.42931  298.49092  180.03006  612.44256  137.08649  645.97640
#> [13]   93.33452  738.43345  190.49178  904.52323 1108.61344  189.61919
#> [19]  304.29314  682.69417  457.10036
```

For the indirect benchmark costs use:

``` r
## indirect benchmark product costs
costs$PC_B_indirect
#>  [1]  7.064298 15.814268 13.860748  6.535262 10.955852 22.339040 19.923788
#>  [8] 15.339007  8.341295 20.216062  9.088077 24.831001  8.398908 21.611419
#> [15] 10.942364 16.149353 22.780298 10.923000 15.523857 12.774971 13.561942
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
#>    DMD_f DMD_r     PC_B_f    PC_B_r  PC_var_f  PC_var_r  PC_fix_f   PC_fix_r
#> 1     97     0  188.63854       NaN  77.01932   0.00000 111.61923        NaN
#> 2    102     0  418.85422       NaN 168.21073   0.00000 250.64349        NaN
#> 3     54     0  554.71796       NaN 149.41288   0.00000 405.30508        NaN
#> 4     50     0  190.17190       NaN  55.08223   0.00000 135.08968        NaN
#> 5     76     0  373.57268       NaN 118.12758   0.00000 255.44511        NaN
#> 6    194     0  457.06087       NaN 263.28660   0.00000 193.77427        NaN
#> 7     39     0 1122.42931       NaN 225.01709   0.00000 897.41222        NaN
#> 8    302     0  298.49092       NaN 198.08093   0.00000 100.41000        NaN
#> 9    322     0  180.03006       NaN 121.77266   0.00000  58.25740        NaN
#> 10    61     0  612.44256       NaN 179.34935   0.00000 433.09321        NaN
#> 11   513   513  137.08649  160.8032 107.93239 107.93239  29.15411   52.87078
#> 12    99    99  645.97640  985.3607 259.64152 259.64152 386.33488  725.71913
#> 13   583   583   93.33452  110.9780  74.21814  74.21814  19.11638   36.75984
#> 14    73    73  738.43345 1168.5480 248.87513 248.87513 489.55832  919.67286
#> 15   361   361  190.49178  239.8642 135.00141 135.00141  55.49037  104.86283
#> 16    34    34  904.52323 1551.9010 171.85581 171.85581 732.66742 1380.04516
#> 17    39    39 1108.61344 1891.1843 231.66248 231.66248 876.95095 1659.52183
#> 18   204   204  189.61919  260.5679 111.07255 111.07255  78.54665  149.49530
#> 19   228   228  304.29314  406.5427 185.42337 185.42337 118.86976  221.11932
#> 20    30    30  682.69417 1147.4733 122.34386 122.34386 560.35031 1025.12941
#> 21    86    86  457.10036  724.6074 161.22652 161.22652 295.87384  563.38084
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

Meyer, M., Meßerschmidt, O., & Mertens, K. G. (2019). *How much does
variety-induced complexity actually cost? Linking axiomatic design with
cost modelling*. In M. Schröder & K. Wegner (Eds.), *Logistik im Wandel
der Zeit – Von der Produktionssteuerung zu vernetzten Supply Chains*
(pp. 813–827). Wiesbaden: Springer Fachmedien Wiesbaden.
<https://doi.org/10.1007/978-3-658-25412-4_39>

Schmidt M, Mertens, K. G., Meyer M (2023) Cost hierarchies and the
pattern of product cost cross-subsidization: Extending a computational
model of costing system design. PLOS ONE 18(9): e0290370.
<https://doi.org/10.1371/journal.pone.0290370>
