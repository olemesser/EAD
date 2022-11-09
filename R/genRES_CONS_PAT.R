#' @name genRES_CONS_PAT
#' @title Generates a Resorce Consumption Pattern Matrix
#' @description Based on the choosen method a RES_CONS_PAT is created by this function.
#' @param ProductionEnvironment A structured list. See 'examples' for further details.  \cr
#' @param DENS number of customers which are generated.
#' @param COR1 Correlation on how the unit size resources correlate with the base. Can be a list in the format list(lb=...,ub=...),
#' where lb stands for the lower bound of COR1 and ub for upper bound. See Anand et al. (2017) for further details.
#' @param COR2 Correlation value of how the batch resources should correlate. Base on the choosen method the base is different.
#' See 'RES_method' for further details."
#' @param RES_METHOD Default: "ANAND". Options are: "KGM" and "BALA". The original implementation
#' from Anand et al. (2017) generates correlated numbers where COR2 ist the correlation between base and batch resources.
#' "KGM" and  "BALA" generates the batch resource correlation based on the unit size. THis is equivalent to Balakrishnan et al. (2011)
#' @return Returns RES_CONS_PAT and some statistics as list. \cr
#' @examples
#' ProductionEnvironment[['NUMB_RES']]<-NUMB_RES #Amount of resources
#' ProductionEnvironment[['NUMB_PRO']]<-NUMB_PRO #Amount of product
#'
#' ProductionEnvironment[['UnitSize']]<-floor(VOL_SHARE_RES*ProductionEnvironment[['NUMB_RES']])
#' ProductionEnvironment[['BatchSize']]<-floor((1-VOL_SHARE_RES)*ProductionEnvironment[['NUMB_RES']])
#'
#' # THe base and unit resources correlate with the coefficient of 0.8
#' # unit size resources and batch res. correlate with a level of 0.33 cause of choosen method
#' genRES_CONS_PAT(ProductionEnvironment,DENS=0.4,COR1=0.8,COR2=0.33,RES_METHOD="BALA")
#'
#' @export
NULL


#' @rdname genRES_CONS_PAT
genRES_CONS_PAT<-function(ProductionEnvironment,DENS,COR1,COR2,RES_METHOD="ANAND"){


  repeat {

    #### ---- STEP 3.1 CALCULATING THE RES_CONS_PAT  ---- ####

      RES_CONS_PAT<-COR_ANAND(ProductionEnvironment,COR1=COR1,COR2=COR2)

            # in Anand et al. Flowchart 5.1 (d)
      # Multiply each resource consumption with 10
      for(j in 1:ProductionEnvironment[['NUMB_RES']]){ # Adapted from Anand et al 2017.
        RES_CONS_PAT[,j] = ceiling(abs(RES_CONS_PAT[,j])*10)
      }


      # in Anand et al. Flowchart 5.1 (e)
      ## DENSITY IMPLEMENTATION ##
      for(j in 2:ProductionEnvironment[['NUMB_RES']]){
        d<-runif(ProductionEnvironment[['NUMB_PRO']])
        RES_CONS_PAT[,j] <- RES_CONS_PAT[,j]*ifelse(d<DENS,1,0)
      }


    ####---- in Anand et al. Flowchart 5.1 (f) ----#####
    # EXPECTION HANDLER  & CHECKS
    # Do any products have no resource and do any resources are not used?
    PRO_ZEROS<-any(rowSums(RES_CONS_PAT[,])==0)
    RES_ZEROS<-any(colSums(RES_CONS_PAT[,])==0)
    BASE_ZEROS<-any(RES_CONS_PAT[,1]==0)



    if(PRO_ZEROS==FALSE & RES_ZEROS==FALSE & BASE_ZEROS==FALSE){
      break
    }

  }

  ####---- Check ---####
  # create list object
  CHECK<-list()

  # check density
  DISP1<-ProductionEnvironment$DISP1
  NUMB_RES<-ProductionEnvironment$NUMB_RES
  CHECK[['DENS_h']]<-sum(ifelse(RES_CONS_PAT>0,1,0))/(dim(RES_CONS_PAT)[1]*dim(RES_CONS_PAT)[2])
  CHECK[['COR1_h']]<-mean(cor(RES_CONS_PAT[,1:DISP1])[1,])
  CHECK[['COR2_h']]<-mean(cor(RES_CONS_PAT[,c(1,(DISP1+1):NUMB_RES)])[1,])

  # Assignment of resources to cost objects
  countNonZero<-colSums(RES_CONS_PAT[,]>0)

  elements<-dim(RES_CONS_PAT)[1]*dim(RES_CONS_PAT)[2]
  CHECK[['DENS_INPUT']]<-DENS
  CHECK[['ZerosinMatrix']]<-sum(RES_CONS_PAT==0)/elements
  CHECK[['AverageNUMBofProductsConsumeResource']]<-mean(countNonZero)

  # Average range in Consumption
  RES_CONS_PATp<- as.data.frame(scale(RES_CONS_PAT, center=FALSE, scale=colSums(RES_CONS_PAT)))
      # range in percentage of resource cost consumed by products for a resource, averaged across resources
    CHECK[['AverageRangeinConsumption']]<-mean(do.call(pmax, RES_CONS_PATp)-do.call(pmin,RES_CONS_PATp))


  return(list(RES_CONS_PAT=RES_CONS_PAT,CHECK=CHECK))

}
