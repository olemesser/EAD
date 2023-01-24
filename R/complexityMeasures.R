#' @title Calculate Design Complexity
#' @description This measure is based on the Shannon entropy. Adapted from \insertCite{Modrak.2015;textual}{EAD}.
#' @param A A domain mapping matrix (DMM) or dependency structure matrix (DSM), reflecting the transition between two domains.
#' @param norm Boolean, default\code{TRUE}. Decides if the normalized design complexity measures should be returned.
#' @return The calculated design complexity measure.
#' @references
#' \insertAllCited{}
#' @examples
#'
#' A<-matrix(c(1,0,0,
#'             1,1,0,
#'             1,1,1),
#'          nrow=3,
#'          ncol=3,
#'          byrow=T)
#'
#' measure_designComplexity(A)
measure_designComplexity<-function(A,norm=T){
  DC <- sum(apply(A,2,function(x){
    z<-sum(x[x>0])*log(sum(x[x>0]))
    if(sum(x)==0) z<-0
    return(z)
  }))

  if(norm){
    rM<-apply(A,1,max)
    A_max<-matrix(rep(rM, dim(A)[2]),nrow = dim(A)[1],ncol = dim(A)[2])
    DC <- DC/sum(apply(A_max,2,function(x) sum(x[x>0])*log(sum(x[x>0]))))
  }
  return(DC)
}



#' @title Calculate Reangularity
#' @description This measure calculates the reangularity for a given domain mapping matrix (DMM) according to \insertCite{Suh.1995;textual}{EAD} and \insertCite{Delas.2018;textual}{EAD}.
#' @param DMM A domain mapping matrix, reflecting the transition between two domains.
#' @return The calculated reangularity measure.
#' @references
#' \insertAllCited{}
#' @examples
#'
#' DMM<-matrix(c(1,0,0,
#'             1,1,0,
#'             1,1,1),
#'          nrow=3,
#'          ncol=3,
#'          byrow=T)
#'
#' measure_reangularity(DMM)
measure_reangularity<-function(A){
  if(dim(A)[1]==dim(A)[2]){
    output<-1
      for(i in 1:(NROW(A)-1)){
        for(j in (1+i):NCOL(A)){
          output <- output * sqrt(1-sum(A[,i]*A[,j])^2/(sum(A[,i]^2)*sum(A[,j]^2)))
        }
      }
  }else{
    output<-NA
  }
  return(output)
}


#' @title Calculate Semangularity
#' @description This measure calculates the semangularity for a given domain mapping matrix (DMM) according to \insertCite{Suh.1995;textual}{EAD} and \insertCite{Delas.2018;textual}{EAD}.
#' @param DMM A domain mapping matrix, reflecting the transition between two domains.
#' @return The calculated semangularity measure.
#' @references
#' \insertAllCited{}
#' @examples
#'
#' DMM<-matrix(c(1,0,0,
#'             1,1,0,
#'             1,1,1),
#'          nrow=3,
#'          ncol=3,
#'          byrow=T)
#'
#' measure_semiangularity(DMM)
measure_semiangularity<-function(A){
  if(dim(A)[1]==dim(A)[2]){
    output<-1
    for(j in 1:NCOL(A)){
      output <- output * abs(A[j,j])/sqrt(sum(A[,j]^2))
    }
  }else{
    output<-NA
  }

  return(output)
}


#' @title Calculates Jung's system design complexity
#' @description This measure calculates the system design complexity  for a given domain mapping matrix (DMM) according to \insertCite{Jung.2022;textual}{EAD}.
#' @param DMM A domain mapping matrix, reflecting the transition between two domains.
#' @return The calculated system design complexity measure.
#' @references
#' \insertAllCited{}
#' @examples
#' ## example as shown in Jung et al. (2022)
#' data('RDU1')
#' measure_JSDC(RDU1) # must be 389.5
#'
#'
#' data('RDU2')
#' measure_JSDC(RDU2) # must be 466.1
measure_JSDC<-function(DMM,norm=F){
  energy<-sum(svd(DMM)$d)
  JSDC <- sum(DMM) * energy/min(dim(DMM))

  if(norm){
    JSDC<-JSDC
    # DMM_max<-matrix(1,nrow=dim(DMM)[1],ncol=dim(DMM)[2])
    # JSDC_max<-measure_JSDC(DMM_max,norm=F)
    # JSDC<-JSDC/JSDC_max
  }
  return(JSDC)
}




#' @title measure_modularity
#' @description Calculates the directed of undirected modularity measure of a network.
#' The directed approach is described by \insertCite{Blondel.2008;textual}{EAD} as well as \insertCite{Newman.2006;textual}{EAD}.
#' The undirected measure is based on \insertCite{Leicht.2008;textual}{EAD}.
#' @param A A data.frame containing the dependency structure matrix with information on element dependencies.
#' @param preMember A numerical vector containing the assignment of vertices to cluster.
#' The default value is \code{preMember=NULL} which means that the communities are identified using the \link[igraph]{fastgreedy.community}
#' @param plot_op Boolean, default \code{plot_op=FALSE}. If set to true the graph is visualized.
#' @return The calculated modularity measure
#' @references
#' \insertAllCited{}
#' @examples
#'
#' ## unsymmetrical matrix as input
#' ## is transformed into a symmetric matrix first
#' ## use the undirected approach by Blondel, Guillaume, Lambiotte, & Lefebvre (2008) as well as Newman (2006)
#'   DSM<-matrix(c(1,0,0,
#'               1,1,0,
#'               1,1,1),
#'            nrow=3,
#'            ncol=3,
#'            byrow=T)
#'   measure_modularity(DSM,symmetric=T)
#'
#' ## unsymmetrical matrix as input
#' ## not transformed into a symmetric matrix first
#' ## use the directed approach by Leicht & Newmann (2008)
#'   measure_modularity(DSM,symmetric=F)
#'
measure_modularity<-function(A,preMember=NULL,plot_op=F){
   require(igraph)

  if(dim(A)[1]==dim(A)[2]){
      ## binarize matrix ##
      A[A>1]<-1
      ## case symmetric input matrix
      diag(A)<-0 #remove diagonal elements
      net<-graph_from_adjacency_matrix(A,mode = "undirected") #if A is non symmetric transform into symmetric matrix
      if(is.null(preMember)){
        net_group<-fastgreedy.community(net,merges=T, modularity=T,membership = T)
        output<-modularity(net_group)
      }else{
        output<-modularity(net,as_membership(preMember))
      }

      if(plot_op & is.null(preMember)){
        plot(net, col = membership(net_group),
             mark.groups = communities(net_group),vertex.size=20,edge.width=5)
      }
  }else{
    output<-NA
  }
  return(output)
}

#' @title Calculate structural complexity
#' @description This measure is defined by \insertCite{Sinha.2018;textual}{EAD} and \insertCite{Sinha.2014;textual}{EAD}.
#' The C1 term is neglected here. Therefore only C2 and C3 are calculated. \code{SC=C2*C3*1/NROW(A)}.
#' Non-symmetric matrices are transformed into a symmetric binary matrix.
#' @param DSM A dependency structure matrix, reflecting the interfaces between elements.
#' @param nrom Boolean, default=F. If set to true, the normalized measure is returned.
#' \insertCite{Sinha.2016;textual}{EAD} proof that the maximum matrix energy is always \eqn{E_max<=n^{3/2}}
#' @return The calculated structural complexity measure.
#' @references
#' \insertAllCited{}
#' @examples
#'
#' DSM<-matrix(c(1,0,0,
#'             1,1,0,
#'             1,1,1),
#'          nrow=3,
#'          ncol=3,
#'          byrow=T)
#'
#' measure_structuralcomplexity(DSM)
measure_structuralcomplexity<-function(A,norm=F){
  if(dim(A)[1]==dim(A)[2]){
    A[A>1]<-1
    A<-makeMatrixsymmetric(A)
    diag(A)<-0
    C_1<-NROW(A)
    C_2<-sum(A)
    C_3<-sum(svd(A)$d)
    output <- C_1+C_2*C_3/NROW(A)
    if(norm){
      C_3_max<-dim(A)[1]*sqrt(dim(A)[1]-1)
      C2_max<-(prod(dim(A))-NROW(A)) # number of elements minus diagonal
      output=C_2*C_3/(NROW(A)*C_3_max*C2_max)
    }
  }else{
    output<-NA
  }
  return(output)
}

#' @title Calculate the von Neumann entropy for graphs.
#' @description This function takes a DSM and calculates the von Neumann entropy. For further details see: \insertCite{Passerini.2009;textual}{EAD}.
#' @param DSM A dependency structure matrix, reflecting the interfaces between elements.
#' @return The calculated von Neumann entropy measure.
#' @references
#' \insertAllCited{}
#' @examples
#'
#' DSM<-matrix(c(1,0,0,
#'             1,1,0,
#'             1,1,1),
#'          nrow=3,
#'          ncol=3,
#'          byrow=T)
#'
#' measure_neumannEntropy(DSM)
measure_neumannEntropy<-function(A,norm=T){
  if(dim(A)[1]==dim(A)[2]){
    d<-A
    d<-makeMatrixsymmetric(d)
    diag(d)<-0
    L<-diag(rowSums(d))-d # create Laplacian matrix
    lambda<-eigen(L)$values
    S<-sum(sapply(lambda,function(x){
      if(x<=0) 0 else  x*log2(x)
      }))
    if(norm){
      lambda_max<-rep(dim(A)[1],dim(A)[1])
      S_max<-sum(sapply(lambda_max,function(x){
        if(x<=0) 0 else  x*log2(x)
      }))
      S<-S/S_max
    }
  }else{
    S<-NA
  }
    return(S)
}





#### Product Measures ####

#' @title Calculate normalized dissimilarity
#' @description Takes a product matrix with products in rows and attributes in columns. This measure is defined by \insertCite{Buchholz.2012;textual}{EAD}.
#' @param P A product matrix with products in rows and attributes in columns.
#' @return The calculated normalized dissimilarity measure.
#' @references
#' \insertAllCited{}
#' @examples
#'
#' P<-matrix(c(1,1,1,1,
#'             0,1,1,0,
#'             0,0,1,1),
#'            nrow = 4,
#'            ncol = 3)
#'
#' measure_DISSnrom(P)
measure_DISS<-function(P,norm=T){
  # if(norm){
  #   s<-sqrt(sum(apply(P,2,max)^2))
  # }else{
  #   s<-1
  # }
  # diss<-measure_DISTsum(P)/((dim(P)[1]*(dim(P)[1]-1))*s)

  P_max<-expand.grid(lapply(1:NCOL(P),function(x) sort(unique(P[,x]))))
  P_max<-P_max[apply(P_max,1,function(x) sum(x)>0),]
  diss<-measure_DISTsum(P)
  if(norm){
    diss<-diss/measure_DISTsum(P_max)
  }
  return(diss)
}


#' @title Calculate sum of distances
#' @description Calculates the absolute sum of euclidean distances for a given matrix with elements in rows and attributes in columns.
#' For further details on this measure see \insertCite{Buchholz.2012;textual}{EAD}.
#' @param P A product matrix with products in rows and attributes in columns.
#' @return The calculated sum of distances.
#' @references
#' \insertAllCited{}
#' @examples
#'
#' P<-matrix(c(1,1,1,1,
#'             0,1,1,0,
#'             0,0,1,1),
#'            nrow = 4,
#'            ncol = 3)
#'
#' measure_DISTsum(P)
measure_DISTsum<-function(P){
  if(NROW(P)==1) m<-0 else m<-sum(as.matrix(dist(P,method = "euclidean")))
  return(m)
}

#' @title product distance mean coefficient of variation.
#' @description Calculates the mean distance of each product. This vector is used to calculate the coefficient of variation (cv) in a second step.
#' A cv of zero means, that all products are equally distributed. Larger cv values imply a more heterogeneous product mix.
#' @param P A product matrix with products in rows and attributes in columns.
#' @return The coefficient of variation of products' mean distances.
#' @examples
#'
#' P<-matrix(c(1,1,1,1,0,
#'             0,1,1,0,1,
#'             0,0,1,1,1),
#'            nrow = 5,
#'            ncol = 3)
#'
#' measure_DISTcv(P)
measure_DISTcv<-function(P){
  dist_P<-as.matrix(dist(P))
  dist_P<-rowMeans(dist_P)
  return(sd(dist_P)/mean(dist_P))
}


#' @title Calculate intra-product heterogeneity
#' @description Calculates the INTRA measure defined by \insertCite{Gupta.1993;textual}{EAD} and adapted by \insertCite{Mertens.2020;textual}{EAD}.
#' @param P A product matrix (e.g, RES_CONS_PAT) with products in rows and features, components, process or resources in columns.
#' @param norm Boolean, default \code{FALSE}. If set to \code{TRUE}, then the normalized measures is used.
#' @return Returns the calculated intra value.
#' @references
#' \insertAllCited{}
#' @examples
#'
#' P<-matrix(c(2,1,1,1,
#'            2,1,1,3,
#'            2,2,2,2,
#'            1,1,1,1),
#'          nrow = 4,
#'          ncol = 4,byrow=T)
#' measure_INTRA(P)
measure_INTRA<-function(P,norm=F){
  if(norm){
    P_scaled<-P/rowSums(P)
  }else{
    P_scaled<-apply(P,2,function(x) x/sum(x))
  }
  rmP<-rowMeans(P_scaled)
  INTRA<-P_scaled-rmP
  INTRA<-t(sapply(1:NROW(P_scaled),function(i){
    INTRA[i,]/rmP[i]
  }))
  INTRA<-rowSums(INTRA^2)
  if(norm) INTRA/mean(P)
  return(INTRA)
}


#' @title Calculate inter-product heterogeneity
#' @description Calculates the INTRA measure defined by \insertCite{Gupta.1993;textual}{EAD} and adapted by \insertCite{Mertens.2020;textual}{EAD}.
#' @param P A product matrix (e.g, RES_CONS_PAT) with products in rows and features, components, process or resources in columns.
#' @return Returns the calculated inter value.
#' @references
#' \insertAllCited{}
#' @examples
#'
#' P<-matrix(c(2,1,1,1,
#'            2,1,1,3,
#'            2,2,2,2,
#'            1,1,1,1),
#'          nrow = 4,
#'          ncol = 4,byrow=T)
#' measure_INTER(P)
measure_INTER<-function(P){
  P_scaled<-apply(P,2,function(x) x/sum(x))
  cm_mat<-matrix(rep(colMeans(P_scaled),NROW(P_scaled)),nrow =NROW(P_scaled),byrow = T)
  INTER<-rowSums(((P_scaled-cm_mat)/cm_mat)^2)
  return(INTER)
}



measure_EXOTICNESS<-function(P){
  dist_mat<-as.matrix(dist(P))
  prod_mean_dist_uw<-rowMeans(dist_mat)
  return(prod_mean_dist_uw/mean(prod_mean_dist_uw))
}


measure_NN<-function(P,prop=0.2){
  NUMB_n<-ceiling(NROW(P)*prop) # number of nearest neighbors
  P_scale<-as.matrix(scale(P))
  P_dist<-as.matrix(dist(P_scale)) # calculate distance matrix
  P_dist_mean<-mean(rowMeans(P_dist))
  diag(P_dist)<-NA # set diag NA
  NN<-apply(P_dist,1,function(x){
    mean(sort(as.numeric(na.omit(x)))[1:NUMB_n])
  })
  return(as.numeric(NN))
}



#' @title Calculate the Product Mix Heterogeneity (PMH)
#' @description The product mix heterogeneity (PMH)as defined by \insertCite{Anderson.1995;textual}{EAD} for a given product matrix is calculated with this function.
#' Since the PMH is calculated for each individual resource the mean and standard deviation are returned.
#' @param P A product matrix (e.g, RES_CONS_PAT) with products in rows and features, components, process or resources in columns.
#' @return Returns the corresponding PMH value
#' @references
#' \insertAllCited{}
#' @examples
#'
#' data('csd_EAD')
#' P<-CSD_EAD$RES_CONS_PAT
#'
#' measure_PMH(P)
measure_PMH<-function(P){
  a_mean<-colMeans(P)
  mach_hour<-rowSums(P)
  PMH<-sapply(1:NCOL(P),function(j){
        sqrt(sum(mach_hour * (P[,j]-a_mean[j])^2) / (sum(mach_hour)-1))
      })
  return(PMH)
}

#' @title Calculate the Product Diversification Index
#' @description Calculates the Herfindahl-based product diversification index (D) as defined by \insertCite{Gollop.1997;textual}{EAD} and \insertCite{Gollop.1991;textual}{EAD}.
#' @param P A product matrix (e.g, RES_CONS_PAT) with products in rows and features, components, process or resources in columns.
#' @param DMD A demand vector. If no demand vector is given, then a normal distribution is assumed where each product has the quantity one.
#' @return Returns the corresponding index value
#' @references
#' \insertAllCited{}
#' @examples
#'
#' data('csd_EAD')
#' P<-CSD_EAD$RES_CONS_PAT
#' DMD<-crt_DEMAND(NROW(P),1000,1)$DEMAND
#'
#' measure_diversificationINDEX(P,DMD)
measure_diversificationINDEX<-function(P,DMD=NULL){
  if(is.null(DMD)) DMD<-rep(1,NROW(P))
  DMD<-DMD/sum(DMD)
  P_dist<-t(apply(P,1,function(x) x/sum(x)))
  sigma<-matrix(NA,ncol = NROW(P),nrow = NROW(P))
  for(i in 1:NROW(P)){
    for(k in 1:NROW(P)){
      sigma[i,k]<-sqrt(sum(abs(P_dist[k,]-P_dist[i,])/2))
      sigma[i,k]<-sigma[i,k]*DMD[i]*DMD[k]
    }
  }
  return(1/2*(1-sum(DMD^2)+sum(sigma)))
}

#' @title Calculate the Product Line Commonality Index
#' @description This index was introduced by \insertCite{Kota2000;textual}{EAD}. The f-values as defined by Kota are set to one.
#' @param P A product matrix representing the components per product (P_DD) with products in rows and components in columns.
#' @return Returns the corresponding index value
#' @references
#' \insertAllCited{}
#' @examples
#'
#' P<-matrix(c(1,0,0,1,
#'               0,1,0,0,
#'               1,0,0,0,
#'               0,1,1,0,
#'               0,1,0,1),
#'        nrow = 4,
#'        ncol = 5)
#'
#' measure_PCI(P)
measure_PCI<-function(P){
  # transform into binary matrix
  P[P>1]<-1
  n_i<-colSums(P)
  MinCCI<-sum(1/n_i^2)
  PCI<-(sum(n_i)-MinCCI)/(NROW(P)*NCOL(P)-MinCCI)
  return(PCI)

}


#' @title Calculate Local Outlier Factor (LOF)
#' @description Calculates the local outlier factor of observation (rows) in \code{P} as defined by \insertCite{Breunig2000;textual}{EAD} and \insertCite{XU2022;textual}{EAD} The input is passed to \link[DescTools]{LOF}
#' @param P A product matrix (e.g, RES_CONS_PAT) with products in rows and features, components, process or resources in columns.
#' @param n A value between \code{2<n<NROW(P)} indicating the number of products used to calculate the LOF's.
#' @return Returns the corresponding index value
#' @references
#' \insertAllCited{}
#' @examples
#'
#' data('csd_EAD')
#' P<-CSD_EAD$RES_CONS_PAT
#'
#' measure_LOF(P,n=2)
measure_LOF<-function(P,perc=0.1){
  require(DescTools)
  LOF<-LOF(P,k = ceiling(NROW(P)*perc))
  # LOF<-LOF(P,k = 2)
  return(LOF)
}


#' @title Measure density of a product matrix
#' @description This function measures the density defined as the proportion of non-zero entries in P.
#' This measure is used in several studies \insertCite{Anand2019,BLH2011}{EAD}.
#' @param P A product matrix (e.g, RES_CONS_PAT) with products in rows and features, components, process or resources in columns.
#' @return Returns the corresponding density value
#' @references
#' \insertAllCited{}
#' @examples
#'
#' data('csd_EAD')
#' P<-CSD_EAD$RES_CONS_PAT
#'
#' measure_DENS(P)
measure_DENS<-function(P){
  return(sum(P>0)/prod(dim(P)))
}

#' @title Measure the skewness of an vector
#' @description This function measures the skewness of an given input vector such as demand or resource costs.
#' It calculates the proportion of top ten percentage of elements on the total sum.
#' \insertCite{Anand2019;textual}{EAD} use this measure to operationalize resource cost distribution.
#' @param x A numeric vector which is strictly positive \code{0<x}.
#' @return Returns the corresponding proportion
#' @references
#' \insertAllCited{}
#' @examples
#'
#' data('csd_EAD')
#' measure_TOP10(CSD_EAD$DMD)
measure_TOP10<-function(x){
  x<-sort(x,decreasing = T)
  x<-sum(x[1:ceiling(length(x)*.1)])/sum(x)
  return(ifelse(is.nan(x),0,x))
}


#' @title Measure the proportion of product variants
#' @description This function measures the proportion of product variants present in \code{P} compared to the free combination.
#' @param P A product matrix (e.g, RES_CONS_PAT) with products in rows and features, components, process or resources in columns.
#' @return Returns the NPV value
#' @examples
#'
#' P<-matrix(c(1,0,0,1,
#'               0,1,0,0,
#'               1,0,0,0,
#'               0,1,1,0,
#'               0,1,0,1),
#'        nrow = 4,
#'        ncol = 5)
#' measure_NPV(P)
measure_NPV<-function(P){
  temp<-lapply(1:NCOL(P),function(x){
    unique(P[,x])
  })
  temp<-prod(sapply(temp,length))-1
  return(NROW(P)/temp)
}

#' @title Calculates the option variability
#' @description This function calculates the option variability measure defined by \insertCite{MacDuffie1996;textual}{EAD}
#' @param P A product matrix containing the features with products in rows and features in columns.
#' @param norm Boolean, default=T normalizes the HVM in a range of 0-1 where HVM=1 is the maximum possible value for a given matrix.
#' @return Returns the OV value
#' @references
#' \insertAllCited{}
#' @examples
#'
#' P<-matrix(c(1,0,0,1,
#'               0,1,0,0,
#'               1,0,0,0,
#'               0,1,1,0,
#'               0,1,0,1),
#'        nrow = 4,
#'        ncol = 5)
#' measure_OV(P)
measure_OV<-function(P,norm=T){
  P[P>0]<-1
  mue<-apply(P,2,function(x) sum(x>0))/NROW(P)
  OV <- sum(sqrt(mue*(1-mue)))
  if(norm){
    OV_max<-sum(sqrt(0.5*(1-0.5)))*length(mue)
    OV <- OV / OV_max
  }
  return(OV)
}

#' @title Calculates the Commonality Index (CI)
#' @description This function calculates the commonality index as defined by \insertCite{MartinIshii1996;textual}{EAD}
#' @param P A product matrix containing the products in rows and and features, components, process or resources in columns.
#' @return Returns the CI value
#' @references
#' \insertAllCited{}
#' @examples
#'
#' P<-matrix(c(1,0,0,1,
#'               0,1,0,0,
#'               1,0,0,0,
#'               0,1,1,0,
#'               0,1,0,1),
#'        nrow = 4,
#'        ncol = 5)
#' measure_CI(P)
measure_CI<-function(P){
  return(NCOL(P)/sum(P[P>0]))
}


#' @title Calculates the Coupling Complexity (CC)
#' @description This function coupling complexity (CC) as defined by \insertCite{Ameri.2008;textual}{EAD} and \insertCite{Summers.2010;textual}{EAD}.
#' Note, this function is not implemented yet since description of authors is not sufficient to reproduce the algorithm.
#' @param DSM A dependency structure matrix, reflecting the interfaces between elements.
#' @return Returns the CC value
#' @references
#' \insertAllCited{}
#' @examples
#' require(EAD)
#'
#' ## water Sprinkler example by Ameri et al. (2008) ##
#' data("waterSprinkler")
#' DSM<-waterSprinkler
#' measure_CoupComp(DSM)
measure_CoupComp<-function(DSM){
  require(data.table)
  require(igraph)
  require(dplyr)

  #### Create Entity Relationship Graph ####
  DSM<-as.data.frame(DSM)
  DSM$src<-colnames(DSM)
  DSM<-setDT(DSM)
  DSM<-melt.data.table(setDT(DSM),id.vars = "src",variable.name="tgt") %>%
        as.data.frame() %>%
        filter(value>0) %>%
        tidyr::uncount(value) %>%
        mutate(node_right=paste0(src,"-",tgt))
  DSM<-rbind(data.frame(node_left=DSM$src,node_right=DSM$node_right),
             data.frame(node_left=DSM$tgt,node_right=DSM$node_right)) %>%
    as.matrix()

    g<-graph_from_edgelist(DSM)
    V(g)$x<-ifelse((V(g)$name) %in% DSM[,2],2,0)
    V(g)$y<-ifelse(V(g)$x==0,
                   match(V(g)$name,sort(unique(DSM[,1]))),
                   match(V(g)$name,sort(unique(DSM[,2]))))


  #### Start Separation Algorithm ####
    level<-1
    total<-0
    ## remove unary constraints
    g<-delete.vertices(g,which(degree(x)==1 & V(x)$name %in% DSM[,2]))

    repeat{
      g<-decompose.graph(g)
      size<-1
      x<-g[[1]]
      lapply(g,function(x){
        on.exit(size=1)
        nodes<-which(degree(x)==size)
        i<-1
        temp_n<-lapply(nodes,function(i)  unlist(lapply(all_simple_paths(x, from = i),function(y) y$name)))
        temp_n<-unique(unlist(temp_n))
        numbSets<-sum(temp_n %in% DSM[,2])
        x_new<-delete.vertices(x,temp_n)
        plot(x_new)
        decompose.graph(x_new)
      })
      size<-1
      nodes<-which(degree(g)==size)
      plot(g_new)
    total<-level*set_size*numbSets+total
      if(cond){
        break
      }else{
        level<-level+1
      }
    }
  return(CC)
}

#' @title Calculates the Halstead-derived volume measure complexity (HVM)
#' @description This function calculates the Halstead-derived volume measure complexity (HVM) as defined by \insertCite{Halstead.1979;textual}{EAD}, \insertCite{Prather.1984;textual}{EAD} and \insertCite{Hennig.2022;textual}{EAD}.
#' @param DSM A dependency structure matrix, reflecting the interfaces between elements.
#' @param norm Boolean, default=T normalizes the HVM in a range of 0-1 where HVM=1 is the maximum possible value for a given matrix.
#' @return Returns the HVM value
#' @references
#' \insertAllCited{}
#' @examples
#' require(EAD)
#'
#' ## water Sprinkler example by Ameri et al. (2008) ##
#' data("waterSprinkler")
#' DSM<-waterSprinkler
#'
#' measure_HVM(DSM)
measure_HVM<-function(DSM,norm=T){
  DSM<-makeMatrixsymmetric(DSM)
  diag(DSM)<-0
  E<-sum(DSM)
  N<-NROW(DSM)
  HVM<-(E+N)*log(E+N)
  if(norm){
    E_max<-N^2-N
    HVM_max<- (E_max+N)*log(E_max+N)
    HVM <- HVM / HVM_max
  }
  return(HVM)
}



#' @title Calculates the Interface Complexity (HIC)
#' @description This function calculates the Interface Complexity (HIC) as defined by \insertCite{Holtta.2005;textual}{EAD}.
#' Note, this function is not implemented yet since description of authors is not sufficient to reproduce the algorithm.
#' @param DSM A dependency structure matrix, reflecting the interfaces between elements.
#' @param norm Boolean, default=T normalizes the HIC in a range of 0-1 where HVM=1 is the maximum possible value for a given matrix.
#' @return Returns the HIC value
#' @references
#' \insertAllCited{}
#' @examples
#' require(EAD)
#'
#' ## water Sprinkler example by Ameri et al. (2008) ##
#' data("waterSprinkler")
#' DSM<-waterSprinkler
#'
#' measure_HIC(DSM)
measure_HIC<-function(DSM,norm=T){
  DSM<-makeMatrixsymmetric(DSM)
  diag(DSM)<-0
  HIC<-sum(DSM)
  if(norm){
    N <- NROW(DSM)
    HIC_max <- N^2-N
    HIC<- HIC / HIC_max
  }
  return(HIC)
}


#' @title Calculates the McCabe's Cyclomatic Complexity (MCC)
#' @description This function calculates the McCabe's cyclomatic complexity (MCC)  \insertCite{McCabe.1976;textual}{EAD}.
#' @param DSM A dependency structure matrix, reflecting the interfaces between elements.
#' @param norm Boolean, default=T normalizes the MCC in a range of 0-1 where HVM=1 is the maximum possible value for a given matrix.
#' @return Returns the MCC value
#' @references
#' \insertAllCited{}
#' @examples
#' require(EAD)
#'
#' ## water Sprinkler example by Ameri et al. (2008) ##
#' data("waterSprinkler")
#' DSM<-waterSprinkler
#'
#' measure_MCC(DSM)
measure_MCC<-function(DSM,norm=T){
  diag(DSM)<-0
  DSM[DSM>1]<-1
  if(sum(DSM)==0){
    MCC<-0
  }else{
    MCC<-sum(DSM)-NROW(DSM)+2
    if(norm){
      DSM_max<-matrix(1,nrow = NROW(DSM),ncol = NROW(DSM))
      diag(DSM)<-1
      MCC_max<-measure_MCC(DSM_max,norm=F)
      MCC <- MCC / MCC_max
    }
  }
  return(MCC)
}
