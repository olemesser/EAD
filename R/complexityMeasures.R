#' @title Calculate Design Complexity
#' @description This measure is based on the Shannon entropy. Adapted from Modrak & Bednar (2015).
#' @param A A domain mapping matrix (DMM) or dependency structure matrix (DSM), reflecting the transition between two domains.
#' @param norm Boolean, default\code{TRUE}. Decides if the normalized design complexity measures should be returned.
#' @return The calculated design complexity measure.
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
#' @description This measure calculates the reangularity for a given domain mapping matrix (DMM) according to Suh (1995) and Delaš et al. (2018).
#' @param DMM A domain mapping matrix, reflecting the transition between two domains.
#' @return The calculated reangularity measure.
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


#' @title Calculate Semiangularity
#' @description This measure calculates the semiangularity for a given domain mapping matrix (DMM) according to Suh (1995) and Delaš et al. (2018).
#' @param DMM A domain mapping matrix, reflecting the transition between two domains.
#' @return The calculated semiangularity measure.
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

#' @title measure_modularity
#' @description Calculates the directed of undirected modularity measure of a network.
#' The directed approach is descriped by Blondel, Guillaume, Lambiotte, & Lefebvre (2008) as well as Newman (2006).
#' The undirected measure is based on Leicht & Newmann (2008).
#' @param A A data.frame containing the dependency structure matrix with information on element dependencies.
#' @param preMember A numerical vector containing the assignment of vertices to cluster.
#' The default value is \code{preMember=NULL} which means that the communities are identified using the \link[igraph]{fastgreedy.community}
#' @param plot_op Boolean, default \code{plot_op=FALSE}. If set to true the graph is visualized.
#' @return The calculated modularity measure
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
             mark.groups = communities(net_group))
      }
  }else{
    output<-NA
  }
  return(output)
}

#' @title Calculate structural complexity
#' @description This measure is defined by Kim et al. (2013), Sinha & de Weck (2014) as well as Sinha & de Weck (2016).
#' The C1 term is neglected here. Thervore only C2 and C3 are calculated. \code{SC=C2*C3*1/NROW(A)}.
#' Non-symmetric matrices are transformed into a symmetric binary matrix.
#' @param DSM A dependency structure matrix, reflecting the interfaces between elements.
#' @param nrom Boolean, default=F. If set to true, the normalized measure is returned.
#' Sinha & de Weck (2016) proof that the maximum matrix energy is always \eqn{E_max<=n^{3/2}}
#' @return The calculated structural complexity measure.
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
    # A<-makeMatrixsymmetric(A)
    A[A>1]<-1
    diag(A)<-0
    C_2<-sum(A)
    C_3<-sum(svd(A)$d)
    output <- C_2*1/NROW(A)*C_3
    if(norm){
      # B<-A
      # B[A>=0]<-1
      # output <- output/measure_structuralcomplexity(B)
      C_3_max<-dim(A)[1]/2*(1+sqrt(dim(A)[1]))
      C2_max<-(prod(dim(A))-NROW(A)) # number of elements minus diagonal
      output=output/(C_3_max*C2_max)
    }
  }else{
    output<-NA
  }
  return(output)
}

#' @title Calculate the von Neumann entropy for graphs.
#' @description This function takes a DSM and calculates the von Neumann entropy. For further details see: Passerini & Severini (2009).
#' @param DSM A dependency structure matrix, reflecting the interfaces between elements.
#' @return The calculated von Neumann entropy measure.
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
      B<-A
      B[A>=0]<-1
      S<-S/measure_neumannEntropy(B,norm=F)
    }
  }else{
    S<-NA
  }
    return(S)
}





#### Product Measures ####

#' @title Calculate normalized dissimilarity
#' @description Takes a product matrix with products in rows and attributes in columns. This measure is defined by Buchholz (2012).
#' @param P A product matrix with products in rows and attributes in columns.
#' @return The calculated normalized dissimilarity measure.
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
#' For further details on this measure see Buchholz (2012).
#' @param P A product matrix with products in rows and attributes in columns.
#' @return The calculated sum of distances.
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
#' @description Calculates the INTRA measure defined by Gupta (1993)
#' @param P A product matrix (e.g, RES_CONS_PAT) with products in rows and features, components, process or resources in columns.
#' @param norm Boolean, default \code{FALSE}. If set to \code{TRUE}, then the normalized measures is used.
#' @return Returns the calculated intra value.
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
#' @description Calculates the INTRA measure defined by Gupta (1993) and adapted by Mertens (2020).
#' @param P A product matrix (e.g, RES_CONS_PAT) with products in rows and features, components, process or resources in columns.
#' @return Returns the calculated inter value.
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
#' @description The product mix heterogeneity (PMH)as defined by Anderson (1995) for a given product matrix is calculated with this function.
#' Since the PMH is calculated for each individual resource the mean and standard deviation are returned.
#' @param P A product matrix (e.g, RES_CONS_PAT) with products in rows and features, components, process or resources in columns.
#' @return Returns the corresponding PMH value
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
#' @description Calculates the Herfindahl-based product diversification index (D) as defined by Gollop (1997)
#' @param P A product matrix (e.g, RES_CONS_PAT) with products in rows and features, components, process or resources in columns.
#' @param DMD A demand vector. If no demand vector is given, then a normal distribution is assumed where each product has the quantity one.
#' @return Returns the corresponding index value
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
#' @description This index was introduced by Kota et al. (2000). The f-values as defined by Kota are set to one.
#' @param P A product matrix representing the components per product (P_DD) with products in rows and components in columns.
#' @return Returns the corresponding index value
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
#' @description Calculates the local outlier factor of observation (rows) in \code{P}. The input is passed to \link[DescTools]{LOF}
#' @param P A product matrix (e.g, RES_CONS_PAT) with products in rows and features, components, process or resources in columns.
#' @param n A value between \code{2<n<NROW(P)} indicating the number of products used to calculate the LOF's.
#' @return Returns the corresponding index value
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
#' @param P A product matrix (e.g, RES_CONS_PAT) with products in rows and features, components, process or resources in columns.
#' @return Returns the corresponding density value
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
#' Anand et al. (2019) use this measure to operationalize resource cost distribution.
#' @param x A numeric vector which is strictly positive \code{0<x}.
#' @return Returns the corresponding proportion
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
#' measure_NPV(CSD_EAD$DMD)
measure_NPV<-function(P){
  temp<-lapply(1:NCOL(P),function(x){
    unique(c(0,P[,x]))
  })
  temp<-prod(sapply(temp,length))-1
  return(NROW(P)/temp)
}
