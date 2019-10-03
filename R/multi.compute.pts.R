#' @title multi.compute.pts
#' @description Computes PathTracer deregulation scores based on gene expressions for a set of pathways. Calls compute.pts
#' @param data Matrix of gene expression data with rows corresponding to genes and columns to samples.
#' @param reference Bolean vector indicating which of the columns in \emph{dat} that are reference samples (e.g normal samples).
#' @param ncomp Number of principal components.
#' @param normalize Boolean. If true, gene expression values are normalized (for each gene, the mean is subtracted and the standard deviation divided upon).
#' @param pathwaydatabase Currently, code only for reactome.db is implemented
#' @param ncores Number of cores used for parallellization. Default value is 1, which implies no parellillization.
#' @param min.n.genes Minimal number of genes in a pathway required to calculate deregulations scores for the pathway
#' @return res List where each element is the result of compute.pts for one pathway.
#' @references Nygard S, Lingjaerde OC, Caldas C, Hovig E, Borresen-Dahle, Helland Aa, Haakensen V.
#' "PathTracer: High-sensitivity detection of differential pathway activity in tumours".
#' @details The main PathTracer function. The input is a matrix of gene expressions
#' for a particular pathway and a set of samples, with rows representing genes, and colums samples.
#' The function first calculates the principal curve for the first \emph{ncomp} principal components and projects all samples onto the principal curve.
#' The PathTracer deregulation score (PTS) for each sample is calculated as the Euclidean distance from the sample's projection onto the principal curve to the
#' reference point (defined as the projection onto the curve for the reference sample with median distance along the curve to the start of the curve). See Nygard et al (2019), for furhter details.

multi.compute.pts = function(data, gene_names=row.names(data), gene.identifier="symbol",reference,
                             pathwaydatabase="reactome.db",min.n.genes=10,
                             normalize=T,ncomp=4, ncores=1){
  ptwy = get.pathways(pathwaydatabase=pathwaydatabase,id=gene.identifier, min.n.genes=min.n.genes)
  p<-length(ptwy$id)
  print(p)
  if (normalize) {
    data<-t(scale(t(data)))
  }
  res<-list()
  if (ncores==1) {#no parallellization
    print("OK2")
    for (i in 1:p){
      print(i)
      keep = gene_names %in% ptwy$genes[[i]]
      print(length(which(keep==TRUE)))
      #print(data[keep,])
      print(max(data[keep,]))
      if(length(which(keep==TRUE))>min.n.genes){
        res[[i]]= compute.pts(data[keep,], reference=reference,ncomp=ncomp,pathwayindex=i)
      }
    }
  }
  else {
    registerDoMC(cores=ncores)
    res<-foreach(i=1:p) %dopar% {
      print(i)
      keep = gene_names %in% ptwy$genes[[i]]
      if(length(which(keep==TRUE))>min.n.genes){
        compute.pts(data[keep,], reference=reference,ncomp=ncomp,pathwayindex=i)
      }
    }
  }
  # Return results
  return(res)
}

