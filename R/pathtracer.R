#' @title pathtracer
#' @description Computes PathTracer deregulations scores based on gene expressions for a set of genes belonging to a particular pathway, and a set of samples, e.g normal and tumor samples.
#' @param dat Matrix of gene expression data with rows corresponding to genes and columns to samples.
#' @param reference Bolean vector indicating which of the columns in \emph{dat} that are reference samples (e.g normal samples).
#' @param ncomp Number of principal components.
#' @param normalize Boolean. If true, gene expression values are normalized (for each gene, the mean is subtracted and the standard deviation divided upon).
#' @param pathwayindex Integer indicating which pathway the set of input genes belongs to.
#' @return pts Vector of PathTracer deregulation scores
#' @return pds Vector of deregulation scores based on distance along the curve (similar to the pathifier algorithm).
#' @return auc.pts Area under the ROC based on pts values, where reference samples are treated as negatives and non-reference samples as positives.
#' @return auc.pds Area under the ROC based on pds values, where reference samples are treated as negatives and non-reference samples as positives.
#' @return data Normalized expression matrix.
#' @return ref Same as input variable \emph{reference}.
#' @return res Output from principal_curve.
#' @return xcen Reference point.
#' @return v Matrix whose colums contain the left singular values of the normalized data matrix.
#' @return d Vector containing the singular values of the normalized data matrix.
#' @return pathwayindex Same as input \emph{pathwayindex}
#' @references Nygard S, Lingjaerde OC, Caldas C, Hovig E, Borresen-Dahle, Helland Aa, Haakensen V.
#' "PathTracer: High-sensitivity detection of differential pathway activity in tumours". Submitted.
#' @details The main PathTracer function. The input is a matrix of gene expressions
#' for a particular pathway and a set of samples, with rows representing genes, and colums samples.
#' The function first calculates the principal curve for the first \emph{ncomp} principal components and projects all samples onto the principal curve.
#' The PathTracer deregulation score (PTS) for each sample is calculated as the Euclidean distance from the sample's projection onto the principal curve to the
#' reference point (defined as the projection onto the curve for the reference sample with median distance along the curve to the start of the curve). See Nygard et al (2019), for furhter details.

pathtracer = function(dat, reference, ncomp=4, normalize=T,pathwayindex) {
  # Check reference argument
  if (mode(reference) == "logical") {
    if (length(reference) != ncol(dat)) {
      stop(paste("Boolean argument has wrong length: reference"))
    } else if (!any(reference)) {
      stop(paste("Reference population must have at least one member"))
    } else if (all(reference)) {
      stop(paste("Non-reference population must have at least one member"))
    }
  } else if (mode(reference) == "numeric") {
    reference = (1:ncol(dat)) %in% reference
  } else {
    stop("Argument must be logical vector or numeric vector: reference")
  }

  # Normalize features if requested
  if (normalize) {
    for (i in 1:nrow(dat)) {
      dat[i,] = (dat[i,]-mean(dat[i,]))/sd(dat[i,])
    }
  }

  # Perform PCA and keep requested number of components
  m = min(ncomp, min(dim(dat)))
  udv = svd(scale(t(dat), scale=F))
  dat.pca = udv$u[,1:m] %*% diag(udv$d[1:m])  # samples x m

  # Compute principal curve and PDS/QDS scores
  res = principal_curve(dat.pca, start=rbind(dat.pca[reference,], dat.pca[!reference,]))
  pds = res$lambda
  opr = order(pds[reference])
  xcen = res$s[which(reference)[opr[1+length(opr)/2]],]
  pts = sqrt(apply(res$s, 1, function(x) {sum((x-xcen)^2)}))

  # Test for differential QDS score in reference samples vs others
  pval = wilcox.test(qds[reference], pts[!reference])$p.value

  #auc
  roc.res = roc(as.numeric(reference),pts)
  auc.res = auc(roc.res)
  roc.res = roc(as.numeric(reference),pds)
  auc.res.2 = auc(roc.res)

  # Return results
  list(pts=pts, pds=pds, p.value=pval, auc.pts=auc.res, auc.pds=auc.res.2, data=dat, ref=reference, res=res,
       xcen=xcen, v=udv$v[,1:m], d=udv$d[1:m],pathwayindex=pathwayindex)
}
