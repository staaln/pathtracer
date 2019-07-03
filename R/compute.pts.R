# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

#' Computes PathTracer scores

compute.pts = function(dat, reference, ncomp=4, normalize=T,pathwayindex) {
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
  qds = sqrt(apply(res$s, 1, function(x) {sum((x-xcen)^2)}))

  # Test for differential QDS score in reference samples vs others
  pval = wilcox.test(qds[reference], qds[!reference])$p.value

  #auc
  roc.res = roc(as.numeric(reference),qds)
  auc.res = auc(roc.res)
  roc.res = roc(as.numeric(reference),pds)
  auc.res.2 = auc(roc.res)

  # Return results
  list(qds=qds, pds=pds, p.value=pval, auc.qds=auc.res, auc.pds=auc.res.2, data=dat, ref=reference, res=res,
       xcen=xcen, v=udv$v[,1:m], d=udv$d[1:m],pathwayindex=pathwayindex)
}
