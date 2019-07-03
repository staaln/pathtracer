#' Method to apply principal curve to new observations
#' @param obj output from compute.pts
#' @param newdata expression data ({genes in pathway} x samples)
#' @param normalize Boolean variable indicating whether normalization (substraction by mean, and division by standard deviation) of the rows (genes) should be applied
#' @return A list of pathtracer (pts), pathifier (pds) and p-value

predict.pts = function(obj, newdata, normalize=T) {
  if (normalize) {
    for (i in 1:nrow(dat)) {
      newdata[i,] = (newdata[i,]-mean(newdata[i,]))/sd(newdata[i,])
    }
  }
  newdata.pca = t(newdata) %*% obj$v   # nsamples x m
  res = project_to_curve(newdata.pca, obj$res$s)
  pds = res$lambda
  pts = sqrt(apply(res$s, 1, function(x) {sum((x-obj$xcen)^2)}))
  pval = wilcox.test(obj$pts[obj$reference], pts)$p.value
  list(pts=pts, pds=pds, p.value=pval)
}
