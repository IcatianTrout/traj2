#'@title Classify the Longitudinal Data Based on the Selected Measures.
#'@description Classifies the trajectories by applying the K-mean clustering
#'  algorithm to the measures selected by \code{Step2Selection}.
#'@param trajSelection Object of class \code{trajSelection} as returned by
#'  \code{Step2Selection}.
#'@param forced.measures Either NULL or a numerical vectors corresponding to the
#'  measures used for clustering. If NULL, the measures selected by
#'  \code{Step2Selection} are used.
#'@param nstart to be passed to the nstart argument of
#'  \code{\link[stats]{kmeans}}.
#'@param iter.max to be passed to the iter.max argument of
#'  \code{\link[stats]{kmeans}}.
#'@param nclusters either NULL or the number of clusters to form. If NULL, the
#'  number of clusters to form will be determined using the GAP criterion as
#'  implemented in the \code{\link[cluster]{clusGap}} function.
#'@param K.max to be passed to the \code{K.max} argument of
#'  \code{\link[cluster]{clusGap}}.
#'@param B to be passed to the \code{SE.factor} argument of
#'  \code{\link[cluster]{clusGap}}.
#'@param d.power to be passed to the \code{B} argument of
#'  \code{\link[cluster]{clusGap}}.
#'@param spaceH0 to be passed to the \code{spaceH0} argument of
#'  \code{\link[cluster]{clusGap}}.
#'@param method to be passed to the \code{method} argument of
#'  \code{\link[cluster]{clusGap}}.
#'@param SE.factor to be passed to the \code{SE.factor} argument of
#'  \code{\link[cluster]{clusGap}}.
#'@return An object of class \code{trajClusters}; a list containing the result
#'  of the clustering as well as curated forms of the data and time matrices.
#'@import cluster
#'
#'@importFrom stats kmeans
#'
#' @examples
#' \dontrun{
#'data = example.data$data
#'
#'m = Step1Measures(data, ID=TRUE)
#'s = Step2Selection(m)
#'
#'s$loadings
#'
#'s2 = Step2Selection(m,select=c(10,12,8,4))
#'
#'groups3Cl <- step3(s2,nclusters=3)$partition
#'groups4Cl <- step3(s2,nclusters=4)$partition
#'groups5Cl <- step3(s2,nclusters=5)$partition
#'
#'}
#'
#'@seealso \code{\link[traj2]{Step2Selection}}
#'
#'@rdname Step3Clusters
#'
#'@export
Step3Clusters <- function (trajSelection, forced.measures = NULL, nstart = 50, iter.max = 20, nclusters = NULL, K.max = 8, B = 500, d.power = 2, spaceH0 = "scaledPCA", method = "Tibs2001SEmax", SE.factor = 1){
  
  if ((!is.null(forced.measures)) && (!is.numeric(forced.measures))) {
    stop("Argument 'forced.measures' must be a numerical vector.")
  }
  if ((!is.null(forced.measures)) & (FALSE %in% (paste("m", forced.measures, sep = "") %in% colnames(trajSelection$measures)))) {
    stop("The measures in the argument 'forced.measures' must be measures present in the argument 'trajMeasures'.")
  }
  
  ID <- trajSelection$selection[, 1]
  
  #standardize the measures to be clustered:
  if (is.null(forced.measures)) {
    data <- data.frame(apply(data.frame(trajSelection$selection[, -1]), 2, scale))
  } else {
    data <- data.frame(apply(data.frame(trajSelection$measures[, paste("m", forced.measures,sep = "")]), 2, scale))
  }
  if (!is.null(nclusters) && (nclusters > nrow(data))) { 
    stop("The number 'nclusters' of requested clusters cannot exceed the number of trajectories.")
  }
  if (is.null(nclusters)) {
    ## The clusGap function takes as argument a function that will perform
    ## the clustering but it does not allow us to set the arguments of that
    ## function, so we have to do so beforehand by defining a new function:
    kmeans.nstart <- function (x, k) {
      return(kmeans(x = x, centers = k, nstart = nstart, iter.max = iter.max))
    }
    GAP <- cluster::clusGap(x = data, FUNcluster = kmeans.nstart, K.max = K.max, B = B, d.power = d.power, spaceH0 = spaceH0)
    nclusters <- maxSE(f = GAP$Tab[, "gap"], SE.f = GAP$Tab[, "SE.sim"], method = method, SE.factor = SE.factor)
    cluster.est2 <- kmeans(x = data, centers = nclusters, iter.max = iter.max, nstart = nstart)
    nb.clust <- nclusters
    partition <- cluster.est2$cluster 
    partition.summary <- summary(factor(partition))
  } else {
    cluster.est2 <- kmeans(x = data, centers = nclusters, iter.max = iter.max, nstart = nstart)
    nb.clust <- nclusters
    partition <- cluster.est2$cluster 
    partition.summary <- summary(factor(partition))
  }
  
  trajClusters <- structure(list(data = trajSelection$data, time = trajSelection$time, nb.clust = nb.clust, nclusters = nclusters, 
                                 partition = partition, partition.summary = partition.summary, class = "trajClusters"))
  
  return(trajClusters)
  
} 
