#' @title Classify the Functions Based on the Selected Measures. 
#' @description Classifies the functions by applying the K-mean clustering algorithm to the measures selected by \code{step2selection}.
#' @param trajSelection Object of class \code{trajSelection} as returned by \code{step2selection}. 
#' @param forced.measures Either NULL or a numerical vectors corresponding to the measures used for clustering. If NULL, the measures selected by \code{step2selection} are used.
#' @param nstart to be passed to the nstart argument of \code{\link[stats]{kmeans}}. 
#' @param iter.max to be passed to the iter.max argument of \code{\link[stats]{kmeans}}. 
#' @param nclusters either NULL or the number of clusters to form. If NULL, the number of clusters to form will be determined using the GAP criterion as implemented in the \code{\link[psych]{clusGap}} function.
#' @param K.max to be passed to the \code{K.max} argument of \code{\link[cluster]{clusGap}}.
#' @param B to be passed to the \code{SE.factor} argument of \code{\link[cluster]{clusGap}}.
#' @param d.power to be passed to the \code{B} argument of \code{\link[cluster]{clusGap}}.
#' @param spaceH0 to be passed to the \code{spaceH0} argument of \code{\link[psych]{clusGap}}.
#' @param method to be passed to the \code{method} argument of \code{\link[cluster]{clusGap}}.
#' @param SE.factor to be passed to the \code{SE.factor} argument of \code{\link[cluster]{clusGap}}.
#' @return An object of class \code{trajClusters}; a list containing the result of the clustering as well as curated forms of the data and time matrices. 
#' @import cluster
#' 
#' @importFrom stats kmeans
#'
#' @examples
#' \dontrun{
##'data = example.data$data
#'
#'m = step1measures(data, ID=TRUE)
#'s = step2selection(m)
#'
#'s$loadings
#'
#'s2 = step2selection(m,select=c(10,12,8,4))
#'
#'groups3Cl <- step3(s2,nclusters=3)$partition
#'groups4Cl <- step3(s2,nclusters=4)$partition
#'groups5Cl <- step3(s2,nclusters=5)$partition
#'
#'}
#'
#'@seealso 
#' \code{\link[traj]{step2selection}}
#'
#' @rdname step3clusters
#'
#' @export
step3clusters <- function (trajSelection, forced.measures = NULL, nstart = 50, iter.max = 20, nclusters = NULL, K.max = 8, B = 500, d.power = 2, spaceH0 = "scaledPCA", method = "Tibs2001SEmax", SE.factor = 1){
  
  if( (!is.null(forced.measures)) && (!is.numeric(forced.measures))){stop("Argument 'forced.measures' must be a numerical vector.")}
  if( (!is.null(forced.measures)) & (FALSE %in% (paste("m",forced.measures,sep="") %in% colnames(trajSelection$measures)))){stop("The measures in the argument 'forced.measures' must be measures present in the argument 'trajMeasures'.")}
  
  ID <- trajSelection$selection[,1]
  
  #standardize the measures to be clustered:
  if (is.null(forced.measures)) {
    data = data.frame(apply(data.frame(trajSelection$selection[, -1]), 2, scale))
  } else {
    data = data.frame(apply(data.frame(trajSelection$measures[, paste("m",forced.measures,sep="")]), 2, scale))
  }
  
  if (!is.null(nclusters) && (nclusters > nrow(data)) ){ stop("The number 'nclusters' of requested clusters cannot exceed the number of trajectories.") }
  
  if (is.null(nclusters)) {
    
    #the clusgap function asks us for a function that'll do the clustering but does not allow us to set the arguments of that function, so we have to do it beforehand, defining a new fct
    kmeans.nstart <- function(x, k){
      return(kmeans(x=x, centers=k, nstart = nstart, iter.max = iter.max))
    }
    
    GAP = cluster::clusGap(x=data, FUNcluster = kmeans.nstart, K.max = K.max, B = B, d.power = d.power, spaceH0 = spaceH0)
    nclusters <- maxSE(f = GAP$Tab[,"gap"], SE.f = GAP$Tab[,"SE.sim"], method = method, SE.factor = SE.factor)
    
    dev.off()
    par(mar=c(5.1, 5.1, 4.1, 2.1),lwd = 4, cex.axis=2,cex.lab=2,cex.main=2,bty="o",fg="black") ###set graphical parameters on the device. lwd is "line width" and controls the main curve of the graph (including distribution lines), cex.axis, cex.lab and cex.main controls width of the axis, labels and title resp. Default margins are mar=c(5.1, 4.1, 4.1, 2.1)
    plot(x=GAP, type = "b", xlab = "k", ylab = expression(Gap[k]),
         main = "", do.arrows = TRUE,
         arrowArgs = list(col="red3", length=1/16, angle=90, code=3))
    
    cluster.est2 <- kmeans(x=data, centers=nclusters, iter.max=iter.max, nstart=nstart)
    
    nb.clust <- nclusters
    
    partition <- cluster.est2$cluster #with k-means
    
    partition.summary <- summary(factor(partition))
    
  } else{
    
    cluster.est2 <- kmeans(x=data, centers=nclusters, iter.max=iter.max, nstart=nstart)
    
    nb.clust <- nclusters
    
    partition <- cluster.est2$cluster #with k-means
    
    partition.summary <- summary(factor(partition))
  }
  
  trajClusters <- structure(list(data=trajSelection$data, time=trajSelection$time, nb.clust = nb.clust, nclusters = nclusters, 
                                 partition = partition, partition.summary = partition.summary, 
                                 class = "trajClusters"))
  
  
  
  return(trajClusters)
  
} 
