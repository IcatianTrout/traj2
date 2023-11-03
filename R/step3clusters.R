#'@title Classify the Longitudinal Data Based on the Selected Measures.
#'@description Classifies the trajectories by applying the K-mean clustering
#'  algorithm to the measures selected by \code{Step2Selection}.
#'@param trajSelection Object of class \code{trajSelection} as returned by
#'  \code{Step2Selection}.
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
#'@return An object of class \code{mytrajClusters}; a list containing the result
#'  of the clustering as well as curated forms of the data and time matrices.
#'@import cluster
#'@importFrom stats kmeans na.omit
#'@importFrom grDevices palette.colors
#'@importFrom graphics legend lines par
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
Step3Clusters <- function (trajSelection, nstart = 50, iter.max = 20, nclusters = NULL, K.max = 8, B = 500, d.power = 2, spaceH0 = "scaledPCA", method = "Tibs2001SEmax", SE.factor = 1){
  
  nclusters.input <- nclusters
  
  ID <- trajSelection$selection[, 1]
  
  #standardize the measures to be clustered:
  
  data <-
    data.frame(apply(data.frame(trajSelection$selection[,-1]), 2, scale))
  
  if (!is.null(nclusters) && (nclusters > nrow(data))) {
    stop(
      "The number 'nclusters' of requested clusters cannot exceed the number of trajectories."
    )
  }
  if (is.null(nclusters)) {
    ## The clusGap function takes as argument a function that will perform
    ## the clustering but it does not allow us to set the arguments of that
    ## function, so we have to do so beforehand by defining a new function:
    kmeans.nstart <- function (x, k) {
      return(kmeans(
        x = x,
        centers = k,
        nstart = nstart,
        iter.max = iter.max
      ))
    }
    GAP <-
      cluster::clusGap(
        x = data,
        FUNcluster = kmeans.nstart,
        K.max = K.max,
        B = B,
        d.power = d.power,
        spaceH0 = spaceH0
      )
    nclusters <-
      maxSE(
        f = GAP$Tab[, "gap"],
        SE.f = GAP$Tab[, "SE.sim"],
        method = method,
        SE.factor = SE.factor
      )
    cluster.est2 <-
      kmeans(
        x = data,
        centers = nclusters,
        iter.max = iter.max,
        nstart = nstart
      )
    partition <- cluster.est2$cluster
    partition.summary <- summary(factor(partition))
  } else {
    GAP <- NULL
    cluster.est2 <-
      kmeans(
        x = data,
        centers = nclusters,
        iter.max = iter.max,
        nstart = nstart
      )
    partition <- cluster.est2$cluster
    partition.summary <- summary(factor(partition))
  }
  
  mytrajClusters <-
    structure(
      list(
        data = trajSelection$data,
        time = trajSelection$time,
        selection = trajSelection$selection,
        GAP = GAP,
        nclusters = nclusters,
        nclusters.input = nclusters.input,
        partition = partition,
        partition.summary = partition.summary,
        class = "mytrajClusters"
      )
    )
  
  return(mytrajClusters)
  
} 

#'@export
print.mytrajClusters <- function(mytrajClusters){
  
  cat(paste("The trajectories were grouped in ", mytrajClusters$nclusters, " clusters labeled ", paste(names(mytrajClusters$partition.summary), collapse = ", ", sep = ""), " of respective size ", paste(mytrajClusters$partition.summary, collapse=", ", sep = ""), ". The exact clustering is as follows.\n\n", sep=""))
  
  clust.by.id <- cbind(mytrajClusters$data[, 1],mytrajClusters$partition)
  colnames(clust.by.id) <- c("ID", "Cluster")
  print(clust.by.id)
}

#'@export
summary.mytrajClusters <- function(mytrajClusters) {
  cat("Distribution among clusters:\n")
  clust.dist <- data.frame(matrix(nrow = 2,ncol = (mytrajClusters$nclusters + 1)))
  clust.dist[1, ] <- signif(c(mytrajClusters$partition.summary, sum(mytrajClusters$partition.summary)))
  clust.dist[2, ] <- signif(c(mytrajClusters$partition.summary/sum(mytrajClusters$partition.summary), sum(mytrajClusters$partition.summary)/sum(mytrajClusters$partition.summary)),2)
  rownames(clust.dist) <- c("Absolute", "Relative")
  colnames(clust.dist) <- c(1:mytrajClusters$nclusters,"Total")
  print(clust.dist)
}

#'@export
plot.mytrajClusters <- function(mytrajClusters, sample.size = 5){
  
  # Plot the gap statistic
  if(!is.null(mytrajClusters$GAP)){
  plot(mytrajClusters$GAP, main="Gap statistic up to one SE")
  }
  
  # Plot sample.size randomly sampled trajectories from each clusters
  colors <- palette.colors(palette = "Okabe-Ito", alpha = 1)
  
  traj.by.clusters <- list()
  for(k in 1:mytrajClusters$nclusters){
    traj.by.clusters[[k]] <- mytrajClusters$data[which(mytrajClusters$partition == k), -c(1)]
  }
  
  time.by.clusters <- list()
  for(k in 1:mytrajClusters$nclusters){
    time.by.clusters[[k]] <- mytrajClusters$time[which(mytrajClusters$partition==k), -c(1)]
  }
  
  # Plot (max) sample.size random trajectories from each group
  smpl.traj.by.clusters <- list()
  smpl.time.by.clusters <- list()
  
  smpl.traj <- matrix(nrow = 0, ncol = ncol(mytrajClusters$data) - 1)
  smpl.time <- matrix(nrow = 0, ncol = ncol(mytrajClusters$time) - 1)
  
  for(k in 1:mytrajClusters$nclusters){
    size <- min(sample.size, nrow(traj.by.clusters[[k]]))
    smpl <- sample(x = c(1:nrow(traj.by.clusters[[k]])), size = size, replace=FALSE)
    smpl <- smpl[order(smpl)]
    
    smpl.traj.by.clusters[[k]] <- traj.by.clusters[[k]][smpl,]
    smpl.time.by.clusters[[k]] <- time.by.clusters[[k]][smpl,]
    
    smpl.traj <- rbind(smpl.traj, smpl.traj.by.clusters[[k]])
    smpl.time <- rbind(smpl.time, smpl.time.by.clusters[[k]])
  }
    
  plot(x = 0, y = 0, xlim = c(min(smpl.time, na.rm = T), max(smpl.time, na.rm = T)), ylim = c(min(smpl.traj, na.rm = T), max(smpl.traj, na.rm = T)), type = "n", xlab = "", ylab = "", main = "Sample trajectories")
  
  for(k in 1:mytrajClusters$nclusters){
    for(i in 1:size){
      lines(x = na.omit(smpl.time.by.clusters[[k]][i, ]), y = na.omit(smpl.traj.by.clusters[[k]][i, ]), type = "l", col = colors[k])
    }
    legend("topright",col=colors[1:k], legend=paste(1:3)[1:k], lty=rep(1,k))
  }
  
  
  # Plot dispersion plots of the selected measures
  if(ncol(mytrajClusters$selection) > 2){
    
    selection <- mytrajClusters$selection[, -c(1)]
    
    selection.by.clusters <- list()
    for(k in 1:mytrajClusters$nclusters){
      selection.by.clusters[[k]] <- selection[which(mytrajClusters$partition == k), ]
    }
  
    nb.measures <- ncol(selection) 
  
    for(m in 1:nb.measures){
      
      X <- sqrt(nb.measures - 1)
      
      if(X - floor(X)==0){
        par(mfrow = c(X,X))
      }
      
      if((X - floor(X) > 0) & (X - floor(X) < 0.5)){
        par(mfrow = c(X,X+1))
      }
      
      if(X - floor(X) >= 0.5){
        par(mfrow = c(X+1,X+1))
      }
      
      for(n in 1:(nb.measures - 1)){
        
        plot(x = 0, y = 0, xlim = c(min(selection[, m]), max(selection[, m])), ylim = c(min(selection[,-c(m)][, n]), max(selection[,-c(m)][, n])), type = "n", xlab = paste(colnames(selection[m])), ylab = paste(colnames(selection[,-c(m)])[n]), main = "")
        
        for(k in 1:mytrajClusters$nclusters){
          
          lines(x = selection.by.clusters[[k]][, m], y = selection.by.clusters[[k]][, -c(m)][, n], type = "p", pch = 20, col = colors[k], bg = colors[k])
          
          legend("topright", lty=rep(0,k), pch = rep(16,k), col=colors[1:k], legend=paste(1:3)[1:k])
        }
      }
    }
  }
}

