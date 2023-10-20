#'@title Select a Subset of the Measures Using Factor Analysis
#'@description Applies the following dimension reduction algorithm to the
#'  measures computed by \code{\link[traj2]{Step1Measures}}:
#' \enumerate{
#'   \item Use principal component analysis (PCA) on the measures;
#'   \item Drop the factors whose variance is less than that of any individual measure;
#'   \item Performs a varimax rotation on the remaining factors;
#'   \item For each rotated factor, select the measure that has the highest correlation (aka factor loading) with it and that hasn't yet been selected;
#'   \item Drop the remaining measures.
#' }
#'@param trajMeasures Object of class \code{trajMeasures} as returned by
#'  \code{\link[traj2]{Step1Measures}}.
#'@param num.select Defaults to NULL. If a numeric positive integer is supplied,
#'  then in step 2 of the algorithm, the \code{num.select} factors with the
#'  highest variance are kept and the others are dropped.
#'@param discard Defaults to NULL. If a numeric vector of positive integers is
#'  supplied, then the corresponding measures will be dropped. See
#'  \code{\link[traj2]{Step1Measures}} for the list of measures.
#'@param select Can be NULL or a numeric vector corresponding to the numerical
#'  identifier of measures present in \code{trajMeasures}. If a numeric vector
#'  is supplied, then four-steps selection algorithm described above is bypassed
#'  and the corresponding measures are selected instead.
#'@return An object of class \code{trajSelection}; a list containing the values
#'  of the selected measures, the loadings of the rotated factors on each
#'  measure as well as curated forms of the data and time matrices.
#'@importFrom psych principal
#'@importFrom stats cor
#'
#'@details Prior to using PCA on the data, if a measure is perfectly or almost
#'perfectly correlated (corr. >0.98) with some other measure that appears before
#'it in the list, then it is discarded. Likewise, if a measure is constant, it
#'is discarded. Quotient measures which turn out to be of the form 0/0 are set
#'to 1.
#'
#'By virtue of being quotients, some measures (4, 7, 8, 15-17, 21-26) can turn
#'out extremely large or even infinite (division by 0). Because the K-means
#'algorithm is sensitive to outliers, it is necessary to act on them. Therefor,
#'if a measure takes values beyond the 0.3% probability threshold as computed
#'from Nishiyama's improved Chebychev bound, those values are imputed to the
#'0.3% probability threshold.
#'
#'The function \code{\link[psych]{principal}} from the \code{psych} package is
#'used for PCA.
#'
#'@references Leffondre K, Abrahamowicz M, Regeasse A, Hawker GA, Badley EM,
#'McCusker J, Belzile E. Statistical measures were proposed for identifying
#'longitudinal patterns of change in quantitative health indicators. J Clin
#'Epidemiol. 2004 Oct;57(10):1049-62. doi: 10.1016/j.jclinepi.2004.02.012. PMID:
#'15528056. Nishiyama T, Improved Chebyshev inequality: new probability bounds
#'with known supremum of PDF, arXiv:1808.10770v2 stat.ME
#'https://doi.org/10.48550/arXiv.1808.10770
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
#'}
#'
#'
#'@seealso \code{\link[psych]{principal}} \code{\link[traj2]{Step1Measures}}
#'
#'@rdname Step2Selection
#'
#'@export
Step2Selection <- function (trajMeasures, num.select = NULL, discard = NULL, select = NULL) {
  
  if ( (!is.null(select)) & (!is.numeric(select))) {
    stop("Argument 'select' must be either NULL or a numerical vector.")
  }
  
  data <- data.frame(trajMeasures$measures)
  ID <- data[, 1]
  data <- data.frame(data[, -1])
  
  if (!is.null(select)) {
    m.select <- paste("m", select, sep = "")
    if (FALSE %in% (m.select %in% colnames(data))) {
      stop("Select 'select' from the measures included in step1measure.")
    } else { 
      output <- cbind(ID,data[, m.select]) 
      colnames(output) <- c("ID", paste("m", select, sep = "")) 
    }
    
    trajSelection <- structure(list(selection = output, loadings = NULL, measures = trajMeasures$measures, 
                                   data = trajMeasures$data, time = trajMeasures$time), 
                              class = "trajSelection")
  } else {
    
    if ((!is.null(discard)) & (!is.numeric(discard))) {
      stop("Argument 'discard' must be either NULL or a numerical vector.")
    }
    if (!is.null(discard)) {
      mes.to.discard <- paste("m", discard, sep = "")
      if ( FALSE %in% (mes.to.discard %in% colnames(data))) {
        stop("Can't discard a measure which was not included in step1measure.")
      }
      w <- which( colnames(data) %in% mes.to.discard)
      data <- data[, -w]
    }
    
    if (!is.null(num.select)) {
      if (!is.numeric(num.select)) {
        stop("Argument 'num.select' must be a numerical vector of length 1.")
      }
      if (!is.vector(num.select)) {
        stop("Argument 'num.select' must be a numerical vector of length 1.")
      }
      if (!(length(num.select)==1) ) {
        stop("Argument 'num.select' must be a numerical vector of length 1.")
      }
      if (!is.null(discard) & (num.select > ncol(data)) ) {
        stop("After discarding the measures specified in 'discard', the requested number 'num.select' of measures to retain exceeds the number of measures available.")
      }
      if (is.null(discard) & (num.select > ncol(data)) ) {
        stop("The requested number 'num.select' of measures to retain exceeds the number of measures included in step1measure.")
      }
    }
    
    corr.vars <- CheckCorrelation(data, verbose = FALSE, is.return = TRUE)
    
    if (!is.null(corr.vars)) {
      corr.vars.pos <- which(names(data) %in% corr.vars[, 1])
      data <- data[, -corr.vars.pos]
      print(paste(corr.vars[, 1], "is discarded because it is perfectly or almost perfectly correlated with", corr.vars[, 2]))
    }
    
    if (num.select > ncol(data) && !is.null(num.select)) {
      stop("After discarding the perfectly or almost perfectly correlated measures, there are ", ncol(data), " measures left, which is less than the number 'num.select' of requested measures to retain.")
    }
    
    Z <- scale(x = data, center = TRUE, scale = TRUE) 
    
    if (is.null(num.select)) {
      eigen.values <- psych::principal(Z, rotate = "none", nfactors = ncol(data))$values
      num.select <- max(1, length(which(eigen.values > 1)))
    }
    
    if (ncol(data) > 1) {
      principal.factors <- psych::principal(Z, rotate = "varimax", nfactors = num.select)
      
      principal.variables <- c()
      
      if (is.null(select)) {
        for (j in 1:num.select) {
          if (j == 1){ 
            aux <- principal.factors$loadings 
          }
          if (j > 1){
            w <- which(row.names(principal.factors$loadings) %in% principal.variables)
            aux <- principal.factors$loadings[-w, ]
          }
          principal.variables[j] <- names(which.max(abs(aux[, j])))
        }
        output <- cbind(ID,data[, principal.variables])
      }
    } else {
      output <- trajMeasures$measures
    }
    
    trajSelection = structure(list(selection = output, loadings = principal.factors$loadings, measures = trajMeasures$measures, 
                                   data = trajMeasures$data, time = trajMeasures$time), class = "trajSelection")
  }
  
  return(trajSelection)
  
}
