#' @title Compute Measures for Identifying Patterns of Change in Functional Data
#' @description \code{step1measures} computes up to 26 measures for each function. See details for the list of measures.
#' @param Data A matrix or data frame in which each row contains the (ordered) observations of a given function.
#' @param Time Either NULL or a vector or a matrix (or data frame) of the same dimension as \code{Data}. If a vector, matrix or data frame is supplied, it is assumed that the entries correspond to the values of the arguments at which the functions in \code{Data} are observed. When set to \code{NULL} (the default), the observations are assumed equidistant.
#' @param ID Logical. Set to \code{TRUE} if the first column of \code{Data} and \code{Time} corresponds to an \code{ID} variable. Defaults to \code{FALSE}.
#' @param measures A numerical vector containing the measures to compute. The default, c(1:23), leaves out the measures which require specifying a midpoint. 
#' @param midpoint Specifies which column of \code{Time} to use as the midpoint in measures 24-26. Can be an integer, an integer vector (of length the number of rows in \code{Time}) or NULL (the default). If NULL, the midpoint for a given function is taken to be the time closest to the median.
#' @return An object of class \code{trajMeasures}; a list containing the values of the measures, a table of outliers which have been imputed, as well as a curated form of the \code{Data} and \code{Time} arguments.
#' @details 
#' Each function must have a minimum of 3 observations otherwise it will be omitted from the analysis. 
#'
#'The 26 measures are:
#'\enumerate{
#'\item  Range\cr
#'\item  Mean of the function\cr
#'\item  Functional standard deviation (SD)\cr
#'\item  Coefficient of variation (ratio of measure 3 to measure 2)\cr
#'\item  Overall change (initial value - final value)\cr
#'\item  Mean change per unit time\cr
#'\item  Overall change relative to initial value\cr
#'\item  Overall change relative to functional mean (ratio of measure 5 to measure 5)\cr
#'\item  Slope of the linear model\cr
#'\item  R^2: Proportion of variance explained by the linear model\cr
#'\item  Maximum value of the speed\cr
#'\item  Functional SD of the speed\cr
#'\item  Mean absolute speed\cr
#'\item  Maximum absolute speed\cr
#'\item  Maximum absolute speed relative to the functional mean (ratio of measure 14 to measure 2)\cr
#'\item  Maximum absolute speed relative to the slope (ratio of measure 14 to measure 9)\cr
#'\item  Functional SD of the speed relative to the slope\cr(ratio of measure 12 to measure 9)
#'\item  Mean acceleration\cr
#'\item  Mean absolute acceleration\cr
#'\item  Maximum of the absolute acceleration\cr
#'\item  Maximum of the absolute acceleration relative to the functional (ratio of measure 20 to measure 2)\cr
#'\item  Maximum of the absolute acceleration relative to the mean absolute speed (ratio of measure 20 to measure 13)\cr
#'\item  Mean absolute acceleration relative to the mean absolute speed (ratio of measure 19 to measure 13)\cr
#'\item  Early change relative to later change\cr
#'\item  Early change relative to overall change\cr
#'\item  Later change relative to overall change\cr
#'}
#'
#' @references 
#' Leffondre K, Abrahamowicz M, Regeasse A, Hawker GA, Badley EM, McCusker J, Belzile E. Statistical measures were proposed for identifying longitudinal patterns of change in quantitative health indicators. J Clin Epidemiol. 2004 Oct;57(10):1049-62. doi: 10.1016/j.jclinepi.2004.02.012. PMID: 15528056.
#'
#'@examples 
#'\dontrun{
#'data = example.data$data
#'
#'m1 = step1measures(data, ID=TRUE, measures=24:26, midpoint=NULL)
#'m2 = step1measures(data, ID=TRUE, measures=24:26, midpoint=as.integer(3))
#'
#'identical(s1$measures,s2$measures)
#'}
#' 
#' @rdname step1measures
#' @export  
step1measures <- function (Data, Time=NULL, ID = FALSE, measures=1:23, midpoint=NULL){
  
  if(is.null(Time)){
    if(ID==TRUE){Time <- c(1:(ncol(Data)-1))}
    if(ID==FALSE){Time <- c(1:(ncol(Data)))}
  }
  
  Time.is.vector <- is.vector(Time)
  
  data <- Data #This line will prompt an error message if the 'Data' argument is missing. 
  data <- data.frame(data)
  names(data) <- NULL
  
  time <- Time #This line will prompt an error message if the 'Time' argument is missing
  time <- data.frame(time) 
  names(time) <- NULL
  
  #prep data:
  
  #remove ID column (if applicable):
  if (ID == TRUE) { 
    IDvector = data[, 1]
    data = data[, -1]
    if(Time.is.vector==FALSE){  
      ID.time = time[, 1]
      time = time[, -1]
      if(identical(IDvector,ID.time)==FALSE){stop("ID vector in Data differs from ID vector in Time.")}
    }
  } else{IDvector <- c(1:nrow(data))} #if ID == FALSE, create an IDvector.
  
  #check that the structure of Data and Time is as expected and, provided Time is supplied as a vector, manipulate data and time so that they are in the same format as if Time is supplied as a data frame.
  if(identical(dim(time),dim(data))){ #this signifies that we're in the generic situation reg. Data and Time
    #...(check that the "NA structures" are compatible and as expected and check that the times are never decreasing
    if(identical(is.na(data),is.na(time))==FALSE){stop("Data has cells with no corresponding time or Time has cells with no corresponding data.")}
    
    #remove any line that contains less than 3 observations and issue a warning that this has been done:
    data2 <- data
    time2 <- time
    for(i in 1:nrow(data)){
      NA.str_i <- is.na(data)[i,]
      w <- unname(which(NA.str_i==FALSE))
      if(length(w)<3){
        data2 <- data2[-i,] ; time2 <- time2[-i,]; warning(paste("Row ",i," of Data contains less than 3 observations; it has been removed.",sep="")); IDvector <- IDvector[-i]
      } else if(!identical(w,c(1:length(w)))){stop(paste("Row ",i," of Data is not formatted correctly. Rows should be of the form X Y ... Z NA ... NA.",sep=""))}
    }
    data <- data2
    time <- time2
    
    for(i in 1:nrow(data)){
      w2 <- unname(which(!is.na(time)[i,]))
      non.NA <- unname(unlist(time[i,w2]))
      if(!length(unique(non.NA))==length(non.NA)){stop(paste("Line ",i," of Time contains duplicates. The rows of Time should be strictly increasing sequences.",sep=""))}
      if(!identical(w2,order(non.NA))){stop(paste("The elements in row ",i," of Time are not ordered chronologically.",sep=""))}
    }
    
  } else if ( is.vector(Time) & (length(Time)==ncol(data)) ) { #this signifies that we're in the equal time measurements situation
    #check that Time contains no NA and that the times are strictly increasing:
    time <- unname(unlist(Time))
    if(TRUE %in% is.na(time)){stop("If Time is supplied as a vector, it must not contain NA.")}
    if(!length(unique(time))==length(time)){stop("The Time vector contains duplicates. The elements of Time should form a strictly increasing sequence.")}
    if(!identical(order(time),c(1:length(time)))){stop("The elements of Time are not ordered chronologically.")}
    if(length(time)<3){stop("Time must contain at least 3 elements.")}
    
    #remove any line of Data that's mono-NA and issue a warning that this has been done:
    data2 <- data
    for(i in 1:nrow(data)){
      NA.str_i <- is.na(data)[i,]
      w <- unname(which(NA.str_i==FALSE))
      if(length(w)<3){
        data2 <- data2[-i,];  warning(paste("Row ",i," of Data contains less than 3 observations; it has been removed.",sep="")); if(ID == TRUE){IDvector <- IDvector[-i]}
      }
    }
    data <- data2
    
    #Change the format of data and time so that it fits the case where Time is supplied as a data frame. 
    data2 <- unname(data.frame(matrix(NA,ncol=ncol(data),nrow=nrow(data))))
    time2 <- unname(data.frame(matrix(NA,ncol=ncol(data),nrow=nrow(data))))
    
    for(i in 1:nrow(data)){
      w <- which(!is.na(data[i,]))
      data2[i,c(1:length(w))] <- data[i,w]
      time2[i,c(1:length(w))] <- Time[w]
    }
    data <- data2
    time <- time2
    
  } else(stop("The dimension of Data and Time are incompatible."))
  
  #From here on out we have data and time, both in the same format.
  
  #if some measures contrasting early vs. later change are included, we must construct a vector of "mid points" based on the argument midpoint
  if(TRUE %in% (c(25:27) %in% measures)){
    
    mid.position <- c()
    
    if(is.null(midpoint)){
      median.time <- (max(time,na.rm=T) + min(time,na.rm=T))/2
      flag <- c()
      for(i in 1:nrow(data)){
        v <- time[i,][!is.na(time[i,])] #remove NAs from time[i,]
        w <- which.min(abs(v - median.time))
        mid.position[i] <- w
        if(w==length(v) | w==1){ flag <- c(flag,i) } #the mid point can't be the first or the last point, so if this occurs, the ith row gets "flagged" 
      }
      if(length(flag)>0){
        mid.position <- mid.position[-flag]
        data <- data[-flag,]
        time <- time[-flag,]
        if(ID==TRUE){Lines <- which(Data[,1]%in%IDvector[flag])} else{Lines <- IDvector[flag]}
        if(length(Lines)==1){ warning(paste("When left blank, the 'midpoint' argument defaults to the observation time closests to the median time 0.5*(max(Time)-min(Time)), but this can't be either the first or last observation time. As a result, row ",Lines," has been removed. To avoid this, consider excluding measures 24-26 from the analysis or providing custom 'midpoint' values.",sep=""))}
        if(length(Lines)>1){ warning(paste("When left blank, the 'midpoint' argument defaults to the observation time closests to the median time 0.5*(max(Time)-min(Time)), but this can't be either the first or last observation time. As a result, rows ",noquote(paste(Lines,collapse=", "))," have been removed. To avoid this, consider excluding measures 24-26 from the analysis or providing custom 'midpoint' values.",sep="")) } 
        IDvector <- IDvector[-flag]
      }
    } else if( !(is.vector(midpoint) & is.integer(midpoint))){ 
      stop("'midpoint' should be either NULL, an integer or an integer vector.")
    } else{
      if(length(midpoint)==nrow(data)){
        mid.position <- midpoint
      } else if(is.vector(Time) & (length(Time)==ncol(data)) & length(midpoint)==1){
        mid.position <- rep(midpoint,nrow(data))
      } else{stop("'midpoint' does not have the correct format.")}
    } 
    #test to see if midpoints are all greater than 1 but lesser than the number of observations
    for(i in 1:nrow(data)){
      v <- time[i,][!is.na(time[i,])] #remove NAs from time[i,]
      if(mid.position[i]<=1){stop(paste("Error in 'midpoint' for subject ",i,"; 'midpoint' must be greater than 1.",sep=""))}
      else if (mid.position[i]>=length(v)){stop(paste("Error in 'midpoint' for subject ",i,"; 'midpoint' must be less than the number of observations.",sep=""))}
    }
  }
  
  
  output = data.frame(matrix(ncol = 1+length(measures), nrow = nrow(data))) 
  colnames(output) <- c("ID", paste("m",measures,sep=""))
  output$ID <- IDvector
  
  data <- as.matrix(data) #this is so that lines are treated as vectors instead of data frames.
  time <- as.matrix(time) #this is so that lines are treated as vectors instead of data frames.
  
  
  #compute the measures of change requested:
  if( (1 %in% measures)){
    for (i in 1:nrow(data)) {
      output$m1[i] = max(data[i, ], na.rm = TRUE) - min(data[i, 
      ], na.rm = TRUE)
    }}
  
  if(sum(c(2,3,4,8,15,21) %in% measures)>0){
    m2 <- c()
    for (i in 1:nrow(data)) {
      y <- data[i,complete.cases(data[i,])]
      x <- time[i,complete.cases(time[i,])]
      m2[i] = fct.mean(x,y)
    }
    if(2 %in% measures){output$m2 <- m2}
  }
  
  if(sum(c(3,4) %in% measures)>0){
    m3 <- c()
    for (i in 1:nrow(data)) {
      y <- data[i,complete.cases(data[i,])]
      x <- time[i,complete.cases(time[i,])]
      m3[i] = sqrt(fct.mean(x,(y-m2[i])^2))
    }
    if(3 %in% measures){output$m3 <- m3}
  }
  
  if((4 %in% measures)){
    for (i in 1:nrow(data)) {
      output$m4[i] = output$m3[i]/output$m2[i]
    }}
  
  if(sum(c(5,6,7,8,25,26) %in% measures)>0){
    m5 <- c()
    for (i in 1:nrow(data)) {
      m5[i] = last(data[i, ]) - first(data[i,]) 
    }
    if(5 %in% measures){output$m5 <- m5}
  }
  
  if((6 %in% measures)){
    for (i in 1:nrow(data)) {
      output$m6[i] = m5[i]/(last(time[i, ]) - first(time[i, ])) #that +1 makes NO SENSE ! I removed it...
    }}
  
  if((7 %in% measures)){
    for (i in 1:nrow(data)) {
      output$m7[i] = m5[i]/first(data[i, ])
    }}
  
  if((8 %in% measures)){
    for (i in 1:nrow(data)) {
      output$m8[i] = m5[i]/m2[i]
    }}
  
  if(sum(c(9,16,17) %in% measures)>0){
    m9 <- c()
    for (i in 1:nrow(data)) {
      b = coefficients(lm(data[i, ] ~ time[i, ]))
      m9[i] = b[2]
    }
    if(9 %in% measures){output$m9 <- m9}
  }
  
  defaultW <- getOption("warn")
  options(warn = -1)  #turn off warnings emited from summary(model) when the residuals are 0, which can happen when the trajectory is a straight line
  if((10 %in% measures)){
    for (i in 1:nrow(data)) {
      model = lm(data[i, ] ~ time[i, ])
      output$m10[i] = summary(model)$r.squared
    }}
  options(warn = defaultW)
  
  if((11 %in% measures)){
    for (i in 1:nrow(data)) {
      y <- data[i,complete.cases(data[i,])]
      x <- time[i,complete.cases(time[i,])]
      D <- der(x,y)
      output$m11[i] = max(D)
    }}
  
  if(sum(c(12,17) %in% measures)>0){
    m12 <- c()
    for (i in 1:nrow(data)) {
      y <- data[i,complete.cases(data[i,])]
      x <- time[i,complete.cases(time[i,])]
      D <- der(x,y)
      m12[i] = sqrt(fct.mean(x,(D-fct.mean(x,D))^2))
    }
    if(12 %in% measures){output$m12 <- m12}
  }
  
  if(sum(c(13,22,23) %in% measures)>0){
    m13 <- c()
    for (i in 1:nrow(data)) {
      y <- data[i,complete.cases(data[i,])]
      x <- time[i,complete.cases(time[i,])]
      D <- der(x,y)
      m13[i] = fct.mean(x,abs(D))
    }
    if(13 %in% measures){output$m13 <- m13}
  }
  
  if(sum(c(14,15,16) %in% measures)>0){
    m14 <- c()
    for (i in 1:nrow(data)) {
      y <- data[i,complete.cases(data[i,])]
      x <- time[i,complete.cases(time[i,])]
      D <- der(x,y)
      m14[i] = max(abs(D))
    }
    if(14 %in% measures){output$m14 <- m14}
  }
  
  if((15 %in% measures)){
    for (i in 1:nrow(data)) {
      output$m15[i] = m14[i]/m2[i]
    }}
  
  if((16 %in% measures)){
    for (i in 1:nrow(data)) {
      output$m16[i] = m14[i]/m9[i]
    }}
  
  if((17 %in% measures)){
    for (i in 1:nrow(data)) {
      output$m17[i] = m12[i]/m9[i]
    }}
  
  if((18 %in% measures)){
    for (i in 1:nrow(data)) {
      y <- data[i,complete.cases(data[i,])]
      x.i <- time[i,complete.cases(time[i,])]
      D <- der(x.i,y)
      A <- der(x.i,D)
      output$m18[i] = fct.mean(x.i,A)
    }}
  
  if(sum(c(19,23) %in% measures)>0){
    m19 <- c()
    for (i in 1:nrow(data)) {
      y <- data[i,complete.cases(data[i,])]
      x.i <- time[i,complete.cases(time[i,])]
      D <- der(x.i,y)
      A <- der(x.i,D)
      m19[i] = fct.mean(x.i,abs(A))
    }
    if(19 %in% measures){output$m19 <- m19}
  }
  
  if(sum(c(20,21,22) %in% measures)>0){
    m20 <- c()
    for (i in 1:nrow(data)) {
      y <- data[i,complete.cases(data[i,])]
      x.i <- time[i,complete.cases(time[i,])]
      D <- der(x.i,y)
      A <- der(x.i,D)
      m20[i] = max(abs(A))
    }
    if(20 %in% measures){output$m20 <- m20}
  }
  
  if((21 %in% measures)){
    for (i in 1:nrow(data)) {
      output$m21[i] = m20[i]/m2[i]
    }}
  
  if((22 %in% measures)){
    for (i in 1:nrow(data)) {
      output$m22[i] = m20[i]/m13[i]
    }}
  
  if((23 %in% measures)){
    for (i in 1:nrow(data)) {
      output$m23[i] = m19[i]/m13[i]
    }}
  
  if(TRUE %in% (c(24:26) %in% measures)){
    EC <- c() #early change
    for (i in 1:nrow(data)) {
      EC[i] <- last(data[i, 1:mid.position[i]]) - first(data[i, 1:mid.position[i]])
    }
    
    LC <- c() #later change
    for (i in 1:nrow(data)) {
      LC[i] <- last(data[i, mid.position[i]:ncol(data)]) - first(data[i, mid.position[i]:ncol(data)
      ])
    }
  }
  
  if((24 %in% measures)){
    for (i in 1:nrow(data)) {
      output$m24[i] = EC[i]/LC[i]
    }}
  
  if((25 %in% measures)){
    for (i in 1:nrow(data)) {
      output$m25[i] = EC[i]/m5[i]
    }}
  
  if((26 %in% measures)){
    for (i in 1:nrow(data)) {
      output$m26[i] = LC[i]/m5[i]
    }}
  
  output[is.na(output)] <- 1 #NAs correspond to quotient measures of the form 0/0
  
  
  #remove the measures that are constant because (1) these are not useful for discriminating between the trajectories (2) they cause problem with the check.correlation.b function below because their SD is 0 so division by 0 occurs when computing correlation.
  flag <- c()
  for(j in 1:ncol(output)){
    col <- output[,j]
    if(max(col)==min(col)){flag <- c(flag,j)}
  }
  if(length(flag)>0){ 
    output <- output[,-flag]
    warning(paste("Measures ",paste(colnames(output)[3:6],collapse=", "), " have been discarded due to being constant.",sep=""))
  }
  
  if(ncol(output)==1){stop("All the measures are constant.")}
  
  #Now let's deal with the outliers.
  
  outliers <- output
  outliers <- data.frame(matrix(NA,nrow=nrow(output),ncol=ncol(output)))
  colnames(outliers) <- colnames(output)
  outliers[, 1] <- output$ID
  
  for(j in 2:ncol(output)){
    y <- output[,j]
    
    y.TRUE <- y
    n <- length(y)
    
    which.inf <- which(is.infinite(y.TRUE))
    if(length(which.inf)>0){y.TRUE <- y.TRUE[-which.inf]}
    
    if(length(y.TRUE)>2){
    top <- rev(order(abs(y.TRUE-median(y.TRUE))))[1:ceiling(n*0.01)] #if n<100, remove 1 element, so this is never empty
    y.TRUE <- y.TRUE[-top]
    }
    
    mu <- mean(y.TRUE)
    sigma <- sd(y.TRUE)
    
    #In analogy with the 3sigma rules for the normal distribution, we wish to consider as outlier an element that, according to (the finite sample version of) Chebychev, has probability <0.3% of being observed
    k_Cheb <- sqrt(100/0.3) #approximately 18.26
    k <- seq(from=0.1,to=18.26,by=0.1)
    M <- c()
    for(i in 1:length(k)){
      max.left <- max(density(y.TRUE,from=(mu-18.3*sigma),to=(mu+18.3*sigma))$y[density(y.TRUE,from=(mu-18.3*sigma),to=(mu+18.3*sigma))$x > mu+k[i]*sigma])
      max.right <- max(density(y.TRUE,from=(mu-18.3*sigma),to=(mu+18.3*sigma))$y[density(y.TRUE,from=(mu-18.3*sigma),to=(mu+18.3*sigma))$x < mu-k[i]*sigma])
      M[i] <- max(max.left,max.right)
    }
    
    p <- 2*pi*k^2*(exp(1)-2*pi/3)
    q <- 2*(2*pi*k)^3/27 - (2*pi)^2*k^3*exp(1)/3 - 2*pi*exp(1)/(sigma*M)
    
    root <- CubeRoot(-q/2+sqrt(q^2/4+p^3/27)) + CubeRoot(-q/2-sqrt(q^2/4+p^3/27))+2*pi*k/3
    
    w <- which(root*sigma*M<.3/100)
    if(length(w)>0){k.opt <- min(k_Cheb,k[w[1]])} else{k.opt <- k_Cheb}
    
    cap <- which(abs(y-mu) > k.opt*sigma)
    outliers[cap,j] <- y[cap]
    
    y[cap] <- mu + sign(y[cap])*k.opt*sigma
    output[,j] <- y
    
    y[cap]
    if(length(cap)==1){warning(paste("In measure ",colnames(output)[j], ", ", length(cap), " outlier has been imputed to ",max(abs(round(mu, 2)), abs(signif(mu, 1))),"±",max(abs(round(sigma, 2)), abs(signif(sigma, 1))),"*",k.opt, " = ",max(abs(round(mu-sigma*k.opt, 2)), abs(signif(mu-sigma*k.opt, 1)))," or ",max(abs(round(mu+sigma*k.opt, 2)), abs(signif(mu+sigma*k.opt, 1))),".",sep=""))}
    if(length(cap)>1){warning(paste("In measure ",colnames(output)[j], ", ", length(cap), " outliers have been imputed to ",max(abs(round(mu, 2)), abs(signif(mu, 1))),"±",max(abs(round(sigma, 2)), abs(signif(sigma, 1))),"*",k.opt, " = ",max(abs(round(mu-sigma*k.opt, 2)), abs(signif(mu-sigma*k.opt, 1)))," or ",max(abs(round(mu+sigma*k.opt, 2)), abs(signif(mu+sigma*k.opt, 1))),".",sep=""))}
  }
  
  row.rm <- which(rowSums(!is.na(outliers[,-c(1)]))==0)
  outliers <- outliers[-row.rm,]
  col.rm <- which(colSums(!is.na(outliers))==0)
  outliers <- outliers[, -col.rm]
  
  ID <- IDvector

  trajMeasures = structure(list(measures = output, imputed = outliers, data = cbind(ID,data), time = cbind(ID,time)), class = "trajMeasures")
  return(trajMeasures)
 
}
