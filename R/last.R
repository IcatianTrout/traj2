last <- function(v){ #returns the last non-NA coordinate of a vector
  if(!(FALSE %in% is.na(v))){stop("Argument must contain at least one non-NA entry.")}
  w <- v[complete.cases(v)]
  m<-length(w)
  return(w[m])
}
