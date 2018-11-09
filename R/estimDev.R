
#' Estimation of the minimal deviations from rational expectations with unconstrained information set g*
#'
#' This function estimates of the minimal deviations from rational expectations with unconstrained information set. Both vectors should have the same length. If not, one can randomly select a subset of the longer vector with length equal to that of the shorter one. The function returns a function via the approxfun of the package stats. This function can then be evaluated directly on a desired grid.
#' @param psi vector of subjective expectations
#' @param y vector of realisations of an individual outcome.
#' @return
#' @export
#'
#' @examples
estimDev <- function(psi,y ){
alpha1 <- mean(y)-mean(psi)
y  <- y - alpha1
psi <- sort(psi,decreasing = TRUE )
y <- sort(y,decreasing = TRUE  )
g <-matrix(0,1, length(psi))
n_p <- length(psi)
length(psi)
g_0 <- g
i_t = 0
i_max = n_p
library(matlab)

t=1
for(t in 2:i_max){
  res <- NULL
  for(j in 1:n_p ){
    res <- c(  res, c_fun(j, i_t, y ,psi) )
  }
  abs1 <- which.min2(res[(i_t+1):n_p], na.rm=TRUE) +i_t
  if(abs(abs1)!=Inf){
    for(j in (i_t+1):abs1 ){
      g[1,j] <- psi[j] + c_fun(abs1, i_t, y ,psi)
    }
    i_t <-  abs1
  }else{break;
  }
  if( i_t==n_p){
    break;
  }
  if(t%%1000==0){
    delta <- sum(sum(abs(g-g_0)))
    g_0 <- g
    #     cat("iteration", t," delta: " ,delta , " \n ")
  }
}

t<- seq(-5,5, length.out=300)
# x11()
g <- as.vector(g)



start <-  sort(g)
for (i in 1:length(g)){
  start[i] <-max(start[1:i])
}
start[is.na(start)]<-max(start, rm.na=T)
# plot(psi,sort(start,decreasing = T),type="l")
# lines(psi,psi,col=2)
g <- sort(start,decreasing = T)

mean(g) - mean(psi)
g_star  <- approxfun(psi, g,method = "linear")

return(g_star)
}
