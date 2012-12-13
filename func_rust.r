## ----------[ lodading necessary libraries ]--------------------
library(stats)
## library(pracma)                       #library for lsqnonlin
## library(Matrix)                       #matrix specific functions
library(matrixStats)                  #use rowSds from library
#library(pls)                            #use stdize
library(BB)                             #use for spg
library(expm)

## ----------[ Functions for least squares fitting ]----------------------------------------
##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title rust.fit.nStt
##' @param gData gene expression data
##' @param tData time points of data
##' @param lambda L1 penalty parameter (default = 0.01)
##' @param n.states number of states in the fitting
##' @param fit.as Fitting lin, log2Dat, logDat, log2Al (default = lin)
##' @return 
##' @author anas ahmad rana
rust.fit.nStt <- function(gData, tData, lambda = 0.01, n.states = 3, fit.as='lin'){
  x0 <- runif( (n.states - 1) + nrow(gData) * n.states)
  p <- nrow(gData)

  fun <- function(x){
    ## initialise w matrix and assign non zero values
    if(n.states==2){
      wFit <- matrix(0, n.states, n.states)
      wFit[2,1] <- x[1]
    } else if(n.states >=3) {
      wFit <- matrix(0, n.states, n.states)
      diag(wFit[-1, ]) <- x[1:(n.states - 1)]
    }

    ## Assign other x values to
    if(n.states==1){
      betaFit <- matrix( x, p, n.states)
      fit <- rust.kStt(wFit=NULL, betaFit, tData)
    } else{
      lnx <- length(x)
      betaFit <- matrix( x[-c(1:(n.states - 1))], {lnx - n.states + 1}/ n.states, n.states)
      fit <- rust.kStt(wFit, betaFit, tData)
    }

    if(fit.as=='lin'){
      rss <- (fit$S - gData)^2
    } else if(fit.as=='log2Dat'){
      rss <- (log2(fit$S) - gData)^2
    } else if(fit.as=='logDat'){
      rss <- (log(fit$S) - gData)^2
    } else if(fit.as=='log2Al'){
      rss <- (log2(fit$S) - log2(gData))^2
    } else {
      stop('String fit.as not recognised')
    }
    penalty <- lambda*sum(abs(betaFit)) 

    ss <- c(rss,penalty)
    sum(ss)
  }

  res <- nlminb(x0, fun, lower = 0, upper=max(tData)*2, control=list(iter.max = 3000,
                                                          eval.max=4000, rel.tol=10^-14))
  ## n.states=n.states, lambda=lambda, gData=gData, tData=tData) ##add if using external function
  rss <- res$objective - lambda*sum(res$par[-c(1:(n.states -1))])
  tmpPar <- res$par
  tmpPar[abs(tmpPar)<10^-3] <- 0
  n.zeros <- sum(tmpPar==0)
  Df <- p*n.states + n.states -1 - n.zeros
  n <- ncol(gData) * nrow(gData)

  bicSc <- n*log(rss/(n-1)) + log(n) * Df
  aicSc <- n*log(rss/(n-1)) + 2 * Df
  par <- rust.par(gData, res$par, n.states)
  
  return(list(fit=res, w=par$w, beta=par$beta, bic=bicSc, aic=aicSc))

}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title rust.par
##' @param gData data matrix used for naming beta
##' @param x parameter vector to reshape into beta and w
##' @param n.states number of states in model
##' @return 
##' @author anas ahmad rana
rust.par <- function(gData, x, n.states){
  ## W matrix from x[1:n-1]
  if(n.states==2){                      
    wFit <- matrix(0, n.states, n.states)
    wFit[2,1] <- x[1]
  } else if(n.states >=3) {
    wFit <- matrix(0, n.states, n.states)
    diag(wFit[-1, ]) <- x[1:(n.states - 1)]
  } else {
    wFit <- NULL
  }

  ## Assign other x values to beta
  if(n.states==1){
    betaFit <- matrix( x, p, n.states)
  } else{
    lnx <- length(x)
    betaFit <- matrix( x[-c(1:(n.states - 1))], {lnx - n.states + 1}/ n.states, n.states)
  }
  rownames(betaFit) <- rownames(gData)
  
  return(list(w=wFit, beta=betaFit))
}

## ----------[ RUST model (general) state indep ]--------------------
##' Function that takes as arguments the w_fit matrix and the beta_fit matrix, it then 
##'
##' .. content for \details{} ..
##' @title 
##' @param wFit W transition matrix
##' @param betaFit beta matrix p x T(n)
##' @param t time points
##' @return S unlogged trajectory for the arguments 
##' @return P state occupation probabilities
##' @author anas ahmad rana
rust.kStt <- function(wFit = NULL, betaFit, t){
  n <- length(t)
  p <- nrow(betaFit)
  ## create new matrix containing -offdiag on diag
  if(!is.null(wFit)){
    odiagElements <- diag(wFit[-1, , drop=FALSE ])
    diag(wFit) <- c(- odiagElements,0)

    k <- nrow(wFit)
    ## initial condition
    p0 <- rep(0,k); p0[1] <- 1
    P <- NULL

    for(i in 1:n){
      pt <- expm::expm(wFit *t[i], method='Ward77') %*% p0
      P <- cbind(P, pt)
      ## mean gene expression
      S <- betaFit %*% P
    }
  } else {
    S <- rep(betaFit, n)
    P <- 1
  }
  return(list(S=S, prob=P))
}

## ----------[ Miscelanious ]--------------------

##' Function to plot heatmap of a given set of names
##'
##' .. content for \details{} ..
##' @title 
##' @param names 
##' @param betas 
##' @return 
##' @author anas ahmad rana
betaHeatmaps <- function(names, betas){
  hBetas <- betas[names,]
  heatmapAnas(hBetas)
}



##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param g 
##' @param nClusters 
##' @param fName 
##' @return 
##' @author anas ahmad rana
rustClusters <- function(g, nClusters, fName){
  ## Perform kmeans clustering
  clKmeans <- kmeans(g, nClusters)
}



##' k-means clustering and plots the centroids
##'
##' 
##' @title 
##' @param time 
##' @param g microarray data columns = timepoints; rows = genes
##' @param k number of clusters
##' @return clst 
##' @author anas ahmad rana
rustClustering <- function(time, g, k){
  clst <- kmeans(g,k, nstart=1000)
  matplot(time, t(clst$centers), pch=1, t='b', ylab=NULL)
  
  return(clst)
  
}

  
plotFit4stt <- function(w,betaFit, tData, gData, gene){
  t <- seq(t[1],t[length(t)],length=300)
  fit4stt <- rust4sttFwd(w,betaFit[gene,],t)
  s <- fit4stt$S
  plot(t,log2(s), t='l', col='darkgreen', ylim=c(min(log2(s),gData[gene,]),max(log2(s),gData[gene,])),
       xlab='t [hoursf]', ylab='log_2 (expression)')
  points(tData, gData[gene,], col='blue', type='b')
  title(gene)
}

