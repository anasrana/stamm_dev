## ----------[ lodading necessary libraries ]--------------------
library(stats)
library(matrixStats)                  #use rowSds from library
library(expm)                           #Matrix exponential also loads library(Matrix)
library(multicore)

##' Function to fit w for centroids
##'
##' .. content for \details{} ..
##' @title 
##' @param gData 
##' @param tData 
##' @param lambda 
##' @param n.states 
##' @param fit.as 
##' @param rSmpls 
##' @return 
##' @author anas ahmad rana
rust.clst.fit <- function(gData, tData, lambda, n.states, fit.as='log2Dat', rSmpls){
  
  cl <- rust.clustering(gData, km.k=6, rSmpl.size=rSmpls, n.randSmpl=10)
  cl.fit <- rust.fit.nStt(gData=cl$clst$centers, tData=tData, lambda=lambda, n.states=n.states, fit.as=fit.as)
  bic <- cl.fit$bic
  aic <- cl.fit$aic
  rss <- cl.fit$rss

  w <- cl.fit$w
  beta.centroid <- cl.fit$beta
  print(paste('w results',w))
  print('now fitting per gene')
  for(i in 1:length(cl$rnd.sample)){
    smpl <- cl$rnd.sample[[i]]
    for(j in 1:ncol(smpl)){
      gVec <- as.list(smpl[,j])
      browser()
      fit.gnes <- mclapply(gVec, function(x)
           cl.fit = rust.fit.nStt(gData=gData[x,], tData=tData, lambda=lambda,
             n.states=n.states, fit.as=fit.as, w=w, fix.w=TRUE))
    }
  }

  return(fit.gnes)
}



##' Clusters genes using k-means and returns centers, cluster rep genes, and random samples from genes
##'
##' .. content for \details{} ..
##' @title rust.clustering
##' @param gData time course of data
##' @param km.init number of intialisations for k-means (default value is 100)
##' @param km.k number of states for clustering
##' @param rSmpl.size sizes for random samples, if NULL (default) nothing returned
##' @param n.randSmpl number of random samples for each size (default is 50)
##' @return clst k-means clustering output important variable is $cluster and $centers
##' @return rep.gns genes closest to the centroids
##' @return rnd.sample random samples from gene list without rep.gns
##' @author anas ahmad rana
rust.clustering <- function(gData, km.init=100, km.k=NULL, rSmpl.size=NULL, n.randSmpl=50){
  if(is.null(km.k)){
    stop('choose number of states km.k') #check in number of clusters is given STOP if not
  }

  ## Initialise random sample list
  km.rnd <- NULL

  ## k-means cluster all the genes from the data set
  tmp.obj <- 10^30
  for(i.init in 1:km.init){
    tmp <- kmeans(x=gData, centers=km.k, iter.max=100)
    tmp.obj <- min(tmp.obj, tmp$tot.withinss)
    if(tmp.obj == tmp$tot.withinss){
      cl.kmns <- tmp
    }
  }
  
  ## Find representative gene clusters
  rep.gns <- NULL
  rep.gns <- sapply(as.list(as.data.frame(t(cl.kmns$centers))),
                    function(x){
                      rep.gns <- cbind(rep.gns, which.min(apply(abs(gData-x), 1,sum)))
                    } )
  rep.gns <- rownames(gData[rep.gns,]) #have gene names instead of index to avoid issues

  ## km.cntrd contains the cluster centers trajectory and rep. gene trajectory
  km.cntrd <- list(centers=cl.kmns$centers, cent.dat=gData[rep.gns,])

  if(!is.null(rSmpl.size)){
    for(iR in 1:length(rSmpl.size)){
      km.rnd <- c(km.rnd, list(replicate(n.randSmpl, sample(x=rownames(gData[!rownames(gData) %in% rep.gns,])
                                             ,size=rSmpl.size[iR]))))
                                      
    }
  }

  return(list(clst=cl.kmns, rep.gns=rep.gns, rnd.sample=km.rnd))
}

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
##' @param fix.w Logical variable, fit with fixed "W" matrix only the beta parameter if true (default FALSE)
##' @param w pass the W matrix if 
##' @return 
##' @author anas ahmad rana
rust.fit.nStt <- function(gData, tData, lambda = 0.01, n.states = 3, fit.as='lin', fix.w=FALSE, w=NULL){
  p <- nrow(gData)
  if(fix.w){
    x0 <- runif( p* n.states)
  } else if(fix.w==FALSE) {
    x0 <- runif( (n.states - 1) + nrow(gData) * n.states)
  }

  ## Function for the fitting procedure
  fun <- function(x){
    tmp <- rust.par(x=x, n.states=n.states, p=p, fix.w=fix.w)
    if(fix.w){
      wFit <- w
    } else {
      wFit <- tmp$w
    }
    betaFit <- tmp$beta

    fit <- rust.kStt(wFit, betaFit, tData)


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

  res <- nlminb(x0, fun, lower = 0, control=list(iter.max = 3000,
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
  
  return(list(fit=res, w=par$w, beta=par$beta,rss=rss, bic=bicSc, aic=aicSc))

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
rust.par <- function(gData=NULL, x, n.states, p=nrow(gData), fix.w=FALSE){
  
  if(fix.w){

    betaFit <- matrix( x, p, n.states)
    return(list(beta=betaFit))
  } else{
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
      p <- {lnx - n.states + 1}/ n.states
      betaFit <- matrix( x[-c(1:(n.states - 1))], p, n.states)
    }
    if(!is.null(gData)){
      rownames(betaFit) <- rownames(gData)
    }
    return(list(w=wFit, beta=betaFit))
  }
     

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

  

