## ----------[ lodading necessary libraries ]--------------------
## source("func_plotting.r")               
library(stats)
## library(pracma)                       #library for lsqnonlin
## library(Matrix)                       #matrix specific functions
library(matrixStats)                  #use rowSds from library
#library(pls)                            #use stdize
library(BB)                             #use for spg
library(expm)

## ----------[ Simulations on single cells]--------------------
##' Forward simulation of model using single cells
##'
##' @param nCells number of cells used in simulation (default 200)
##' @param nGenes number of genes simulated (default 9)
##' @param tau average jumpTime per state
##' @param nStates number of states simulated (default 2)
##' @param betaVec predefined beta values for simulation (default randomly generated)
##' @param dt timestep used throughout simulation
##' @param endTime final timepoint (in "real" time)
##' @param avNoise standard deviation of Gaussian noise added to the average trajectory (default 0.001)
##' @param sttNoise vector of standard deviations for Gaussian noise per state (default 0)
singelCellSML <- function(nCells = 200,nGenes = 9, tau = c(3.5,5,14.5), nStates = 2,
                          betaVec=NULL, dt = 0.01, endTime = 30, avNoise=0.001, sttNoise=rep(0,nStates),
                          nonzero=FALSE){
  ## Checking if some of the arguments passed are the right format and adjusting if possible
  if(!is.vector(betaVec) & !is.null(betaVec)){
    stop('beta must be passed as vector')
  }
  if(length(sttNoise)!=nStates & is.vector(sttNoise)){
    sttNoise <- rep(sttNoise, nStates)
  }
  jumpTime <- NULL
  for(i in 1:(nStates-1)){
    jumpTime <- rbind(jumpTime, rexp(nCells, 1/tau[i]))
  }
  ## Checks if beta values are passed in the argument (as vector)
  ## Assigns randomvalues to beta if not passed as argument
  if(is.null(betaVec)){
    betaVec <- rnorm(nGenes*nStates, sd=5)
    if(nonzero)
      betaVec[betaVec<0] <- abs(betaVec[betaVec<0])
    else
      betaVec[betaVec<0]= 10^(-40)
  }
  ## Reshape betavector as matrix/
  betaVals <- matrix(betaVec, nGenes, nStates)
  ## Initialise results matrix with zeros
  gSim <- matrix( rep(0, nGenes * endTime / dt), nGenes, endTime / dt)
  nPoints <- endTime / dt
  for(iCells in 1:nCells){
    gSimC <- NULL
    for(jStates in 1:(nStates-1)){
      nSentries <- ceiling(jumpTime[jStates, iCells]/dt)
      val_tmp <- betaVals[,rep.int(jStates,nSentries)]
      gSimC <- cbind(gSimC,  val_tmp + rnorm(length(val_tmp),sd=sttNoise[jStates]) )
      rm(val_tmp)
    }
    nSentries <- (nPoints - ncol(gSimC))
    if(nSentries>0){
      val_fin <- betaVals[,rep.int(nStates,nSentries)]
      gSimC <- cbind(gSimC, val_fin + rnorm(length(val_fin),sd=sttNoise[nStates]))
      rm(val_fin)
    }
    gSim <- gSim + gSimC[,1:nPoints]/nCells
  }
  ## Add gaussian noise to the data
  datasim <- log(gSim) + matrix(rnorm(length(gSim), sd=avNoise), dim(gSim))
  dataSim <- exp(datasim)
  ## Return values are full simulated data all time points, beta and if t given gSim with t-pts
  return(list(gsim=gSim, beta=betaVals, dataSim=dataSim)) 
}


## ----------[ Simulations on single cells]--------------------
##' Forward simulation of model using single cells
##'
##' @param nCells number of cells used in simulation (default 200)
##' @param nGenes number of genes simulated (default 9)
##' @param tau average jumpTime per state
##' @param nStates number of states simulated (default 2)
##' @param betaVec predefined beta values for simulation (default randomly generated)
##' @param dt timestep used throughout simulation
##' @param endTime final timepoint (in "real" time)
##' @param avNoise standard deviation of Gaussian noise added to the average trajectory (default 0.001)
##' @param sttNoise vector of standard deviations for Gaussian noise per state (default 0)
rust.simSnglCl <- function(nCells = 200,nGenes = 9, tau = c(3.5,5,14.5), nStates = 2,
                          betaVec=NULL, dt = 0.01, endTime = 30, avNoise=0.001, sttNoise=rep(0,nStates)){
  ## Checking if some of the arguments passed are the right format and adjusting if possible
  if(!is.vector(betaVec) & !is.null(betaVec)){
    stop('beta must be passed as vector')
  }
  if(length(sttNoise)!=nStates & is.vector(sttNoise)){
    sttNoise <- rep(sttNoise, nStates)
  }
  ## get jump time from the state occupation time $/tau$
  jumpTime <- NULL
  for(i in 1:(nStates-1)){
    jumpTime <- rbind(jumpTime, rexp(nCells, 1/tau[i]))
  }
  ## Checks if beta values are passed in the argument (as vector)
  ## Assigns randomvalues to beta if not passed as argument
  if(is.null(betaVec)){
    betaVec <- rnorm(nGenes*nStates, sd=5)
    betaVec[betaVec<0]= 10^(-40)
  }
  ## Reshape betavector as matrix/
  betaVals <- matrix(betaVec, nGenes, nStates)
  ## Initialise results matrix with zeros
  gSim <- matrix( rep(0, nGenes * endTime / dt), nGenes, endTime / dt)
  nPoints <- endTime / dt
  for(iCells in 1:nCells){
    gSimC <- NULL
    for(jStates in 1:(nStates-1)){
      nSentries <- ceiling(jumpTime[jStates, iCells]/dt)
      val_tmp <- betaVals[,rep.int(jStates,nSentries)]
      gSimC <- cbind(gSimC,  val_tmp + rnorm(length(val_tmp),sd=sttNoise[jStates]) )
      rm(val_tmp)
    }
    nSentries <- (nPoints - ncol(gSimC))
    if(nSentries>0){
      val_fin <- betaVals[,rep.int(nStates,nSentries)]
      gSimC <- cbind(gSimC, val_fin + rnorm(length(val_fin),sd=sttNoise[nStates]))
      rm(val_fin)
    }
    gSim <- gSim + gSimC[,1:nPoints]/nCells
  }
  ## Add gaussian noise to the data
  datasim <- log(gSim) + matrix(rnorm(length(gSim), sd=avNoise), dim(gSim))
  dataSim <- exp(datasim)
  ## Return values are full simulated data all time points, beta and if t given gSim with t-pts
  return(list(gsim=gSim, beta=betaVals, dataSim=dataSim)) 
}


##'  Function that add noise to normal data
##' 
##' @title addNoise
##' @param gData 
##' @param nsSd 
##' @param nsTyp 
##' @return datasim
##' @author anas ahmad rana
addNoise <- function(gData, nsSd, nsTyp){
  if(nsTyp==1){
    datsim <- log(gData) + abs(rowMeans(log(gData)))*matrix(rnorm(length(gData), sd=nsSd), dim(gData))
    datasim <- exp(datsim)
    return(datasim)
  } else if(nsTyp==2) {
    datsim <- log(gData) + matrix(rnorm(length(gData), sd=nsSd), dim(gData))
    datasim <- exp(datsim)
    return(datasim)
  } else if(nsTyp==3) {
    datasim <- gData + matrix(rnorm(length(gData), sd=nsSd), dim(gData))
    return(datasim)
  } 
}

## ----------[ Functions for least squares fitting ]----------------------------------------
##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param gData 
##' @param tData 
##' @param tuning 
##' @param n.states 
##' @return 
##' @author anas ahmad rana
fitnSttModel <- function(gData, tData, tuning = 0.01, n.states = 3){
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


    rss <- (log(fit$S) - log(gData))^2
    penalty <- tuning*sum(abs(betaFit)) ##why did I do the sqrt
    ## penalty <- tuning*sum(abs(x))
    ss <- c(rss,penalty)
    sum(ss)
  }


  res <- nlminb(x0, fun, lower = 0, upper=max(tData)*2, control=list(iter.max = 3000,
                                                          eval.max=4000, rel.tol=10^-14))
  ## n.states=n.states, tuning=tuning, gData=gData, tData=tData) ##add if using external function
  rss <- res$objective - tuning*sum(res$par[-c(1:(n.states -1))])
  tmpPar <- res$par
  tmpPar[abs(tmpPar)<10^-3] <- 0
  n.zeros <- sum(tmpPar==0)
  Df <- p*n.states + n.states -1 - n.zeros
  n <- ncol(gData) * nrow(gData)

  bicSc <- n*log(rss/(n-1)) + log(n) * Df
  aicSc <- n*log(rss/(n-1)) + 2 * Df

  return(list(fit=res, bic=bicSc, aic=aicSc))

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
