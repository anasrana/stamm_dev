## ----------[ lodading necessary libraries ]--------------------
library(stats)          #fitting procedure
library(expm)           #Matrix exponential also loads library(Matrix)
library(multicore)      #Parallelisation of code

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
##' @param fit.pll 
##' @param n.smpl 
##' @param reps 
##' @return 
##' @author anas ahmad rana
rust.clst.fit <- function(gData, tData, lambda, n.states, fit.as='log2Dat', rSmpls,
                          fit.pll=FALSE, n.smpl=10,reps=2){
  
  ## Cluster gene trajectories
  cl <- rust.clustering(gData, km.k=6)
  rnd.sample <- rust.sampl(gData, rSmpl.size=rSmpls, n.randSmpl=10) 

  if(length(n.states)==1){

    ## Fit the cluster centroids to fix w
    ct.fit <- rust.centroid.fit(dat.cntrd=cl$clst$centers, t.dat=tData, lambda,
                                n.stt=n.states, fit.as, rSmpls, rand.sample=rnd.sample,
                                reps=reps, n.smpl=n.smpl, fit.pll=fit.pll, gData=gData)
    
    ms.stts <- list(rss=ct.fit$rss, bic=ct.fit$bic, aic=ct.fit$aic)
  } else if(length(n.states)>1) {
    ## Fit the cluster centroids to fix w
    ms.stts <- vector('list', length(n.states))
    names(ms.stts) <- paste('stt',n.states, sep='')
    for(i.s in 1:length(n.states)){
      n.stt <- n.states[i.s]
      ms.rss <-  10^10

      ct.fit <- rust.centroid.fit(dat.cntrd=cl$clst$centers, t.dat=tData, lambda, n.stt, fit.as, reps=reps)
      ms.stts[[i.s]] <- list(rss=ct.fit$rss, bic=ct.fit$bic, aic=ct.fit$aic)
    }
  }
  
  return(list(ms.cl=ms.stts, states=n.states, sample.N=rSmpls))
}

 
##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param dat.cntrd 
##' @param t.dat 
##' @param lambda 
##' @param n.stt 
##' @param fit.as 
##' @param rand.sample 
##' @param gData 
##' @param reps 
##' @param n.smpl 
##' @param fit.pll 
##' @return 
##' @author anas ahmad rana
rust.centroid.fit <- function(dat.cntrd, t.dat, lambda, n.stt, fit.as, rSmpls, rand.sample,
                              gData, reps, n.smpl, fit.pll){
  ms.rss <-  10^10
  for(is in 1:reps){
    tmp.fit <- rust.fit.kStt(gData=dat.cntrd, tData=t.dat, lambda=lambda, n.states=n.stt, fit.as=fit.as)
    if(tmp.fit$ms[1]<ms.rss){
      cl.fit <- tmp.fit
      ms.rss <- tmp.fit$ms[1]
    }
  }
  ms <- cl.fit$ms
  beta.centroid <- cl.fit$beta
  w <- cl.fit$w

  cat('\ntransition_rates_fit= \n')
  if(n.stt>2)
    print(diag(w[-1,]))
  else
    print(w)
  cat('\nrss = \n')
  print(ms[1])
  cat('\n now fitting genes \n')

  aic <- matrix(0,length( rand.sample), n.smpl)
  bic <- matrix(0,length( rand.sample), n.smpl)
  rss <- matrix(0,length( rand.sample), n.smpl)

  for(i in 1:length(rand.sample)){
    smpl <- rand.sample[[i]]
    for(j in 1:ncol(smpl)){
      gVec <- smpl[,j]
      fit.v <- rust.fit.gnlst(gVec, gData=gData, tData=t.dat, lambda, n.states=n.stt, fit.as, w, fit.pll=fit.pll)
      rss[i,j] <- fit.v$ms[1]
      bic[i,j] <- fit.v$ms[2]
      aic[i,j] <- fit.v$ms[3]
    }
  }
  rownames(rss) <- paste('samples',rSmpls,sep='')
  rownames(bic) <- paste('samples',rSmpls,sep='')
  rownames(aic) <- paste('samples',rSmpls,sep='')
  return(list(rss=rss, bic=bic, aic=aic, w=w))
}

                              

##' Fits betas for a list of genes given w
##'
##' 
##' @title rust.fit.gnlst
##' @param gName vector of gene names to be fitted
##' @param gData data matrix for fitting
##' @param tData time vector of data points
##' @param lambda L1 penalty parameter
##' @param n.states number of states to be fitted
##' @param fit.as string determining log manipulation of fit
##' @param w transition matrix
##' @param fit.pll logical for parallel fitting
##' @return 
##' @author anas ahmad rana
rust.fit.gnlst <- function(gName, gData, tData, lambda, n.states, fit.as, w, fit.pll=FALSE){

  gName <- as.list(gName)

  if(fit.pll){
    fit.gnes <- mclapply(gName, function(x)
                     cl.fit = rust.fit.kStt(gData=gData[x,], tData=tData, lambda=lambda,
                       n.states=n.states, fit.as=fit.as, w=w, fix.w=TRUE))
  } else {
    fit.gnes <- lapply(gName, function(x)
                         cl.fit = rust.fit.kStt(gData=gData[x,], tData=tData, lambda=lambda,
                           n.states=n.states, fit.as=fit.as, w=w, fix.w=TRUE))
  }

  ## take out important information from the fitting data
  betas <- NULL
  rss.v <- NULL
  for(il in 1:length(fit.gnes)){
    betas <- rbind(betas, fit.gnes[[il]]$beta)
    rss.v <- c(rss.v, fit.gnes[[il]]$ms)
  }

  rownames(betas) <- unlist(gName)
  names(rss.v) <- unlist(gName)
  rss <- sum(rss.v)
  n <- ncol(gData)*nrow(betas)
  ms <- rust.bic(rss, n, betas)
  
  return(list(beta=betas, rss.v=rss.v, w=w, ms=ms))
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
rust.clustering <- function(gData, km.init=100, km.k=NULL, rSmpl.size=NULL){
  if(is.null(km.k)){
    stop('choose number of states km.k') #check in number of clusters is given STOP if not
  }

  ## Initialise random sample list


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

  return(list(clst=cl.kmns, rep.gns=rep.gns))
 }

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param gData 
##' @param rSmpl.size 
##' @return 
##' @author anas ahmad rana
rust.sampl <- function(gData, rSmpl.size, n.randSmpl, rep.gns=NULL){
  km.rnd <- NULL
  if(!is.null(rep.gns)){
    for(iR in 1:length(rSmpl.size)){
      km.rnd <- c(km.rnd, list(replicate(n.randSmpl, sample(x=rownames(gData[!rownames(gData) %in% rep.gns,])
                                                            ,size=rSmpl.size[iR]))))
    }
  } else {
   for(iR in 1:length(rSmpl.size)){
     km.rnd <- c(km.rnd, list(replicate(n.randSmpl, sample(x=rownames(gData),size=rSmpl.size[iR]))))
   }
  }
  
  return(km.rnd)
}

                       
## ****************************************************************************************
## ----------[ Functions for least squares fitting ]---------------------------------------
## ****************************************************************************************
##' Fits data to aggregate Markov Chain model
##'
##' .. content for \details{} ..
##' @title rust.fit.kStt
##' @param gData gene expression data
##' @param tData time points of data
##' @param lambda L1 penalty parameter (default = 0.01)
##' @param n.states number of states in the fitting
##' @param fit.as Fitting lin, log2Dat, logDat, log2Al (default = 'lin')
##' @param fix.w Logical variable, fit with fixed "W" matrix only the beta parameter
##' if true (default FALSE)
##' @param w pass the W matrix if 
##' @return The function returns fit which is returned from the fitting function.
##' It returns the fitted $w$ matrix and the $/beta$ matrix. It also returns a  obj vector
##' that contains the rss, bic and aic scores for the fit.
##' @author anas ahmad rana
rust.fit.kStt <- function(gData, tData, lambda = 0.01, n.states = 3, fit.as='lin',
                          fix.w=FALSE, w=NULL){
  p <- nrow(gData)
  if(fix.w){
    p <- 1
    x0 <- runif( p* n.states)
    wFit <- w
  } else if(fix.w==FALSE) {
    x0 <- runif( (n.states - 1) + nrow(gData) * n.states)
    wFit <- NULL
  }

  if(0 %in% gData & fit.as %in% c('log2Dat', 'logDat', 'log2Al'))
    gData = gData+1

  ## Functions for the fitting procedure depending on how to fit
  if(fit.as=='lin'){
    fun <- function(x){
      tmp <- rust.par(x=x, n.states=n.states, p=p, fix.w=fix.w, wFit=wFit)
      wFit <- tmp$w
      betaFit <- tmp$beta
      fit <- rust.kStt(wFit, betaFit, tData)
      rss <- (fit$y - gData)^2
      penalty <- lambda*sum(abs(betaFit)) 
      ss <- c(rss,penalty)
      sum(ss)
    }
  } else if(fit.as=='log2Dat'){
    fun <- function(x){
      tmp <- rust.par(x=x, n.states=n.states, p=p, fix.w=fix.w, wFit=wFit)
      wFit <- tmp$w
      betaFit <- tmp$beta
      fit <- rust.kStt(wFit, betaFit, tData)
      rss <- (log2(fit$y) - gData)^2
      penalty <- lambda*sum(abs(betaFit)) 
      ss <- c(rss,penalty)
      sum(ss)
    }
  } else if(fit.as=='logDat'){
    fun <- function(x){
      tmp <- rust.par(x=x, n.states=n.states, p=p, fix.w=fix.w, wFit=wFit)
      wFit <- tmp$w
      betaFit <- tmp$beta
      fit <- rust.kStt(wFit, betaFit, tData)
      rss <- (log(fit$y) - gData)^2
      penalty <- lambda*sum(abs(betaFit)) 
      ss <- c(rss,penalty)
      sum(ss)
    }
  } else if(fit.as=='log2Al'){
    fun <- function(x){
      tmp <- rust.par(x=x, n.states=n.states, p=p, fix.w=fix.w, wFit=wFit)
      wFit <- tmp$w
      betaFit <- tmp$beta
      fit <- rust.kStt(wFit, betaFit, tData)
      rss <- (log2(fit$y) - log2(gData))^2
      penalty <- lambda*sum(abs(betaFit)) 
      ss <- c(rss,penalty)
      sum(ss)
    }

  } else {
    stop('String fit.as not recognised')
  }
  ##

  res <- nlminb(x0, fun, lower=0, control=list(iter.max = 10000, eval.max=7000, rel.tol=10^-14, sing.tol=10.^-7))
  par <- rust.par(gData, res$par, n.states, fix.w=fix.w, wFit=w)
  
  

  n <- ncol(gData)
  
  if(fix.w){
    rss <- res$objective - lambda*sum(res$par)
    obj <- rss
    names(obj) <- 'rss'
  } else {
    rss <- res$objective - lambda*sum(res$par[-c(1:(n.states -1))])
    obj <- rust.bic(rss, n, par$beta)
    
  }

  return(list(fit=res, w=par$w, beta=par$beta, ms=obj))
}

##' Calculates BIC and AIC 
##'
##' BIC and AIC calculated for least squared fit
##' @title rust.bic
##' @param rss The residual sum of squares of the fit
##' @param n numer of data points used when fitting
##' @param beta the number of parameters used to calculate Df
##' @return vector of rss bic and aic
##' @author anas ahmad rana
rust.bic <- function(rss, n, beta){
  Df <- sum(beta>10^-3)
  bicSc <- n*log(rss/(n-1)) + log(n) * Df
  aicSc <- n*log(rss/(n-1)) + 2 * Df
  ms <- c(rss, bicSc, aicSc)
  names(ms) <- c('rss', 'bic', 'aic')
  return(ms)
}

##' Reshapes parameter vecot, x, into beta matrix and w matrix (or only beta matrix)
##'
##' .. content for \details{} ..
##' @title rust.par
##' @param gData data matrix used for naming beta
##' @param x parameter vector to reshape into beta and w
##' @param n.states number of states in model
##' @param p no of genes to be fitted (default # rows in data) 
##' @param fix.w logical if w is kept fixed or not
##' @return The function returns w only if fix.w=F and it returns a beta matrix rearranged from the x vector.
##' @author anas ahmad rana
rust.par <- function(gData=NULL, x, n.states, p=nrow(gData), fix.w=FALSE, wFit=NULL){
  ## only rearrange 
  if(fix.w){
    if(length(x)!=n.states)
      stop('No of parameters has to equal no of states when fitting per gene')
    betaFit <- matrix( x, 1, n.states)
    return(list(w=wFit, beta=betaFit))
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
##' Function that takes as arguments the w_fit matrix, the beta_fit matrix and time points, it then calculates a
##' trajectory
##'
##' 
##' @title rust.kStt
##' @param wFit W transition matrix default=NULL
##' @param betaFit beta matrix p x T(n)
##' @param t time points
##' @return S unlogged trajectory for the arguments. P state occupation probabilities.
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
  return(list(y=S, prob=P))
}

