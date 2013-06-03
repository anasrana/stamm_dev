## ----------[ lodading necessary libraries ]--------------------
library(stats)          #fitting procedure
library(expm)           #Matrix exponential also loads library(Matrix)
library(multicore)      #Parallelisation of code


## ****************************************************************************************
## ----------[ Functions for least squares fitting ]---------------------------------------
## ****************************************************************************************
##' Fits data to aggregate Markov Chain model
##'
##' .. content for \details{} ..
##' @title RustFitKstt
##' @param g.dat gene expression data
##' @param t.dat time points of data
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
RustFitKstt <- function(g.dat, t.dat, lambda = 0.01, n.states = 3, fix.w=FALSE, w=NULL) {
    p <- nrow(g.dat)
    if (fix.w) {
        p <- 1
        x0 <- runif (p * n.states)
        wFit <- w
    } else if (fix.w == FALSE)  {
        x0 <- runif ((n.states - 1) + nrow(g.dat) * n.states )
        wFit <- NULL
    }

    g.dat.l <- asinh(g.dat)
    if ( !is.vector(g.dat) )
        g.nl <- apply(g.dat, 1, sd)
    else
        g.nl <- sd(g.dat)

    fun <- function(x){
        tmp <- RustPar(x = x, n.states = n.states, p = p, fix.w = fix.w, wFit = wFit)
        wFit <- tmp$w
        betaFit <- tmp$beta
        fit <- RustKstt(wFit, betaFit, t.dat)
        rss <- {asinh(fit$y) - g.dat.l}^2
        penalty <- lambda * sum(betaFit / g.nl)
        ss <- c(rss, penalty)
        obj <- sum(ss)
        obj
    }
    if (is.vector(g.dat)) {
        par.scale <- c(rep(1, n.states - 1), rep(max(g.dat), n.states))
    } else {
        par.scale <- c(rep(1, n.states - 1), rep(apply(g.dat, 1, max), n.states))
    }
    res <- nlminb(x0, fun, lower = 0, upper = max(g.dat), scale = 1 / par.scale,
                  control=list(iter.max=10000, eval.max=7000, rel.tol=10^-14, sing.tol=10^-14))

    par <- RustPar(g.dat, res$par, n.states, fix.w=fix.w, wFit=w)

    n <- ncol(g.dat)
    if (fix.w) {
        rss <- res$objective - lambda * sum(par$beta / g.nl)
        obj <- rss
        names(obj) <- 'rss'
    } else {
        rss <- res$objective - lambda * sum(par$beta / g.nl)
        obj <- RustBic(rss, n, par$beta)
    }

    return(list(fit = res, w = par$w, beta = par$beta, ms = obj))
}

##' Calculates BIC and AIC
##'
##' BIC and AIC calculated for least squared fit
##' @title RustBic
##' @param rss The residual sum of squares of the fit
##' @param n numer of time points used when fitting
##' @param beta the number of parameters used to calculate Df
##' @return vector of rss bic and aic
##' @author anas ahmad rana
RustBic <- function(rss, n, beta, b.thresh = 10^-4) {
    Df <- sum(beta > b.thresh)
    bicSc <- n * log(rss / (n - 1)) + log(n) * Df
    aicSc <- n * log(rss / (n - 1)) + 2 * Df
    ms <- c(rss, bicSc, aicSc)
    names(ms) <- c('rss', 'bic', 'aic')
    return(ms)
}

##' Reshapes parameter vecot, x, into beta matrix and w matrix (or only beta matrix)
##'
##' .. content for \details{} ..
##' @title RustPar
##' @param g.dat data matrix used for naming beta
##' @param x parameter vector to reshape into beta and w
##' @param n.states number of states in model
##' @param p no of genes to be fitted (default # rows in data)
##' @param fix.w logical if w is kept fixed or not
##' @return The function returns w only if fix.w=F and it returns beta matrix rearranged from the x vector.
##' @author anas ahmad rana
RustPar <- function(g.dat=NULL, x, n.states, p=nrow(g.dat), fix.w=FALSE, wFit=NULL, fit.cent=FALSE) {
    ## only rearrange
    if (fit.cent) {
        p <- nrow(g.dat)
        betaFit <- matrix(x, p, n.states)
        return(list(w = wFit, beta = betaFit))
    } else if (fix.w && !fit.cent) {
        if (length(x) != n.states )
            stop('No of parameters has to equal no of states when fitting per gene')
        betaFit <- matrix(x, 1, n.states)
        return(list(w = wFit, beta = betaFit))
    } else {
        ## W matrix from x[1:n-1]
        if (n.states == 2) {
            wFit <- matrix(0, n.states, n.states)
            wFit[2, 1] <- x[1]
        } else if (n.states >= 3)  {
            wFit <- matrix(0, n.states, n.states)
            diag(wFit[-1, ]) <- x[1:(n.states - 1)]
        } else {
            wFit <- NULL
        }
        ## Assign other x values to beta
        if (n.states == 1) {
            betaFit <- matrix(x, p, n.states)
        } else {
            lnx <- length(x)
            p <- {lnx - n.states + 1} / n.states
            betaFit <- matrix( x[-c(1:(n.states - 1))], p, n.states)
        }
        if (!is.null(g.dat)) {
            rownames(betaFit) <- rownames(g.dat)
        }
        return(list(w = wFit, beta = betaFit))
    }
}

##' Function that takes as arguments the w_fit matrix, the beta_fit matrix and time
##' points, it then calculates a trajectory
##'
##'
##' @title RustKstt
##' @param wFit W transition matrix default=NULL
##' @param betaFit beta matrix p x T(n)
##' @param t time points
##' @return S unlogged trajectory for the arguments. P state occupation probabilities.
##' @author anas ahmad rana
RustKstt <- function(wFit=NULL, betaFit, t){
    n <- length(t)
    p <- nrow(betaFit)
    ## create new matrix containing -offdiag on diag
    if (!is.null(wFit)) {
        odiagElements <- diag(wFit[-1, , drop=FALSE ])
        diag(wFit) <- c(-odiagElements, 0)
        k <- nrow(wFit)
        ## initial condition
        p0 <- rep(0, k)
        p0[1] <- 1
        P <- NULL
        for (i in 1:n) {
            pt <- expm::expm(wFit * t[i], method='Ward77') %*% p0
            P <- cbind(P, pt)
            ## mean gene expression
            S <- betaFit %*% P
        }
    } else {
        S <- rep(betaFit, n)
        P <- 1
    }
    return(list(y = S, prob = P))
}

## ******************************************************************************************
## Clustering and fitting centroids
## ******************************************************************************************

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title
##' @param g.dat
##' @param t.dat
##' @param stt.v
##' @param n.cl
##' @param n.core
##' @return
##' @author anas ahmad rana
RustCrossValClust.p <- function(g.dat, t.dat, stt.v, n.cl = c(2, seq(5, 30, 5)), n.core = 20) {
    ## cluster genes
    cl.dat <- RustCluster(g.dat, n.cl = n.cl)
    g.cen <- cl.dat[['cent.dat']]
    ## define matrix for parameters in mclapply
    cs.mat <- cbind(rep(stt.v, length(n.cl)), rep(n.cl, length(stt.v)))
    ## start mcl apply for all time-p knockout
    n.ml <- length(stt.v) * length(n.cl)
    fit.k.m <- mclapply(1:n.ml, function(x) {
        x.name <- paste('m.', cs.mat[x, 2], sep='')
        i.v <- 1
        ## First fit all the data for different m and k
        fit.cl <- vector('list', length(t.dat))
        fit.conv <- 10^10
        fit.iter <- 1
        while (fit.conv != 0 | fit.iter > 5) {
            fit <- RustFitKstt(g.cen[[x.name]], t.dat, lambda=0, n.states = cs.mat[x, 1])
            fit.conv <- fit$fit$convergence
            fit.iter <- fit.iter + 1
        }
        print(paste('done all:', cs.mat[x, 1]))
        fit.cl[[i.v]]  <- fit
        ## Fit 2:end t-ko for different m and k
        for (i.v in 2:length(t.dat)) {
            fit.conv <- 10^10
            fit.iter <- 1
            while (fit.conv != 0 | fit.iter > 5) {
                fit <- RustFitKstt(g.cen[[x.name]][, -i.v], t.dat[-i.v], lambda=0,
                                   n.states = cs.mat[x, 1])
                fit.conv <- fit$fit$convergence
                fit.iter <- fit.iter + 1
            }
            fit.cl[[i.v]]  <- fit
        }
        names(fit.cl) <- paste('tpt.', 1:length(t.dat), sep='')
        print(paste('done t.ko:', cs.mat[x, 1]))
    }, #end of function in mclapply
                        mc.preschedule = TRUE, mc.cores = n.core)
    names(fit.k.m) <- paste('k.', cs.mat[, 1],'_m.', cs.mat[, 2], sep='')
    return(list(cl.dat = cl.dat, fit = fit.k.m))
}

##' Fits centroid for given number of cluster sizes and fixed number of states fits
##'
##' .. content for \details{} ..
##' @title
##' @param g.dat
##' @param t.dat
##' @param g.cent
##' @param n.cl
##' @param n.core
##' @param pen
##' @param n.stt
##' @return
##' @author anas ahmad rana
RustFitCentroid.p <- function(g.dat, t.dat, g.cent = NULL, n.cl = c(2, seq(5, 30, 5)),
                              n.core = 20, pen = 0, n.stt = 4) {
    ## First cluster all genes
    if (is.null(g.cent)) {
        cl.dat <- RustCluster(g.dat, n.cl)
        g.cent <- cl.dat[["cent.dat"]]
        rep.gns <- cl.dat[["rep.gns"]]
    }
    ## Fit all cluster sizes in paralell
    options(cores=n.core)
    fit <- mclapply(1:length(g.cent), function(x){
        fit.conv <- 10^10
        fit.iter  <- 1
        max.fit.iter <- 5
        while (fit.conv != 0 | fit.iter > max.fit.iter) {
            fit <- RustFitKstt(g.cent[[x]], t.dat , lambda = pen, n.states = n.stt)
            fit.conv <- fit$fit$convergence
            fit.iter  <- fit.iter + 1
        }
        if (fit.conv != 0){
            print('ERROR: did not converge')
        }
        print(paste('fit: cl = ', n.cl[x], 'with RSS =',fit$ms[1], '... done!'))
        return(fit)
    }, #end of function in mclapply
                    mc.preschedule = TRUE)

    return(list(fit = fit, cl.dat = cl.dat))
}

##' Scale data and calculate k-means cluster for all vector elements
##' of n.cl
##'
##' .. content for \details{} ..
##' @title RustCluster
##' @param g.dat data to be clustered
##' @param n.cl vector of arbitrary length
##' @return
##' @author anas ahmad rana
RustCluster <- function(g.dat, n.cl) {
    ## Scale input data to unit variance
    g.sd <- apply(g.dat, 1, sd)
    g.norm <- (g.dat)/g.sd
    ## initialise empty lists
    g.cl <- vector('list', length(n.cl))
    g.cent <- vector('list', length(n.cl))
    g.cl.names <- vector('list', length(n.cl))
    ## perform a k-means clustering
    for (i.n in 1:length(n.cl)) {
        i.m  <- n.cl[i.n]
        g.kn.cl <- kmeans(g.norm,i.m, iter.max=100, nstart=100)
        g.cl[[i.n]] <- g.kn.cl
        g.cent[[i.n]] <- g.kn.cl$centers
        g.rep.names <- rep(NA, i.m)
        for (i.clg in 1:nrow(g.kn.cl$centers)) {
            d.g <- apply((abs(g.norm - g.kn.cl$centers[i.clg, ])), 1, sum)
            tmp <- which(d.g == min(d.g))
            g.rep.names[i.clg] <- rownames(g.norm)[tmp]
        }
        g.cl.names[[i.n]] <- g.rep.names
    }
    names(g.cl) <- paste('m.', n.cl, sep='')
    names(g.cent) <- paste('m.', n.cl, sep='')
    names(g.cl.names)  <- paste('m.', n.cl, sep='')
    ## return variables k-means centroid, fit and names of
    ## representative genes
    return(list(cent.dat=g.cent, cl.kmean=g.cl, rep.gns=g.cl.names))
}


##' Cluster (kmeans) and fit under timepoint deletions
##'
##' .. content for \details{} ..
##' @title
##' @param g.dat
##' @param t.dat
##' @param m.cl
##' @param t.ko
##' @param k.stt
##' @return
##' @author anas ahmad rana
ParClusterCV <- function(g.dat, t.dat=NULL, m.cl=seq(5, 20, 2), t.ko=2:ncol(g.dat), k.stt=4) {
    ## some parameters determined from the input
    n.genes <- nrow(g.dat)
    g.norm <- g.dat / apply(g.dat, 1, sd)
    vec.mt <- cbind(rep(1:length(m.cl), length(t.ko)), rep(t.ko, each=length(m.cl)))

    fit <- mclapply(1:nrow(vec.mt), function(x){
        t.dl <- vec.mt[x, 2]
        ## Clustered after time point deletion using k-means
        g.km <- kmeans(g.norm[, -t.dl], vec.mt[x, 1], iter.max=100, nstart=20)
        g.ct <- g.km$centers
        ## fit cluster centroids with fixed states k
        fit.m <- RustFitKstt(g.ct, t.dat[-t.dl], lambda=0, n.states=k.stt)
        ## Fix w and fit all genes used in clustering
        w <- fit.m$w
        beta <- matrix(0, n.genes, k.stt)
        rownames(beta) <- rownames(g.dat[1:n.genes, ])
        for (i.j in 1:n.genes) {
            fit.tmp <- RustFitKstt(g.dat[i.j, -t.dl], t.dat[-t.dl], lambda=0,
                                   n.states=k.stt, fix.w=TRUE, w=w)
            beta[i.j, ] <- fit.tmp$beta
        }
        return(list(cl=g.km, cl.fit=fit.m, w=w, beta=beta))
    }, mc.preschedule = TRUE)

    names(fit) <- paste('m', vec.mt[, 1], 'td', vec.mt[, 2], sep='.')
    ## output from the function
    return(list(g.norm=g.norm, fit=fit))
}

FitClGns <- function(g.dat, t.dat, l.pen=0, k.stt, m, pll=FALSE, w=NULL) {
    if (is.null(w)) {
        ## Normalise data to be univariate and fit clusters
        g.norm <- g.dat / apply(g.dat, 1, sd)
        g.km <- kmeans(g.norm, m, iter.max=100, nstart=100)
        g.ct <- g.km$centers
        ## fit cluster centroids with fixed states k
        fit.m <- RustFitKstt(g.ct, t.dat, lambda=0, n.states=k.stt)
        w <- fit.m$w
    }
    ## Fit all genes in g.dat with w from clustering and single penalty
    if (is.numeric(l.pen)) {
        fit.g <- RustFitGns(g.dat=g.dat, t.dat=t.dat, lambda=l.pen, n.states=k.stt, w=w, pll=pll)
    } else if (is.vector(l.pen)) {
        fit.g <- vector('list', length(l.pen))
        for (i.l in 1:length(l.pen)) {
            fit.g[[i.l]] <- RustFitGns(g.dat=g.dat, t.dat=t.dat, lambda=l.pen[i.l],
                                       n.states=k.stt, w=w, pll=pll)
        }
    }
    return(list(fit.g=fit.g, fit.m=fit.m, cl.km=g.km))
}

##' Fits betas for a list of genes
##'
##'
##' @title rust.fit.gnlst
##' @param g.names
##' @param g.dat data matrix for fitting
##' @param t.dat time vector of data points
##' @param lambda L1 penalty parameter
##' @param n.states number of states to be fitted
##' @param w transition matrix
##' @param pll logical run in parallel if set to true
##' @return
##' @author anas ahmad rana
RustFitGns <- function(g.names=NULL, g.dat, t.dat, lambda=0, n.states, w, pll=FALSE){
    if (!is.null(g.names)) {
        g.name.idx <- which(rownames(g.dat)==g.names)
    } else {
        g.name.idx <- 1:nrow(g.dat)
    }

    if (pll) {
        fit.g <- mclapply(g.name.idx, function(x)
                          cl.fit = RustFitKstt(g.dat=g.dat[x, ], t.dat=t.dat, lambda=lambda,
                              n.states=n.states, w=w, fix.w=TRUE)
                          )
    } else {
        fit.g <- lapply(g.name.idx, function(x)
                        cl.fit = RustFitKstt(g.dat=g.dat[x, ], t.dat=t.dat, lambda=lambda,
                            n.states=n.states, w=w, fix.w=TRUE))
    }

    ## take out important information from the fitting data
    betas <- NULL
    rss.v <- NULL
    for (il in 1:length(fit.g)) {
        betas <- rbind(betas, fit.g[[il]]$beta)
        rss.v <- c(rss.v, fit.g[[il]]$ms)
    }

    rownames(betas) <- rownames(g.dat[g.name.idx, ])
    names(rss.v) <- rownames(g.dat[g.name.idx, ])
    rss <- sum(rss.v)
    n <- length(t.dat)
    ms <- RustBic(rss, n, betas)

    return(list(beta=betas, rss.v=rss.v, w=w, ms=ms))
}
