StammKmeans <- function(g.dat, m.v=2:30, iter.max=100, nstart=100) {
    ## Transform data and standardize for kmeans
    g.sdA <- apply(asinh(g.dat), 1, sd)
    g.normA <- asinh(g.dat) / g.sdA

    
    within.ss <- rep(NA, length(m.v))
    for (i.m in 1:length(m.v)) {
        g.km <- kmeans(g.normA, m.v[i.m], iter.max=iter.max, nstart=nstart)
        within.ss[i.m] <- g.km$tot.withinss
    }

    J <- within.ss[-1] / within.ss[1:28]
    return(list(DeltaJ=(1-J), ))
}
