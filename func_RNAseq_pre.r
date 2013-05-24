library('ggplot2')
library('edgeR')


LoadMcf10aDat <- function(f.name = NULL, dat.col, col.order) {
  if (is.null(f.name)) {
    f.name  <- "/Users/anas/Dropbox/complexity/phd/data/Sandra2012/gene_exp.txt"
    dat.col  <- 2:21
    dat.order <- c(2, 12, 13:20, 1, 3:7)
    dat.null  <-  8:11
  }

  raw.dat <- read.table(f.name, header = TRUE)
  g.raw <- as.matrix(raw.dat[dat.col])
  g.raw <- g.raw[, dat.order]
  rownames(g.raw)  <- raw.dat$external_gene_id
  group.name  <- unlist(lapply(data.frame(strsplit(colnames(g.raw), '_'))[2, ], as.character))

  y  <- DGEList(counts = g.raw, group = group.name )

  keep  <- rowSums(cpm(y) > 1) >= 3
  y <- y[keep, ]
  y$sample$lib.size  <- colSums(y$counts)

  norm.factor <- NULL
  lib.size <- NULL
  cut.off  <- seq(0.1, 0.4, 0.02)
  y  <- calcNormFactors(y, method = 'TMM')

  g.norm <- cpm(y, normalized.lib.sizes=TRUE)
  colnames(g.norm) <- group.name
  return(g.norm)
}
