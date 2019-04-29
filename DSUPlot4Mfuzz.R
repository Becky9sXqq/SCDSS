#' SJ plot for DSU results
#'
#' This function are used to plot the results of \code{BetaBinomialDSU} functions.
#' \code{\link[DSUPlot.heatmap]{DSUPlot.heatmap}},
#' \code{\link[DSUPlot.Mfuzz]{DSUPlot.Mfuzz}},
#' \code{\link[DSUPlot.GO]{DSUPlot.GO}},
#' \code{\link[DSUPlot.PCA]{DSUPlot.PCA}},
#' \code{\link[DSUPlot.distance]{DSUPlot.distance}},
#' \code{\link[DSUPlot.River]{DSUPlot.River}}
#' are useful for visulizating Betabinomial results for  further detection for cell specific SJs, BetaBinomialDSU only gives the important SJs
#' which can used to classify different groups. In addition, the function also can output a data.table of all or given group-specific (or group-classify-improtant) 
#' SJs that have been calculated.
#' @rdname DSUPlot.Mfuzz
#' @name DSUPlot.Mfuzz
#' @importFrom data.table setDT
#' @importFrom data.table setkey
#' @importFrom data.table as.data.table
#' @import Mfuzz
#' @param cell.psi an \code{SCASDataSet-class} object obtained from BetaBinomialDSU.
#' @param sampleInfo a table for sample information
#' @param used_celltype for \code{DSU.plot}: the cells that you want to plot or output, only used for Betabinomial.
#' @param del_psi the cutoff of deltaPSI for delta_PSI, only consider the SJ which deltaPSI
#' @param MEAN_psi the mean value of psi per sample, used for cutoff
#' @param NA_frac the NA fraction cutoff for the Mfuzz clustering object
#' @param NA_fill_meth the method for the NA vlaue filling in Mfuzz modified PSI table
#' @param min_std the standard deviance for variance filter in Mfuzz clustering object
#' @param clu_num the number for clustering groups, generally same as the celltype numbers
#' @param pdf_path the pathway for pdf output, parameter show defines the plot show out or save in pdf file
#' @param adjp the cutoff of adj_pval in heatmapPlot, which means the adjusted p-values after betabinomial test.
#' @param ... futher parameters of pheatmap.
#' @return a heatmap or group-specific SJs' data.table
#' @seealso \code{\link[Mfuzz]{Mfuzz}}
#' @example 
#' DSUplot.Mfuzz <- function(PSI_table,
#' used_celltype = unique(SampleInfo$cell.type),
#' NA_frac = 0.25,NA_fill_meth = 'knnw',min_std = 0.3,
#' clu_num = number_of_celltypes,show = F,pdf_path)
#' @export
DSUPlot.Mfuzz <- function(PSI_table,
                          used_celltype = cell_order,
                          NA_frac = 0.05,
                          NA_fill_meth = 'knnw',
                          min_std = 0.2,
                          clu_num = number_of_celltypes,
                          show = F,
                          pdf_path,
                          ...
){
  for(pack in c('dplyr','reshape2','ggplot2',
                'splines','data.table','pheatmap','RColorBrewer')){
    usePackage(pack)}
  for(pack in c('Mfuzz','Biobase','clusterProfiler','DOSE','org.Mm.eg.db','org.Hs.eg.db','AnnotationDbi')){
    useBiocPackage(pack)}

  options(stringsAsFactors = F)
  set.seed(201904)
  branch = PSI_table[,as.vector(used_celltype)]

  print( '----------- data processing  ----------')
  eset <- new('ExpressionSet', exprs=as.matrix(branch))
  ### cutoff based on NA fraction
  eset <- filter.NA(eset, thres = NA_frac)
  ### impute NA based on knn (default)
  eset <- fill.NA(eset, mode = NA_fill_meth)
  ### filtering and standadization
  eset <- filter.std(eset,min.std = min_std,visu=F)
  eset <- standardise(eset)
  ### running cluster --------
  cl <- mfuzz(eset,c = clu_num, m = mestimate(eset)[1])
  clu <- as.data.frame(table(as.data.frame(cl$cluster)[,1]))
  clu$Var1 <- NULL
  cluster=cl$cluster
  cluster=as.data.table(cluster,keep.rownames='junction')

  print('------into plotting of mfuzz result-----')
  if(show == FALSE){
    pdf(pdf_path, width = 10, height = 8)
    par(mar = c(8,6,3,0.5))
    cos_fuc(eset, cex.lab = 2.5, cex.axis = 1, cex.main = 3, cl=cl, clu = clu,
            ylab ='PSI', mfrow=c(3,3), x11 = F, time.labels = used_celltype)
    cos_bar(seq(0,1,0.01),k=12,horizontal = F,cex.axis = 1)
    dev.off()
  }else{
    par(mar = c(8,6,3,0.5))
    cos_fuc(eset, cex.lab = 2.5, cex.axis = 1, cex.main = 3, cl=cl, clu = clu,
            ylab ='PSI', mfrow=c(3,3), x11 = F, time.labels = used_celltype)
    cos_bar(seq(0,1,0.01),k=12,horizontal = TRUE,cex.axis = 1)
  }
  return(cl)
}


#' @rdname DSUPlot.Mfuzz
#' @name DSUPlot.Mfuzz
#' @title DSUPlot.Mfuzz
cos_fuc <- function (eset, cl, mfrow = c(1, 1), colo, min.mem = 0, time.labels,
                     time.points, ylim.set = c(0, 0), xlab = "", ylab = "PSI",
                     x11 = TRUE, ax.col = "black", bg = "white", col.axis = "black",
                     col.lab = "black", col.main = "black", col.sub = "black",
                     col = "black", centre = FALSE, centre.col = "black", centre.lwd = 2,
                     Xwidth = 5, Xheight = 5, single = FALSE,clu = clu, ...)
{
  clusterindex <- cl[[3]]
  memship <- cl[[4]]
  memship[memship < min.mem] <- -1
  colorindex <- integer(dim(exprs(eset))[[1]])
  if (missing(colo)) {
    colo <- colorRampPalette(rev(brewer.pal(n = 8, name = "Spectral")))(50)
  }
  else {
    if (colo == "fancy") {
      fancy.blue <- c(c(255:0), rep(0, length(c(255:0))),
                      rep(0, length(c(255:150))))
      fancy.green <- c(c(0:255), c(255:0), rep(0, length(c(255:150))))
      fancy.red <- c(c(0:255), rep(255, length(c(255:0))),
                     c(255:150))
      colo <- rgb(b = fancy.blue/255, g = fancy.green/255,
                  r = fancy.red/255)
    }
  }
  colorseq <- seq(0, 1, length = length(colo))
  for (j in 1:dim(cl[[1]])[[1]]) {
    if (single)
      j <- single
    tmp <- exprs(eset)[clusterindex == j, , drop = FALSE]
    tmpmem <- memship[clusterindex == j, j]
    if (((j - 1)%%(mfrow[1] * mfrow[2])) == 0 | single) {
      if (x11)
        X11(width = Xwidth, height = Xheight)
      if (sum(clusterindex == j) == 0) {
        ymin <- -1
        ymax <- +1
      }
      else {
        ymin <- min(tmp)
        ymax <- max(tmp)
      }
      if (sum(ylim.set == c(0, 0)) == 2) {
        ylim <- c(ymin, ymax)
      }
      else {
        ylim <- ylim.set
      }
      if (!is.na(sum(mfrow))) {
        par(mfrow = mfrow, bg = bg, col.axis = col.axis,
            col.lab = col.lab, col.main = col.main, col.sub = col.sub,
            col = col)
      }
      else {
        par(bg = bg, col.axis = col.axis, col.lab = col.lab,
            col.main = col.main, col.sub = col.sub, col = col)
      }
      xlim.tmp <- c(1, dim(exprs(eset))[[2]])
      if (!(missing(time.points)))
        xlim.tmp <- c(min(time.points), max(time.points))
      plot.default(x = NA, xlim = xlim.tmp, ylim = ylim,
                   xlab = xlab, ylab = ylab, main = paste("Cluster",j,'(',clu[j,],')'), axes = FALSE, ...)
      if (missing(time.labels) && missing(time.points)) {
        axis(1, 1:dim(exprs(eset))[[2]], c(1:dim(exprs(eset))[[2]]),
             col = ax.col, ...)
        text(1:dim(exprs(eset))[[2]],par("usr")[3] - 0.2,labels=1:dim(exprs(eset))[[2]],srt = 45, pos = 1, xpd = TRUE)
        axis(2, col = ax.col, ...)
      }
      if (missing(time.labels) && !(missing(time.points))) {
        axis(1, time.points, 1:length(time.points), time.points,
             col = ax.col, ...)
        text(1:dim(exprs(eset))[[2]],par("usr")[3] - 0.2,labels=1:dim(exprs(eset))[[2]],srt = 45, pos = 1, xpd = TRUE)
        axis(2, col = ax.col, ...)
      }
      if (missing(time.points) & !(missing(time.labels))) {
        axis(1, 1:dim(exprs(eset))[[2]],labels = FALSE,col = ax.col, ...)
        text(1:dim(exprs(eset))[[2]],par("usr")[3] - 0.2,labels=time.labels,srt = 45, pos = 1, xpd = TRUE)
        axis(2, col = ax.col, ...)
      }
      if (!(missing(time.points)) & !(missing(time.labels))) {
        axis(1, time.points, time.labels, col = ax.col,
             ...)
        text(1:dim(exprs(eset))[[2]],par("usr")[3] - 0.2,labels=time.labels,srt = 45, pos = 1, xpd = TRUE)
        axis(2, col = ax.col, ...)
      }
    }
    else {
      if (sum(clusterindex == j) == 0) {
        ymin <- -1
        ymax <- +1
      }
      else {
        ymin <- min(tmp)
        ymax <- max(tmp)
      }
      if (sum(ylim.set == c(0, 0)) == 2) {
        ylim <- c(ymin, ymax)
      }
      else {
        ylim <- ylim.set
      }
      xlim.tmp <- c(1, dim(exprs(eset))[[2]])
      if (!(missing(time.points)))
        xlim.tmp <- c(min(time.points), max(time.points))
      plot.default(x = NA, xlim = xlim.tmp, ylim = ylim,
                   xlab = xlab, ylab = ylab, main = paste("Cluster",
                                                          j,'(',clu[j,],')'), axes = FALSE, ...)
      if (missing(time.labels) && missing(time.points)) {
        axis(1, 1:dim(exprs(eset))[[2]], c(1:dim(exprs(eset))[[2]]),
             col = ax.col, ...)
        text(1:dim(exprs(eset))[[2]],par("usr")[3] - 0.2,labels=time.labels,srt = 45, pos = 1, xpd = TRUE)
        axis(2, col = ax.col, ...)
      }
      if (missing(time.labels) && !(missing(time.points))) {
        axis(1, time.points, 1:length(time.points), time.points,
             col = ax.col, ...)
        text(1:dim(exprs(eset))[[2]],par("usr")[3] - 0.2,labels=1:dim(exprs(eset))[[2]],srt = 45, pos = 1, xpd = TRUE)
        axis(2, col = ax.col, ...)
      }
      if (missing(time.points) & !(missing(time.labels))) {
        axis(1, 1:dim(exprs(eset))[[2]],labels = FALSE,
             col = ax.col, ...)
        text(1:dim(exprs(eset))[[2]],par("usr")[3] - 0.2,labels=time.labels,srt = 45, pos = 1, xpd = TRUE)
        axis(2, col = ax.col, ...)
      }
      if (!(missing(time.points)) & !(missing(time.labels))) {
        axis(1, time.points, time.labels, col = ax.col,
             ...)
        text(1:dim(exprs(eset))[[2]],par("usr")[3] - 0.2,labels=time.labels,srt = 45, pos = 1, xpd = TRUE)
        axis(2, col = ax.col, ...)
      }
    }
    if (length(tmpmem) > 0) {
      for (jj in 1:(length(colorseq) - 1)) {
        tmpcol <- (tmpmem >= colorseq[jj] & tmpmem <=
                     colorseq[jj + 1])
        if (sum(tmpcol) > 0) {
          tmpind <- which(tmpcol)
          for (k in 1:length(tmpind)) {
            if (missing(time.points)) {
              lines(tmp[tmpind[k], ], col = colo[jj])
            }
            else lines(time.points, tmp[tmpind[k], ],
                       col = colo[jj])
          }
        }
      }
      1}
    if (centre) {
      lines(cl[[1]][j, ], col = centre.col, lwd = centre.lwd)
    }
    if (single)
      return()
  }
}


#' @rdname DSUPlot.Mfuzz
#' @name DSUPlot.Mfuzz
#' @title DSUPlot.Mfuzz
cos_bar <- function (x, horizontal = TRUE, cex.axis = cex.axis,col = NULL, scale = 1:length(x),
                     k = 10, ...)
{
  if (is.null(col) == TRUE) {
    col <- colorRampPalette(rev(brewer.pal(n = 8, name = "Spectral")))(50)
  }
  if (is.numeric(x)) {
    x <- x
    colmap <- col
  }
  else {
    colmap <- x
    low <- range(scale)[1]
    high <- range(scale)[2]
    x <- seq(low, high, length = length(x))
  }
  if (length(x) > k)
    x.small <- seq(x[1], x[length(x)], length = k)
  else x.small <- x
  if (horizontal) {
    image(x, 1, matrix(x, length(x), 1), axes = FALSE, xlab = "",
          ylab = "", col = colmap, ...)
    axis(1, at = rev(x.small),cex.axis = cex.axis, labels = signif(rev(x.small),
                                                                   2), srt = 270)
  }
  if (!horizontal) {
    image(1, x, matrix(x, 1, length(x)), axes = FALSE, xlab = "",
          ylab = "", col = colmap, ...)
    par(las = 1)
    axis(4, at = rev(x.small),cex.axis = cex.axis, labels = signif(rev(x.small),
                                                                   2))
    par(las = 0)
  }
  box()
}


#' @rdname DSUPlot.Mfuzz
#' @name DSUPlot.Mfuzz
#' @title DSUPlot.Mfuzz
standardise <- function (eset)
{
  data <- exprs(eset)
  for (i in 1:dim(data)[[1]]) {
    data[i, ] <- (data[i, ] - mean(data[i, ], na.rm = TRUE))/sd(data[i,], na.rm = TRUE)
  }
  eset <- new('ExpressionSet', exprs=data)
  eset
}

