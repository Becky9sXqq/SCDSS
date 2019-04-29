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
#' @rdname DSUPlot.distance
#' @name DSUPlot.distance
#' @importFrom data.table setDT
#' @importFrom data.table setkey
#' @importFrom  data.table as.data.table
#' @import pheatmap grDevices RColorBrewer
#' @param raw.psi an \code{RawPsiTable} result for raw counts PSI table 
#' @param pdf_path the pathway for pdf output
#' @param canvas_width the width of canvas
#' @param canvas_length the length of canvas
#' @param fontsize_row by default is 5
#' @param ... futher parameters of sample distance plot
#' @example 
#' DSUplot.distance(raw.psi = raw.psi,canvas_length=15,
#' fontsize_row=8,canvas_width =15,
#' pdf_path="./distance_plot.pdf")
#' @seealso \code{\link[ggplot2]{ggplot2}}
#' @export
DSUPlot.distance=function(raw.psi = raw.psi,canvas_length,fontsize_row=8,canvas_width ,pdf_path="./distance_plot.pdf"){
  require(pheatmap)
  require(RColorBrewer)
  mat=as.matrix(raw.psi)
  sampleDists <- dist(t(mat))
  sampleDistMatrix <- as.matrix(sampleDists)
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "YlOrRd")) )(255)
  print('------starting for pheatmap plotting------------')
  if(is.null(canvas_length)){
    canvas_length = 15
  }
  if(is.null(canvas_width)){
    canvas_width = 15
  }
  pdf(pdf_path, canvas_length , canvas_width)
  p = pheatmap(sampleDistMatrix,main='sample_distance_heatmap_measuredByPSI',
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,labels_row = SampleInfo$cell.type,
           fontsize_row = fontsize_row,
           col=colors)
  p
  dev.off()
}
