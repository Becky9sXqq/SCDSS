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
#' @rdname DSUPlot.PCA
#' @name DSUPlot.PCA
#' @importFrom data.table setDT
#' @importFrom data.table setkey
#' @importFrom data.table as.data.table
#' @import RColorBrewer factoextra FactoMineR
#' @param raw.psi an \code{RawPsiTable} result for raw counts PSI table 
#' @param used_cell The cell types user names to plot PCA
#' @param canvas_width the width of canvas
#' @param canvas_length the length of canvas
#' @param pdf_path the pathway for pdf output
#' @param ... futher parameters for PCA plotting
#' @seealso \code{\link[FactoMineR]{PCA}}  \code{\link[DESeq2]{plotPCA}}
#' @example DSUplot.PCA(raw.psi = raw.psi , canvas_width, canvas_length, used_cell = unique(SampleInfo$cell.type), pdf_path)
#' @export
DSUPlot.PCA = function(raw.psi = raw.psi , canvas_width, canvas_length, used_cell , pdf_path)
{
  for(pack in c("RColorBrewer","factoextra","FactoMineR")){
    usePackage(pack)
  }
  if(missing(used_cell)){
    used_cell= unique(SampleInfo$cell.type)
  }
  
  colnames(raw.psi) <- SampleInfo$cell.type
  raw.psi[,colnames(raw.psi) %in% as.vector(used_cell)] -> raw.psi
  colnames(raw.psi) <- paste0(colnames(raw.psi),'_',SampleInfo[SampleInfo$cell.type %in% colnames(raw.psi),]$SampleInfo)
  res_PCA=PCA(t(raw.psi),ncp=5,graph=F)
  pdf(pdf_path,canvas_length,canvas_width)
  p = fviz_pca_ind(res_PCA,xlab=paste0('PC1',':',round(res_PCA$eig[1,2],2),'%'),ylab=paste0('PC2',':',round(res_PCA$eig[2,2],2),'%'),
                    pointsize = "contrib",#Contributions of the Observations to the Component
                    col.ind = "cos2",# Squared Cosines of the Observations shows the importance of a principal component for a given observation
               pointshape = 21, fill = "#B0E0E6",gradient.cols = c("#B0E0E6", "#66B2FF", "blue"),
               repel = TRUE, # Avoid text overlapping (slow if many points)
               subtitle='--by loadings--',
               title='PCA scatter Plot for PSI'
  )
  message(" Plotting ...")
  p
  message(" Plotting finished")
  dev.off()
}
