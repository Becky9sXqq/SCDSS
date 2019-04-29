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
#' @rdname DSUPlot.heatmap
#' @name DSUPlot.heatmap
#' @importFrom data.table setDT
#' @importFrom data.table setkey
#' @importFrom data.table as.data.table
#' @import pheatmap 
#' @param PSI_table an \code{SCASDataSet-class} cell specific object obtained from BetaBinomialDSU.
#' @param annotation_col The annotation for column in heatmap 
#' @param annotation_row The annotation for row in heatmap 
#' @param main The main title for heatmap output description
#' @param pdf_path the pathway for pdf output
#' @param ... futher parameters of pheatmap.
#' @return a heatmap or group-specific SJs' data.table
#' @example DSUplot.heatmap(PSI_table, pdf_path='~/DSUhp.pdf')
#' @seealso \code{\link[pheatmap]{pheatmap}}
#' @export 
DSUPlot.heatmap <- function(PSI_table = PSI_table_list$cell_specific_PSItable,
                            annotation_col = PSI_table_list$annotation_col,
                            annotation_row = PSI_table_list$annotation_row, 
                            main = "Mouse--AS events  853 , p.value = 0.05, del_psi = 0.5 , sd = 0.3, MEAN_psi = 0.3",
                            pdf_path='~/hp.pdf',...){
  for(pack in c('dplyr','reshape2','ggplot2',
                'splines','dendsort','grid','data.table','pheatmap','RColorBrewer')){
    usePackage(pack)}
  for(pack in c('Mfuzz','Biobase','clusterProfiler','DOSE','org.Mm.eg.db','org.Hs.eg.db','AnnotationDbi')){
    useBiocPackage(pack)}
  PSI_table[is.na(PSI_table)] = 0
  PSI_table <- round(PSI_table,2)
  annotation_row <- annotation_row
  rownames(annotation_col) <- colnames(PSI_table)
  hm_palette = colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlGn"))[c(-1,-2,-10,-11)])(100)
  pdf(pdf_path, width = 7, height = 10)
  p = pheatmap(PSI_table, 
               annotation_col = annotation_col,
               annotation_row = annotation_row,
               show_rownames = F, angle_col = "45",
               cluster_row = F,cluster_col = F,
               gaps_row = cumsum(as.numeric( table(annotation_row$cell.type)[as.vector(cell.orders)])),
               breaks = NA,
               color = hm_palette,
               fontsize = 6,
               legend_breaks =NA,
               cellwidth = 15,
               fontsize_col= 6,
               fontsize_row = 3, 
               main = main)
  p
  dev.off()
  raw.psi[rownames(raw.psi) %in% rownames(PSI_table_list$cell_specific_PSItable),] -> PSI_tableFromRawPSi
  colnames(PSI_tableFromRawPSi) <- SampleInfo[colnames(PSI_tableFromRawPSi),'cell.type']
  k <- list()
  for (i in seq_along(colnames(PSI_tableFromRawPSi))) {
    k[[i]] = PSI_tableFromRawPSi[,colnames(PSI_tableFromRawPSi) %in% cell.orders[i]]
  }
  k = do.call(cbind,k)
  k[is.na(k)] = 0
  PSI_tableFromRawPSi_col = data.frame(cell.type = do.call(rbind,strsplit(colnames(k),split = '[.]'))[,1] ,row.names = colnames(k))
  PSI_tableFromRawPSi = k
 pheatmap(PSI_tableFromRawPSi, labels_col = PSI_tableFromRawPSi_col$cell.type,
           annotation_col = PSI_tableFromRawPSi_col,
           annotation_row = annotation_row,
           show_rownames = F, angle_col = "45",
           cluster_row = F,cluster_col = F,
           gaps_row = cumsum(as.numeric( table(annotation_row$cell.type)[as.vector(cell.orders)])),
           breaks = NA,
           color = hm_palette,
           fontsize = 6,
           legend_breaks =NA,
           cellwidth = 3,
           fontsize_col= 3,
           fontsize_row = 3, 
           main = main)
  
  print('---finish plotting pheatmap!!-----')
}



