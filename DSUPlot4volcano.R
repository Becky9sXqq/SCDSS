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
#' @rdname DSUPlot.volcano
#' @name DSUPlot.volcano
#' @importFrom data.table setDT
#' @importFrom data.table setkey
#' @importFrom data.table as.data.table
#' @import ggplot2
#' @param raw.psi an \code{RawPsiTable} result for raw counts PSI table 
#' @param contrast The assigned cell types to be compared in volcano plot
#' @param varCutoff Cutoff for the variance value of PSI table
#' @param pdf_path the pathway for pdf output
#' @param ... futher parameters of volcano
#' @return a volcano plot graph
#' @seealso \code{\link[DESeq2]{plotMA}}
#' @example DSUplot.volcano(raw.psi , contrast = c('EB','MK') , varCutoff = 0.3, pdf_path)
#' @export
DSUPlot.volcano <- function(raw.psi,contrast,varCutoff,pdf_path){
    require(genefilter)
    colnames(raw.psi) <- SampleInfo$cell.type
    mat <- data.table::as.data.table(raw.psi)
    tab <- mat[,contrast,with=F]
    tab$Delta_PSI <- tab[,contrast[1],with=F]-tab[,contrast[2],with=F]
    tab$Var <- rowVars(tab[,contrast,with=F])
    tab$sig <- as.factor(ifelse(tab$Var > varCutoff & abs(tab$Delta_PSI) >= 0.1,ifelse(tab$Delta_PSI > 0.15 ,'Up','Down'),'NOT_diff'))
    pdf(pdf_path)
    volcano <- ggplot(tab,aes(x= Delta_PSI ,y= 10*Var,color=sig)
    )+geom_point(alpha=0.8,shape=16
    )+scale_color_manual(values=c("#6B3FAE","grey","#D02090")
    )+labs(x= " Delta_PSI ",y= " Variance ",title='Volcano plot for AS events',subtitle = paste0('---  ',contrast[1],'  VS  ',contrast[2],'  ---')
    )+geom_hline(yintercept=1, cex=1, colour="darkgrey", linetype="dashed"
    )+geom_vline(xintercept=c(-0.3,0.3),cex=1, colour="darkgrey", linetype="dashed"
    )+theme(plot.title=element_text(colour='maroon'))+theme(plot.subtitle=element_text(colour='grey'))
    print(volcano)
    dev.off()
}
