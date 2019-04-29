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
#' @rdname DSUPlot.GO
#' @name DSUPlot.GO
#' @importFrom data.table setDT
#' @importFrom data.table setkey
#' @importFrom data.table as.data.table
#' @import Biobase clusterProfiler DOSE org.Mm.eg.db org.Hs.eg.db AnnotationDbi
#' @param host_gene an \code{HostGeneforID} result for AS events host genes mapping 
#' @param species.db the genome database for human or mouse geneID switch
#' @param pvaluecutoff the cutoff of pvalue for GO enrichment
#' @param showcategory Numbers of GO pathway terms show out 
#' @param pdf_path the pathway for pdf output
#' @param canvas_len the length of enrichGO graph, for layout adjustment
#' @param canvas_wid the width of enrichGO graph, for layout adjustment
#' @param qvaluecutoff the cutoff of adj_pval in heatmapPlot, which means the adjusted p-values after betabinomial test.
#' @param ... futher parameters of pheatmap.
#' @example DSUplot.GO(host_gene,org.Mm.eg.db,
#' pvaluecutoff=0.05,qvaluecutoff=0.01,showcategory=25,
#' pdf_path='./DSUGO.pdf',canvas_len=10,canvas_wid=10)
#' @return a heatmap or group-specific SJs' data.table
#' @seealso \code{\link[enrichGO]{enrichGO}}
#' @export


myGOslim <- function(x){
  requireNamespace("topGO") || stop("package topGO is required")
  groupGOTerms <- rvcheck::get_fun_from_pkg("topGO", "groupGOTerms")
  annFUN.gene2GO <- rvcheck::get_fun_from_pkg("topGO", "annFUN.gene2GO")
  if (!class(x) %in% c("gseaResult", "enrichResult")) {
    stop("x should be output of gseGO or enrichGO...")
  }
  gene2GO <- inverseList(x@geneSets)
  if (is(x, "gseaResult")) {
    ont <- x@setType
    allgenes <- x@geneList
    core_genes <- unique(unlist(geneInCategory(x)))
    allgenes[!names(allgenes) %in% core_genes] <- -1
    allgenes[core_genes] <- 1
  }
  else {
    ont <- x@ontology
    universe <- x@universe
    allgenes <- numeric(length(universe))
    names(allgenes) <- universe
    allgenes[x@gene] <- 1
  }
  selector <- function(scores) return(scores == 1)
  if (!ont %in% c("BP", "MF", "CC")) {
    stop("ontology should be one of 'BP', 'MF' or 'CC'...")
  }
  pvalue <- x@result$p.adjust
  names(pvalue) <- x@result$ID
  groupGOTerms()
  GOdata <- new("topGOdata", description = "clusterProfiler enrichment results", 
                ontology = ont, 
                allGenes = allgenes, 
                geneSel = selector, 
                annot = annFUN.gene2GO, 
                gene2GO = gene2GO)
  result <- buildLevels(graph(GOdata))
  return(result)
}

DSUPlot.GO= function(host_gene,species.db,pvaluecutoff,qvaluecutoff,showcategory,pdf_path,canvas_len,canvas_wid){
  print('-----loading necessary packages--------')
  for(pack in c('Mfuzz','Biobase','clusterProfiler','DOSE','org.Mm.eg.db','org.Hs.eg.db','AnnotationDbi')){
    useBiocPackage(pack)}
  entreze_gene=mapIds(species.db,host_gene,keytype="ENSEMBL", column="ENTREZID")
  print('--------starting to enrich GO--------')
  GO_result=enrichGO(gene = unique(entreze_gene),
                     OrgDb  = species.db,
                     ont  = "ALL", pAdjustMethod = 'BH',
                     pvalueCutoff  = pvaluecutoff,qvalueCutoff  = qvaluecutoff)
  print('-------staring to plot------')
  pdf(pdf_path,canvas_len,canvas_wid)
  p=dotplot(GO_result, showCategory=showcategory)
  print(p)
  dev.off()
  return(GO_result)
}

###GOslim by AS type respectively
###Circle plot for the merged GO terms cts table
SE_BP <- enrichGO(gene = SE,
                  universe      = bkg,
                  OrgDb         = org.Mm.eg.db,
                  keyType       = 'ENSEMBL',
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.05)
AFE_BP <- enrichGO(gene = AFE,
                   universe      = bkg,
                   OrgDb         = org.Mm.eg.db,
                   keyType       = 'ENSEMBL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05)
ALE_BP <- enrichGO(gene = ALE,
                   universe      = bkg,
                   OrgDb         = org.Mm.eg.db,
                   keyType       = 'ENSEMBL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05)
A3SS_BP <- enrichGO(gene = A3SS,
                    universe      = bkg,
                    OrgDb         = org.Mm.eg.db,
                    keyType       = 'ENSEMBL',
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.05)
A5SS_BP <- enrichGO(gene = A5SS,
                    universe      = bkg,
                    OrgDb         = org.Mm.eg.db,
                    keyType       = 'ENSEMBL',
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.05)
MXE_BP <- enrichGO(gene = MXE,
                   universe      = bkg,
                   OrgDb         = org.Mm.eg.db,
                   keyType       = 'ENSEMBL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05)
############
SE_CC <- enrichGO(gene = SE,
                  universe      = bkg,
                  OrgDb         = org.Mm.eg.db,
                  keyType       = 'ENSEMBL',
                  ont           = "CC",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.05)
AFE_CC <- enrichGO(gene = AFE,
                   universe      = bkg,
                   OrgDb         = org.Mm.eg.db,
                   keyType       = 'ENSEMBL',
                   ont           = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05)
ALE_CC <- enrichGO(gene = ALE,
                   universe      = bkg,
                   OrgDb         = org.Mm.eg.db,
                   keyType       = 'ENSEMBL',
                   ont           = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05)
A3SS_CC <- enrichGO(gene = A3SS,
                    universe      = bkg,
                    OrgDb         = org.Mm.eg.db,
                    keyType       = 'ENSEMBL',
                    ont           = "CC",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.05)
A5SS_CC <- enrichGO(gene = A5SS,
                    universe      = bkg,
                    OrgDb         = org.Mm.eg.db,
                    keyType       = 'ENSEMBL',
                    ont           = "CC",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.05)
MXE_CC <- enrichGO(gene = MXE,
                   universe      = bkg,
                   OrgDb         = org.Mm.eg.db,
                   keyType       = 'ENSEMBL',
                   ont           = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05)
#############

SE_CC_GOslim <- myGOslim(SE_CC)
AFE_CC_GOslim <- myGOslim(AFE_CC)
ALE_CC_GOslim <- myGOslim(ALE_CC)
MXE_CC_GOslim <- myGOslim(MXE_CC)
A3SS_CC_GOslim <- myGOslim(A3SS_CC)
A5SS_CC_GOslim <- myGOslim(A5SS_CC)

SE_BP_GOslim <- myGOslim(SE_BP)
AFE_BP_GOslim <- myGOslim(AFE_BP)
ALE_BP_GOslim <- myGOslim(ALE_BP)
MXE_BP_GOslim <- myGOslim(MXE_BP)
A3SS_BP_GOslim <- myGOslim(A3SS_BP)
A5SS_BP_GOslim <- myGOslim(A5SS_BP)

SE_MF_GOslim <- myGOslim(SE_MF)
AFE_MF_GOslim <- myGOslim(AFE_MF)
ALE_MF_GOslim <- myGOslim(ALE_MF)
MXE_MF_GOslim <- myGOslim(MXE_MF)
A3SS_MF_GOslim <- myGOslim(A3SS_MF)
A5SS_MF_GOslim <- myGOslim(A5SS_MF)

all_CC <- Reduce(function(x,y) {merge(x,y,all = T,by="Description")},list(SE_CC_L3 , AFE_CC_L3 , ALE_CC_L3,MXE_CC_L3, A3SS_CC_L3, A5SS_CC_L3))
all_MF <- Reduce(function(x,y) {merge(x,y,all = T,by="Description")},list(SE_MF_L3 , AFE_MF_L3 , ALE_MF_L3,MXE_MF_L3, A3SS_MF_L3, A5SS_MF_L3))
all <- Reduce(function(x,y) {merge(x,y,all = T,by="Description")},list(SE_L3 , AFE_L3 , ALE_L3,MXE_L3, A3SS_L3, A5SS_L3))
all_MF[is.na(all_MF)] = 0
all_MF$cate = rep('MF',dim(all_MF)[1])

all_CC[is.na(all_CC)] = 0
all_CC$cate = rep('CC',dim(all_CC)[1])

all[is.na(all)] = 0
all$cate = rep('BP',dim(all)[1])

all2p = rbind(all,all_CC,all_MF)
all2p[rowSums(all2p[,2:7])>50,] -> all2p
as.matrix(all2p) -> all2pm
all2pm = t(apply(all2pm[,2:7],1,as.integer))
colnames(all2pm) = colnames(all2p)[2:7]
rownames(all2pm) = seq(22)

library(circlize)
library(RColorBrewer)
set.seed(201904)
col_fun = colorRamp2(range(all2pm), c("#FFEEEE", "#FF0000"), transparency = 0.5)
grid.col = c(brewer.pal(12, "Set3"),
             brewer.pal(8, "Set2"),
             brewer.pal(9, "Set1"))[1:28]
names(grid.col) = c(colnames(all2pm), seq(1:22))
circos.par(start.degree = 180, gap.after = c(rep(2, nrow(all2pm)-1), 60, rep(2, ncol(all2pm)-1), 60))
chordDiagram(all2pm, annotationTrack = c("name", "grid"),
             annotationTrackHeight = c(0.03, 0.035),
             transparency = 0.5,  
             grid.col =  grid.col,
             preAllocateTracks = list(track.height = 0.1),
             col = col_fun)

chordDiagram(all2pm,directional=1,
             grid.col=grid.col,
             direction.type=c('diffHeight','arrows'),
             link.arr.type='big.arrow',annotationTrack = c("name", "grid"),
             annotationTrackHeight = c(0.03, 0.035),
             preAllocateTracks=list(track.height=0.1))

lgd_points = ComplexHeatmap::Legend(at = rownames(all2pm),
                    labels = paste0(seq(1:22),':',all2p$Description,'from',all2p$cate),
                    type = "points", 
                    legend_gp = gpar(col = grid.col[7:28]), 
                    background = grid.col[7:28],
                    title_position = "topleft", 
                    title = "GO Term", 
                    nrow = 22)

lgd_list_vertical = packLegend(lgd_points, direction = "vertical")
pushViewport(viewport(x = unit(45, "mm"), y = unit(50, "mm"), 
                      width = grobWidth(lgd_list_vertical), 
                      height = grobHeight(lgd_list_vertical)))
grid.draw(lgd_list_vertical)
upViewport()

#######BP + CC + MF
all2 = rbind(all,all_CC,all_MF)
all3 = as.data.frame(t(apply(all2[,2:7], 1, function(x) (x-mean(x))/sd(x))))
all3$Description = all2$Description
all3$cate = all2$cate
all3_melt = reshape2::melt(all3)  
colnames(all3_melt) = c('Description','Category','Species','Cts')

ggplot(all3_melt, aes(x = Species, y = Description, fill = Cts)) + 
  geom_tile(width = .9, height = .9) +
  scale_fill_gradient2(low = "darkblue", mid = "darkblue", high = "red", space = "Lab", name = "Gene Counts") +
  facet_grid(Category ~ ., scales = "free_y") +
  scale_x_discrete(position = "top") +
  labs( title = "Top 70 enriched GO terms", x = "", y = "")

ggplot(all3_melt, aes(x = Cts, y = Description,color = Category,shape = Species  ,size = Cts)) +
  geom_point(cex= 2) +
  labs(
    title = "Common enriched GO terms",
    x = "Gene Counts",
    y = ""
  ) +
  facet_grid(Category ~ Species, scales = "free_y")


