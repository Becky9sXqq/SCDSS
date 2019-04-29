#for all cell type AS events homologous analysis
##LIFTOVER
library(rtracklayer)
chain <- import.chain("/Users/becky9s/Downloads/mm10ToHg38.over.chain")
lapply(PSI_table_list$cell_specific_PSItable,function(x)
unique(rownames(x))) -> cp_GRanges

cell.orders -> names(cp_GRanges)

lapply(cp_GRanges,function(cp_GRanges){
  lapply(strsplit(stringr::str_replace_all(cp_GRanges, "-(?=[[:digit:]])", ":"), "[@:]"), function(x){
as.data.frame(matrix(x, ncol = 4, byrow = T, dimnames = list(NULL, c("seqnames", "start", "end", "strand")))) }) -> tmp
lapply(tmp, function(tmp) {
  tmp$seqnames = paste0('chr',tmp$seqnames)
  return(tmp)
})-> tmp
tmp = do.call(rbind,tmp)
return(tmp)
}) -> cell_specific_GRanges
names(cell_specific_GRanges) = cell.orders

tmp = list()
liftover2human = list()
for(i in 1:17){
lapply(cell_specific_GRanges[[i]],function(x) {
  GenomicRanges::GRanges(seqnames = paste0('chr',as.character(x$seqnames)),
                         ranges = IRanges::IRanges(start = as.integer(as.character(x$start)), end = as.integer(as.character(x$end))),
                         strand = as.character(x$strand))}) -> tmp[[i]]
  lapply(tmp[[i]],function(h) liftOver(h, chain)) -> liftover2human[[i]]
}
##for cell specific liftover in human and mouse
liftover2human

