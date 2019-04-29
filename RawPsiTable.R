#' @title RawPsiTable
#' @name RawPsiTable
#' @rdname RawPsiTable
#' @export
RawPsiTable <- function(object , contrast= colData(object)$cell.type , NT = 5,...
) {
for(p in c("data.table",'parallel','plyr','dplyr','VGAM')){
  usePackage(p)
}

if(!is(object, "SCASDataSet")){
  stop(paste0("object must be a SCASDataSet class."))
}

junc <- assay(object,1)
paste(rowRanges(object),sep=':') -> names
junc <- as.data.table(cbind(junctions = names,junc))
setkey(junc,junctions)

SameStartSJ <- object@SameStartSJ
SameEndSJ <- object@SameEndSJ
junc.as <- as.list(c(unique(SameStartSJ$ID),unique(SameEndSJ$ID)))
names(junc.as) <- c(unique(SameStartSJ$ID),unique(SameEndSJ$ID))
junc.as <- lapply(junc.as, function(x) unlist(strsplit(x,split = '@')))
junc.as <- lapply(junc.as, function(sjs) do.call(rbind,strsplit(sjs,split = '[-@:]')))
for(a in seq_along(junc.as)){
  colnames(junc.as[[a]]) = c('chr','start','end','strand')
  junc.as[[a]] = as.data.frame(junc.as[[a]])
  
  if(unique(junc.as[[a]]$strand) =='+'){
    #paste2 <- function(x, y, sep = ":") paste(x, y, sep = sep)
    junc.as[[a]]$names = paste0(junc.as[[a]]$chr,':',junc.as[[a]]$start,'-',junc.as[[a]]$end,':',junc.as[[a]]$strand)
  }else{
    junc.as[[a]]$strand = '-'
    junc.as[[a]]$names = paste0(junc.as[[a]]$chr,':',junc.as[[a]]$start,'-',junc.as[[a]]$end,':',junc.as[[a]]$strand)
  }
}

design <- design(object)
SampleInfo <- as.data.table(colData(object))

options(warn=1)

raw.psi = lapply(junc.as,function(sjs){
  sjs.names = sjs$names[order(as.numeric(sjs$end)-as.numeric(sjs$start))]
  contrast = contrast
  tab0 <- as.matrix(t(junc[sjs.names, SampleInfo[cell.type %in% contrast, SampleInfo], with = F]))
  tab0 <- t(apply(tab0,1,as.numeric))
  colnames(tab0)<- sjs.names
  tab0 <- tab0/rowSums(tab0)
  tab0[is.nan(tab0)] <- 0
  if(dim(tab0)[2]==2){
  res <- as.numeric(tab0[,1])
  colnames(tab0)[1] -> names(res)
  } else{
  res <- as.numeric(as.character(tab0[,which.max(apply(tab0, 2, function(z) mean(as.numeric(z))))]))
  names(which.max(apply(tab0, 2, function(z) mean(as.numeric(z))))) -> names(res)
  }
  return(res)
})

target_junc = vector()
for (i in seq_along(raw.psi)) {
  names(raw.psi[[i]][1]) -> target_junc[i]
  names(raw.psi[[i]]) <- as.character(SampleInfo[cell.type %in% contrast, cell.type])
}
target_junc = unlist(target_junc)
raw.psi = do.call(rbind,raw.psi)
apply(raw.psi,2,function(x) round(as.numeric(x),2)) -> tmp
raw.psi <- data.frame(tmp,row.names = rownames(raw.psi))
colnames(raw.psi) = as.character(SampleInfo[cell.type %in% contrast, cell.type])
raw.psi$target_junc <- target_junc
SCDSS::RawPsiTable(object) <- raw.psi
return(object)
}
  
  
  