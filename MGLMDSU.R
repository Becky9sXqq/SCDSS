#' @title MGLMDSU
#' @name MGLMDSU
#' @rdname MGLMDSU
#' @export
MGLMDSU <- function(object , contrast= colData(object)$cell.type , NT = 5,...) 
  {
  
  for(p in c("data.table",'parallel','aod','plyr','dplyr','MGLM')){
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
  
  design <- object@design
  SampleInfo <- as.data.table(colData(object))
  options(warn=1)
  print(paste("Starting PSI fitting at ", print(Sys.time()),sep=""))
cell.psi_MGLM = parallel::mclapply(junc.as,function(sjs){
  print(sjs)
  sjs.names = sjs$names[order(as.numeric(sjs$end)-as.numeric(sjs$start))]
  contrast = contrast
  tab0 <- as.matrix(t(junc[sjs.names, SampleInfo[cell.type %in% contrast,SampleInfo], with = F]))
  tab0 <- t(apply(tab0,1,as.numeric))
  colnames(tab0) <- sjs.names
  tab0<-tab0[rowSums(tab0)>=10, ,drop = FALSE]
  if (sum(tab0)>100 && sum(colSums(tab0)/sum(tab0)>0.05)==2) {
    tab <- tab0[,(colSums(tab0)/sum(tab0))>0.05,drop = FALSE]
    event <- sum((colSums(tab)/sum(tab))>0.05)
    sjs.use <- sjs[sjs$names%in%colnames(tab),]
    sjs.use$chr <- as.numeric(as.character(sjs.use$chr))
    sjs.use$start <- as.numeric(as.character(sjs.use$start))
    sjs.use$end <- as.numeric(as.character(sjs.use$end))
    dist.all<- c(diff(as.numeric(sjs.use$start)),diff(as.numeric(sjs.use$end)))
    dist <- max(c(diff(as.numeric(sjs.use$start)),diff(as.numeric(sjs.use$end))))
    for (i in seq_along(all.vars(design))) {
      assign(paste0("fac", i), tryCatch(as.factor(unlist(SampleInfo[SampleInfo %in% rownames(tab),all.vars(design)[i],with=F])),error=function(e) NA))
      tmp <- paste("fac", seq_along(all.vars(design)), collapse = "+",sep='')
    }
    tab <- as.data.frame(tab)
    cell_PSI_withreplis <- tab/rowSums(tab)
    all_cell_PSI_withreplis = paste0(SampleInfo[SampleInfo %in% rownames(cell_PSI_withreplis),cell.type],':',apply(cell_PSI_withreplis,1,function(x) paste(round(x,3),collapse=':')),collapse = "|")
    tab[,length(colnames(tab))] = rowSums(tab)
    scrip <- paste("fit1=aod::betabin(cbind(tab[,1],tab[,2]-tab[,1]) ~ ", tmp, ", ~ 1 , data = tab , warnings = F , na.action = na.omit)", sep = "")
    tryCatch(eval(parse(text = scrip)), error=function(e) NA)
    fit0 <- tryCatch(aod::betabin(cbind(tab[,1],tab[,2]-tab[,1]) ~ 1, ~ 1, data = tab , warnings = F , na.action = na.omit), error=function(e) NA)
    beta.model <- tryCatch(fitted(fit1), error=function(e) NA)
    cell_PSI <- tryCatch(paste(round(beta.model,3),round(1-beta.model,3), sep= ":"), error=function(e) NA)
    all_cell_PSI <- tryCatch(paste(apply(unique(cbind(as.character(SampleInfo[SampleInfo %in% rownames(tab),cell.type]),cell_PSI)),1,function(x) paste(x,collapse = ":")),collapse = "|"), error=function(e) NA)
    pvalue <- tryCatch(anova(fit1,fit0)@anova.table$`P(> Chi2)`, error=function(e) NA)[2]
  } else if (sum(tab0)>100 && sum((colSums(tab0)/sum(tab0))>0.05)>2) {
    tab <- tab0[,(colSums(tab0)/sum(tab0))>0.05,drop = FALSE]
    event <- sum((colSums(tab)/sum(tab))>0.05)
    sjs.use <- sjs[sjs$names%in%colnames(tab),]
    dist.all <- c(diff(as.numeric(sjs.use$start)),diff(as.numeric(sjs.use$end)))
    dist = paste0(abs(dist.all[abs(dist.all)>0]),collapse=";")
    for (i in seq_along(all.vars(design))) {
      n <- as.character(unlist(SampleInfo[SampleInfo %in% rownames(tab),all.vars(design)[i],with=F]))
      assign_value <- data.frame(runif(length(unique(n)),1,1.01) , row.names = unique(n))
      assign(paste0("fac", i), assign_value[n,])
      tmp <- paste("fac", seq_along(all.vars(design)), collapse = "+",sep='')
    }
    tab <- as.data.frame(tab)
    cell_PSI_withreplis=tab/rowSums(tab)
    all_cell_PSI_withreplis = paste0(SampleInfo[SampleInfo %in% rownames(cell_PSI_withreplis),cell.type],':',apply(cell_PSI_withreplis,1,function(x) paste(round(x,3),collapse=':')),collapse = "|")
    tab$totalcounts <- rowSums(tab)
    for (i in seq_along(colnames(tab)[-1])) {
      assign(paste0("tab", i), tab[,seq_along(colnames(tab)[-1])[i],])
      for4mula <- paste("tab", seq_along(colnames(tab)[-1]), collapse = ",",sep='')
    }
    scrip <- paste("fit1=tryCatch(MGLM::MGLMreg(formula = cbind(" , for4mula , ") ~ totalcounts +", tmp, ",data = tab, dist = 'DM'),error=function(e) NA)", sep = "")
    tryCatch(eval(parse(text = scrip)), error=function(e) NA)
    scrip0 <- paste("fit0=tryCatch(MGLM::MGLMreg(formula = cbind(" , for4mula , ") ~ totalcounts , data = tab, dist = 'DM'),error=function(e) NA)", sep = "")
    tryCatch(eval(parse(text = scrip0)), error=function(e) NA)
    dir.model = tryCatch(fit1@fitted, error=function(e) NA)
    cell_PSI=tryCatch(cbind(SampleInfo[SampleInfo%in%rownames(dir.model),"cell.type"], round(dir.model,3)), error=function(e) NA)
    all_cell_PSI=tryCatch(paste(apply(unique(cell_PSI),1,function(x) paste(x,collapse = ":")),collapse = "|"), error=function(e) NA)
    pvalue <- tryCatch(1-pchisq(2*(fit1@logL-fit0@logL),0.01), error=function(e) NA)
    tab <- rev(rev(tab)[,-1])
  } else {
    tab <- tab0
    event <- ncol(tab)
    sjs.use <- sjs[sjs$names%in%colnames(tab),]
    dist.all <- c(diff(as.numeric(sjs.use$start)),diff(as.numeric(sjs.use$end)))
    dist=paste0(abs(dist.all[abs(dist.all)>0]),collapse=",")
    dist
    fit1=NULL
    fit0=NULL
  }
  if(!is.null(fit1) & !is.null(fit0)){
    return(
      data.frame(     id = paste(colnames(tab),collapse="@"),
                      event = event,
                      pvalue = pvalue,
                      all_cell_PSI = all_cell_PSI,
                      all_cell_PSI_withreplis = all_cell_PSI_withreplis,
                      all_cell_sum = paste0(colSums(tab),collapse = ";"),
                      dist = dist,
                      cell="All_Cell",stringsAsFactors=FALSE)
    )
  } else {
    return(data.frame(
      id=paste(colnames(tab0),collapse="@"),
      event=event,
      pvalue=NA,
      all_cell_PSI=NA,
      all_cell_PSI_withreplis=NA,
      all_cell_sum=tryCatch(paste0(colSums(tab0),collapse = ";"), error=function(e) NA),
      dist=dist,
      cell="All_Cell",stringsAsFactors=FALSE))
  }
}, mc.cores = NT)



cell.psi_MGLM <- do.call(rbind,cell.psi_MGLM)
cell.psi_MGLM <- na.omit(cell.psi_MGLM)
SCDSS::MGLMDSU(object) <- cell.psi_MGLM
return(object)
}
