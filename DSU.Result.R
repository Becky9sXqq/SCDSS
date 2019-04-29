#' DSU result obtain
#'
#' For \code{DSU.Results} we obtain DSU results from cell.psi, the output is the PSI value table calculated by BetaBinomial Method.
#'
#' @title DSU.Results
#' @name DSU.Results
#' @rdname DSU.Results
#' @importFrom data.table as.data.table
#' @importFrom parallel mclapply
#' @param cell.psi function \code{BetaBinomialDSU} abstracts out the SCASDataSet object result.
#' @param p.value The signifcance of difference in AS events among all cell types
#' @param del_psi  The minimum cutoff  of delta psi for PSI table, by default is 0.1
#' @param sd The minimum cutoff  of standard  deviation for PSI table, by default is 0.1
#' @param MEAN_psi The minimum cutoff  of mean value for every PSI table row, by default is 0.1
#' @param cell.orders The biological order for cell type in heatmap
#' @param Method The different processing method for PSI table，by default is
#' 'BBfit', which is the acronym of betabinomial fitting method; another one is Mean
#' @export
DSU.Results <- function(cell.psi ,cell.orders, Method = "BBfit" , 
                        p.value = 0.05 , 
                        del_psi = 0.3 , sd = 0.2 , 
                        MEAN_psi = 0.1){
  options(stringsAsFactors = F)
  PSI_table <-  na.omit(cell.psi)
  PSI_table <-  PSI_table[PSI_table$pvalue < p.value , ]
  PSI_table$event_num <- unlist(lapply(PSI_table$event,function(x) ifelse(x > 2,'more_junc','two_junc')))

  print('--------constructing the psi table--------')
  if(Method == "BBfit"){
  tcl = lapply(seq_len(nrow(PSI_table)), function(i) {
    if(PSI_table[i, 2] == 2){
      res <- as.numeric(mapply(unlist(strsplit(as.character(PSI_table[i, 4]), "\\|")), FUN = function(y) unlist(strsplit(y, ":"))[2]))
      names(res) <- as.character(mapply(unlist(strsplit(as.character(PSI_table[i, 4]), "\\|")), FUN = function(y) unlist(strsplit(y, ":"))[1]))
      data.frame(res=res,names=names(res)) %>% group_by(names) %>% summarise(mean=mean(res)) -> tmp
      res <- as.vector(tmp$mean)
      names(res) <- tmp$names
    } else {
      tmp <- do.call(rbind, lapply(unlist(strsplit(as.character(PSI_table[i, 4]), "\\|")), function(y) unlist(strsplit(y, ":"))))
      tmpo <- data.frame(tmp,  stringsAsFactors=F)
      for (i in seq_along(tmpo)[-1]) {
        bb <- paste0("mean(as.numeric(X",seq_along(tmpo)[-1],'))',collapse = ",")
      }
      scrip <- paste("tmpo <- tmpo %>% group_by(X1,.drop=FALSE) %>%  summarise(", bb, ")", sep = "")
      eval(parse(text = scrip))
      tmp = as.data.frame(tmpo)[,-1]
      res <- as.numeric(as.character(tmp[ ,which.max(apply(tmp, 2, function(z) sd(as.numeric(z))))]))
      names(res) <- as.character(unlist(tmpo[,1]))
    }
    return(res)
  })
    }else if(Method == "Mean"){
    tcl = lapply(seq_len(nrow(PSI_table)), function(i) {
      if(PSI_table[i, 2] == 2){
        res <- as.numeric(mapply(unlist(strsplit(as.character(PSI_table[i, 5]), "\\|")), FUN = function(y) unlist(strsplit(y, ":"))[2]))
        names(res) <- as.character(mapply(unlist(strsplit(as.character(PSI_table[i, 5]), "\\|")), FUN = function(y) unlist(strsplit(y, ":"))[1]))
        data.frame(res=res,names=names(res)) %>% group_by(names) %>% summarise(mean=mean(res)) -> tmp
        res <- as.vector(tmp$mean)
        names(res) <- tmp$names
        } else {
        tmp <- do.call(rbind, lapply(unlist(strsplit(as.character(PSI_table[i, 5]), "\\|")), function(y) unlist(strsplit(y, ":"))))
        tmpo <- data.frame(tmp,  stringsAsFactors=F)
        for (i in seq_along(tmpo)[-1]) {
          assign(paste0("X", i), paste0('mean(as.numeric(X',i,'))'))
          bb <- paste("X",seq_along(tmpo)[-1], collapse = ",",sep='')
        }
        scrip <- paste("tmpo <- tmpo %>% group_by(X1) %>%  summarise(", bb, ")", sep = "") 
        eval(parse(text = scrip))
        tmp = as.data.frame(tmpo)[,-1]
        res <- as.numeric(as.character(tmp[,which.max(apply(tmp, 2, function(z) sd(as.numeric(z))))]))
        names(res) <- as.character(unlist(tmpo[,1]))
         }
      return(res)
    })
    }
  

  PSI_table_toplot <- data.frame(row.names = names(tcl[[which.max(mapply(length, tcl))]]))
  for (i in seq_along(tcl)) {
    PSI_table_toplot[names(tcl[[i]]), i] <- tcl[[i]]
  }
  PSI_table_toplot <- t(PSI_table_toplot)
  rownames(PSI_table_toplot) <- rownames(PSI_table)
  PSI_table_toplot[is.na(PSI_table_toplot)] <- 0
  PSI_table_toplot  <- PSI_table_toplot[,as.vector(cell.orders)]
  
cell_specific_in17 = lapply(cell.orders, function(select_cell) {
    as.factor(ifelse(colnames(PSI_table_toplot) %in% select_cell,1,0)) -> fac_index
    delta_psi = apply(PSI_table_toplot, 1, function(x) diff(c(mean(x[fac_index == 0]),mean(x[fac_index == 1]))))
    SD = apply(PSI_table_toplot, 1, function(x) sd(c(mean(x[fac_index == 0]),mean(x[fac_index == 1]))))
    meanpsi = rowMeans(PSI_table_toplot)
    tsd <- as.data.frame(cbind(PSI_table_toplot , meanpsi= meanpsi , SD = SD , delta_psi = delta_psi))
    MAT <- subset(tsd,meanpsi > MEAN_psi & delta_psi > del_psi &  SD > sd)
    return(MAT[,1:length(cell.orders)])
  })

  annotation_row = data.frame(cell.type = rep(cell.orders,c(unlist(lapply(cell_specific_in17, function(x) dim(x)[1])))) , row.names = rownames(do.call(rbind,cell_specific_in17)))
  cell_specific_PSItable = do.call(rbind,cell_specific_in17)
  annotation_col = data.frame(cell.type = cell.orders, row.names = cell.orders)
  if(missing(cell.orders)){
     cell.orders=unique(SampleInfo$cell.type)
  }
  round(cell_specific_PSItable,2)
  return(list(fitout_PSI_table = PSI_table_toplot,cell_specific_PSItable =cell_specific_in17 ,annotation_col = annotation_col,annotation_row = annotation_row))
}


