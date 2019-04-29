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
#' @rdname DSUPlot.River 
#' @name DSUPlot.River
#' @import ggplot2 
#' @param PSI_table Table contians meta info from BetaBinomialDSU
#' @param cells.order The order for hematopoiesis cells river plot, which must be accurately corresponded with the colnames of PSI table
#' @param species The different layout of riverplot for different species,Human OR Mouse
#' @param genes The selected AS events ID for Riverplot, which is constrained less than 8
#' @param white.line.points The parameter for PSI value adjustment 
#' @param white.space The parameter for river plot ribbons gap
#' @param multiplier The parameter for PSI value adjustment 
#' @param ... futher parameters of river plot
#' @example 
#' cell.orders <- c('HSC','MPP','CLP','CMP','GMP','Mono','Neutro','activated_macrophage','macrophage_inf','macrophage','MEP','EB','MK','DC','B','CD4T','CD8T','NK')
#' cell.orders <- c('LT-HSC','ST-HSC','MPP','CLP','CMP','GMP','MEP','EryA','EryB','MK','Gn','Ma','Mo','B','CD4T','CD8T','NK')
#' genes <- rownames(PSI_table)[1:4]
#' DSUPlot.River4humanpsi <- function( PSI_table = PSI_table ,white.line.points =c(),white.space = 0.001,multiplier = 1, cells.order = cells.order, genes)
#' DSUPlot.River4mousepsi <- function( PSI_table = PSI_table ,white.line.points =c(),white.space = 0.001,multiplier = 1, cells.order = cells.order, genes)
#' DSUPlot.River4human <- function( PSI_table = PSI_table ,white.line.points =c(),white.space = 0.001,multiplier = 1, cells.order = cells.order, genes)
#' DSUPlot.River4mouse <- function( PSI_table = PSI_table ,white.line.points =c(),white.space = 0.001,multiplier = 1, cells.order = cells.order, genes)
#' @seealso \code{\link[ggplot2]{ggplot2}}
#' @export
DSUPlot.River4human <- function( PSI_table = PSI_table ,white.line.points =c(),white.space = 0.001,multiplier = 1, cells.order , genes){
  if(length(genes) > 24){
    message('the numbers of genes must be less than 24')
  }else{
    genes=genes
  }
  
  if( isTRUE(PSI_table < 0) || isTRUE(PSI_table > 1)){
    message('the values in PSI table must between 0 and 1')
  }else{
    PSI_table=PSI_table[rownames(PSI_table) %in% genes,]
    PSI_table=PSI_table[,colnames(PSI_table) %in% cells.order]
  }
  cell.type.expression <- as.data.frame(PSI_table)
  cell.type.expression[is.na(cell.type.expression)] <- 0.1
  cell.type.expression$mean <- rowMeans(cell.type.expression)
  cell.type.expression$gene.name <- rownames(cell.type.expression)
  
  cell.type.expression$col <- rev(brewer.pal(ifelse(length(rownames(cell.type.expression))>3,length(rownames(cell.type.expression)), 3), "Paired"))
  
  cell.type.expression <- melt(cell.type.expression, variable_name = "cell.type", id=c("gene.name", "mean", "col"))
  colnames(cell.type.expression)[4] <- 'cell.type'
  
  cell.type.expression$value <- log10(cell.type.expression$value+1)
  
  cell.type.expression$cell.type <- factor(cell.type.expression$cell.type, levels=cells.order)
  
  print(cell.type.expression)
  
  cell.type.expression$gene.name <- factor(cell.type.expression$gene.name, levels = genes)
  cell.type.expression <- cell.type.expression[order(cell.type.expression$gene.name), ]
  
  cum.sum <- data.frame()
  
  for (ct in levels(cell.type.expression$cell.type))
  {
    max.fpkm.sum <- cumsum(rev(cell.type.expression$value[cell.type.expression$cell.type==ct]))
    min.fpkm.sum <- c(0, max.fpkm.sum[1:(length(max.fpkm.sum)-1)])
    
    ### If I want to have space between all of them
    if (length(white.line.points) == 0)
    {
      white.line <- c(0, cumsum(rep(white.space, (length(max.fpkm.sum)-1))))
      max.fpkm.sum <- max.fpkm.sum + white.line
      min.fpkm.sum <- min.fpkm.sum + white.line
    } else
    {
      for (ind in cumsum(rev(white.line.points[-1])))
      {
        max.fpkm.sum[(ind + 1):length(max.fpkm.sum)] <- max.fpkm.sum[(ind + 1):length(max.fpkm.sum)] + white.space
        min.fpkm.sum[(ind + 1):length(min.fpkm.sum)] <- min.fpkm.sum[(ind + 1):length(min.fpkm.sum)] + white.space
      }
    }
    if(length(levels(cell.type.expression$gene.name))>1)
    {
      cum.sum <- rbind(cum.sum, cbind(gene.name = rev(levels(cell.type.expression$gene.name)), cell.type = rep(ct, length(levels(cell.type.expression$gene.name))), cum.fpkm = max.fpkm.sum, min.fpkm = min.fpkm.sum))
    } else
    {
      cum.sum <- rbind(cum.sum, cbind(gene.name = rev(levels(cell.type.expression$gene.name)), cell.type = rep(ct, length(levels(cell.type.expression$gene.name))), cum.fpkm = max.fpkm.sum, min.fpkm = 0))
    }
  }
  
  cell.type.expression <- merge(cell.type.expression, cum.sum, by = c("cell.type", "gene.name"))
  cell.type.expression <- cell.type.expression[order(cell.type.expression$gene.name), ]
  cell.type.expression$cum.fpkm <- multiplier*(as.numeric(as.character(cell.type.expression$cum.fpkm)))
  cell.type.expression$min.fpkm <- multiplier*(as.numeric(as.character(cell.type.expression$min.fpkm)))
  cells.order <- cells.order
  if(!identical(cells.order,c('HSC','MPP','CLP','CMP','GMP','Mono','Neutro','activated_macrophage','macrophage_inf','macrophage','MEP','EB','MK','DC','B','CD4T','CD8T','NK'))){
    stop('cells.order elements must in HSC MPP CLP CMP GMP Mono Neutro activated_macrophage macrophage_inf macrophage MEP EB MK B CD4T CD8T DC NK')
  }
  cell.type.coordinates <- data.frame(stringsAsFactors = F,cbind(cell.type = c('HSC','MPP','MPP','CMP','CMP','CLP','CLP','CLP','CLP','CLP','MEP','MEP','GMP','GMP','GMP','EB','MK','Mono','Mono','Mono','Neutro','DC','DC','CD4T','CD8T','B','NK','activated_macrophage','macrophage_inf','macrophage'),
                                                                 branch=c('A','A','B','A','C','B','N','X','Y','Z','A','D','C','E','F','A','D','F','P','Q','C','E','N','Z','Y','X','B','F','P','Q'),
                                                                 x = c(0,1,1,2,2,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,5,5,5),
                                                                 y = c(0,0,0,-200,-200,300,300,300,300,300,-270,-270,-10,-10,-10,-430,-300,0,0,0,100,200,200,300,400,500,600,-150,-50,50)))
  
  
  text <- data.frame(stringsAsFactors = F,cbind(labels = c('HSC','MPP','CMP','CLP','MEP','GMP','EB','MK','Mono','Neutro','DC','CD4T','CD8T','B','NK','activated_macrophage','macrophage_inf','macrophage'),
                                                x = c(0,1,2,2,3,3,4.1,4.1,4.2,4.2,4.2,4.1,4.1,4.1,4.1,4.8,5,5),
                                                y = c(-10,-10,-220,280,-275,-10,-300,-430,-100,100,200,330,450,520,600,-200,-100,0)))
  text$x <- as.numeric(text$x)
  text$y <- as.numeric(text$y)
  
  for (ct in unique(cell.type.expression$cell.type))
  {
    mean.norm <- (max(cell.type.expression$cum.fpkm[cell.type.expression$cell.type==ct])/2)
    print(paste("For cell type ", ct, "substract: ", mean.norm, sep=""))
    cell.type.expression$cum.fpkm[cell.type.expression$cell.type==ct] <- cell.type.expression$cum.fpkm[cell.type.expression$cell.type==ct] - mean.norm
    cell.type.expression$min.fpkm[cell.type.expression$cell.type==ct] <- cell.type.expression$min.fpkm[cell.type.expression$cell.type==ct] - mean.norm
    text$y[text$labels==ct] <- text$y[text$labels==ct] - mean.norm
  }
  
  data2plot <- merge(cell.type.coordinates, cell.type.expression, by = "cell.type")
  
  for(clmn in c("x", "y", "value", "cum.fpkm",'min.fpkm')) data2plot[, clmn] <- as.numeric(data2plot[, clmn])
  data2plot <- data2plot[order(data2plot$gene.name), ]
  
  h <- ggplot(data=data2plot, aes(x=x, group=gene.name))
  for (i in unique(data2plot$branch))
  {
    h <- h + geom_ribbon(data = data2plot[data2plot$branch==i, ], aes(x = x, ymin= y + min.fpkm, ymax= y + cum.fpkm,fill = gene.name))
  }
  h <- h + annotate('text',x = text$x,y = text$y,label = text$labels, col = 'black', size=5)+
    labs(title = 'Gene Expr River Plot For Human Hematopoiesis')+theme_bw() + labs(x="", y="")+theme(axis.text.x=element_blank())+
    theme(axis.text.y=element_blank())+guides(fill = guide_legend(title = "Gene Name : ",nrow = ceiling(length(levels(data2plot$gene.name))/2)))+
    theme(legend.position='bottom', legend.text = element_text(size = 10, face = "bold"))
  return(h)
}

#' @rdname DSUPlot.River 
#' @name DSUPlot.River
#' @export
DSUPlot.River4humanpsi <- function( PSI_table = PSI_table ,white.line.points =c(),white.space = 0.001,multiplier = 1, cells.order , genes){
  if(length(genes) > 8){
    message('the numbers of genes must be less than 8')
  }else{
    genes=genes
  }
  
  if( isTRUE(PSI_table < 0) || isTRUE(PSI_table > 1)){
    message('the values in PSI table must between 0 and 1')
  }else{
    PSI_table=PSI_table[rownames(PSI_table) %in% genes,]
  }
  cell.type.expression <- as.data.frame(PSI_table)
  cell.type.expression[is.na(cell.type.expression)] <- 0.1
  cell.type.expression$mean <- rowMeans(cell.type.expression)
  cell.type.expression$gene.name <- rownames(cell.type.expression)
  
  cell.type.expression$col <- rev(brewer.pal(ifelse(length(rownames(cell.type.expression))>3,length(rownames(cell.type.expression)), 3), "Paired"))
  
  cell.type.expression <- melt(cell.type.expression, variable_name = "cell.type", id=c("gene.name", "mean", "col"))
  colnames(cell.type.expression)[4] <- 'cell.type'
  
  cell.type.expression$value <- 10*(cell.type.expression$value)
  
  cell.type.expression$cell.type <- factor(cell.type.expression$cell.type, levels=cells.order)
  
  print(cell.type.expression)
  
  cell.type.expression$gene.name <- factor(cell.type.expression$gene.name, levels = genes)
  cell.type.expression <- cell.type.expression[order(cell.type.expression$gene.name), ]
  
  cum.sum <- data.frame()
  
  for (ct in levels(cell.type.expression$cell.type))
  {
    max.fpkm.sum <- cumsum(rev(cell.type.expression$value[cell.type.expression$cell.type==ct]))
    min.fpkm.sum <- c(0, max.fpkm.sum[1:(length(max.fpkm.sum)-1)])
    
    ### If I want to have space between all of them
    if (length(white.line.points) == 0)
    {
      white.line <- c(0, cumsum(rep(white.space, (length(max.fpkm.sum)-1))))
      max.fpkm.sum <- max.fpkm.sum + white.line
      min.fpkm.sum <- min.fpkm.sum + white.line
    } else
    {
      for (ind in cumsum(rev(white.line.points[-1])))
      {
        max.fpkm.sum[(ind + 1):length(max.fpkm.sum)] <- max.fpkm.sum[(ind + 1):length(max.fpkm.sum)] + white.space
        min.fpkm.sum[(ind + 1):length(min.fpkm.sum)] <- min.fpkm.sum[(ind + 1):length(min.fpkm.sum)] + white.space
      }
    }
    if(length(levels(cell.type.expression$gene.name))>1)
    {
      cum.sum <- rbind(cum.sum, cbind(gene.name = rev(levels(cell.type.expression$gene.name)), cell.type = rep(ct, length(levels(cell.type.expression$gene.name))), cum.fpkm = max.fpkm.sum, min.fpkm = min.fpkm.sum))
    } else
    {
      cum.sum <- rbind(cum.sum, cbind(gene.name = rev(levels(cell.type.expression$gene.name)), cell.type = rep(ct, length(levels(cell.type.expression$gene.name))), cum.fpkm = max.fpkm.sum, min.fpkm = 0))
    }
  }
  
  cell.type.expression <- merge(cell.type.expression, cum.sum, by = c("cell.type", "gene.name"))
  cell.type.expression <- cell.type.expression[order(cell.type.expression$gene.name), ]
  cell.type.expression$cum.fpkm <- multiplier*(as.numeric(as.character(cell.type.expression$cum.fpkm)))
  cell.type.expression$min.fpkm <- multiplier*(as.numeric(as.character(cell.type.expression$min.fpkm)))
  cells.order <- cells.order
  if(!identical(cells.order,c('HSC','MPP','CLP','CMP','GMP','Mono','Neutro','activated_macrophage','macrophage_inf','macrophage','MEP','EB','MK','DC','B','CD4T','CD8T','NK'))){
    stop('cells.order elements must in HSC MPP CLP CMP GMP Mono Neutro activated_macrophage macrophage_inf macrophage MEP EB MK B CD4T CD8T DC NK')
  }
  cell.type.coordinates <- data.frame(stringsAsFactors = F,cbind(cell.type = c('HSC','MPP','MPP','CMP','CMP','CLP','CLP','CLP','CLP','CLP','MEP','MEP','GMP','GMP','GMP','EB','MK','Mono','Mono','Mono','Neutro','DC','DC','CD4T','CD8T','B','NK','activated_macrophage','macrophage_inf','macrophage'),
                                                                 branch=c('A','A','B','A','C','B','N','X','Y','Z','A','D','C','E','F','A','D','F','P','Q','C','E','N','Z','Y','X','B','F','P','Q'),
                                                                 x = c(0,1,1,2,2,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,5,5,5),
                                                                 y = c(0,0,0,-200,-200,300,300,300,300,300,-270,-270,-10,-10,-10,-430,-300,0,0,0,100,200,200,300,400,500,600,-150,-50,50)))
  
  
  text <- data.frame(stringsAsFactors = F,cbind(labels = c('HSC','MPP','CMP','CLP','MEP','GMP','EB','MK','Mono','Neutro','DC','CD4T','CD8T','B','NK','activated_macrophage','macrophage_inf','macrophage'),
                                                x = c(0,1,2,2,3,3,4.1,4.1,4,4.2,4.2,4.1,4.1,4.1,4.1,4.5,4.7,4.7),
                                                y = c(-10,-10,-220,280,-275,-10,-300,-430,0,100,200,330,450,520,600,-140,0,150)))
  text$x <- as.numeric(text$x)
  text$y <- as.numeric(text$y)
  
  for (ct in unique(cell.type.expression$cell.type))
  {
    mean.norm <- (max(cell.type.expression$cum.fpkm[cell.type.expression$cell.type==ct])/2)
    print(paste("For cell type ", ct, "substract: ", mean.norm, sep=""))
    cell.type.expression$cum.fpkm[cell.type.expression$cell.type==ct] <- cell.type.expression$cum.fpkm[cell.type.expression$cell.type==ct] - mean.norm
    cell.type.expression$min.fpkm[cell.type.expression$cell.type==ct] <- cell.type.expression$min.fpkm[cell.type.expression$cell.type==ct] - mean.norm
    text$y[text$labels==ct] <- text$y[text$labels==ct] - mean.norm
  }
  
  data2plot <- merge(cell.type.coordinates, cell.type.expression, by = "cell.type")
  
  for(clmn in c("x", "y", "value", "cum.fpkm",'min.fpkm')) data2plot[, clmn] <- as.numeric(data2plot[, clmn])
  data2plot <- data2plot[order(data2plot$gene.name), ]
  
  h <- ggplot(data=data2plot, aes(x=x, group=gene.name))
  for (i in unique(data2plot$branch))
  {
    h <- h + geom_ribbon(data = data2plot[data2plot$branch==i, ], aes(x = x, ymin= y + min.fpkm, ymax= y + cum.fpkm,fill = gene.name))
  }
  h <- h + annotate('text',x = text$x,y = text$y,label = text$labels, col = 'black', size=5)+
    labs(title = 'PSI  River Plot For Human Hematopoiesis')+theme_bw() + labs(x="", y="")+theme(axis.text.x=element_blank())+
    theme(axis.text.y=element_blank())+guides(fill = guide_legend(title = "Gene Name : ",nrow = ceiling(length(levels(data2plot$gene.name))/2)))+
    theme(legend.position='bottom', legend.text = element_text(size = 10, face = "bold"))
  return(h)
}

#' @rdname DSUPlot.River 
#' @name DSUPlot.River
#' @export
DSUPlot.River4mousepsi <- function( PSI_table = PSI_table ,white.line.points =c(),white.space = 0.001,multiplier = 1, cells.order , genes){
  if(length(genes) > 24){
    message('the numbers of genes must be less than 24')
  }else{
    genes=genes
  }
  
  if( isTRUE(PSI_table < 0) || isTRUE(PSI_table > 1)){
    message('the values in PSI table must between 0 and 1')
  }else{
    PSI_table=PSI_table[rownames(PSI_table) %in% genes,]
  }
  cell.type.expression <- as.data.frame(PSI_table)
  cell.type.expression[is.na(cell.type.expression)] <- 0.1
  cell.type.expression$mean <- rowMeans(cell.type.expression)
  cell.type.expression$gene.name <- rownames(cell.type.expression)
  
  cell.type.expression$col <- rev(brewer.pal(ifelse(length(rownames(cell.type.expression))>3,length(rownames(cell.type.expression)), 3), "Paired"))
  
  cell.type.expression <- melt(cell.type.expression, variable_name = "cell.type", id=c("gene.name", "mean", "col"))
  colnames(cell.type.expression)[4] <- 'cell.type'
  
  cell.type.expression$value <- 10*(cell.type.expression$value)
  
  cell.type.expression$cell.type <- factor(cell.type.expression$cell.type, levels=cells.order)
  
  #print(cell.type.expression)
  
  cell.type.expression$gene.name <- factor(cell.type.expression$gene.name, levels = genes)
  cell.type.expression <- cell.type.expression[order(cell.type.expression$gene.name), ]
  
  cum.sum <- data.frame()
  
  for (ct in levels(cell.type.expression$cell.type))
  {
    max.fpkm.sum <- cumsum(rev(cell.type.expression$value[cell.type.expression$cell.type==ct]))
    min.fpkm.sum <- c(0, max.fpkm.sum[1:(length(max.fpkm.sum)-1)])
    
    ### If I want to have space between all of them
    if (length(white.line.points) == 0)
    {
      white.line <- c(0, cumsum(rep(white.space, (length(max.fpkm.sum)-1))))
      max.fpkm.sum <- max.fpkm.sum + white.line
      min.fpkm.sum <- min.fpkm.sum + white.line
    } else
    {
      for (ind in cumsum(rev(white.line.points[-1])))
      {
        max.fpkm.sum[(ind + 1):length(max.fpkm.sum)] <- max.fpkm.sum[(ind + 1):length(max.fpkm.sum)] + white.space
        min.fpkm.sum[(ind + 1):length(min.fpkm.sum)] <- min.fpkm.sum[(ind + 1):length(min.fpkm.sum)] + white.space
      }
    }
    if(length(levels(cell.type.expression$gene.name))>1)
    {
      cum.sum <- rbind(cum.sum, cbind(gene.name = rev(levels(cell.type.expression$gene.name)), cell.type = rep(ct, length(levels(cell.type.expression$gene.name))), cum.fpkm = max.fpkm.sum, min.fpkm = min.fpkm.sum))
    } else
    {
      cum.sum <- rbind(cum.sum, cbind(gene.name = rev(levels(cell.type.expression$gene.name)), cell.type = rep(ct, length(levels(cell.type.expression$gene.name))), cum.fpkm = max.fpkm.sum, min.fpkm = 0))
    }
  }
  
  cell.type.expression <- merge(cell.type.expression, cum.sum, by = c("cell.type", "gene.name"))
  cell.type.expression <- cell.type.expression[order(cell.type.expression$gene.name), ]
  cell.type.expression$cum.fpkm <- multiplier*(as.numeric(as.character(cell.type.expression$cum.fpkm)))
  cell.type.expression$min.fpkm <- multiplier*(as.numeric(as.character(cell.type.expression$min.fpkm)))
  
  
  if(!identical(cells.order,c('LT-HSC','ST-HSC','MPP','CLP','CMP','GMP','MEP','EryA','EryB','MK','Gn','Ma','Mo','B','CD4T','CD8T','NK'))){
    stop('cells.order elements must in LT-HSC ST-HSC MPP CLP CMP GMP MEP EryA EryB MK Gn Ma Mo B CD4T CD8T NK')
  }
  cell.type.coordinates <- data.frame(stringsAsFactors = F,cbind(cell.type = c('LT-HSC','ST-HSC','MPP','MPP','CMP','CMP','CLP','CLP','CLP','CLP','MEP','MEP','GMP','GMP','GMP','EryA','MK','EryB','Gn','Ma','Mo','B','CD4T','CD8T','NK'),
                                                                 branch=c('A','A','A','B','A','C','B','M','N','P','A','D','C','E','F','D','A','D','F','E','C','P','N','M','B'),
                                                                 x = c(0,0.5,1,1,2,2,2,2,2,2,3,3,3,3,3,3.5,4,4,4,4,4,4,4,4,4),
                                                                 y = c(0,0,0,0,-200,-200,250,250,250,250,-250,-250,0,0,0,-200,-300,-150,-50,0,50,200,250,300,350)))
  
  text <- data.frame(stringsAsFactors = F,cbind(labels = c('LT-HSC','ST-HSC','MPP','CLP','CMP','MEP','GMP','EryB','MK','EryA','Gn','Ma','Mo','B','CD4T','CD8T','NK'),
                                                x = c(0, 0.5,1, 2, 2, 3, 3,3.5,4.1,4.1,4.1,4.1,4.1,4.1,4.1,4.1,4.1),
                                                y= c(-40,-40,-40,245,-225,-265,-40,-200,-300,-150,-50,0,50,200,250,300,350)))
  text$x <- as.numeric(text$x)
  text$y <- as.numeric(text$y)
  for (ct in unique(cell.type.expression$cell.type))
  {
    mean.norm <- (max(cell.type.expression$cum.fpkm[cell.type.expression$cell.type==ct])/2)
    #print(paste("For cell type ", ct, "substract: ", mean.norm, sep=""))
    cell.type.expression$cum.fpkm[cell.type.expression$cell.type==ct] <- cell.type.expression$cum.fpkm[cell.type.expression$cell.type==ct] - mean.norm
    cell.type.expression$min.fpkm[cell.type.expression$cell.type==ct] <- cell.type.expression$min.fpkm[cell.type.expression$cell.type==ct] - mean.norm
    text$y[text$labels==ct] <- text$y[text$labels==ct] - mean.norm
  }
  
  data2plot <- merge(cell.type.coordinates, cell.type.expression, by = "cell.type")
  
  for(clmn in c("x", "y", "value", "cum.fpkm",'min.fpkm')) data2plot[, clmn] <- as.numeric(data2plot[, clmn])
  data2plot <- data2plot[order(data2plot$gene.name), ]
  
  m <- ggplot(data=data2plot, aes(x=x, group=gene.name))
  for (i in unique(data2plot$branch))
  {
    m <- m + geom_ribbon(data = data2plot[data2plot$branch==i, ], aes(x = x, ymin= y + min.fpkm, ymax= y + cum.fpkm,fill = gene.name))
  }
  m <- m +
    annotate('text',x = text$x,y = text$y,label = text$labels, col = 'black', size = 3)+
    labs(title = 'PSI Plot For Mouse Hematopoiesis')+theme_bw() + labs(x="", y="")+theme(axis.text.x=element_blank())+
    theme(axis.text.y=element_blank())+guides(fill = guide_legend(title = "Event Name : ",nrow = ceiling(length(levels(data2plot$gene.name))/1)))+
    theme(legend.position='bottom', legend.text = element_text(size = 5, face = "bold"))
  return(m)
}


#' @rdname DSUPlot.River 
#' @name DSUPlot.River
#' @export
DSUPlot.River4mouse <- function( PSI_table = PSI_table ,white.line.points =c(),white.space = 0.001,multiplier = 1, cells.order , genes){
  if(length(genes) > 24){
    message('the numbers of genes must be less than 24')
  }else{
    genes=genes
  }
  
  if( isTRUE(PSI_table < 0) || isTRUE(PSI_table > 1)){
    message('the values in PSI table must between 0 and 1')
  }else{
    PSI_table=PSI_table[rownames(PSI_table) %in% genes,]
    PSI_table=PSI_table[,colnames(PSI_table) %in% cells.order]
  }
  cell.type.expression <- as.data.frame(PSI_table)
  cell.type.expression[is.na(cell.type.expression)] <- 0.1
  cell.type.expression$mean <- rowMeans(cell.type.expression)
  cell.type.expression$gene.name <- rownames(cell.type.expression)
  
  cell.type.expression$col <- rev(brewer.pal(ifelse(length(rownames(cell.type.expression))>3,length(rownames(cell.type.expression)), 3), "Paired"))
  
  cell.type.expression <- melt(cell.type.expression, variable_name = "cell.type", id=c("gene.name", "mean", "col"))
  colnames(cell.type.expression)[4] <- 'cell.type'
  
  cell.type.expression$value <- log10(cell.type.expression$value+1)
  
  cell.type.expression$cell.type <- factor(cell.type.expression$cell.type, levels=cells.order)
  
  #print(cell.type.expression)
  
  cell.type.expression$gene.name <- factor(cell.type.expression$gene.name, levels = genes)
  cell.type.expression <- cell.type.expression[order(cell.type.expression$gene.name), ]
  
  cum.sum <- data.frame()
  
  for (ct in levels(cell.type.expression$cell.type))
  {
    max.fpkm.sum <- cumsum(rev(cell.type.expression$value[cell.type.expression$cell.type==ct]))
    min.fpkm.sum <- c(0, max.fpkm.sum[1:(length(max.fpkm.sum)-1)])
    
    ### If I want to have space between all of them
    if (length(white.line.points) == 0)
    {
      white.line <- c(0, cumsum(rep(white.space, (length(max.fpkm.sum)-1))))
      max.fpkm.sum <- max.fpkm.sum + white.line
      min.fpkm.sum <- min.fpkm.sum + white.line
    } else
    {
      for (ind in cumsum(rev(white.line.points[-1])))
      {
        max.fpkm.sum[(ind + 1):length(max.fpkm.sum)] <- max.fpkm.sum[(ind + 1):length(max.fpkm.sum)] + white.space
        min.fpkm.sum[(ind + 1):length(min.fpkm.sum)] <- min.fpkm.sum[(ind + 1):length(min.fpkm.sum)] + white.space
      }
    }
    if(length(levels(cell.type.expression$gene.name))>1)
    {
      cum.sum <- rbind(cum.sum, cbind(gene.name = rev(levels(cell.type.expression$gene.name)), cell.type = rep(ct, length(levels(cell.type.expression$gene.name))), cum.fpkm = max.fpkm.sum, min.fpkm = min.fpkm.sum))
    } else
    {
      cum.sum <- rbind(cum.sum, cbind(gene.name = rev(levels(cell.type.expression$gene.name)), cell.type = rep(ct, length(levels(cell.type.expression$gene.name))), cum.fpkm = max.fpkm.sum, min.fpkm = 0))
    }
  }
  
  cell.type.expression <- merge(cell.type.expression, cum.sum, by = c("cell.type", "gene.name"))
  cell.type.expression <- cell.type.expression[order(cell.type.expression$gene.name), ]
  cell.type.expression$cum.fpkm <- multiplier*(as.numeric(as.character(cell.type.expression$cum.fpkm)))
  cell.type.expression$min.fpkm <- multiplier*(as.numeric(as.character(cell.type.expression$min.fpkm)))
  
  
  if(!identical(cells.order,c('LT-HSC','ST-HSC','MPP','CLP','CMP','GMP','MEP','EryA','EryB','MK','Gn','Ma','Mo','B','CD4T','CD8T','NK'))){
    stop('cells.order elements must in LT-HSC ST-HSC MPP CLP CMP GMP MEP EryA EryB MK Gn Ma Mo B CD4T CD8T NK')
  }
  cell.type.coordinates <- data.frame(stringsAsFactors = F,cbind(cell.type = c('LT-HSC','ST-HSC','MPP','MPP','CMP','CMP','CLP','CLP','CLP','CLP','MEP','MEP','GMP','GMP','GMP','EryA','MK','EryB','Gn','Ma','Mo','B','CD4T','CD8T','NK'),
                                                                 branch=c('A','A','A','B','A','C','B','M','N','P','A','D','C','E','F','D','A','D','F','E','C','P','N','M','B'),
                                                                 x = c(0,0.5,1,1,2,2,2,2,2,2,3,3,3,3,3,3.5,4,4,4,4,4,4,4,4,4),
                                                                 y = c(0,0,0,0,-200,-200,250,250,250,250,-250,-250,0,0,0,-200,-300,-150,-50,0,50,200,250,300,350)))
  
  text <- data.frame(stringsAsFactors = F,cbind(labels = c('LT-HSC','ST-HSC','MPP','CLP','CMP','MEP','GMP','EryB','MK','EryA','Gn','Ma','Mo','B','CD4T','CD8T','NK'),
                                                x = c(0, 0.5,1, 2, 2, 3, 3,3.5,4.1,4.1,4.1,4.1,4.1,4.1,4.1,4.1,4.1),
                                                y= c(-40,-40,-40,245,-225,-265,-40,-200,-300,-150,-50,0,50,200,250,300,350)))
  
  text$x <- as.numeric(text$x)
  text$y <- as.numeric(text$y)
  for (ct in unique(cell.type.expression$cell.type))
  {
    mean.norm <- (max(cell.type.expression$cum.fpkm[cell.type.expression$cell.type==ct])/2)
    #print(paste("For cell type ", ct, "substract: ", mean.norm, sep=""))
    cell.type.expression$cum.fpkm[cell.type.expression$cell.type==ct] <- cell.type.expression$cum.fpkm[cell.type.expression$cell.type==ct] - mean.norm
    cell.type.expression$min.fpkm[cell.type.expression$cell.type==ct] <- cell.type.expression$min.fpkm[cell.type.expression$cell.type==ct] - mean.norm
    text$y[text$labels==ct] <- text$y[text$labels==ct] - mean.norm
  }
  
  data2plot <- merge(cell.type.coordinates, cell.type.expression, by = "cell.type")
  
  for(clmn in c("x", "y", "value", "cum.fpkm",'min.fpkm')) data2plot[, clmn] <- as.numeric(data2plot[, clmn])
  data2plot <- data2plot[order(data2plot$gene.name), ]
  
  m <- ggplot(data=data2plot, aes(x=x, group=gene.name))
  for (i in unique(data2plot$branch))
  {
    m <- m + geom_ribbon(data = data2plot[data2plot$branch==i, ], aes(x = x, ymin= y + min.fpkm, ymax= y + cum.fpkm,fill = gene.name))
  }
  m <- m +
    annotate('text',x = text$x,y = text$y,label = text$labels, col = 'black', size=3)+
    labs(title = 'Gene Expr Plot For Mouse Hematopoiesis')+theme_bw() + labs(x="", y="")+theme(axis.text.x=element_blank())+
    theme(axis.text.y=element_blank())+guides(fill = guide_legend(title = "Gene Name : ",nrow = ceiling(length(levels(data2plot$gene.name))/1)))+
    theme(legend.position='bottom', legend.text = element_text(size = 10, face = "bold"))
  return(m)
}


DSUPlot.River4cellline<- function( PSI_table = PSI_table ,white.line.points =c(),white.space = 0.001,multiplier = 1, cells.order , genes){
  if(length(genes) > 8){
    message('the numbers of genes must be less than 8')
  }else{
    genes=genes
  }
  
  if( isTRUE(PSI_table < 0) || isTRUE(PSI_table > 1)){
    message('the values in PSI table must between 0 and 1')
  }else{
    PSI_table=PSI_table[rownames(PSI_table) %in% genes,]
  }
  cell.type.expression <- as.data.frame(PSI_table)
  cell.type.expression[is.na(cell.type.expression)] <- 0.1
  cell.type.expression$mean <- rowMeans(cell.type.expression)
  cell.type.expression$gene.name <- rownames(cell.type.expression)
  
  cell.type.expression$col <- rev(brewer.pal(ifelse(length(rownames(cell.type.expression))>3,length(rownames(cell.type.expression)), 3), "Paired"))
  
  cell.type.expression <- melt(cell.type.expression, variable_name = "cell.type", id=c("gene.name", "mean", "col"))
  colnames(cell.type.expression)[4] <- 'cell.type'
  
  cell.type.expression$value <- log10(cell.type.expression$value+1)
  
  cell.type.expression$cell.type <- factor(cell.type.expression$cell.type, levels=cells.order)
  
  #print(cell.type.expression)
  
  cell.type.expression$gene.name <- factor(cell.type.expression$gene.name, levels = genes)
  cell.type.expression <- cell.type.expression[order(cell.type.expression$gene.name), ]
  
  cum.sum <- data.frame()
  
  for (ct in levels(cell.type.expression$cell.type))
  {
    max.fpkm.sum <- cumsum(rev(cell.type.expression$value[cell.type.expression$cell.type==ct]))
    min.fpkm.sum <- c(0, max.fpkm.sum[1:(length(max.fpkm.sum)-1)])
    
    ### If I want to have space between all of them
    if (length(white.line.points) == 0)
    {
      white.line <- c(0, cumsum(rep(white.space, (length(max.fpkm.sum)-1))))
      max.fpkm.sum <- max.fpkm.sum + white.line
      min.fpkm.sum <- min.fpkm.sum + white.line
    } else
    {
      for (ind in cumsum(rev(white.line.points[-1])))
      {
        max.fpkm.sum[(ind + 1):length(max.fpkm.sum)] <- max.fpkm.sum[(ind + 1):length(max.fpkm.sum)] + white.space
        min.fpkm.sum[(ind + 1):length(min.fpkm.sum)] <- min.fpkm.sum[(ind + 1):length(min.fpkm.sum)] + white.space
      }
    }
    if(length(levels(cell.type.expression$gene.name))>1)
    {
      cum.sum <- rbind(cum.sum, cbind(gene.name = rev(levels(cell.type.expression$gene.name)), cell.type = rep(ct, length(levels(cell.type.expression$gene.name))), cum.fpkm = max.fpkm.sum, min.fpkm = min.fpkm.sum))
    } else
    {
      cum.sum <- rbind(cum.sum, cbind(gene.name = rev(levels(cell.type.expression$gene.name)), cell.type = rep(ct, length(levels(cell.type.expression$gene.name))), cum.fpkm = max.fpkm.sum, min.fpkm = 0))
    }
  }
  
  cell.type.expression <- merge(cell.type.expression, cum.sum, by = c("cell.type", "gene.name"))
  cell.type.expression <- cell.type.expression[order(cell.type.expression$gene.name), ]
  cell.type.expression$cum.fpkm <- multiplier*(as.numeric(as.character(cell.type.expression$cum.fpkm)))
  cell.type.expression$min.fpkm <- multiplier*(as.numeric(as.character(cell.type.expression$min.fpkm)))
  
  cell.type.coordinates <- data.frame(stringsAsFactors = F,cbind(cell.type = c('K562','K562','EB','MK'),
                                                                 branch=c('A','B','A','B'),
                                                                 x = c(1,1,2,2),
                                                                 y = c(0,0,-200,200)))
  
  
  text <- data.frame(stringsAsFactors = F,cbind(labels = c('K562','EB','MK'),
                                                x = c(1,2,2),
                                                y = c(-10,-210,180)))
  text$x <- as.numeric(text$x)
  text$y <- as.numeric(text$y)
  
  for (ct in unique(cell.type.expression$cell.type))
  {
    mean.norm <- (max(cell.type.expression$cum.fpkm[cell.type.expression$cell.type==ct])/2)
    #print(paste("For cell type ", ct, "substract: ", mean.norm, sep=""))
    cell.type.expression$cum.fpkm[cell.type.expression$cell.type==ct] <- cell.type.expression$cum.fpkm[cell.type.expression$cell.type==ct] - mean.norm
    cell.type.expression$min.fpkm[cell.type.expression$cell.type==ct] <- cell.type.expression$min.fpkm[cell.type.expression$cell.type==ct] - mean.norm
    text$y[text$labels==ct] <- text$y[text$labels==ct] - mean.norm
  }
  
  data2plot <- merge(cell.type.coordinates, cell.type.expression, by = "cell.type")
  
  for(clmn in c("x", "y", "value", "cum.fpkm",'min.fpkm')) data2plot[, clmn] <- as.numeric(data2plot[, clmn])
  data2plot <- data2plot[order(data2plot$gene.name), ]
  
  h <- ggplot(data=data2plot, aes(x=x, group=gene.name))
  for (i in unique(data2plot$branch))
  {
    h <- h + geom_ribbon(data = data2plot[data2plot$branch==i, ], aes(x = x, ymin= y + min.fpkm, ymax= y + cum.fpkm,fill = gene.name))
  }
  h <- h + annotate('text',x = text$x,y = text$y,label = text$labels, col = 'black', size=5)+
    labs(title = 'PSI  River Plot For Human Hematopoiesis')+theme_bw() + labs(x="", y="")+theme(axis.text.x=element_blank())+
    theme(axis.text.y=element_blank())+guides(fill = guide_legend(title = "Gene Name : ",nrow = ceiling(length(levels(data2plot$gene.name))/2)))+
    theme(legend.position='bottom', legend.text = element_text(size = 10, face = "bold"))
  return(h)
}


DSUPlot.River4cellline_psi<- function( PSI_table = PSI_table ,white.line.points =c(),white.space = 0.001,multiplier = 1, cells.order , genes){
  if(length(genes) > 8){
    message('the numbers of genes must be less than 8')
  }else{
    genes=genes
  }
  
  if( isTRUE(PSI_table < 0) || isTRUE(PSI_table > 1)){
    message('the values in PSI table must between 0 and 1')
  }else{
    PSI_table=PSI_table[rownames(PSI_table) %in% genes,]
  }
  cell.type.expression <- as.data.frame(PSI_table)
  cell.type.expression[is.na(cell.type.expression)] <- 0.1
  cell.type.expression$mean <- rowMeans(cell.type.expression)
  cell.type.expression$gene.name <- rownames(cell.type.expression)
  
  cell.type.expression$col <- rev(brewer.pal(ifelse(length(rownames(cell.type.expression))>3,length(rownames(cell.type.expression)), 3), "Paired"))
  
  cell.type.expression <- melt(cell.type.expression, variable_name = "cell.type", id=c("gene.name", "mean", "col"))
  colnames(cell.type.expression)[4] <- 'cell.type'
  
  cell.type.expression$value <- 10*(cell.type.expression$value)
  
  cell.type.expression$cell.type <- factor(cell.type.expression$cell.type, levels=cells.order)
  
  print(cell.type.expression)
  
  cell.type.expression$gene.name <- factor(cell.type.expression$gene.name, levels = genes)
  cell.type.expression <- cell.type.expression[order(cell.type.expression$gene.name), ]
  
  cum.sum <- data.frame()
  
  for (ct in levels(cell.type.expression$cell.type))
  {
    max.fpkm.sum <- cumsum(rev(cell.type.expression$value[cell.type.expression$cell.type==ct]))
    min.fpkm.sum <- c(0, max.fpkm.sum[1:(length(max.fpkm.sum)-1)])
    
    ### If I want to have space between all of them
    if (length(white.line.points) == 0)
    {
      white.line <- c(0, cumsum(rep(white.space, (length(max.fpkm.sum)-1))))
      max.fpkm.sum <- max.fpkm.sum + white.line
      min.fpkm.sum <- min.fpkm.sum + white.line
    } else
    {
      for (ind in cumsum(rev(white.line.points[-1])))
      {
        max.fpkm.sum[(ind + 1):length(max.fpkm.sum)] <- max.fpkm.sum[(ind + 1):length(max.fpkm.sum)] + white.space
        min.fpkm.sum[(ind + 1):length(min.fpkm.sum)] <- min.fpkm.sum[(ind + 1):length(min.fpkm.sum)] + white.space
      }
    }
    if(length(levels(cell.type.expression$gene.name))>1)
    {
      cum.sum <- rbind(cum.sum, cbind(gene.name = rev(levels(cell.type.expression$gene.name)), cell.type = rep(ct, length(levels(cell.type.expression$gene.name))), cum.fpkm = max.fpkm.sum, min.fpkm = min.fpkm.sum))
    } else
    {
      cum.sum <- rbind(cum.sum, cbind(gene.name = rev(levels(cell.type.expression$gene.name)), cell.type = rep(ct, length(levels(cell.type.expression$gene.name))), cum.fpkm = max.fpkm.sum, min.fpkm = 0))
    }
  }
  
  cell.type.expression <- merge(cell.type.expression, cum.sum, by = c("cell.type", "gene.name"))
  cell.type.expression <- cell.type.expression[order(cell.type.expression$gene.name), ]
  cell.type.expression$cum.fpkm <- multiplier*(as.numeric(as.character(cell.type.expression$cum.fpkm)))
  cell.type.expression$min.fpkm <- multiplier*(as.numeric(as.character(cell.type.expression$min.fpkm)))
  
  cell.type.coordinates <- data.frame(stringsAsFactors = F,cbind(cell.type = c('K562','K562','EB','MK'),
                                                                 branch=c('A','B','A','B'),
                                                                 x = c(1,1,2,2),
                                                                 y = c(0,0,-200,200)))
  
  
  text <- data.frame(stringsAsFactors = F,cbind(labels = c('K562','EB','MK'),
                                                x = c(1,2,2),
                                                y = c(-10,-210,180)))
  text$x <- as.numeric(text$x)
  text$y <- as.numeric(text$y)
  
  for (ct in unique(cell.type.expression$cell.type))
  {
    mean.norm <- (max(cell.type.expression$cum.fpkm[cell.type.expression$cell.type==ct])/2)
    print(paste("For cell type ", ct, "substract: ", mean.norm, sep=""))
    cell.type.expression$cum.fpkm[cell.type.expression$cell.type==ct] <- cell.type.expression$cum.fpkm[cell.type.expression$cell.type==ct] - mean.norm
    cell.type.expression$min.fpkm[cell.type.expression$cell.type==ct] <- cell.type.expression$min.fpkm[cell.type.expression$cell.type==ct] - mean.norm
    text$y[text$labels==ct] <- text$y[text$labels==ct] - mean.norm
  }
  
  data2plot <- merge(cell.type.coordinates, cell.type.expression, by = "cell.type")
  
  for(clmn in c("x", "y", "value", "cum.fpkm",'min.fpkm')) data2plot[, clmn] <- as.numeric(data2plot[, clmn])
  data2plot <- data2plot[order(data2plot$gene.name), ]
  
  h <- ggplot(data=data2plot, aes(x=x, group=gene.name))
  for (i in unique(data2plot$branch))
  {
    h <- h + geom_ribbon(data = data2plot[data2plot$branch==i, ], aes(x = x, ymin= y + min.fpkm, ymax= y + cum.fpkm,fill = gene.name))
  }
  h <- h + annotate('text',x = text$x,y = text$y,label = text$labels, col = 'black', size = 3)+
    labs(title = 'PSI  River Plot For Human Hematopoiesis')+theme_bw() + labs(x="", y="")+theme(axis.text.x=element_blank())+
    theme(axis.text.y=element_blank())+guides(fill = guide_legend(title = "Gene Name : ",nrow = ceiling(length(levels(data2plot$gene.name))/2)))+
    theme(legend.position='bottom', legend.text = element_text(size = 10, face = "bold"))
  return(h)
}

