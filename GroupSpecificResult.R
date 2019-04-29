#' Group-specific SJ plot
#'
#' This function are used to plot the results of \code{GroupSpecific} functions.
#' \code{GroupSepPlot} only useful for SSMD and Hypergeometric results but not
#' BetaBinomial or RandomForest, because only SSMD and Hypergeometric can output
#' gruop-specific SJs, BetaBinomial or RandomForest only give the important SJs
#' which can used to classify different groups. In addition, the function
#' also can output a data.table of all or given group-specific (or
#' group-classify-improtant) SJs that have been calculated.
#'
#' @rdname GroupSpecificResult
#' @name GroupSepPlot
#' @importFrom data.table setDT
#' @importFrom data.table setkey
#' @importFrom data.table as.data.table
#' @import pheatmap
#' @param object an \code{SCAS-class} object after group-specific selection.
#' @param method for \code{GroupSepPlot}: the group-specific result you want to
#' plot or output, only used for SSMD and Hypergeometric, or one of 'SSMD' or
#' 'Hypergeometric'.
#' @param Pvalue the cutoff of SSMD 'P', only consider the SJ which P (wilcox.test) of
#' SSMD result less than Pvalue.
#' @param ssmd the cutoff of SSMD 'SSMD', only consider the SJ which SSMD of
#' SSMD result more than ssmd.
#' @param delta.psi the cutoff of deltaPSI of SSMD, only consider the SJ which deltaPSI
#' more than delta.psi.
#' @param Cell.obse the cutoff of Cell_obse in Hypergeometric \code{\link{HypergeometricForGroupSpecific}} of
#' result, Cell_obse mean that the fraction of observation which PSI more than
#' PSICutoff in Hypergeometric test.
#' @param adj.pval the cutoff of adj_pval in Hypergeometric result, which
#' mean the adjusted p-values after Hypergeometric test.
#' @param fc the cutoff of FC (fold change) in Hypergeometric result, which
#' mean the fold change of tested group's expressed SJ over other groups.
#' @param select the result selection of different methods, if the
#' method you given 'SSMD' and 'Hypergeometric', the \code{select} mean
#' that you want to select the 'intersect' or 'union' of two methods.
#' @param onlyASSJ a logical value, if you just want to consider the
#' alternative splicing SJs or not.
#' @param groupBy The group you have chosen to find specific SJs. We
#' recommend the groupBy you used in \code{GroupSpecific}.
#' @param show_rownames,show_colnames,NA.col The parameters of \link{pheatmap}.
#' @param annotation_col The columns you want to annotate in pheatmap,
#' if NULL, we chose the groupBy.
#' @param fittedORrawPSI for BetaBinomialPlot, plot the 'fitted' PSI or 'raw' PSI
#' @param ... futher parameters of pheatmap.
#'
#' @return a heatmap or group-specific SJs' data.table
#'
#' @seealso \code{\link[pheatmap]{pheatmap}}
#' @export
GroupSepPlot <- function(object,
                         method = c("SSMD", "Hypergeometric"),
                         Pvalue = NULL,
                         ssmd = NULL,
                         delta.psi = NULL,
                         Cell.obse = NULL,
                         adj.pval = NULL,
                         fc = NULL,
                         select = "union",
                         onlyASSJ = FALSE,
                         plot = TRUE,
                         groupBy = NULL,
                         show_rownames = FALSE,
                         show_colnames = FALSE,
                         NA.col = "white",
                         annotation_col = NULL,
                         ...){
  if(!is(object, "SCASDataSet")){
    stop(paste0("object must be a SCASDataSet class."))
  }

  if(is.null(method) || !all(is.element(method, c("SSMD", "Hypergeometric", "BetaBinomial", "RandomForest")))){
    stop("method must be one or more of 'SSMD', 'Hypergeometric', 'BetaBinomial' and 'RandomForest'")
  }

  if( is.element("SSMD", method) ) {
    if( is.null(Pvalue) || is.null(ssmd) || is.null(delta.psi) ) {
      stop("If you want to select SSMD result, all values of 'Pvalue', 'ssmd' and 'deltaPSI' can not be NULL.")
    } else if ( Pvalue > 1 || Pvalue < 0 ) {
      stop( paste0("invalid Pvalue", Pvalue) )
    } else if ( ssmd <= 0 ) {
      warning("ssmd should a positive value.")
    } else if ( delta.psi <= 0 ) {
      warning("delta.psi should a positive value.")
    }
  }

  if( is.element("Hypergeometric", method) ) {
    if( is.null(adj.pval) || is.null(Cell.obse) || is.null(fc) ) {
      stop("If you want to select Hypergeometric result, all values of 'adj.pval', 'Cell.obse' and 'fc' can not be NULL.")
    } else if ( adj.pval > 1 || adj.pval < 0 ) {
      stop(paste0("invalid adj.pval", adj.pval))
    } else if ( Cell.obse > 1 || Cell.obse < 0 ) {
      stop(paste0("invalid Cell.obse", Cell.obse))
    } else if ( fc < 0 ) {
      stop("fc should a positive value.")
    }
  }

  if( is.element("SSMD", method) & is.null(object@SSMD) ) {
    stop("The result of SSMDForGroupSpecific are NULL.")
  }
  if( is.element("Hypergeometric", method) & is.null(object@Hypergeometric) ) {
    stop("The result of HypergeometricForGroupSpecific are NULL.")
  }

  if( is.null(select) & length(method) > 1 ) {
    stop( "select should be one of 'intersect' or 'union'" )
  }

  if( !is.null(select) ) {
    if( !is.element(select, c("union", "intersect")) ) {
      stop( paste0("invalid select", select) )
    }
  }

  if( !is.logical(plot) ) {
    stop("The plot parameter mustbe a logical value.")
  }

  if( !is.null(annotation_col) & all(!is.element(annotation_col, colnames(colData(object)))) ) {
    stop( "invalid annotation_col, and your annotation_col must one or more of colnames of samplein formation (colData)" )
  }

  if( is.element("SSMD", method) & !is.null(object@SSMD) ) {
    SSMD.tab <- object@SSMD
    SSMD.tab <- dplyr::filter(SSMD.tab, P < Pvalue & SSMD > ssmd & deltaPSI > delta.psi)
  } else {
    SSMD.tab <- NULL
  }
  if( is.element("Hypergeometric", method) & !is.null(object@Hypergeometric) ) {
    Hyper.tab <- object@Hypergeometric
    Hyper.tab <- dplyr::filter(Hyper.tab, adj_pval < adj.pval & Cell_obse > Cell.obse & Cell_p/Other_p > fc)
  } else {
    Hyper.tab <- NULL
  }

  if ( plot ) {

    if ( !is.null(SSMD.tab) & !is.null(Hyper.tab) ) {
      SSMD.tab.tu <- SSMD.tab[, c("SJ", "Group", "direction")]
      SSMD.tab.tu$SJ <- as.character(SSMD.tab.tu$SJ)
      Hyper.tab.tu <- Hyper.tab[, c("SJ", "Group", "direction")]
      Hyper.tab.tu$SJ <- as.character(Hyper.tab.tu$SJ)
      if( select == "intersect" ) {
        sep.tab <- SSMD.tab.tu[SSMD.tab.tu$SJ %in% Hyper.tab.tu$SJ, ]
      } else {
        sep.tab <- unique(rbind(SSMD.tab.tu, Hyper.tab.tu))
      }

    } else if ( !is.null(SSMD.tab) & is.null(Hyper.tab) ) {
      sep.tab <- SSMD.tab[, c("SJ", "Group", "direction")]
      sep.tab$SJ <- as.character(sep.tab$SJ)
    } else {
      sep.tab <- Hyper.tab[, c("SJ", "Group", "direction")]
      sep.tab$SJ <- as.character(sep.tab$SJ)
    }

    data.table::setDT(sep.tab)
    sep.tab <- sep.tab[, .SD[.N == 1, ], by = "SJ"]

    if ( onlyASSJ ) {
      as_sj <- union(with(data.frame(object@SameStartSJ), paste0(seqnames, ":", start, "-", end, ":", strand)),
                     with(data.frame(object@SameEndSJ), paste0(seqnames, ":", start, "-", end, ":", strand)))
      message(paste0(round(mean(sep.tab$SJ %in% as_sj)*100, 2), "% of all ", nrow(sep.tab) , " SJs are AS SJ"))
      sep.tab <- sep.tab[SJ %in% as_sj, ]
    }

    same_start <- assay(object, 2)
    row.names(same_start) <- as.character(rowRanges(object))
    same_start <- as.data.frame(same_start)
    psi_same_start <- same_start[sep.tab[direction == "SameStart", SJ], ]

    same_end <- assay(object, 3)
    row.names(same_end) <- as.character(rowRanges(object))
    same_end <- as.data.frame(same_end)
    psi_same_end <- same_end[sep.tab[direction == "SameEnd", SJ], ]
    psi_tab <- rbind.data.frame(psi_same_start, psi_same_end)

    coldata <- as.data.table(as.data.frame(colData(object)), keep.rownames = "ID")

    # group ----
    if(is.null(groupBy)){
      groupBy <- all.vars(design(object))
    } else if (!all(is.element(groupBy, colnames(colData(object))))) {
      stop(paste("groupBy", groupBy, "are not in the colnames of sample information."))
    }

    if (length(groupBy) == 1) {
      eval(parse(text = paste0("data.table::setkey(coldata, ", groupBy, ")")))
    } else {
      eval(parse(text = paste0("data.table::setkey(coldata, ", paste0(groupBy, collapse = ", "), ")")))
    }

    data.table::setkey(sep.tab, Group)
    psi_tab <- psi_tab[sep.tab$SJ, coldata$ID]

    if(is.null(annotation_col)) annotation_col = groupBy
    coldata <- as.data.frame(colData(object))
    col_nano <- subset.data.frame(as.data.frame(colData(object)), select = annotation_col)

    pheatmap::pheatmap(psi_tab,
                       cluster_rows = FALSE,
                       cluster_cols = FALSE,
                       show_rownames = show_rownames,
                       show_colnames = show_colnames,
                       na_col = NA.col,
                       # breaks = breaksList,
                       annotation_col = col_nano,
                       ...)

  } else {
    # return a merged table

    if ( !is.null(SSMD.tab) & !is.null(Hyper.tab) ) {
      colnames(SSMD.tab)[-1] <- paste("ssmd", colnames(SSMD.tab), sep = "_")[-1]
      SSMD.tab$SJ <- as.character(SSMD.tab$SJ)

      colnames(Hyper.tab)[-1] <- paste("Hyper", colnames(Hyper.tab), sep = "_")[-1]
      Hyper.tab$SJ <- as.character(Hyper.tab$SJ)

      if( select == "intersect" ) {
        sep.tab <- merge(SSMD.tab, Hyper.tab, all = FALSE)
      } else {
        sep.tab <- merge(SSMD.tab, Hyper.tab, all = TRUE)
      }

    } else if ( !is.null(SSMD.tab) & is.null(Hyper.tab) ) {
      sep.tab <- SSMD.tab
      sep.tab$SJ <- as.character(sep.tab$SJ)
    } else {
      sep.tab <- Hyper.tab
      sep.tab$SJ <- as.character(sep.tab$SJ)
    }
    data.table::setDT(sep.tab)
    if ( onlyASSJ ) {
      as_sj <- union(with(data.frame(object@SameStartSJ), paste0(seqnames, ":", start, "-", end, ":", strand)),
                     with(data.frame(object@SameEndSJ), paste0(seqnames, ":", start, "-", end, ":", strand)))
      message(paste0(nrow(sep.tab), "(", round(mean(sep.tab$SJ %in% as_sj)*100, 2), "%) of SJs are AS SJ"))
      sep.tab <- sep.tab[SJ %in% as_sj, ]
    }
    return(sep.tab)
  }
}

#' @rdname GroupSpecificResult
#' @export
BetaBinomialPlot <- function(object,
                             Pvalue = 0.05,
                             delta.psi = .1,
                             groupBy = NULL,
                             onlyASSJ = TRUE,
                             fittedORrawPSI = "fitted",
                             plot = TRUE,
                             NA.col = "white",
                             show_rownames = FALSE,
                             show_colnames = TRUE,
                             cluster_cols = FALSE,
                             annotation_col = NULL,
                             ...) {

  if(!is(object, "SCASDataSet")){
    stop(paste0("object must be a SCASDataSet class."))
  }

  if( is.null(Pvalue) | !is.numeric(Pvalue) ) {
    stop( paste0("invalid Pvalue", Pvalue) )
  } else if ( Pvalue > 1 || Pvalue < 0 ) {
    stop( paste0("invalid Pvalue", Pvalue) )
  }

  if( is.null(delta.psi) | !is.numeric(delta.psi) ) {
    stop( paste0("invalid delta.psi", delta.psi) )
  } else if ( delta.psi > 1 || delta.psi <= 0 ) {
    stop( paste0("invalid delta.psi", delta.psi) )
  }

  if( is.null(fittedORrawPSI) | length(fittedORrawPSI) != 1 ) {
    stop( paste0("invalid fittedORrawPSI", fittedORrawPSI) )
  } else if( !is.element(fittedORrawPSI, c("fitted", "raw")) ) {
    stop( paste0("invalid fittedORrawPSI", fittedORrawPSI, "must be one of 'fitted' or 'raw'") )
  }

  if( !is.logical(onlyASSJ) ) {
    stop( "onlyASSJ must be a logical value" )
  }

  if( !is.logical(cluster_cols) ) {
    stop( "cluster_cols must be a logical value" )
  }

  isColor2 <- function(x){
    return(x%in%colors() | grepl("^#(\\d|[a-f]){6,8}$", x, ignore.case=TRUE))
  }
  if( !isColor2(NA.col) | length(NA.col) != 1 ) {
    warning( paste0("invalid NA.col", NA.col) )
    message( "Using default 'white' for NA.col")
    NA.col <- "white"
  }

  if( is.null(object@BetaBinomial) ) {
    stop("The result of BetaBinomialForGroupSpecific are NULL.")
  }

  Sep.tab <- object@BetaBinomial
  Sep.tab <- dplyr::filter(Sep.tab, P < Pvalue & deltaPSI > delta.psi)
  Sep.tab$SJ <- as.character(Sep.tab$SJ)

  if( onlyASSJ ) {
    as_sj <- union(with(data.frame(object@SameStartSJ), paste0(seqnames, ":", start, "-", end, ":", strand)),
                   with(data.frame(object@SameEndSJ), paste0(seqnames, ":", start, "-", end, ":", strand)))
    message(paste0(round(mean(Sep.tab$SJ %in% as_sj)*100, 2), "% of all ", nrow(Sep.tab) , " SJs are AS SJ"))
    Sep.tab <- dplyr::filter(Sep.tab, SJ %in% as_sj)
  }


  if( plot ) {

    fit.psi <- Sep.tab[, !colnames(Sep.tab) %in% c("SJ", "P", "deltaPSI", "direction")]
    row.names(fit.psi) <- Sep.tab$SJ

    fit.psi <- Sep.tab[, !colnames(Sep.tab) %in% c("SJ", "P", "deltaPSI", "direction")]
    row.names(fit.psi) <- Sep.tab$SJ

    if(is.null(groupBy)){
      groupBy <- all.vars(design(object))
    } else if (!all(is.element(groupBy, colnames(colData(object))))) {
      stop(paste("groupBy", groupBy, "are not in the colnames of sample information."))
    }

    # group ----
    if (length(groupBy) == 1) {
      groups <- eval(parse(text = paste0("colData(object)$", groupBy)))
      names(groups) <- row.names(colData(object))
      group_type <- as.character(unique(eval(parse(text = paste0("colData(object)$", groupBy)))))
    } else {
      colData <- colData(object)
      groups <- apply(data.frame(colData(object)[, groupBy], stringsAsFactors = F), 1, paste, collapse = "_")
      group_type <- as.character(unique(groups))
      colData$groups <- groups
    }
    colnames(fit.psi) <- group_type
    tmp <- fit.psi
    tmp[is.na(tmp)] <- -0.01

    hc <- stats::hclust(dist(tmp), "complete")
    hc_c <- stats::hclust(dist(t(tmp)), "complete")

    if( fittedORrawPSI == "fitted" ) {
      ## plot fitted PSI

      if( cluster_cols ) {
        fit.psi_tu <- fit.psi[hc$labels[hc$order], hc_c$labels[hc_c$order]]
      } else {
        fit.psi_tu <- fit.psi[hc$labels[hc$order], order(colnames(fit.psi))]
      }

      pheatmap::pheatmap(fit.psi_tu,
                         show_rownames = show_rownames,
                         show_colnames = show_colnames,
                         cluster_rows = FALSE,
                         na_col = NA.col,
                         cluster_cols = FALSE,
                         ...)
    } else {
      ## plot raw PSI
      data.table::setDT(Sep.tab)
      same_start <- assay(object, 2)
      row.names(same_start) <- as.character(rowRanges(object))
      same_start <- as.data.frame(same_start)
      psi_same_start <- same_start[Sep.tab[direction == "SameStart", SJ], ]

      same_end <- assay(object, 3)
      row.names(same_end) <- as.character(rowRanges(object))
      same_end <- as.data.frame(same_end)
      psi_same_end <- same_end[Sep.tab[direction == "SameEnd", SJ], ]
      psi_tab <- rbind.data.frame(psi_same_start, psi_same_end)

      coldata <- as.data.table(as.data.frame(colData(object)), keep.rownames = "ID")

      # group ----
      if(is.null(groupBy)){
        groupBy <- all.vars(design(object))
      } else if (!all(is.element(groupBy, colnames(colData(object))))) {
        stop(paste("groupBy", groupBy, "are not in the colnames of sample information."))
      }

      if (length(groupBy) == 1) {
        eval(parse(text = paste0("data.table::setkey(coldata, ", groupBy, ")")))
      } else {
        eval(parse(text = paste0("data.table::setkey(coldata, ", paste0(groupBy, collapse = ", "), ")")))
      }


      psi_tab <- psi_tab[hc$labels[hc$order], coldata$ID]

      if(is.null(annotation_col)) annotation_col = groupBy
      coldata <- as.data.frame(colData(object))
      col_nano <- subset.data.frame(as.data.frame(colData(object)), select = annotation_col)

      pheatmap::pheatmap(psi_tab,
                         cluster_rows = FALSE,
                         cluster_cols = FALSE,
                         show_rownames = show_rownames,
                         show_colnames = show_colnames,
                         na_col = NA.col,
                         annotation_col = col_nano,
                         ...)
    }
  } else {
    ## return Sep.tab
    return(Sep.tab)
  }
}

#' @rdname GroupSpecificResult
#' @export
RandomForestPlot <- function(object,
                             groupBy = NULL,
                             onlyASSJ = TRUE,
                             plot = TRUE,
                             NA.col = "white",
                             show_rownames = FALSE,
                             show_colnames = TRUE,
                             cluster_cols = FALSE,
                             annotation_col = NULL,
                             ...) {

  if(!is(object, "SCASDataSet")){
    stop(paste0("object must be a SCASDataSet class."))
  }

  if( !is.logical(onlyASSJ) ) {
    stop( "onlyASSJ must be a logical value" )
  }

  if( !is.logical(cluster_cols) ) {
    stop( "cluster_cols must be a logical value" )
  }

  isColor2 <- function(x){
    return(x%in%colors() | grepl("^#(\\d|[a-f]){6,8}$", x, ignore.case=TRUE))
  }
  if( !isColor2(NA.col) | length(NA.col) != 1 ) {
    warning( paste0("invalid NA.col", NA.col) )
    message( "Using default 'white' for NA.col")
    NA.col <- "white"
  }

  if( is.null(object@RandomForest) ) {
    stop("The result of RandomForestForGroupSpecific are NULL.")
  }

  Sep.tab <- object@RandomForest
  Sep.tab$SJ <- as.character(Sep.tab$SJ)

  if( onlyASSJ ) {
    as_sj <- union(with(data.frame(object@SameStartSJ), paste0(seqnames, ":", start, "-", end, ":", strand)),
                   with(data.frame(object@SameEndSJ), paste0(seqnames, ":", start, "-", end, ":", strand)))
    message(paste0(round(mean(Sep.tab$SJ %in% as_sj)*100, 2), "% of all ", nrow(Sep.tab) , " SJs are AS SJ"))
    Sep.tab <- dplyr::filter(Sep.tab, SJ %in% as_sj)
  }


  if( plot ) {
    data.table::setDT(Sep.tab)
    same_start <- assay(object, 2)
    row.names(same_start) <- as.character(rowRanges(object))
    same_start <- as.data.frame(same_start)
    psi_same_start <- same_start[row.names(same_start) %in% Sep.tab[direction == "SameStart", SJ], ]

    same_end <- assay(object, 3)
    row.names(same_end) <- as.character(rowRanges(object))
    same_end <- as.data.frame(same_end)
    psi_same_end <- same_end[row.names(same_end) %in% Sep.tab[direction == "SameEnd", SJ], ] # Maybe some problem

    psi_tab <- rbind.data.frame(psi_same_start, psi_same_end)

    psi_tab <- psi_tab[apply(psi_tab, 1, sd, na.rm = T) != 0, ]

    tmp <- psi_tab
    tmp[is.na(tmp)] <- -1
    hc <- stats::hclust(dist(tmp), "complete")

    # group ----
    if(is.null(groupBy)){
      groupBy <- all.vars(design(object))
    } else if (!all(is.element(groupBy, colnames(colData(object))))) {
      stop(paste("groupBy", groupBy, "are not in the colnames of sample information."))
    }

    coldata <- as.data.table(as.data.frame(colData(object)), keep.rownames = "ID")
    if (length(groupBy) == 1) {
      eval(parse(text = paste0("data.table::setkey(coldata, ", groupBy, ")")))
    } else {
      eval(parse(text = paste0("data.table::setkey(coldata, ", paste0(groupBy, collapse = ", "), ")")))
    }

    psi_tab <- psi_tab[hc$labels[hc$order], coldata$ID]

    if(is.null(annotation_col)) annotation_col = groupBy
    coldata <- as.data.frame(colData(object))
    col_nano <- subset.data.frame(as.data.frame(colData(object)), select = annotation_col)

    pheatmap::pheatmap(psi_tab,
                       cluster_rows = FALSE,
                       cluster_cols = FALSE,
                       show_rownames = show_rownames,
                       show_colnames = show_colnames,
                       na_col = NA.col,
                       annotation_col = col_nano,
                       ...)

  } else {
    ## return Sep.tab
    return(Sep.tab)
  }
}

