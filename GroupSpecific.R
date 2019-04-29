#' group-specific using SSMD
#'
#' For grouped samples, group by cell type, disease status or any other way,
#' we can use these functions to select the group-specific SJs based their PSI.
#'
#' @title GroupSpecific
#' @name GroupSpecific
#' @rdname GroupSpecific
#' @importFrom dplyr filter
#' @importFrom data.table as.data.table
#' @importFrom parallel mclapply
#' @import VGAM
#' @importFrom varSelRF varSelRF
#'
#' @param object an \code{SCAS-class} object.
#' @param EffecObs The minimum fraction of effective observation (non-NA) in
#' every SJ. The default value is 0.1.
#' @param groupBy The group you choose to find specific SJs. And your groupBy
#' must be one of colnames of sample information. Default is the Names of \code{\link{design}}.
#' @param lamda1,lamda2 penalty coefficient of SSMD and deltaPSI.
#' @param NT The number of cores to use, i.e. at most how many child processes will be run simultaneously in
#' parallel.
#' @param maxit for BetaBinomial: the max iteration, default is 10.
#' @param PSICutoff in Hypergeometric, the PSI more than this cutoff was treated as effective SJ.
#' @param ... other
#' @export
SSMDForGroupSpecific <- function(object,
                           EffecObs = .1,
                           groupBy = NULL,
                           lamda1 = 10,
                           lamda2 = 5,
                           NT = 1,
                           ...
){
  if(!is(object, "SCASDataSet")){
    stop(paste0("object must be a SCASDataSet class."))
  }

  if(is.null(groupBy)){
    groupBy <- all.vars(design(object))
  } else if (!all(is.element(groupBy, colnames(colData(object))))) {
    stop(paste("groupBy", groupBy, "are not in the colnames of sample information."))
  }

  if(!is.numeric(lamda1) || !is.numeric(lamda2)) {
    stop("The penalty coefficient of lamda1 and lamda2 must be numeric value.")
  }

  if(!is.numeric(NT) || round(NT) != NT) {
    stop("The cores number of NT must be a integer value.")
  }

  # group ----
  if (length(groupBy) == 1) {
    groups <- eval(parse(text = paste0("colData(object)$", groupBy)))
    names(groups) <- row.names(colData(object))
    group_type <- unique(eval(parse(text = paste0("colData(object)$", groupBy))))
  } else {
    colData <- colData(object)
    groups <- apply(data.frame(colData(object)[, groupBy], stringsAsFactors = F), 1, paste, collapse = "_")
    group_type <- unique(groups)
    colData$groups <- groups
  }

  # same start PSI ----
  count <- assay(object, 1)
  row.names(count) <- as.character(rowRanges(object))
  count <- as.data.frame(count)
  same_start <- assay(object, 2)
  row.names(same_start) <- as.character(rowRanges(object))
  same_start <- as.data.frame(same_start)
  same_start_sub <- same_start[rowSums(!is.na(same_start)) > EffecObs*ncol(same_start), ]

  sep_list_start <- list()
  for(i in seq_along(group_type)){
    group <- group_type[i]
    te <- same_start_sub[rowSums(!is.na(same_start_sub[, groups == group])) >= EffecObs*sum(groups == group), ]
    parallel::mclapply(seq_len(nrow(te)), function(j){
      x <- te[j, ]
      tab <- data.frame(Prior = as.numeric(x),
                        group = as.numeric(groups == group),
                        row.names = names(x))
      SSMD <- tryCatch(diff(as.vector(by(as.numeric(tab$Prior), as.numeric(tab$group), mean, na.rm = T)))/sqrt(sum(as.vector(by(as.numeric(na.omit(tab)$Prior), as.numeric(na.omit(tab)$group), stats::var, na.rm = T)))), error = function(e) NA)
      SSMD <- SSMD*(lamda1**mean(!is.na(as.numeric(x))))
      deltaPSI <- tryCatch(diff(as.vector(by(as.numeric(tab$Prior), as.numeric(tab$group), mean, na.rm = T))), error = function(e) NA)
      deltaPSI <- deltaPSI*((lamda2**mean(!is.na(as.numeric(x))))/lamda2)
      pvalue <- tryCatch(wilcox.test(x = with(tab, Prior[group == 1]), y = with(tab, Prior[group == 0]), alternative = "greater")$p.value, error = function(e) NA)
      res <- tryCatch(data.frame(SJ = row.names(x), Group = group, P = pvalue, SSMD = SSMD, deltaPSI = deltaPSI), error = function(e) NA)
      return(res)
    }, mc.cores = NT) -> tmp
    sep_list_start[[i]] <- data.table::as.data.table(do.call(rbind, tmp))
    names(sep_list_start)[i] <- group
  }

  same_start_sep_list <- lapply(sep_list_start, function(x) {
    dplyr::filter(x, !is.na(SSMD))
  })


  # same end PSI ----
  same_end <- assay(object, 3)
  row.names(same_end) <- as.character(rowRanges(object))
  same_end <- as.data.frame(same_end)
  same_end_sub <- same_end[rowSums(!is.na(same_end)) > EffecObs*ncol(same_end), ]

  sep_list_end <- list()
  for(i in seq_along(group_type)){
    group <- group_type[i]
    te <- same_end_sub[rowSums(!is.na(same_end_sub[, groups == group])) >= EffecObs*sum(groups == group), ]
    parallel::mclapply(seq_len(nrow(te)), function(j){
      x <- te[j, ]
      tab <- data.frame(Prior = as.numeric(x),
                        group = as.numeric(groups == group),
                        row.names = names(x))
      SSMD <- tryCatch(diff(as.vector(by(as.numeric(tab$Prior), as.numeric(tab$group), mean, na.rm = T)))/sqrt(sum(as.vector(by(as.numeric(tab$Prior), as.numeric(tab$group), stats::var, na.rm = T)))), error = function(e) NA)
      SSMD <- SSMD*(lamda1**mean(!is.na(as.numeric(x))))
      deltaPSI <- tryCatch(diff(as.vector(by(as.numeric(tab$Prior), as.numeric(tab$group), mean, na.rm = T))), error = function(e) NA)
      deltaPSI <- deltaPSI*((lamda2**mean(!is.na(as.numeric(x))))/lamda2)
      pvalue <- tryCatch(wilcox.test(x = with(tab, Prior[group == 1]), y = with(tab, Prior[group == 0]), alternative = "greater")$p.value, error = function(e) NA)
      res <- tryCatch(data.frame(SJ = row.names(x), Group = group, P = pvalue, SSMD = SSMD, deltaPSI = deltaPSI), error = function(e) NA)
      return(res)
    }, mc.cores = NT) -> tmp
    sep_list_end[[i]] <- data.table::as.data.table(do.call(rbind, tmp))
    names(sep_list_end)[i] <- group
  }

  same_end_sep_list <- lapply(sep_list_end, function(x) {
    dplyr::filter(x, !is.na(SSMD))
  })


  lapply(seq_along(group_type), function(i) {
    ssi <- same_start_sep_list[[i]]
    ssi$direction <- "SameStart"
    sei <- same_end_sep_list[[i]]
    sei$direction <- "SameEnd"
    tmp <- data.table::as.data.table(rbind(ssi, sei))
    tmp <- tmp[, .SD[which.min(P), ], by = "SJ"]
  }) -> sep_list
  sep_tab <- do.call(rbind, sep_list)
  SCDSS::SSMD(object) <- sep_tab
  return(object)
}

#' @rdname GroupSpecific
#' @export
HypergeometricForGroupSpecific <- function(object,
                                           EffecObs = .1,
                                           PSICutoff = .5,
                                           groupBy = NULL,
                                           ...
) {
  if(!is(object, "SCASDataSet")){
    stop(paste0("object must be a SCASDataSet class."))
  }

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

  # same start PSI ----
  same_start <- assay(object, 2)
  row.names(same_start) <- as.character(rowRanges(object))
  same_start_sub <- same_start[rowSums(!is.na(same_start)) > EffecObs*ncol(same_start), ]
  te <- as.data.frame(same_start_sub)

  sep_list_start <- list()
  for(i in seq_along(group_type)){
    group <- group_type[i]
    tmp <- t(apply(te, 1, function(x){
      k = sum(x[which(groups == group)] >= PSICutoff, na.rm = T)
      D = sum(x >= PSICutoff, na.rm = T)
      n = sum(!is.na(x)) - D
      N = sum(groups[!is.na(x)] == group)
      pval = tryCatch(phyper(k, D, n, N, lower.tail = FALSE), error = function(e) NA)
      Cell_obse = k/sum(groups == group)
      if(k==0) {
        adj_pval <- pval
      } else {
        adj_pval <- pval * k
      }
      enrichment <- c(as.integer(k), D, n, N, pval, adj_pval, Cell_obse)
      names(enrichment) <- c("k", "D", "n", "N", "pval", "adj_pval", "Cell_obse")
      return(enrichment)
    }))
    tmp <- as.data.frame(tmp)
    tmp$SJ <- row.names(tmp)
    tmp$Group <- group
    tmp$Cell_p <- with(tmp, k/N)
    tmp$Bg_p <- with(tmp, D/(D + n))
    tmp$Other_p <- with(tmp, (D - k)/(D + n - N))
    sep_list_start[[i]] <- tmp
    names(sep_list_start)[i] <- group
  }

  # same end PSI ----
  same_end <- assay(object, 3)
  row.names(same_end) <- as.character(rowRanges(object))
  same_end_sub <- same_end[rowSums(!is.na(same_end)) > EffecObs*ncol(same_end), ]
  te <- as.data.frame(same_end_sub)

  sep_list_end <- list()
  for(i in seq_along(group_type)){
    group <- group_type[i]
    tmp <- t(apply(te, 1, function(x){
      k = sum(x[which(groups == group)] >= PSICutoff, na.rm = T)
      D = sum(x >= PSICutoff, na.rm = T)
      n = sum(!is.na(x)) - D
      N = sum(groups[!is.na(x)] == group)
      pval = tryCatch(phyper(k, D, n, N, lower.tail=FALSE), error = function(e) NA)
      Cell_obse = k/sum(groups == group)
      if(k==0) {
        adj_pval <- pval
      } else {
        adj_pval <- pval * k
      }
      enrichment <- c(as.integer(k), D, n, N, pval, adj_pval, Cell_obse)
      names(enrichment) <- c("k", "D", "n", "N", "pval", "adj_pval", "Cell_obse")
      return(enrichment)
    }))
    tmp <- as.data.frame(tmp)
    tmp$SJ <- row.names(tmp)
    tmp$Group <- group
    tmp$Cell_p <- with(tmp, k/N)
    tmp$Bg_p <- with(tmp, D/(D + n))
    tmp$Other_p <- with(tmp, (D - k)/(D + n - N))
    sep_list_end[[i]] <- tmp
    names(sep_list_end)[i] <- group
  }

  lapply(seq_along(group_type), function(i) {
    ssi <- sep_list_start[[i]]
    ssi$direction <- "SameStart"
    sei <- sep_list_end[[i]]
    sei$direction <- "SameEnd"
    tmp <- data.table::as.data.table(rbind(ssi, sei))
    tmp <- tmp[, .SD[which.min(pval), ], by = "SJ"]
  }) -> sep_list

  sep_tab <- do.call(rbind, sep_list)
  SCDSS::Hypergeometric(object) <- sep_tab
  return(object)
}


#' @rdname GroupSpecific
#' @export
BetaBinomialForGroupSpecific <- function(object,
                                         EffecObs = .1,
                                         maxit = 10,
                                         groupBy = NULL,
                                         NT = 1,
                                         ...
) {
  library(VGAM)
  if(!is(object, "SCASDataSet")){
    stop(paste0("object must be a SCASDataSet class."))
  }

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

  # same start PSI ----
  count <- assay(object, 1)
  row.names(count) <- as.character(rowRanges(object))
  count <- as.data.frame(count)
  same_start <- assay(object, 2)
  row.names(same_start) <- as.character(rowRanges(object))
  same_start <- as.data.frame(same_start)
  same_start_sub <- same_start[rowSums(!is.na(same_start)) >= EffecObs*ncol(same_start), ]
  same_start_sub <- same_start_sub[apply(same_start_sub, 1, sd, na.rm = T) != 0, ]
  as_start_sj <- with(data.frame(object@SameStartSJ), paste0(seqnames, ":", start, "-", end, ":", strand))
  same_start_sub <- same_start_sub[row.names(same_start_sub) %in% as_start_sj, ] # only fit alternative splicing sjunction for betabinomial

  parallel::mclapply(seq_len(nrow(same_start_sub)), function(j){
    x <- same_start_sub[j, ]
    sjs <- data.table::as.data.table(object@SameStartSJ[grep(row.names(x), object@SameStartSJ$ID), 1:4])[, paste0(seqnames, ":", start, "-", end, ":", strand)]
    sjsc <- colSums(count[sjs, ])
    tab <- data.frame(N = sjsc,
                      y = as.numeric(count[row.names(x),]),
                      Prior = as.numeric(x),
                      group = groups,
                      row.names = names(x))
    fit <- tryCatch(VGAM::vglm(cbind(y, N-y) ~ group, betabinomial, data = na.omit(tab), trace = FALSE, maxit = maxit), error = function(e) NULL)

    if(!is.null(fit)){
      fitted.psi <- VGAM::fitted(fit)[,1][!duplicated(as.character(na.omit(tab)$group))];
      names(fitted.psi) <- as.character(na.omit(tab)$group)[!duplicated(as.character(na.omit(tab)$group))]
      all.psi <- fitted.psi[match(group_type, names(fitted.psi))]
      names(all.psi) <- group_type
      deltaPSI <- diff(range(all.psi, na.rm = TRUE))
    }else{
      all.psi <- rep(NA, length(group_type))
      names(all.psi) <- group_type
      deltaPSI <- NA
    }
    pvalue <- tryCatch(VGAM::anova.vglm(fit, type = "I")$`Pr(>Chi)`[2], error = function(e) NA)
    # SSMD <- tryCatch(diff(as.vector(by(tab$Prior, tab$group, mean, na.rm = T)))/sqrt(sum(as.vector(by(tab$Prior, tab$group, var, na.rm = T)))), error = function(e) NA)
    # deltaPSI <- tryCatch(diff(as.vector(by(VGAM::fitted(fit)[, 1], na.omit(tab)$group, mean, na.rm = T))), error = function(e) NA)
    return(data.frame(SJ = row.names(x), P = pvalue, t(as.data.frame(all.psi)), deltaPSI = deltaPSI))
  }, mc.cores = NT) -> tmp
  sep_list_start <- data.table::as.data.table(do.call(rbind, tmp))

  # same end PSI ----
  same_end <- assay(object, 3)
  row.names(same_end) <- as.character(rowRanges(object))
  same_end <- as.data.frame(same_end)
  same_end_sub <- same_end[rowSums(!is.na(same_end)) > EffecObs*ncol(same_end), ]
  same_end_sub <- same_end_sub[apply(same_end_sub, 1, sd, na.rm = T) != 0, ] # only fit alternative splicing sjunction for betabinomial
  as_end_sj <- with(data.frame(object@SameEndSJ), paste0(seqnames, ":", start, "-", end, ":", strand))
  same_end_sub <- same_end_sub[row.names(same_end_sub) %in% as_end_sj, ] # only fit alternative splicing sjunction for betabinomial

  parallel::mclapply(seq_len(nrow(same_end_sub)), function(j){
    x <- same_end_sub[j, ]
    sjs <- data.table::as.data.table(object@SameEndSJ[grep(row.names(x), object@SameEndSJ$ID), 1:4])[, paste0(seqnames, ":", start, "-", end, ":", strand)]
    sjsc <- colSums(count[sjs, ])
    tab <- data.frame(N = sjsc,
                      y = as.numeric(count[row.names(x),]),
                      Prior = as.numeric(x),
                      group = groups,
                      row.names = names(x))
    fit <- tryCatch(VGAM::vglm(cbind(y, N-y) ~ group, betabinomial, data = na.omit(tab), trace = FALSE, maxit = maxit), error = function(e) NULL)
    if(!is.null(fit)){
      fitted.psi <- VGAM::fitted(fit)[,1][!duplicated(as.character(na.omit(tab)$group))];
      names(fitted.psi) <- as.character(na.omit(tab)$group)[!duplicated(as.character(na.omit(tab)$group))]
      all.psi <- fitted.psi[match(group_type, names(fitted.psi))]
      names(all.psi) <- group_type
      deltaPSI <- diff(range(all.psi, na.rm = TRUE))
    }else{
      all.psi <- rep(NA, length(group_type))
      names(all.psi) <- group_type
      deltaPSI <- NA
    }

    pvalue <- tryCatch(VGAM::anova.vglm(fit, type = "I")$`Pr(>Chi)`[2], error = function(e) NA)
    # SSMD <- tryCatch(diff(as.vector(by(tab$Prior, tab$group, mean, na.rm = T)))/sqrt(sum(as vector(by(tab$Prior, tab$group, var, na.rm = T)))), error = function(e) NA)
    # deltaPSI <- tryCatch(diff(as.vector(by(VGAM::fitted(fit)[, 1], na.omit(tab)$group, mean, na.rm = T))), error = function(e) NA)
    return(data.frame(SJ = row.names(x), P = pvalue, t(as.data.frame(all.psi)), deltaPSI = deltaPSI))
  }, mc.cores = NT) -> tmp
  sep_list_end <- data.table::as.data.table(do.call(rbind, tmp))

  sep_list_start$direction <- "SameStart"
  sep_list_end$direction <- "SameEnd"

  sep_tab <- data.table::as.data.table(rbind(sep_list_start, sep_list_end))
  sep_tab <- sep_tab[, .SD[which.min(P), ], by = "SJ"]

  SCDSS::BetaBinomial(object) <- sep_tab
  return(object)
}

#' @rdname GroupSpecific
#' @export
RandomForestForGroupSpecific <- function(object,
                                         EffecObs = .1,
                                         groupBy = NULL,
                                         NA.replace = -1,
                                         NT = 1,
                                         ...
) {
  library(VGAM)
  if(!is(object, "SCASDataSet")){
    stop(paste0("object must be a SCASDataSet class."))
  }

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

  # same start PSI ----
  same_start <- assay(object, 2)
  row.names(same_start) <- as.character(rowRanges(object))
  same_start_sub <- same_start[rowSums(!is.na(same_start)) > EffecObs*ncol(same_start), ]
  te <- as.data.frame(same_start_sub)
  DT.replace.NA(te, NA.replace)
  rft <- varSelRF::varSelRF(xdata = t(te), Class = as.factor(groups))
  same_start <- subset.data.frame(data.frame(rft$initialImportances),
                                  row.names(rft$initialImportances) %in% rft$selected.vars)
  same_start$direction <- "SameStart"

  # same end PSI ----
  same_end <- assay(object, 3)
  row.names(same_end) <- as.character(rowRanges(object))
  same_end_sub <- same_end[rowSums(!is.na(same_end)) > EffecObs*ncol(same_end), ]
  te <- as.data.frame(same_end_sub)
  DT.replace.NA(te, NA.replace)
  rft <- varSelRF::varSelRF(xdata = t(te), Class = as.factor(groups))
  same_end <- subset.data.frame(data.frame(rft$initialImportances),
                                row.names(rft$initialImportances) %in% rft$selected.vars)
  same_end$direction <- "SameEnd"

  sep_tab <- data.table::data.table(rbind(same_start, same_end), keep.rownames = "SJ")
  sep_tab <- sep_tab[, .SD[which.max(MeanDecreaseAccuracy), ], by = "SJ"]

  SCDSS::RandomForest(object) <- sep_tab
  return(object)
}

