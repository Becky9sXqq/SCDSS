#' The SCASDataSet Class
#'
#' The SCASDataSet object is the center of single cell different splice site
#' usage analysis. It stores all information associated with the dataset,
#' including raw data, PSI, annotations, analyes, etc. All that is needed
#' to construct a SCASDataSet object is an Splice Junction (SJ) expression
#' matrix (rows are SJ, columns are cells), and GTF file which have been used
#' in STAR (which we strongly recommend) alignment.
#'
#' Each SCASDataSet object has a number of slots which store information. Key
#' slots to access are listed below.
#'
#' @rdname SCASDataSet
#' @export
#' @import methods
#' @import SummarizedExperiment
#' @import BiocGenerics
#' @importFrom utils packageVersion

#' @slot counts The raw SJ counts are stored in \code{RangedSummarizedExperiment} Container.
#' @slot start.psi The PSI of same start of all SJs are stored in \code{RangedSummarizedExperiment} Container.
#' @slot end.psi The PSI of same end of all SJs are stored in \code{RangedSummarizedExperiment} Container.
#' @slot rowRanges The \code{GRanges} format information of all SJs are stored in \code{RangedSummarizedExperiment} Container.
#' @slot colData The sample information are stored in \code{RangedSummarizedExperiment} Container.
#' @slot design A \code{formula} about your experimental design.
#' @slot HostGene The host gene of all SJs. For novel (unannotated in GTF) SJs
#' we use range overlap method and we give all the covered genes and the genes anchored by the cleavage site at the same time.
#' @slot ASType Alternative splicing classification of all SJs.
#' @slot SameStartSJ All same strat different end SJs.
#' @slot SameEndSJ All same ens different strat SJs.
#' @slot Hypergeometric Hypergeometric test for group-specific SJs.
#' @slot SSMD SSMD for group-specific SJs.
#' @slot BetaBinomial BetaBinomial for group-specific SJs.
#' @slot RandomForest RandomForest for group-specific SJs.
#'
#' @name SCASDataSet
#' @aliases SCASDataSet-class, SCAS-class
#' @exportClass SCASDataSet
.SCASDataSet <- setClass("SCASDataSet",
                         contains = "RangedSummarizedExperiment",
                         slots = representation(
                           design = "ANY",
                           HostGene = "ANY",
                           ASType = "ANY",
                           SameStartSJ = "ANY",
                           SameEndSJ = "ANY",
                           Hypergeometric = "ANY",
                           SSMD = "ANY",
                           BetaBinomial = "ANY",
                           BetaBinomialDSU = 'ANY',
                           MGLMDSU = 'ANY',
                           RawPsiTable = 'ANY',
                           RandomForest = "ANY",
                           project.name = "character"
                         ))

setValidity("SCASDataSet", function(object) {
  if (! ("counts" %in% assayNames(object)) )
    return( "the assays slot must contain a matrix named 'counts'" )
  if (! ("start.psi" %in% assayNames(object)) )
    return( "the assays slot must contain a matrix named 'start.psi'" )
  if (! ("end.psi" %in% assayNames(object)) )
    return( "the assays slot must contain a matrix named 'end.psi'" )

  if ( !is.numeric( assay(object, 1) ) )
    return( "the count data is not numeric" )
  if ( !is.numeric( assay(object, 2) ) )
    return( "the start.psi data is not numeric" )
  if ( !is.numeric( assay(object, 3) ) )
    return( "the end.psi data is not numeric" )

  if ( any( is.na( assay(object, 1) ) ) )
    return( "NA values are not allowed in the count matrix" )
  if ( !is.integer( assay(object, 1) ) )
    return( "the count data is not in integer mode" )
  if ( any( assay(object, 1) < 0 ) )
    return( "the count data contains negative values" )

  if ( any( assay(object, 2) < 0, na.rm = TRUE ) )
    return( "the start.psi data contains negative values" )
  if ( any( assay(object, 3) < 0, na.rm = TRUE ) )
    return( "the end.psi data contains negative values" )


  design <- design(object)
  # 'design' is either a formula
  stopifnot(is(design, "formula"))

  if (is(design, "formula")) {
    designVars <- all.vars(design)
    if (!all(designVars %in% names(colData(object)))) {
      return("all variables in design formula must be columns in colData")
    }
    designVarsClass <- sapply(designVars, function(v) class(colData(object)[[v]]))
    if (any(designVarsClass == "character")) {
      return("variables in design formula are character vectors.
             convert these columns of colData(object) to factors before including in the design formula")
    }
    designFactors <- designVars[designVarsClass == "factor"]
    # levels would duplicate after make.names()
    if (any(sapply(designFactors, function(v) {
      factor.lvls <- levels(colData(object)[[v]])
      factor.nms <- make.names(factor.lvls)
      any(duplicated(factor.nms))
    }))) {
      return("levels of factors in the design have non-unique level names after make.names() is applied.
             best to only use letters and numbers for levels of factors in the design")
    }
    # levels contain characters other than letters, numbers, and underscore
    if (any(sapply(designFactors, function(v) {
      factor.lvls <- levels(colData(object)[[v]])
      any(!grepl("^[A-Za-z0-9_.]+$", factor.lvls))
    }))) {
      # just a warning for now
      message("  Note: levels of factors in the design contain characters other than
  letters, numbers, '_' and '.'. It is recommended (but not required) to use
  only letters, numbers, and delimiters '_' or '.', as these are safe characters
  for column names in R. [This is a message, not an warning or error]")
    }
    } else {
      stop("design is not a formula")
    }
  TRUE
  })


#' SCASDataSet object and constructors
#'
#' \code{SCASDataSet} is a subclass of \code{RangedSummarizedExperiment},
#' used to store the input values, intermediate calculations and results of an
#' analysis of single cell different splice site usage analysis.
#' The \code{SCASDataSet} class
#' enforces non-negative integer values in the "counts" matrix of all SJs stored
#' as the first element in the assay list. And the same start and same end PSI as
#' the second and third element in the assay respectively. It is worth noting that
#' the column names and row names of the three assays must be identical.
#' In addition, a formula which specifies the design of the experiment must be provided.
#' The constructor functions create a SCASDataSet object from various types of input:
#' a \code{RangedSummarizedExperiment}, a matrix, or SJ count files generated by
#' the STAR aligner.
#' See the vignette for examples of construction from different types.
#'
#' Note on the error message "assay colnames() must be NULL or equal colData rownames()":
#' this means that the colnames of countData are different than the rownames of colData.
#' Fix this with: \code{colnames(countData) <- NULL}
#'
#' @param se a \code{RangedSummarizedExperiment} with columns of variables
#' indicating sample information in \code{colData}, \code{GRanges} format data of all SJs
#' in \code{rowRanges} and the counts as the first element in the assays list, which will be renamed "counts".
#' And the same start and same end PSI as the second and third element in the assay respectively.
#' A \code{RangedSummarizedExperiment} object can be
#' generated by the function \code{SummarizedExperiment} in the SummarizedExperiment package.
#' @param design a \code{formula}.
#' the \code{formula} expresses how the counts for each gene
#' depend on the variables in \code{colData}. Many R \code{formula} are valid, e.g., \code{~ group},
#' and designs with multiple variables, e.g., \code{~ group + condition}.
#' See \code{\link{results}} for a variety of designs and how to extract results tables.
#' By default, the functions in this package will use the last variable in the formula for building
#' results tables and plotting. \code{~ 1} can be used for no design, although users need to remember
#' to switch to another design for differential testing.
#' @param gtfFile The gtf file used in alignment.
#' @param NT The number of cores to use, i.e. at most how many child processes will be run simultaneously in
#' parallel.
#' @param path for STAR *SJ.out.tab input: the directory of *SJ.out.tab files.
#' @param pattern for STAR *SJ.out.tab input: the input file names' suffix.
#' @param minSJ for STAR *SJ.out.tab input: the minimum mapped reads of every SJs.
#' @param SSMinSJ for STAR *SJ.out.tab input: the minimum junction reads of every splice site.
#' @param minSJs for STAR *SJ.out.tab input: the minimum of sums of every SJs over all samples.
#' @param uniqueMapOnly for STAR *SJ.out.tab input: logical value indicates whether only the
#' unique mapped reads is counted for subsequent operations.
#' @param SampleInfo for STAR *SJ.out.tab input and Matrix input: a \code{data.frame} with the sample information
#' about each sample. Including group variable(s) used in \code{\link{design}}. Will be stored in \code{colData}
#' @param counts for Matrix input: read counts for all SJs.
#' @param start.psi for Matrix input: all same start PSI of all SJs,
#' the dimension and dimnames are identical with \code{counts}.
#' @param end.psi for Matrix input: all same end PSI of all SJs,
#' the dimension and dimnames are identical with \code{counts}.
#' @param ... arguments provided to \code{SummarizedExperiment} including metadata.
#' If a user wants to store metadata columns about the rows of the countData, but
#' does not have GRanges or GRangesList information, first construct the SCASDataSet
#' without rowRanges and then add the DataFrame with \code{mcols(SCASDataSet)}.
#'
#'
#' @return A SCASDataSet object.
#'
#' @aliases SCASDataSet SCASDataSet-class SCASDataSetFromse SCASDataSetFromSJFiles SCASDataSetFromMatrix
#'
#' @docType class
#'
#' @rdname SCASDataSet
#' @export
SCASDataSetFromse <- function( se,
                               design,
                               gtfFile,
                               NT = 1,
                               ...) {
  if (!is(se, "RangedSummarizedExperiment")) {
    if (is(se, "SummarizedExperiment")) {
      se <- methods::as(se, "RangedSummarizedExperiment")
    } else {
      stop("'se' must be a RangedSummarizedExperiment object")
    }
  }

  if(length(assays(se)) != 3 || !identical(c("counts", "start.psi", "end.psi"), assayNames(se))){
    stop("The lenth of assay list must be 3, and names of assays must counts, start.psi and end.psi respectively.")
  }

  # special tximeta processing
  if ("tximetaInfo" %in% names(metadata(se))) {
    se <- processTximeta(se)
  }
  # before validity check, try to convert assay to integer mode
  if (any(is.na(assay(se, 1))))
    stop("NA values are not allowed in the count matrix")
  if (any(assay(se, 1) < 0)) {
    stop("some values in assay are negative")
  }
  if (!is.integer(assay(se, 1))) {
    if (!is.numeric(assay(se, 1))) {
      stop(paste("counts matrix should be numeric, currently it has mode:", mode(assay(se, 1))))
    }
    if (any(round(assay(se, 1)) != assay(se, 1))) {
      stop("some values in assay are not integers")
    }
    message("converting counts to integer mode")
    mode(assay(se, 1)) <- "integer"
  }

  if (all(assay(se, 1) == 0)) {
    stop("all samples have 0 counts for all SJs. check the counting script.")
  }

  if (all(rowSums(assay(se, 1) == assay(se, 1)[,1]) == ncol(se))) {
    stop("All SJs have equal values for all samples. will not be able to perform following analysis")
  }

  if (all(is.na(assay(se, 2))) || all(is.na(assay(se, 3)))) {
    stop("all samples' PSI is NA for all SJs. check the counting script.")
  }

  if (!any(is.na(assay(se, 2))) || !any(is.na(assay(se, 3)))) {
    warning("There are no missing value (NA) in your PSI table, that seems to be unreasonable.
            check the counting script.")
  }

  if(any(assay(se, 2) > 1, na.rm = T) || any(assay(se, 3) > 1, na.rm = T)){
    stop("There are some PSI more than 1, that seems to be unreasonable.")
  }

  if(any(assay(se, 2) < 0, na.rm = T) || any(assay(se, 3) < 0, na.rm = T)){
    stop("There are some PSI are negative, that seems to be unreasonable.")
  }

  if (any(duplicated(rownames(se)))) {
    warning(sum(duplicated(rownames(se)))," duplicate rownames were renamed by adding numbers")
    rnms <- rownames(se)
    dups <- unique(rnms[duplicated(rnms)])
    for (rn in dups) {
      idx <- which(rnms == rn)
      rnms[idx[-1]] <- paste(rnms[idx[-1]], c(seq_len(length(idx) - 1)), sep=".")
    }
    rownames(se) <- rnms
  }

  if (is(design, "formula")) {
    designVars <- all.vars(design)
    if (!all(designVars %in% names(colData(se)))) {
      stop("all variables in design formula must be columns in colData")
    }
    for (v in designVars) {
      if (any(is.na(colData(se)[[v]])))
        stop(paste("variables in design formula cannot contain NA:", v))
    }

    designVarsClass <- sapply(designVars, function(v) class(colData(se)[[v]])[1])

    if (any(designVarsClass == "character")) {
      warning("some variables in design formula are characters, converting to factors")
      for (v in designVars[designVarsClass == "character"]) {
        colData(se)[[v]] <- factor(colData(se)[[v]])
      }
    }

    if (length(designVars) == 1) {
      var <- colData(se)[[designVars]]
      if (all(var == var[1])) {
        stop("design has a single variable, with all samples having the same value.
             use instead a design of '~ 1'.")
      }
      }

    designVarsNumeric <- sapply(designVars, function(v) is.numeric(colData(se)[[v]]))
    if (any(designVarsNumeric)) {
      warnIntVars <- FALSE
      for (v in designVars[designVarsNumeric]) {
        if (all(colData(se)[[v]] == round(colData(se)[[v]]))) {
          warnIntVars <- TRUE
        }
      }
      if (warnIntVars) {
        message("  the design formula contains a numeric variable with integer values,
                specifying a model with increasing fold change for higher values.
                did you mean for this to be a factor? if so, first convert
                this variable to a factor using the factor() function")
      }
      }

    if (any(designVarsClass == "ordered")) {
      stop("the design formula contains an ordered factor. The internal steps
           do not work on ordered factors as a formula.")
    }

    designFactors <- designVars[designVarsClass == "factor"]
    missingLevels <- sapply(designFactors, function(v) any(table(colData(se)[[v]]) == 0))
    if (any(missingLevels)) {
      message("factor levels were dropped which had no samples")
      for (v in designFactors[missingLevels]) {
        colData(se)[[v]] <- droplevels(colData(se)[[v]])
      }
    }

    singleLevel <- sapply(designFactors, function(v) all(colData(se)[[v]] == colData(se)[[v]][1]))
    if (any(singleLevel)) {
      stop("design contains one or more variables with all samples having the same value,
           remove these variables from the design")
    }
          } else {
            stop("'design' should be a formula")
          }

  if( !file.exists(gtfFile) ) {
    stop(paste0("Your GTF file", gtfFile, "does not exist."))
  }

  if( !is.integer(NT) ) {
    stop("Your NT parameter is invalid. See also parallel:mc.cores")
  }

  # Add columns on the columns
  mcolsCols <- DataFrame(type =  rep("input", ncol(colData(se))),
                         description = rep("", ncol(colData(se))))
  mcols(colData(se)) <- if (is.null(mcols(colData(se)))) {
    mcolsCols
  } else if (all(names(mcols(colData(se))) == c("type","description"))) {
    mcolsCols
  } else {
    cbind(mcols(colData(se)), mcolsCols)
  }

  sjInfo <- data.table::as.data.table(rowRanges(se))
  same_start_sj <- sjInfo[, .SD[.N > 1, ], by = list(seqnames, start, strand)]
  same_start_sj.id <- same_start_sj[, .SD[, paste(seqnames, ":", start, "-", end, ":", strand, collapse = "@", sep = "")], by = c("seqnames", "start", "strand")][, V1]
  same_start_sj$ID <- rep(same_start_sj.id, same_start_sj[, .N, by = c("seqnames", "start", "strand")][, N])
  SameStartSJ <- S4Vectors::DataFrame(same_start_sj[, .(seqnames, start, end, strand, width, motif, annotation, ID)])

  same_end_sj <- sjInfo[, .SD[.N > 1, ], by = list(seqnames, end, strand)]
  same_end_sj.id <- same_end_sj[, .SD[, paste(seqnames, ":", start, "-", end, ":", strand, collapse = "@", sep = "")], by = c("seqnames", "end", "strand")][, V1]
  same_end_sj$ID <- rep(same_end_sj.id, same_end_sj[, .N, by = c("seqnames", "end", "strand")][, N])
  SameEndSJ <- S4Vectors::DataFrame(same_end_sj[, .(seqnames, start, end, strand, width, motif, annotation, ID)])
  
  message(" Starting find AStype and HostGenes ...")
  astypeandhostgene <- AStypeAndHostGene(sjInfo = sjInfo, gtfFile = gtfFile, NT = NT)
  HostGene <- astypeandhostgene$HostGene
  ASType <- astypeandhostgene$AS_Type

  object <- .SCASDataSet(se, design = design, HostGene = HostGene, ASType = ASType, SameStartSJ = SameStartSJ, SameEndSJ = SameEndSJ, ...)

  # now we know we have at least an empty GRanges or GRangesList for rowRanges
  # so we can create a metadata column 'type' for the mcols
  # and we label any incoming columns as 'input'

  # this is metadata columns on the rows
  mcolsRows <- DataFrame(type = rep("input", ncol(mcols(object))),
                         description = rep("", ncol(mcols(object))))
  mcols(mcols(object)) <- if (is.null(mcols(mcols(object)))) {
    mcolsRows
  } else if (all(names(mcols(mcols(object))) == c("type","description"))) {
    mcolsRows
  } else {
    cbind(mcols(mcols(object)), mcolsRows)
  }

  # stash the package version
  metadata(object)[["version"]] <- tryCatch(packageVersion("SCDSS"), error = function(e) NA)
  validObject(object)
  return(object)
}




#' @rdname SCASDataSet
#' @export

SCASDataSetFromSJFiles <- function( path,
                                    pattern = "SJ.out.tab",
                                    SampleInfo,
                                    design,
                                    MinReadsPerSJ = 1L,
                                    MinSjReadsInSpliSite = 1L,
                                    MinSJReadsRowSum = 100L,
                                    minSamps = 2L,
                                    uniqueMapOnly = TRUE,
                                    gtfFile,
                                    NumberThreads = 1,
                                    ... ) {
  if( length(list.files(path, pattern, full.names = TRUE) ) == 0){
    stop(paste0("There are no SJ.out.tab files in your directory."))
  }

  if( !file.exists(gtfFile) ) {
    stop(paste0("Your GTF file", gtfFile, "does not exist."))
  }

  if( !is.logical(uniqueMapOnly) ) {
    stop("Your uniqueMapOnly parameter is invalid, and must be a logical value.")
  }

  tmp <- MakeFromSJs(path = path, pattern = pattern, 
                     minSJ = MinReadsPerSJ, SSMinSJ = MinSjReadsInSpliSite, 
                     minSJs = MinSJReadsRowSum, minSamps = minSamps, 
                     uniqueMapOnly = uniqueMapOnly, NT = NumberThreads)
  sj.info <- tmp$sj.info
  rowData <- GenomicRanges::GRanges(seqnames = as.character(sj.info[, seqname]),
                                    ranges = IRanges::IRanges(start = as.integer(sj.info[, start]), end = as.integer(sj.info[, end])),
                                    strand = as.character(sj.info[, strand]),
                                    motif = as.character(sj.info[, motif]),
                                    annotation = as.character(sj.info[, annotation]))

  count <- as.matrix(tmp$count[, -c(1:4)])
  colnames(count) <- NULL
  start.psi <- as.matrix(tmp$start.psi[, -c(1:4)])
  colnames(start.psi) <- NULL
  end.psi <- as.matrix(tmp$end.psi[, -c(1:4)])
  colnames(end.psi) <- NULL

  if(!all(is.element(colnames(as.matrix(tmp$start.psi[, -c(1:4)])), row.names(SampleInfo)))) {
    stop("your colnames of SampleInfo are not contain all samples or names are not same.")
  }
  colData <- S4Vectors::DataFrame(SampleInfo[colnames(as.matrix(tmp$start.psi[, -c(1:4)])),])

  # check that these agree in number
  stopifnot(ncol(count) == nrow(colData))

  se = SummarizedExperiment::SummarizedExperiment(assays = list(counts = count,
                                                                start.psi = start.psi,
                                                                end.psi = end.psi),
                                                  rowRanges = rowData,
                                                  colData = colData)

  if (is(design, "formula")) {
    designVars <- all.vars(design)
    if (!all(designVars %in% names(colData(se)))) {
      stop("all variables in design formula must be columns in colData")
    }
    for (v in designVars) {
      if (any(is.na(colData(se)[[v]])))
        stop(paste("variables in design formula cannot contain NA:", v))
    }

    designVarsClass <- sapply(designVars, function(v) class(colData(se)[[v]])[1])

    if (any(designVarsClass == "character")) {
      warning("some variables in design formula are characters, converting to factors")
      for (v in designVars[designVarsClass == "character"]) {
        colData(se)[[v]] <- factor(colData(se)[[v]])
      }
    }

    if (length(designVars) == 1) {
      var <- colData(se)[[designVars]]
      if (all(var == var[1])) {
        stop("design has a single variable, with all samples having the same value.
             use instead a design of '~ 1'.")
      }
      }

    designVarsNumeric <- sapply(designVars, function(v) is.numeric(colData(se)[[v]]))
    if (any(designVarsNumeric)) {
      warnIntVars <- FALSE
      for (v in designVars[designVarsNumeric]) {
        if (all(colData(se)[[v]] == round(colData(se)[[v]]))) {
          warnIntVars <- TRUE
        }
      }
      if (warnIntVars) {
        message("  the design formula contains a numeric variable with integer values,
                specifying a model with increasing fold change for higher values.
                did you mean for this to be a factor? if so, first convert
                this variable to a factor using the factor() function")
      }
      }

    if (any(designVarsClass == "ordered")) {
      stop("the design formula contains an ordered factor. The internal steps
           do not work on ordered factors as a formula.")
    }

    designFactors <- designVars[designVarsClass == "factor"]
    missingLevels <- sapply(designFactors, function(v) any(table(colData(se)[[v]]) == 0))
    if (any(missingLevels)) {
      message("factor levels were dropped which had no samples")
      for (v in designFactors[missingLevels]) {
        colData(se)[[v]] <- droplevels(colData(se)[[v]])
      }
    }

    singleLevel <- sapply(designFactors, function(v) all(colData(se)[[v]] == colData(se)[[v]][1]))
    if (any(singleLevel)) {
      stop("design contains one or more variables with all samples having the same value,
           remove these variables from the design")
    }
    } else {
      stop("'design' should be a formula")
    }

  # Add columns on the columns
  mcolsCols <- DataFrame(type =  rep("input", ncol(colData(se))),
                         description = rep("", ncol(colData(se))))
  mcols(colData(se)) <- if (is.null(mcols(colData(se)))) {
    mcolsCols
  } else if (all(names(mcols(colData(se))) == c("type","description"))) {
    mcolsCols
  } else {
    cbind(mcols(colData(se)), mcolsCols)
  }

  sjInfo <- data.table::as.data.table(rowRanges(se))
  same_start_sj <- sjInfo[, .SD[.N > 1, ], by = list(seqnames, start, strand)]
  same_start_sj.id <- same_start_sj[, .SD[, paste(seqnames, ":", start, "-", end, ":", strand, collapse = "@", sep = "")], by = c("seqnames", "start", "strand")][, V1]
  same_start_sj$ID <- rep(same_start_sj.id, same_start_sj[, .N, by = c("seqnames", "start", "strand")][, N])
  SameStartSJ <- S4Vectors::DataFrame(same_start_sj[, .(seqnames, start, end, strand, width, motif, annotation, ID)])

  same_end_sj <- sjInfo[, .SD[.N > 1, ], by = list(seqnames, end, strand)]
  same_end_sj.id <- same_end_sj[, .SD[, paste(seqnames, ":", start, "-", end, ":", strand, collapse = "@", sep = "")], by = c("seqnames", "end", "strand")][, V1]
  same_end_sj$ID <- rep(same_end_sj.id, same_end_sj[, .N, by = c("seqnames", "end", "strand")][, N])
  SameEndSJ <- S4Vectors::DataFrame(same_end_sj[, .(seqnames, start, end, strand, width, motif, annotation, ID)])
  
  astypeandhostgene <- AStypeAndHostGene(sjInfo = sjInfo, gtfFile = gtfFile, NT = NumberThreads)
  HostGene <- astypeandhostgene$HostGene
  ASType <- astypeandhostgene$AS_Type

  object <- .SCASDataSet(se, design = design, HostGene = HostGene, ASType = ASType, SameStartSJ = SameStartSJ, SameEndSJ = SameEndSJ, ...)

  # now we know we have at least an empty GRanges or GRangesList for rowRanges
  # so we can create a metadata column 'type' for the mcols
  # and we label any incoming columns as 'input'

  # this is metadata columns on the rows
  mcolsRows <- DataFrame(type = rep("input", ncol(mcols(object))),
                         description = rep("", ncol(mcols(object))))
  mcols(mcols(object)) <- if (is.null(mcols(mcols(object)))) {
    mcolsRows
  } else if (all(names(mcols(mcols(object))) == c("type","description"))) {
    mcolsRows
  } else {
    cbind(mcols(mcols(object)), mcolsRows)
  }

  # stash the package version
  metadata(object)[["version"]] <- tryCatch(packageVersion("SCDSS"), error = function(e) NA)
  validObject(object)
  return(object)
}


#' @rdname SCASDataSet
#' @export
SCASDataSetFromMatrix <- function( counts,
                                   start.psi,
                                   end.psi,
                                   sjInfo,
                                   SampleInfo,
                                   design,
                                   gtfFile,
                                   NT = 1,
                                   ... ) {

  # check that these agree in number
  if(!identical(row.names(counts), row.names(start.psi)) || !identical(row.names(counts), row.names(end.psi))){
    stop("row.names of your counts, start.psi and end.psi tables must be consistent.")
  }

  if(!identical(colnames(counts), colnames(start.psi)) || !identical(colnames(counts), colnames(end.psi))){
    stop("colnames of your counts, start.psi and end.psi tables must be consistent.")
  }

  if (any(is.na(counts)))
    stop("NA values are not allowed in the count matrix")
  if (any(counts < 0)) {
    stop("some values in assay are negative")
  }
  if (!is.integer(counts)) {
    if (!is.numeric(counts)) {
      stop(paste("counts matrix should be numeric, currently it has mode:", mode(counts)))
    }
    if (any(round(counts) != counts)) {
      stop("some values in assay are not integers")
    }
    message("converting counts to integer mode")
    mode(counts) <- "integer"
  }

  if (all(counts == 0)) {
    stop("all samples have 0 counts for all SJs. check the counting script.")
  }

  if (all(rowSums(counts == counts[,1]) == ncol(counts))) {
    warning("all SJs have equal values for all samples. will not be able to perform differential analysis")
  }

  if (all(is.na(start.psi)) || all(is.na(end.psi))) {
    stop("all samples' PSI is NA for all SJs. check the counting script.")
  }

  if (!any(is.na(start.psi)) || !any(is.na(end.psi))) {
    warning("There are no missing value (NA) in your PSI table, that seems to be unreasonable.
            check the counting script.")
  }

  if(any(start.psi > 1, na.rm = T) || any(end.psi > 1, na.rm = T)){
    stop("There are some PSI more than 1, that seems to be unreasonable.")
  }

  if(any(start.psi < 0, na.rm = T) || any(end.psi < 0, na.rm = T)){
    stop("There are some PSI are negative, that seems to be unreasonable.")
  }

  if( !is(sjInfo, "GRanges") ){
    stop("sjInfo must be a GRanges format data related with assayList.")
  } else if ( !identical(as.character(sjInfo), row.names(counts)) ) {
    stop("The content and order of sjInfo must identical with counts, start.psi and end.psi.")
  }

  if(!all(is.element(colnames(start.psi), row.names(SampleInfo)))) {
    stop("your colnames of SampleInfo are not contain all samples or names are not same.")
  }
  colData <- S4Vectors::DataFrame(SampleInfo[colnames(start.psi),])

  se = SummarizedExperiment::SummarizedExperiment(assays = list(counts = counts,
                                                                start.psi = start.psi,
                                                                end.psi = end.psi),
                                                  rowRanges = sjInfo,
                                                  colData = colData)

  if (is(design, "formula")) {
    designVars <- all.vars(design)
    if (!all(designVars %in% names(colData(se)))) {
      stop("all variables in design formula must be columns in colData")
    }
    for (v in designVars) {
      if (any(is.na(colData(se)[[v]])))
        stop(paste("variables in design formula cannot contain NA:", v))
    }

    designVarsClass <- sapply(designVars, function(v) class(colData(se)[[v]])[1])

    if (any(designVarsClass == "character")) {
      warning("some variables in design formula are characters, converting to factors")
      for (v in designVars[designVarsClass == "character"]) {
        colData(se)[[v]] <- factor(colData(se)[[v]])
      }
    }

    if (length(designVars) == 1) {
      var <- colData(se)[[designVars]]
      if (all(var == var[1])) {
        stop("design has a single variable, with all samples having the same value.
             use instead a design of '~ 1'.")
      }
      }

    designVarsNumeric <- sapply(designVars, function(v) is.numeric(colData(se)[[v]]))
    if (any(designVarsNumeric)) {
      warnIntVars <- FALSE
      for (v in designVars[designVarsNumeric]) {
        if (all(colData(se)[[v]] == round(colData(se)[[v]]))) {
          warnIntVars <- TRUE
        }
      }
      if (warnIntVars) {
        message("  the design formula contains a numeric variable with integer values,
                specifying a model with increasing fold change for higher values.
                did you mean for this to be a factor? if so, first convert
                this variable to a factor using the factor() function")
      }
      }

    if (any(designVarsClass == "ordered")) {
      stop("the design formula contains an ordered factor. The internal steps
           do not work on ordered factors as a formula.")
    }

    designFactors <- designVars[designVarsClass == "factor"]
    missingLevels <- sapply(designFactors, function(v) any(table(colData(se)[[v]]) == 0))
    if (any(missingLevels)) {
      message("factor levels were dropped which had no samples")
      for (v in designFactors[missingLevels]) {
        colData(se)[[v]] <- droplevels(colData(se)[[v]])
      }
    }

    singleLevel <- sapply(designFactors, function(v) all(colData(se)[[v]] == colData(se)[[v]][1]))
    if (any(singleLevel)) {
      stop("design contains one or more variables with all samples having the same value,
           remove these variables from the design")
    }
    } else {
      stop("'design' should be a formula")
    }

  # Add columns on the columns
  mcolsCols <- DataFrame(type =  rep("input", ncol(colData(se))),
                         description = rep("", ncol(colData(se))))
  mcols(colData(se)) <- if (is.null(mcols(colData(se)))) {
    mcolsCols
  } else if (all(names(mcols(colData(se))) == c("type","description"))) {
    mcolsCols
  } else {
    cbind(mcols(colData(se)), mcolsCols)
  }

  sjInfo <- data.table::as.data.table(rowRanges(se))
  same_start_sj <- sjInfo[, .SD[.N > 1, ], by = list(seqnames, start, strand)]
  same_start_sj.id <- same_start_sj[, .SD[, paste(seqnames, ":", start, "-", end, ":", strand, collapse = "@", sep = "")], by = c("seqnames", "start", "strand")][, V1]
  same_start_sj$ID <- rep(same_start_sj.id, same_start_sj[, .N, by = c("seqnames", "start", "strand")][, N])
  SameStartSJ <- S4Vectors::DataFrame(same_start_sj[, .(seqnames, start, end, strand, width, motif, annotation, ID)])

  same_end_sj <- sjInfo[, .SD[.N > 1, ], by = list(seqnames, end, strand)]
  same_end_sj.id <- same_end_sj[, .SD[, paste(seqnames, ":", start, "-", end, ":", strand, collapse = "@", sep = "")], by = c("seqnames", "end", "strand")][, V1]
  same_end_sj$ID <- rep(same_end_sj.id, same_end_sj[, .N, by = c("seqnames", "end", "strand")][, N])
  SameEndSJ <- S4Vectors::DataFrame(same_end_sj[, .(seqnames, start, end, strand, width, motif, annotation, ID)])

  astypeandhostgene <- AStypeAndHostGene(sjInfo = sjInfo, gtfFile = gtfFile, NT = NT)
  HostGene <- astypeandhostgene$HostGene
  ASType <- astypeandhostgene$AS_Type

  object <- .SCASDataSet(se, design = design, HostGene = HostGene, ASType = ASType, SameStartSJ = SameStartSJ, SameEndSJ = SameEndSJ, ...)

  # now we know we have at least an empty GRanges or GRangesList for rowRanges
  # so we can create a metadata column 'type' for the mcols
  # and we label any incoming columns as 'input'

  # this is metadata columns on the rows
  mcolsRows <- DataFrame(type = rep("input", ncol(mcols(object))),
                         description = rep("", ncol(mcols(object))))
  mcols(mcols(object)) <- if (is.null(mcols(mcols(object)))) {
    mcolsRows
  } else if (all(names(mcols(mcols(object))) == c("type","description"))) {
    mcolsRows
  } else {
    cbind(mcols(mcols(object)), mcolsRows)
  }

  # stash the package version
  metadata(object)[["version"]] <- tryCatch(packageVersion("SCDSS"), error = function(e) NA)
  validObject(object)
  return(object)
}


