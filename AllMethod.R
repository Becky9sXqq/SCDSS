#' design for SCASDataSet
#'
#' @param object A SCASDataSet object
#' @name design
#' @export
setMethod(f = "design",
          signature = "SCASDataSet",
          definition = function(object) {
            object@design
})

#' show method for SCASDataSet
#'
#' @param object A SCASDataSet object
#' @name show
#' @aliases show
#' @docType methods
#' @rdname show-methods
#'
setMethod(
  f = "show",
  signature = "SCASDataSet",
  definition = function(object) {
    cat(
      "An object of class",
      class(object),
      "in project",
      object@project.name,
      "\n",
      S4Vectors::nrow(object),
      "SJs across",
      S4Vectors::ncol(object),
      "cells.\n"
    )
    invisible(x = NULL)
  }
)

#' @name SSMD
#' @rdname groupSpecific
#' @export
setMethod("SSMD", "SCASDataSet",
                 function(object) {
                   object@SSMD
                 })

#' @name SSMD
#' @rdname groupSpecific
#' @exportMethod "SSMD<-"
setReplaceMethod("SSMD", "SCASDataSet", function(x, value) {
  x@SSMD <- value
  validObject(x)
  x
})

#' @name Hypergeometric
#' @rdname groupSpecific
#' @export
setMethod("Hypergeometric", "SCASDataSet",
          function(object) {
            object@Hypergeometric
          })

#' @name Hypergeometric
#' @rdname groupSpecific
#' @exportMethod "Hypergeometric<-"
setReplaceMethod("Hypergeometric", "SCASDataSet", function(x, value) {
  x@Hypergeometric <- value
  validObject(x)
  x
})

#' @name BetaBinomial
#' @rdname groupSpecific
#' @export
setMethod("BetaBinomial", "SCASDataSet",
          function(object) {
            object@BetaBinomial
          })

#' @name BetaBinomial
#' @rdname groupSpecific
#' @exportMethod "BetaBinomial<-"
setReplaceMethod("BetaBinomial", "SCASDataSet", function(x, value) {
  x@BetaBinomial <- value
  validObject(x)
  x
})

#' @name RandomForest
#' @rdname groupSpecific
#' @export
setMethod("RandomForest", "SCASDataSet",
          function(object) {
            object@RandomForest
          })

#' @name RandomForest
#' @rdname groupSpecific
#' @exportMethod "RandomForest<-"
setReplaceMethod("RandomForest", "SCASDataSet", function(x, value) {
  x@RandomForest <- value
  validObject(x)
  x
})

#' select candidate
#'
#' select a given SJs from a SCASDataSet object
#'
#' @param object A SCASDataSet object
#' @param sj name of splice junction, eg: 1:26460145-26468894:+
#' @name select
#' @export
#'
setMethod("selectSJ",
          signature = "SCASDataSet",
          definition = function(object, sj) {
            if( !is(object, "SCASDataSet") ) {
              stop("object must be a SCASDataSet class.")
            }
            if( !is.element(sj, as.character(rowRanges(object))) ) {
              stop( paste0("invalid sj", sj, "and your sj must be in rowRanges of your object") )
            }
            irow <- which(as.character(rowRanges(object)) == sj)
            tab <- cbind(assay(object, i = 1)[irow, ], assay(object, 2)[irow, ], assay(object, 3)[irow, ])
            colnames(tab) <- assayNames(object)
            tab <- merge(tab, colData(object), by = 0)
            data.table::setDT(tab)
            data.table::setnames(tab, "Row.names", "ID")
            return(tab)
          })


#' Obtain BetaBinomialDSU result from SCASDataSet
#' @name BetaBinomialDSU
#' @title BetaBinomialDSU
#' @rdname groupSpecific
#' @export
setMethod("BetaBinomialDSU", "SCASDataSet",
          function(object) {
            object@BetaBinomialDSU
          })

#' Obtain BetaBinomialDSU result from SCASDataSet
#' @name BetaBinomialDSU
#' @title BetaBinomialDSU
#' @rdname groupSpecific
#' @exportMethod "BetaBinomialDSU<-"
setReplaceMethod("BetaBinomialDSU", "SCASDataSet", function(x, value) {
  x@BetaBinomialDSU <- value
  validObject(x)
  x
})

#' Obtain MGLMDSU result from SCASDataSet
#' @name MGLMDSU
#' @title MGLMDSU
#' @rdname groupSpecific
#' @export
setMethod("MGLMDSU", "SCASDataSet",
          function(object) {
            object@MGLMDSU
          })

#' Obtain MGLMDSU result from SCASDataSet
#' @name MGLMDSU
#' @title MGLMDSU
#' @rdname groupSpecific
#' @exportMethod "MGLMDSU<-"
setReplaceMethod("MGLMDSU", "SCASDataSet", function(x, value) {
  x@MGLMDSU <- value
  validObject(x)
  x
})



#' Obtain RawPsiTable result from SCASDataSet
#' @name RawPsiTable
#' @title RawPsiTable
#' @rdname groupSpecific
#' @export
setMethod("RawPsiTable", "SCASDataSet",
          function(object) {
            object@RawPsiTable
          })

#' Obtain RawPsiTable result from SCASDataSet
#' @name RawPsiTable
#' @title RawPsiTable
#' @rdname groupSpecific
#' @exportMethod "RawPsiTable<-"
setReplaceMethod("RawPsiTable", "SCASDataSet", function(x, value) {
  x@RawPsiTable <- value
  validObject(x)
  x
})

