#' DSUplot.ridge
#'
#' DSUplot.ridge of SJs candidates 
#' @rdname DSUPlot.ridge
#' @name DSUPlot.ridge
#' @import ggplot2
#' @import ggridges
#' @import SummarizedExperiment
#' @importFrom cowplot plot_grid
#'
#' @param object an \code{SCAS-class} object.
#' @param sj name of splice junction, eg: "1:26460145-26468894:+"
#' @param direct the direction of SJs, same 'start' or 'end'
#' @param groupBy the group of plot, same as group-specific analysis.
#' @param point_size point size of rug
#' @param jitter_width,jitter_height jitter width and height of rug point
#' @param xlim xlim of plot
#' @example 
#' ridgeplot <- function(object,sj,direct = NULL,groupBy = NULL,
#' point_size = 2,jitter_width = 0,jitter_height = 0.1,xlim = c(0, 1))
#' @seealso  \code{\link[ggridges]{ggridges}}
#' @export
DSUPlot.ridge <- function(object,
                       sj,
                       direct = NULL,
                       groupBy = NULL,
                       point_size = 2,
                       jitter_width = 0,
                       jitter_height = 0.1,
                       xlim = c(0, 1)) {
  tab <- selectSJ(object, sj)

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
  }

  tab <- merge(tab, data.table::as.data.table(groups, keep.rownames = "ID"), by = "ID")

  if( !is.null(direct) ) {
    if ( length(direct) != 1 ) {
      stop("direct must be one of 'start' or 'end'")
    } else if ( !is.element( direct, c("start", "end") ) ) {
      stop("direct must be one of 'start' or 'end'")
    }
  }

  if( is.null(direct) ) {
    p1 = ggplot(tab, aes(x = start.psi, y = groups, fill = groups)) +
      geom_density_ridges(aes(point_color = groups),
                          scale = .5,
                          color = "white",
                          jittered_points = TRUE,
                          point_shape = 19,
                          alpha = 0.7,
                          point_size = point_size,
                          position = position_raincloud(width = jitter_width,
                                                        height = jitter_height,
                                                        ygap = 0.05,
                                                        adjust_vlines = FALSE,
                                                        seed = NULL)) +
      theme_ridges(center_axis_labels = TRUE)+
      guides(fill = FALSE, point_color = FALSE)+
      labs(x = "PSI", y = paste0(groupBy, collapse = "_"), tag = "Same Start")+
      lims(x = xlim)

    p2 = ggplot(tab, aes(x = end.psi, y = groups, fill = groups)) +
      geom_density_ridges(aes(point_color = groups),
                          scale = .5,
                          color = "white",
                          jittered_points = TRUE,
                          point_shape = 19,
                          alpha = 0.5,
                          point_size = point_size,
                          position = position_raincloud(width = jitter_width,
                                                        height = jitter_height,
                                                        ygap = 0.05,
                                                        adjust_vlines = FALSE,
                                                        seed = NULL)) +
      theme_ridges(center_axis_labels = TRUE)+
      guides(fill = FALSE, point_color = FALSE)+
      labs(x = "PSI", y = paste0(groupBy, collapse = "_"), tag = "Same End")+
      lims(x = xlim)
    print(cowplot::plot_grid(p1, p2))
  } else if( direct == "start" ) {
    p = ggplot(tab, aes(x = start.psi, y = groups, fill = groups)) +
      geom_density_ridges(aes(point_color = groups),
                          scale = .5,
                          color = "white",
                          jittered_points = TRUE,
                          point_shape = 19,
                          alpha = 0.5,
                          point_size = point_size,
                          position = position_raincloud(width = jitter_width,
                                                        height = jitter_height,
                                                        ygap = 0.05,
                                                        adjust_vlines = FALSE,
                                                        seed = NULL)) +
      theme_ridges(center_axis_labels = TRUE)+
      guides(fill = FALSE, point_color = FALSE)+
      labs(x = "PSI", y = paste0(groupBy, collapse = "_"))+
      lims(x = xlim)
    print(p)
  } else if( direct == "end" ) {
    p = ggplot(tab, aes(x = end.psi, y = groups, fill = groups)) +
      geom_density_ridges(aes(point_color = groups),
                          scale = .5,
                          color = "white",
                          jittered_points = TRUE,
                          point_shape = 19,
                          alpha = 0.5,
                          point_size = point_size,
                          position = position_raincloud(width = jitter_width,
                                                        height = jitter_height,
                                                        ygap = 0.05,
                                                        adjust_vlines = FALSE,
                                                        seed = NULL)) +
      theme_ridges(center_axis_labels = TRUE)+
      guides(fill = FALSE, point_color = FALSE)+
      labs(x = "PSI", y = paste0(groupBy, collapse = "_"))+
      lims(x = xlim)
    print(p)
  }
}
