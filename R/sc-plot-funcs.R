#' getlegend
#'
#' Using the ggpubr package, get the legend of a ggplot2 based plot.
#' Useful for Seurat::AugmentPlot plots because the legends are lost.
#'
#' @param plot a ggplot2 plot object
#'
#' @return Returns a data-frame containing the the 10X cell ID and the lintraceR DNA barcode ID.
#'
#' @import ggpubr
#' @export
#'

getlegend <- function(plot){
  # Extract the legend. Returns a gtable
  legend <- ggpubr::get_legend(plot)
  # Convert to a ggplot and print
  ggpubr::as_ggplot(legend)
}


#' plotAugmentDim
#'
#' Wrapper around Seurat::AugmentPlot and Seurat::Dimplot
#' Use in conjunction with pdf() statements
#'
#' @param obj a Seurat object
#' @param reduction a dimensional reduction within obj
#' @param label Logical. Label active identity classes
#' @param label.size pt size of label text
#' @param group.by identity class to group dataset by prior to plotting
#' @param title plot title
#' @param dpi resolution (dots per inch) of raster image
#' @param pt.size size of points on plot
#' @param cells.highlight vector or named list of cell labels to highlight on plot
#' @param cols.highlight vector of colors to highlight cells
#' @param sizes.highlight pt size of highlighted cells
#' @param split.by an ident class to split plots by
#'
#' @return Returns a rasterised ggplot2 object containing the plot (axes are vectors) and on a second page the plot legend in vector format.
#'
#' @import Seurat
#' @export
#'

plotAugmentDim <- function(obj, reduction = "umap",
                           label = F, label.size = 6,
                           group.by = NULL, title = NULL,
                           dpi = 600, pt.size = 1,
                           cells.highlight = NULL, cols.highlight = "red",
                           sizes.highlight = 2, split.by = NULL) {

  # generate plot object
  p <- Seurat::DimPlot(obj, reduction = reduction,
                       label = label, label.size = label.size,
                       group.by = group.by, pt.size = pt.size,
                       cells.highlight = cells.highlight, cols.highlight = cols.highlight,
                       sizes.highlight = sizes.highlight, split.by = split.by) +
    ggplot2::ggtitle(title)

  # return plot and legend on separate page
  print(Seurat::AugmentPlot(p, dpi = dpi))
  print(getlegend(p))
}


#' plotAugmentFeature
#'
#' Wrapper around Seurat::AugmentPlot and Seurat::Dimplot
#' Use in conjunction with pdf() statements
#'
#' @param obj a Seurat object
#' @param reduction a dimensional reduction within obj
#' @param features a vector of features to plot
#' @param label Logical. Label active identity classes
#' @param label.size pt size of label text
#' @param group.by identity class to group dataset by prior to plotting
#' @param title plot title
#' @param dpi resolution (dots per inch) of raster image
#' @param pt.size size of points on plot
#' @param order order by expression value, largest values are plotted last (on top)
#' @param ncol number of columns
#'
#' @return Returns a rasterised ggplot2 object containing the plot (axes are vectors) and on a second page the plot legend in vector format.
#'
#' @import Seurat
#' @export
#'
plotAugmentFeature <- function(obj, reduction = "umap", features = NULL,
                               label = F, label.size = 6, group.by = NULL,
                               title = NULL, dpi = 300, pt.size = 1,
                               order = T, ncol = NULL) {

  # generate plot object
  p <- Seurat::FeaturePlot(obj, reduction = reduction, features = features,
                           label = label, label.size = label.size,
                           pt.size = pt.size, order = order, ncol = ncol) +
    ggplot2::ggtitle(title)

  # return plot and legend on separate page
  print(AugmentPlot(p, dpi = dpi))
  print(getlegend(p))
}


#' plotElbow
#'
#' Wrapper around Seurat::ElbowPlot and PCAtools::findElbowPoint
#'
#' @param obj a Seurat object
#' @param outdir path to output directory
#' @param samplename name of Seurat object
#'
#' @return Returns a rasterised ggplot2 object containing the plot (axes are vectors) and on a second page the plot legend in vector format.
#'
#' @import Seurat
#' @export
#'
plotElbow <- function(obj, outdir = NULL, samplename = "sample"){
    percent.var <- PCAtools::Stdev(obj)
    elbow.dim <- PCAtools::findElbowPoint(percent.var)
    p <- Seurat::ElbowPlot(obj, ndims = 50) +
      ggplot2::geom_vline(ggplot2::aes(xintercept = elbow.dim), color = "red") +
      ggplot2::labs(title = paste("Elbow plot -", samplename), subtitle = paste("Predicted elbow PC =", elbow.dim))

    if(is.null(outdir)){
      print(p)
    } else {
      ggplot2::ggsave(file.path(outdir,paste(prefix,".pdf", sep = '')), width = 8, height = 6)
    }
  }
