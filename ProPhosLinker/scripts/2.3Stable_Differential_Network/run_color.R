#' Generate a qualitative color palette for annotation classes using RColorBrewer.
#'
#' This function automatically assigns distinct colors to unique categories (classes) 
#' in the \code{Class} column of an annotation table. It draws from all available 
#' qualitative palettes in \code{RColorBrewer}, and if more colors are needed than 
#' the total number of qualitative colors (~40), it interpolates a smooth expanded 
#' palette using \code{grDevices::colorRampPalette}.
#'
#' @param annotation_table A data frame containing a column named \code{Class} with 
#'        categorical annotations (e.g., "Kinase", "Transcription Factor", etc.).
#'
#' @return A named character vector where names are unique class labels and values 
#'         are hex color codes (e.g., \code{"Kinase" = "#E69F00"}).
#'
#' @importFrom RColorBrewer brewer.pal brewer.pal.info
#' @importFrom grDevices colorRampPalette
#'
#' @examples
#' # Example annotation table
#' anno <- data.frame(Class = c("A", "B", "C", "A"))
#' colors <- run_color(anno)
#' print(colors)
#'
#' @export

run_color<-function(annotation_table=NULL){
  unique_classes<-unique(annotation_table$Class)
  qual_col_pals <- RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
  cluster_Palette <- unique(unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))
  if(length(unique_classes)>length(cluster_Palette)){
    cluster_Palette <-grDevices::colorRampPalette(cluster_Palette)(length(unique_classes))#避免颜色不够
  }
  color_mapping <- setNames(cluster_Palette[1:length(unique_classes)], unique_classes)
  return(color_mapping)
}
