#RColorBrewer
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
