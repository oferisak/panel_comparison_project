plot_number_of_genes_per_panel_dist<-function(){
  g<-panels_list%>%ggplot(aes(x=number_of_genes,fill=source))+
    geom_density(alpha=0.3)+
    theme_minimal()+scale_fill_nejm()+
    scale_x_log10()+
    theme(legend.position = 'bottom')+
    labs(title = 'Number of genes per panel',fill='',x='Number of genes')
  g
}

plot_panel_gene_heatmap<-function(binary_panel_gene_table,sort_by_similarity_to_selected_panel=T){
  # first remove all the genes that are not found in any of the panels
  bin_for_heatmap<-binary_panel_gene_table%>%select_if(~ !is.numeric(.) || sum(.) != 0)
  # heat_plot<-heatmap(as.matrix(bin_for_heatmap),scale = 'none',
  #         Colv = T,Rowv = T)
  library(pheatmap)
  heat_plot <- pheatmap::pheatmap(
    as.matrix(bin_for_heatmap),
    border_color = 'darkgray',
    cluster_cols = T,
    cluster_rows = !sort_by_similarity_to_selected_panel,
    #cellheight = 20,cellwidth = 10,
    treeheight_row = 0,treeheight_col = 0,
    legend = F,
    color = c('lightgray', 'darkred'),
    scale = 'none'
  )
  heat_plot
}

plot_dm_as_dendogram<-function(dm,selected_panel=NA){
  hc<-hclust(dm)
  #plot(hc,labels = NULL, hang = 0)
  hcd <- as.dendrogram(hc)
  if (!is.na(selected_panel)){
    label_col<-c(rep('black',nrow(selected_bin_data)))
    label_col[which(selected_panel==hcd%>%labels())]<-'darkred'
    hcd<-hcd%>%dendextend::set("labels_col", label_col)
  }
  # Default plot
  plot(hcd, type = "rectangle", ylab = "Height",horiz=T,xlim = c(1,-1))
  # to generate a radial plot
  # colors = pal_npg(palette = c("nrc"), alpha = 1)(10)
  # clus4 = cutree(hc, 10)
  # plot(as.phylo(hc), type = "fan", tip.color = colors[clus4],
  #      label.offset = 0.1, cex = 0.7)
}
