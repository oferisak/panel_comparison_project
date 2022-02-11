
generate_panel_from_panels<-function(joined_panel_names,gen_source_db,gen_dm_df,max_distance=1,with_mitochondrial=F){

  # first tier similar panels
  #similar_panels_df<-get_similar_panels_from_dm_df_by_prop(gen_dm_df,joined_panel_names,top_most_similar = 0.1,maximal_distance = 0.99)
  
  similar_panels_df<-gen_dm_df%>%filter(row %in% joined_panel_names)%>%filter(value<max_distance)
  
  gene_rank_source<-NULL
  message('collect all similar panels and add their genes..')
  for (joined_panel_name in unique(similar_panels_df$col)){
    panel_dist<-similar_panels_df%>%filter(col==joined_panel_name)%>%pull(value)%>%min()# the slice part is in case the panels are found to be related
    #gene_score_by_dist<-1-panel_dist
    gene_score_by_dist<-1-panel_dist
    panel_source<-gen_source_db%>%filter(panel_joined_name==joined_panel_name)%>%select(panel_joined_name,gene_symbol)%>%
      mutate(gene_score=gene_score_by_dist)
    if (!with_mitochondrial){panel_source<-panel_source%>%filter(!grepl('MT-',gene_symbol))}
    gene_rank_source<-gene_rank_source%>%bind_rows(panel_source)
  }
  
  gene_rank<-gene_rank_source%>%group_by(gene_symbol)%>%summarize(num_o_panels=n(),
                                                                  gene_score=sum(gene_score))%>%
    mutate(rank=rank(gene_score,tie='max'))
  gene_rank$rank_percentile=gene_rank$rank/nrow(gene_rank)
  return(gene_rank)
  
}
