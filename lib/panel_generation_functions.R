
generate_panel_from_panels<-function(joined_panel_names,gen_source_db,gen_dm_df,max_distance=1,with_mitochondrial=F){
  panels_per_gene<-gen_source_db%>%group_by(gene_symbol)%>%summarize(number_of_panels_with_gene=n())
  
  # first tier similar panels
  #similar_panels_df<-get_similar_panels_from_dm_df_by_prop(gen_dm_df,joined_panel_names,top_most_similar = 0.1,maximal_distance = 0.99)


  similar_panels_df<-gen_dm_df%>%filter(row %in% joined_panel_names)%>%filter(value<=max_distance)%>%
    group_by(col)%>%summarize(dist=min(value))#take the panel with the smallest distance to one of the given panels
  similar_panels_df<-similar_panels_df%>%left_join(
    gtr_panels_list%>%
      select(panel_joined_name,condition_identifiers),
    by=c('col'='panel_joined_name'))
  
  message('collect all similar panels and add their genes..')
  gene_rank<-gen_source_db%>%mutate(panel_joined_name=as.character(panel_joined_name))%>%
    filter(panel_joined_name%in%similar_panels_df$col)%>%select(panel_joined_name,gene_symbol)%>%
    left_join(similar_panels_df%>%select(col,dist,condition_identifiers)%>%rename('panel_joined_name'=col))%>%
    mutate(gene_score=1-dist)%>%
    group_by(gene_symbol)%>%
    summarize(num_o_panels =n(), 
              #closest_dist=min(gene_score),
              #closest_panel=panel_joined_name[which(gene_score==closest_dist)][1],
              gene_score=sum(gene_score))%>%
    mutate(rank=rank(gene_score,tie='max'))
  
  genes_in_original_panels<-gen_source_db%>%filter(panel_joined_name%in%joined_panel_names)%>%pull(gene_symbol)%>%unique()
  gene_rank$in_original_panels<-ifelse(gene_rank$gene_symbol %in% genes_in_original_panels,1,0)
  
  # now adjust for number of panels per gene
  gene_rank<-gene_rank%>%left_join(panels_per_gene)%>%
    mutate(adjusted_gene_score=gene_score*num_o_panels/number_of_panels_with_gene,
           adjusted_rank=rank(adjusted_gene_score,tie='max'))
    
  
  gene_rank<-gene_rank%>%mutate(rank_percentile=rank/nrow(gene_rank),
                                adjusted_rank_percentile=adjusted_rank/nrow(gene_rank),
                                scaled_gene_score=as.numeric(scale(gene_score,center=F)),# disabled center, so the score is limited by 0 at the bottom
                                scaled_adjusted_gene_score=as.numeric(scale(adjusted_gene_score,center=F)))
  
  # gene_rank$rank_percentile=gene_rank$rank/nrow(gene_rank)
  # gene_rank$adjusted_rank_percentile=gene_rank$adjusted_rank/nrow(gene_rank)
  # 
  return(gene_rank)
  
}
