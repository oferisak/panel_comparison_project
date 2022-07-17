
generate_panel_from_panels<-function(joined_panel_names,gen_source_db,gen_dm_df,max_distance=1,with_mitochondrial=F,
                                     original_panel_max_score=T){
  panels_per_gene<-gen_source_db%>%group_by(gene_symbol)%>%summarize(number_of_panels_with_gene=n())%>%ungroup()
  number_of_genes_per_panel<-gen_source_db%>%group_by(panel_joined_name)%>%summarize(number_of_genes=n())%>%ungroup()
  # first tier similar panels
  #similar_panels_df<-get_similar_panels_from_dm_df_by_prop(gen_dm_df,joined_panel_names,top_most_similar = 0.1,maximal_distance = 0.99)

  similar_panels_df<-gen_dm_df%>%
    filter(row %in% joined_panel_names,value<=max_distance)%>%
    group_by(col)%>%
    summarize(dist=min(value))#take the panel with the smallest distance to one of the given panels
  similar_panels_df <- similar_panels_df %>% 
    left_join(gtr_panels_list %>%
                select(panel_joined_name),
              by = c('col' = 'panel_joined_name')) %>%
    left_join(number_of_genes_per_panel, by = c('col' = 'panel_joined_name'))
  
  message('collect all similar panels and add their genes..')
  gene_rank<-gen_source_db%>%
    mutate(panel_joined_name=as.character(panel_joined_name))%>%
    filter(panel_joined_name%in%similar_panels_df$col)%>%
    select(panel_joined_name,gene_symbol)%>%
    left_join(similar_panels_df%>%
                select(col,dist,number_of_genes)%>%
                rename('panel_joined_name'=col))%>%
    mutate(gene_score=1-dist^2,
           panel_dist=dist)%>%
    group_by(gene_symbol)%>%
    summarize(num_o_panels =n(),
              #closest_dist=min(gene_score),
              #closest_panel=panel_joined_name[which(gene_score==closest_dist)][1],
              gene_score=sum(gene_score),
              min_panel_dist=min(panel_dist),
              max_panel_dist=max(panel_dist),
              median_panel_dist=median(panel_dist),
              mean_panel_dist=mean(panel_dist))%>%
    mutate(rank=rank(gene_score,tie='max'))
  
  genes_in_original_panels<-gen_source_db%>%filter(panel_joined_name%in%joined_panel_names)%>%pull(gene_symbol)%>%unique()
  gene_rank$in_original_panels<-ifelse(gene_rank$gene_symbol %in% genes_in_original_panels,1,0)
  
  # if the gene is in the original panel and the original_panel_max_score flag is T change its score to be the max score
  gene_rank<-gene_rank%>%
    mutate(gene_score=ifelse(original_panel_max_score & in_original_panels,max(gene_score),gene_score))
  
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

# A function that takes in a phenotype name, the results for the autogenerator for that phenotype and 
# the results for the pubmed search for all the genes against that phenotype 
# and calculates the pr curve for the scaled gene score 
calculate_auto_generator_pr_for_phenotype<-function(phenotype_name,
                                                    auto_generator_panel_res,
                                                    all_genes_vs_phenotype,
                                                    gen_source_db,
                                                    score_column='scaled_gene_score'){
  panels_per_gene<-gen_source_db%>%group_by(gene_symbol)%>%summarize(number_of_panels_with_gene=n())%>%ungroup()
  table_for_roc_panelapp_only<-all_genes_vs_phenotype%>%
    left_join(panels_per_gene)%>%
    select(gene_symbol,has_pub,number_of_panels_with_gene)%>%#add the number of panels the gene is in for ALL the genes
    left_join(auto_generator_panel_res%>%select(-number_of_panels_with_gene))%>%
    mutate(scaled_gene_score=ifelse(is.na(scaled_gene_score),0,scaled_gene_score),
           scaled_adjusted_gene_score=ifelse(is.na(scaled_adjusted_gene_score),0,scaled_adjusted_gene_score),
           rank_percentile=ifelse(is.na(rank_percentile),0,rank_percentile),
           num_o_panels=ifelse(is.na(num_o_panels),0,num_o_panels),
           min_panel_dist=ifelse(is.na(min_panel_dist),1,min_panel_dist),
           max_panel_dist=ifelse(is.na(max_panel_dist),1,max_panel_dist),
           median_panel_dist=ifelse(is.na(median_panel_dist),1,median_panel_dist),
           mean_panel_dist=ifelse(is.na(mean_panel_dist),1,mean_panel_dist),
           )

    # add a modified gene score - if you want to test other options
    # mutate(gene_score_modified=ifelse(!is.na(in_original_panels)&(in_original_panels),scaled_adjusted_gene_score+5,scaled_adjusted_gene_score))
  
  if (!(score_column%in%colnames(table_for_roc_panelapp_only))){stop(glue('{score_column} is not a valid column produced by the auto panel generator..'))}
  
  table_for_roc_panelapp_only['score_for_roc']<-table_for_roc_panelapp_only[score_column]
  
  phenotype_auto_generator_pr_curve <-
    data.frame(
      phenotype_name = selected_phenotype_name,
      group = 'gene_score',
      table_for_roc_panelapp_only %>%
        pr_curve(
          truth = has_pub,
          estimate = score_for_roc,
          event_level = 'second'
        )
    )%>%
    # if you want to add a modified score as well - if so MIND THAT YOU NEED TO CHANGE THE GROUP NAME IN THE MAIN SCRIPT
    # bind_rows(data.frame(phenotype_name=selected_phenotype_name,
    #                      group='gene_score_modified',
    #                      table_for_roc_panelapp_only%>%
    #                        pr_curve(truth=has_pub,estimate=gene_score_modified,event_level = 'second')))%>%
    mutate(f1=(2 * precision * recall) / (precision + recall),
           pr_auc=table_for_roc_panelapp_only%>%
             pr_auc(truth = has_pub,
                    estimate = score_for_roc,
                    event_level = 'second')%>%pull(.estimate))
  return(phenotype_auto_generator_pr_curve)
}
