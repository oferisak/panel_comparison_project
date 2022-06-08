calculate_panel_similarity_as_df<-function(panel_a,panel_b){
  sum_ab<-length(panel_a)+length(panel_b)
  common_genes<-length(intersect(panel_a,panel_b))
  similarity<-common_genes/sum_ab
  return(data.frame(similarity,
                    selected_panel_genes=length(panel_a),
                    comparison_panel_genes=length(panel_b),
                    common_genes,
                    genes_not_in_selected_panel=length(setdiff(panel_b,panel_a))))
}

calculate_panel_similarity<-function(panel_a,panel_b){
  sum_ab<-length(panel_a)+length(panel_b)
  common_genes<-length(intersect(panel_a,panel_b))
  similarity<-common_genes/sum_ab
  return(similarity)
}

# A function that takes in a vector of gene panel names (joined panel names) and returns the subset of the dm_df with
# all the pairs of panels
subset_dm_df_by_gene_panel<-function(gene_panels,gtr_dm_df){
  panels_dist<-gtr_dm_df%>%filter(row%in%gene_panels,col%in%gene_panels)%>%filter(row!=col)
  panels_dist<-panels_dist%>%rowwise()%>%mutate(joined_panels=paste0(sort(c(row,col)),collapse='_'))%>%ungroup()
  # remove duplicates
  panels_dist<-panels_dist%>%dplyr::filter(!duplicated(joined_panels))
  return(panels_dist)
}

# A function that takes in a panel (source + id) and calculates similarity with ALL other panels
get_similarity_df_table_for_panel<-function(selected_source,selected_panel_id,panels_db){
  message(glue('Analyzing {selected_source}:{selected_panel_id}..'))
  selected_panel_genes<-panels_db%>%filter(source==selected_source & panel_id==selected_panel_id)%>%pull(gene_symbol)
  message(glue('There are {length(selected_panel_genes)} genes in panel'))
  similar_panels_table<-panels_db%>%group_by(source,panel_name,panel_id)%>%
    filter(!(source==selected_source & panel_id==selected_panel_id))%>%
    summarize(similarity=calculate_panel_similarity_as_df(selected_panel_genes,gene_symbol))
  message(glue('There are {nrow(similar_panels_table)} panels in output'))
  return(similar_panels_table)
}

get_similarity_table_for_panel<-function(selected_source,selected_panel_id,panels_db){
  message(glue('calculating similarity table for {selected_source}: {selected_panel_id}'))
  selected_panel_genes<-panels_db%>%filter(source==selected_source & panel_id==selected_panel_id)%>%pull(gene_symbol)
  similar_panels_table<-panels_db%>%group_by(source,panel_name,panel_id)%>%
    filter(!(source==selected_source & panel_id==selected_panel_id))%>%
    summarize('{selected_source}|{selected_panel_id}':=calculate_panel_similarity(selected_panel_genes,gene_symbol))
  return(similar_panels_table)
}

get_similarity_matrix<-function(source_db){
  clingen_similarity<-panels_list%>%
    group_by(source,panel_id,panel_name)%>%
    summarize(get_similarity_table_for_panel(selected_source = source,selected_panel_id = panel_id,panels_db = source_db))
  similartiy_matrix<-clingen_similarity%>%dplyr::select(c('source','panel_id','panel_name'))%>%distinct()
  for (colname in setdiff(colnames(clingen_similarity),c('source','panel_id','panel_name'))){
    print(colname)
    similartiy_matrix<-similartiy_matrix%>%
      left_join(clingen_similarity%>%
                  dplyr::select('source','panel_id','panel_name',colname)%>%
                  drop_na())
  }
  return(similartiy_matrix)
}

create_binary_gene_table<-function(source_db){
  # remove panelapp low confidence genes from source db
  binary_panel_gene_table<-source_db
  binary_panel_gene_table <-
    binary_panel_gene_table %>%
    unite(col = panel, sep = '|', source, panel_id, panel_name) %>%
    dplyr::select(panel, gene_symbol) %>%
    mutate(val = 1) %>%
    distinct %>%
    spread(key = gene_symbol, value = val, fill = 0)
  binary_panel_gene_table<-data.frame(binary_panel_gene_table)
  rownames(binary_panel_gene_table)<-binary_panel_gene_table$panel
  binary_panel_gene_table$panel<-NULL
  return(binary_panel_gene_table)
}

subset_binary_gene_table<-function(binary_panel_gene_table,selected_panel_names){
  panels_num<-match(selected_panel_names,row.names(binary_panel_gene_table))
  return(binary_panel_gene_table[panels_num,])
}

generate_distance_matrix<-function(binary_panel_gene_table,method=1){
  # one option is to use the base distance metric with the binary method but im not sure its the most appropriate
  #dm<-dist(binary_panel_gene_table,upper = F,diag=F,method='binary')
  # the more appropriate option is probably to use the jaccard method (recommended across the web)
  dm<-ade4::dist.binary(binary_panel_gene_table,method=method)
  return(dm)
}

convert_dm_to_df<-function(dm){
  df <- reshape2::melt(as.matrix(dm), varnames = c("row", "col"))
  return(df)
}

find_panel_in_dm_df<-function(dm_df,source,panel_id=NA,panel_name=NA){
  search_term<-glue('{source}')
  search_term<-ifelse(!is.na(panel_id),glue('{search_term}\\|{panel_id}'),glue('{search_term}\\|.+'))
  search_term<-ifelse(!is.na(panel_name),glue('{search_term}\\|.*{panel_name}.*'),glue('{search_term}\\|.+'))
  selected_panel_name<-dm_df%>%filter(grepl(search_term,dm_df$row,ignore.case = T))%>%pull(row)%>%as.character()%>%unique()
  if (length(selected_panel_name)!=1){stop(glue('find_panel_in_dm_df: using this search term: {search_term} resulted in {length(selected_panel_name)} panels:\n{paste0(selected_panel_name,collapse="\n")}'))}
  return(selected_panel_name)
}

# get similar panels from the distance matrix dataframe. 
# if the distance is based on jaccard, then a value of 0.95 means that the panels are only 5% similar 
# this function gets the top X% most similar panels, with a maximal distance of Y
# if names_only = TRUE - will return only the panel names
get_similar_panels_from_dm_df_by_prop<-function(dm_df,selected_panel_name,top_most_similar=0.05,maximal_distance=0.99,names_only=F){
  similar_panels_df<-dm_df%>%filter(row%in%selected_panel_name)%>%slice_min(prop = top_most_similar,order_by = value)%>%filter(value<maximal_distance)
  
  if (names_only){return(similar_panels_df%>%pull(col))}
  
  return(similar_panels_df)
}

get_similar_panels_from_dm_df_by_n<-function(dm_df,selected_panel_name,top_most_similar=6,maximal_distance=0.99){
  similar_panels_df<-dm_df%>%filter(row==selected_panel_name)%>%slice_min(n = top_most_similar,order_by = value)%>%filter(value<maximal_distance)
  return(similar_panels_df)
}

# split joined panel name
split_joined_panel_name<-function(joined_panel_names){
  components<-data.frame(stringr::str_split(joined_panel_names,pattern = '\\|',simplify = T))
  colnames(components)<-c('source','panel_id','panel_name')
  components$panel_joined_name<-joined_panel_names
  return(components)
}

# Get top similar panels
# Given a list of (joined) panel names (source|panel_id|panel_name), get the top similar panels 
# same_source - whether you want the similar panels to be selected from the same source as well
get_top_similar_panels_table<-function(joined_panel_names,
                                 source_db,
                                 dm_df,
                                 number_of_top_panels=6,
                                 same_source=F){
  
  split_names<-split_joined_panel_name(joined_panel_names)
  
  top_similars_per_selected<-dm_df%>%
    filter(row %in% joined_panel_names)%>%
    group_by(row)%>%left_join(split_names%>%dplyr::select(source,panel_joined_name),by=c('row'='panel_joined_name'))
  
  # if specified - remove the source from the comparison sources
  if (!same_source){
    top_similars_per_selected<-top_similars_per_selected%>%rowwise()%>%filter(!grepl(source,col))
  }  
  
  top_similars_per_selected<-top_similars_per_selected%>%group_by(row)%>%
    slice_min(n=number_of_top_panels,order_by=value)%>%
    filter(value!=1)
  
  comparison_table<-NULL
  for (selected_panel_name in joined_panel_names){
    message(glue('Analyzing {selected_panel_name}..'))
    panel_genes<-source_db%>%filter(panel_joined_name==selected_panel_name)%>%pull(gene_symbol)
    selected_num_o_genes<-length(panel_genes)
    most_similar_panels<-as.character(top_similars_per_selected%>%filter(row==selected_panel_name)%>%pull(col))
    
    for (similar_panel in most_similar_panels){
      similar_panel_genes<-source_db%>%filter(panel_joined_name==similar_panel)%>%pull(gene_symbol)
      similar_panel_num_o_genes<-length(similar_panel_genes)
      common_genes<-sum(panel_genes %in% similar_panel_genes)
      genes_not_in_similar_panel<-selected_num_o_genes-common_genes
      distance<-top_similars_per_selected%>%filter(row==selected_panel_name &col==similar_panel)%>%pull(value)
      comparison_table<-comparison_table%>%
        bind_rows(data.frame(selected_panel_name,
                             similar_panel,
                             selected_num_o_genes,
                             similar_panel_num_o_genes,
                             common_genes,
                             genes_not_in_similar_panel,
                             distance))
    }
  }
  return(comparison_table)
}

get_panel_comparison_metrics<-function(panelA_genes,panelB_genes){
  common_genes<-sum(panelA_genes %in% panelB_genes)
  Agenes_not_in_B<-paste0(setdiff(panelA_genes,panelB_genes),collapse='|')
  A_not_in_B<-length(panelA_genes)-common_genes
  B_not_in_A<-length(panelB_genes)-common_genes
  panel_comp_table<-data.frame(genes_in_a=length(panelA_genes),
                               genes_in_b=length(panelB_genes),
                               common_genes,A_not_in_B,B_not_in_A,
                               Agenes_not_in_B)
  return(panel_comp_table)
}

compare_panels_using_joined_names<-function(panel_a_joined,panel_b_joined,source_db){
  #message(glue('comparing: {panel_a_joined} vs {panel_b_joined}'))
  panelA_genes<-source_db%>%filter(panel_joined_name==panel_a_joined)%>%pull(gene_symbol)
  panelB_genes<-source_db%>%filter(panel_joined_name==panel_b_joined)%>%pull(gene_symbol)
  panel_comp_metrics<-get_panel_comparison_metrics(panelA_genes,panelB_genes)
  return(panel_comp_metrics)
}

compare_panel_vs_many_panels<-function(panel_a_joined,comparison_panel_joined_names,source_db,min_confidence_level=3){
  panel_genes<-source_db%>%filter(panel_joined_name==panel_a_joined)%>%filter(confidence_level>=min_confidence_level)%>%pull(gene_symbol)
  comparison_genes<-unique(source_db%>%filter(panel_joined_name %in% comparison_panel_joined_names)%>%pull(gene_symbol))
  panel_comp_metrics<-get_panel_comparison_metrics(panel_genes,comparison_genes)
  return(panel_comp_metrics)
}

# Compare panels list -
# given a list of panel joined names - create a binary matrix and output metrics regarding their relatedness
# also outputs a plot describing the number of genes that are found in each percentile of panels (for example how many genes are found in >90% of the samples)
# if panel to compare is given, check how the genes in it are represented in the other panels
get_panel_relatedness_stats <- function(panel_names,gtr_with_expert_db,panel_to_compare=NA) {
  panels_db<-gtr_with_expert_db%>%filter(panel_joined_name%in%panel_names)
  panels_bin<-create_binary_gene_table(panels_db)
  relatedness_stats<-list()
  relatedness_stats$gene_list<-colnames(panels_bin)
  relatedness_stats$num_of_panels<-length(panel_names)
  # number of genes in all panels
  relatedness_stats$num_of_genes<-ncol(panels_bin)
  # per gene number of panels
  relatedness_stats$panels_per_gene<-panels_db%>%group_by(gene_symbol)%>%summarize(npanels=n())%>%mutate(npanels_rate=npanels/relatedness_stats$num_of_panels)
  
  # gene content in panels stats 
  relatedness_stats$genes_per_panel<-gtr_with_expert_db%>%
    filter(panel_joined_name%in%panel_names)%>%
    group_by(panel_joined_name)%>%
    summarize(ngenes=n())%>%pull(ngenes)%>%skimr::skim_without_charts()
  per_gene_stats<-data.frame(genes=colnames(panels_bin),
                             percent_of_panels=panels_bin%>%colMeans())%>%
    mutate(percent_of_panels_cat=cut(percent_of_panels,breaks=seq(0,1,0.1),labels=c('0-10%','10-20%','20-30%','30-40%','40-50%','50-60%','60-70%','70-80%','80-90%','90-100%')))
  relatedness_stats$percent_of_panels<-per_gene_stats%>%group_by(percent_of_panels_cat)%>%
    summarize(n=n(),rate=n/relatedness_stats$num_of_genes)
  # generate a plot for the number of genes vs percent of panels
  relatedness_stats$gene_percent_plot<-
    per_gene_stats%>%
    ggplot(aes(x=percent_of_panels_cat))+
    geom_bar()+
    labs(y='Number of genes',x='Percent of panels')+
    theme_minimal()
  if (!is.na(panel_to_compare)){
    # check what is the status only for the genes in the panel to compare
    relatedness_stats$panel_to_compare<-panel_to_compare
    panel_to_compare_genes<-gtr_with_expert_db%>%filter(panel_joined_name==panel_to_compare)%>%
      dplyr::select(gene_symbol,confidence_level)
    relatedness_stats$panel_to_compare_genes<-panel_to_compare_genes%>%pull(gene_symbol)
    panel_to_compare_genes<-panel_to_compare_genes%>%pull(gene_symbol)
    relatedness_stats$panel_to_compare_number_of_genes<-length(panel_to_compare_genes)
    relatedness_stats$panel_to_compare_percent_of_panels<-per_gene_stats%>%
      filter(genes%in%panel_to_compare_genes)%>%
      group_by(percent_of_panels_cat)%>%
      summarize(n=n(),rate=n/relatedness_stats$panel_to_compare_number_of_genes)
    relatedness_stats$panel_to_compare_unique_genes<-relatedness_stats$panels_per_gene%>%filter(gene_symbol%in%panel_to_compare_genes,npanels==1)%>%pull(gene_symbol)
    
  }
  
  return(relatedness_stats)
}

# A function that takes in the output of the get_panel_relatedness_stats function and (optinally) the 
# summarize_panel_vs_phenotype_pubmed function output and produces a text summary
get_relatedness_summary_text<-function(phenotype_name,
                                       phenotype_relatedness,
                                       phenotype_distance_summary,
                                       panel_vs_pubmed=NA){
  genes_in_top_0.9_panels<-phenotype_relatedness$panels_per_gene%>%filter(npanels_rate>=0.9)%>%pull(gene_symbol)
  #genes_in_top_0.9_panels<-phenotype_relatedness$percent_of_panels%>%filter(percent_of_panels_cat=="90-100%")%>%pull(n)
  num_genes_in_top_0.9_panels<-length(genes_in_top_0.9_panels)
  perc_genes_in_top_0.9_panels<-round(num_genes_in_top_0.9_panels/phenotype_relatedness$num_of_genes*100,1)
  relatedness_text<-paste0(glue('For {phenotype_name}, we found {phenotype_relatedness$num_of_panels} panels in the GTR with {phenotype_relatedness$num_of_genes} genes (median {phenotype_relatedness$genes_per_panel$numeric.p50} genes, range [{phenotype_relatedness$genes_per_panel$numeric.p0},{phenotype_relatedness$genes_per_panel$numeric.p100}]).'),
                           glue('Out of the list of genes, only {num_genes_in_top_0.9_panels} genes ({perc_genes_in_top_0.9_panels}%) were found in more than 90% of the panels ({paste0(genes_in_top_0.9_panels,collapse=", ")}).'),
                           glue('There were {phenotype_relatedness$percent_of_panels%>%filter(percent_of_panels_cat=="0-10%")%>%pull(n)} genes ({round(phenotype_relatedness$percent_of_panels%>%filter(percent_of_panels_cat=="0-10%")%>%pull(rate)*100,1)}%) that were found in less than 10% of the panels.'),
                           glue('The average distance between the panels was {phenotype_distance_summary$numeric.mean} (IQR [{phenotype_distance_summary$numeric.p25},{phenotype_distance_summary$numeric.p75}])'),
                           collapse=' ')
  
  # if the panels were compared against pubmed
  if (!is.na(panel_vs_pubmed)){
    ngenes_no_pub<-sum(!panel_vs_pubmed$has_pub)
    ngenes_with_pub_0.1_panels<-nrow(panel_vs_pubmed%>%filter(n_panels_rate<0.1,has_pub))
    relatedness_text<-glue('{relatedness_text} Among those genes, {ngenes_with_pub_0.1_panels} appeared in at least one publication regarding {phenotype_name}. For {ngenes_no_pub} genes ({round(100*ngenes_no_pub/phenotype_relatedness$num_of_genes,1)}%), a pubmed search did not identify any paper that associates them with {phenotype_name}.',
                           collapse=' ')
  }
  
  # if the panels were compared against a specific panel
  if (!is.null(phenotype_relatedness)){
    
    relatedness_text<-paste0(relatedness_text,
                             glue('The {phenotype_relatedness$panel_to_compare} panel contains {phenotype_relatedness$panel_to_compare_number_of_genes} genes.'),
                             glue('There were {length(phenotype_relatedness$panel_to_compare_unique_genes)} genes that were unique to this panel and were not suggested by other panels ({paste0(phenotype_relatedness$panel_to_compare_unique_genes,collapse=", ")}).'),
                             collapse=' ')
    if (!is.na(panel_vs_pubmed)){
      panel_genes_without_pub<-panel_vs_pubmed%>%filter(gene_symbol %in% phenotype_relatedness$panel_to_compare_genes,!has_pub)%>%pull(gene_symbol)
      relatedness_text<-paste0(relatedness_text,
                               glue('For {length(panel_genes_without_pub)} genes in the panel ({paste0(panel_genes_without_pub,collapse=", ")}), we could not find any publication associating them with {phenotype_name}.'),
                               collapse=' ')
    }
  }
  return(relatedness_text)
}


distance_summary_in_panels_group<-function(panels_to_summarize,gtr_dm_df){
  message(glue('calculating distance summary for {length(panels_to_summarize)} panels'))
  panels_dm_df<-subset_dm_df_by_gene_panel(panels_to_summarize,gtr_dm_df)
  distance_summary<-panels_dm_df%>%select(value)%>%skimr::skim_without_charts()
  return(distance_summary)
}

# given a phenotype name and selected panels produce:
# -- a relatedness analysis 
# -- a distance summary
# -- a pubmed analysis (if the all genes vs pubmed file is given)
# -- a summary text
panels_discrepancy_analysis<-function(phenotype_name,
                                      phenotype_panels,
                                      gtr_with_expert_db,
                                      panel_to_compare,
                                      gtr_dm_df,
                                      all_genes_vs_pubmed=NA){
  discrepancy_analysis_res<-list()
  discrepancy_analysis_res$relatedness<-get_panel_relatedness_stats(phenotype_panels,
                                                                    gtr_with_expert_db,
                                                                    panel_to_compare = panel_to_compare)
  
  discrepancy_analysis_res$distance_summary<-distance_summary_in_panels_group(phenotype_panels,gtr_dm_df)
  
  # run pubmed search for each gene in the panels 
  if (is.na(all_genes_vs_pubmed)){
    discrepancy_analysis_res$panel_genes_pubmed_search<-pubmed_gene_list_vs_phenotype(phenotype_name,discrepancy_analysis_res$relatedness$gene_list)
  }else{
    discrepancy_analysis_res$panel_genes_pubmed_search<-all_genes_vs_pubmed%>%filter(gene_symbol%in%discrepancy_analysis_res$relatedness$gene_list)
  }
  discrepancy_analysis_res$panel_vs_pubmed <-
    summarize_panel_vs_phenotype_pubmed(
      selected_panels = phenotype_panels,
      gtr_with_expert_db = gtr_with_expert_db,
      panel_genes_pubmed_search = discrepancy_analysis_res$panel_genes_pubmed_search,
      panelapp_panel = panel_to_compare
    )
  
  # discrepancy_analysis_res$num_publications_per_gene_vs_in_panelapp_plot<-
  #   plot_num_of_publications_per_gene_vs_in_panelapp(discrepancy_analysis_res$panel_vs_pubmed)
  # discrepancy_analysis_res$num_of_publications_vs_rate_of_panels_plot<-
  #   plot_num_of_publications_vs_rate_of_panels(discrepancy_analysis_res$panel_vs_pubmed)
  # discrepancy_analysis_res$num_of_genes_with_publications_vs_rate_of_panels_plot<-
  #   plot_num_of_genes_with_publications_vs_rate_of_panels(discrepancy_analysis_res$panel_vs_pubmed)
  # discrepancy_analysis_res$num_of_genes_with_publications_vs_rate_of_panels_plot<-
  #   plot_num_of_genes_with_publications_vs_rate_of_panels_with_panelapp(discrepancy_analysis_res$panel_vs_pubmed)
  
  discrepancy_analysis_res$summary_text<-get_relatedness_summary_text(phenotype_name,
                                                                      discrepancy_analysis_res$relatedness,
                                                                      discrepancy_analysis_res$distance_summary,
                                                                      discrepancy_analysis_res$panel_vs_pubmed)
  return(discrepancy_analysis_res)
}

# a function that takes in a panelapp panel, a phenotype name, a search text to identify panels for the naive approach, 
# and calculates sensitivity, ppv and f1 and discrepancy analysis for:
# -- the panelapp panel
# -- every panel in the naive search
# -- a unification of all the naive search panels, plus panelapp and without
# -- for a panel corresponding to all the panels in the same cluster as the panelapp panel
compare_panelapp_naive_and_cluster<-function(panelapp_panel,
                             phenotype_name,
                             naive_search,
                             max_dists_for_cluster=c(0.75,0.85,0.95),
                             all_panel_names,
                             all_genes,
                             gtr_with_expert_db,
                             gtr_dm_df,
                             pre_generated_gene_publications_list){
  
  test_res<-list()
  
  all_genes_vs_phenotype<-pubmed_gene_list_vs_phenotype(phenotype_name = phenotype_name,
                                                        gene_list = all_genes,
                                                        pre_generated_genes_pubs = pre_generated_gene_publications_list)
  
  all_genes_with_pheno_pub<-all_genes_vs_phenotype%>%filter(pub_num>1)%>%pull(gene_symbol)
  message(glue('Out of {length(all_genes)} genes, {length(all_genes_with_pheno_pub)} had a publication mentioning {phenotype_name}'))
  
  message(glue('Collecting metrics for panelapp..'))
  # panelapp metrics
  test_res$all_metrics<-NULL
  test_res$panelapp_genes<-gtr_with_expert_db%>%filter(panel_joined_name==panelapp_panel)%>%pull(gene_symbol)
  metrics_panelapp<-calculate_sens_and_ppv_vs_pubmed(test_res$panelapp_genes,all_genes_with_pheno_pub)
  metrics_panelapp$group<-panelapp_panel
  test_res$all_metrics<-test_res$all_metrics%>%bind_rows(as.data.frame(metrics_panelapp))
  
  message(glue('Collecting metrics for naive search (with panelapp)..'))
  # naive search - with panelapp
  pheno_naive_panels<-grep(naive_search,all_panel_names,ignore.case = T,value=T)
  
  if (length(pheno_naive_panels)==0){
    stop(glue('Error: could not find any panel corresponding to the submitted naive search: {naive_search}'))}
  
  # first calculate performance per panel
  for (pheno_panel in pheno_naive_panels){
    if (grepl('panelapp',pheno_panel)){next}#no need to calculate again
    metrics_pheno_panel<-calculate_sens_and_ppv_vs_pubmed(gtr_with_expert_db%>%
                                                            filter(panel_joined_name==pheno_panel)%>%
                                                            pull(gene_symbol),
                                                          all_genes_with_pheno_pub)
    metrics_pheno_panel$group<-pheno_panel
    test_res$all_metrics<-test_res$all_metrics%>%bind_rows(as.data.frame(metrics_pheno_panel))
  }
  
  # then claculate performance for all panels
  
  test_res$naive_discrepancy<-panels_discrepancy_analysis(phenotype_name,
                                                          pheno_naive_panels,
                                                          gtr_with_expert_db,gtr_dm_df = gtr_dm_df,
                                                          panelapp_panel,
                                                          all_genes_vs_pubmed = all_genes_vs_phenotype)
  
  test_res$naive_with_panelapp_genes<-test_res$naive_discrepancy$relatedness$gene_list
  metrics_naive<-calculate_sens_and_ppv_vs_pubmed(test_res$naive_with_panelapp_genes,
                                                  all_genes_with_pheno_pub)
  metrics_naive$group='naive_search - with panelapp'
  test_res$all_metrics<-test_res$all_metrics%>%bind_rows(as.data.frame(metrics_naive))
  message(glue('Collecting metrics for naive search (no panelapp)..'))
  # naive search - no panelapp
  pheno_naive_panels_no_panelapp<-grep('panelapp',pheno_naive_panels,invert = T,value=T)
  
  if (length(pheno_naive_panels_no_panelapp)==0){
    stop(glue('Error: could not find any panel corresponding to the submitted naive search: {naive_search}'))}
  
  # then claculate performance for all panels
  
  test_res$naive_no_panelapp_discrepancy<-panels_discrepancy_analysis(phenotype_name,
                                                                      pheno_naive_panels_no_panelapp,
                                                                      gtr_with_expert_db,gtr_dm_df = gtr_dm_df,
                                                                      panelapp_panel,
                                                                      all_genes_vs_pubmed = all_genes_vs_phenotype)
  
  test_res$naive_no_panelapp_genes<-test_res$naive_no_panelapp_discrepancy$relatedness$gene_list
  metrics_naive<-calculate_sens_and_ppv_vs_pubmed(test_res$naive_no_panelapp_genes,
                                                  all_genes_with_pheno_pub)
  metrics_naive$group='naive_search - no panelapp'
  test_res$all_metrics<-test_res$all_metrics%>%bind_rows(as.data.frame(metrics_naive))
  
  # Cluster based panel generation
  if (!is.na(max_dists_for_cluster)){
    message(glue('Collecting metrics for cluster based panel generation..'))
    for (max_dist in max_dists_for_cluster){
      dist_panels<-get_similar_panels_from_dm_df_by_prop(gtr_dm_df,
                                                         selected_panel_name = panelapp_panel,
                                                         top_most_similar = 1,
                                                         maximal_distance = max_dist,
                                                         names_only = T)
      
      #setdiff(dist_panels,cluster_panels)
      
      discrepancy_name<-glue('dist_discrepancy_{max_dist}')
      test_res[[discrepancy_name]]<-panels_discrepancy_analysis(phenotype_name,
                                                                dist_panels,
                                                                gtr_with_expert_db,gtr_dm_df= gtr_dm_df,
                                                                panelapp_panel,
                                                                all_genes_vs_pubmed = all_genes_vs_phenotype)
      
      test_res[[glue('{discrepancy_name}_genes')]]<-test_res[[discrepancy_name]]$relatedness$gene_list
      metrics_dist<-calculate_sens_and_ppv_vs_pubmed(test_res[[glue('{discrepancy_name}_genes')]],all_genes_with_pheno_pub)
      metrics_dist$group=glue('max_dist={max_dist}')
      
      test_res$all_metrics<-test_res$all_metrics%>%bind_rows(as.data.frame(metrics_dist))
    }
  }
  return(test_res)
}
