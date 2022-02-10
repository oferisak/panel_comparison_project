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
  similartiy_matrix<-clingen_similarity%>%select(c('source','panel_id','panel_name'))%>%distinct()
  for (colname in setdiff(colnames(clingen_similarity),c('source','panel_id','panel_name'))){
    print(colname)
    similartiy_matrix<-similartiy_matrix%>%
      left_join(clingen_similarity%>%
                  select('source','panel_id','panel_name',colname)%>%
                  drop_na())
  }
  return(similartiy_matrix)
}

create_binary_gene_table<-function(source_db){
  # remove panelapp low confidence genes from source db
  binary_panel_gene_table <-
    source_db %>% filter(!(source == 'panelapp' &
                             confidence_level == 1)) %>%
    unite(col = panel, sep = '|', source, panel_id, panel_name) %>%
    select(panel, gene_symbol) %>%
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
get_similar_panels_from_dm_df_by_prop<-function(dm_df,selected_panel_name,top_most_similar=0.05,maximal_distance=0.99){
  similar_panels_df<-dm_df%>%filter(row==selected_panel_name)%>%slice_min(prop = top_most_similar,order_by = value)%>%filter(value<maximal_distance)
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
    group_by(row)%>%left_join(split_names%>%select(source,panel_joined_name),by=c('row'='panel_joined_name'))
  
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
