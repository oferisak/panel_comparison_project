get_source_stats<-function(source_db){
  source_stats <- source_db %>% group_by(source) %>%
    summarize(
      number_of_panels = length(unique(panel_name)),
      number_of_genes = length(unique(gene_symbol))
    ) %>%
    left_join(source_db%>%
                group_by(source,panel_name)%>%
                summarize(number_of_genes=n())%>%
                ungroup()%>%group_by(source)%>%
                summarize(skimr::skim_without_charts(number_of_genes),
                          panels_with_only_one_gene=sum(number_of_genes==1)))%>%
    select(-contains(c('skim_','complete','missing','p75','p25','.sd')))%>%
    mutate(source=stringr::str_replace_all(source,'_',' '))%>%
    rename(
      'Panel source' = source,
      'Panels' = number_of_panels,
      'Genes' = number_of_genes,
      'Genes per panel' = numeric.mean,
      'Largest panel' = numeric.p100,
      'Smallest panel' = numeric.p0,
      'Median number of genes' = numeric.p50,
      'One gene panels' = panels_with_only_one_gene
    )
  return(source_stats)
}

generate_gene_panel_cat_number_table<-function(panels_list){
  gene_num_cat_table<-panels_list%>%group_by(source)%>%
    count(number_of_genes_cat)%>%
    arrange(number_of_genes_cat)%>%
    pivot_wider(names_from = number_of_genes_cat,values_from = n,values_fill = 0)%>%
    ungroup()
  return(gene_num_cat_table)
}

generate_gtr_gene_panel_cat_number_table<-function(gtr_panels_list){
  gene_num_cat_table<-gtr_panels_list%>%group_by(name_of_laboratory)%>%
    count(number_of_genes_cat)%>%
    arrange(number_of_genes_cat)%>%
    pivot_wider(names_from = number_of_genes_cat,values_from = n,values_fill = 0)%>%
    ungroup()
  return(gene_num_cat_table)
}
