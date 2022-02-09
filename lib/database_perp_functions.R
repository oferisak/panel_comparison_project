fix_source_db<-function(source_db){
  fixed_source_db <- source_db %>%
    mutate(gene_symbol = trimws(gene_symbol)) %>%
    mutate(gene_symbol = stringr::str_replace(string = gene_symbol, '\\(.+', '')) %>%
    mutate(gene_symbol = stringr::str_replace(string = gene_symbol, '^.+;', '')) %>%
    mutate(gene_symbol = stringr::str_replace(string = gene_symbol, 'c\\..+', '')) %>%
    filter(stringr::str_length(gene_symbol)<15)%>%
    separate_rows(gene_symbol,sep='\\.')%>%
    filter(!grepl('SMN1[/,]', gene_symbol)) %>%
    mutate(gene_symbol = stringr::str_replace_all(string = gene_symbol, '\\/+|!', '')) %>%
    bind_rows(
      source_db %>% filter(grepl('SMN1[/,]', gene_symbol)) %>%
        separate_rows(gene_symbol, sep = '[/,]')
    )%>%
    separate_rows(gene_symbol,sep=',')%>%
    mutate(gene_symbol=stringr::str_replace_all(gene_symbol,'\\s',''))%>%
    filter(gene_symbol != '')%>%filter(!grepl('\\)',gene_symbol))
  return(fixed_source_db)
}

remove_red_genes_from_panelapp<-function(source_db){
  source_db<-source_db%>%filter(!(source=='panelapp' & confidence_level==1))
  return(source_db)
}

# refuted and no known gene-disease association genes were already removed from the list
remove_disputed_genes_from_clingen<-function(source_db){
  source_db<-source_db%>%filter(!(source=='clingen' & confidence_level==1))
  return(source_db)
}

