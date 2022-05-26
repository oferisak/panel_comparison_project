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
  source_db<-source_db%>%filter(!(source=='panelapp' & confidence_level<=1))
  return(source_db)
}


# refuted and no known gene-disease association genes were already removed from the list
remove_disputed_genes_from_clingen<-function(source_db){
  source_db<-source_db%>%filter(!(source=='clingen' & confidence_level==1))
  return(source_db)
}

# fix gene names to be the approved ones
fix_gene_names<-function(source_db){
  standard_genenames_table<-readr::read_delim('./data/accessory_files/standard_genenames_db_2022-04-07.csv')
  double_prev<-standard_genenames_table%>%filter(!is.na(prev_symbol))%>%count(prev_symbol)%>%filter(n>1)
  # manually fix them 
  manually_fixed_genes<-c('QARS')
  source_db<-source_db%>%mutate(gene_symbol=case_when(
    gene_symbol=='QARS' & grepl('Leukodystrophy',panel_name)~'EPRS1',
    gene_symbol=='QARS' & !grepl('Leukodystrophy',panel_name)~'QARS1',
    TRUE~gene_symbol))
  source_db<-source_db%>%left_join(standard_genenames_table%>%filter(!is.na(prev_symbol)),by=c('gene_symbol'='prev_symbol'))
  # verify there are no cases of panel db with a gene symbol that have more than one approved gene name
  genes_with_more_than_one_name<-fixed_source_db%>%filter(gene_symbol%in%setdiff(double_prev$prev_symbol,manually_fixed_genes))%>%pull(gene_symbol)%>%unique()
  if (length(genes_with_more_than_one_name)>0){stop(glue('There are genes with more than one approved gene symbol ({paste0(genes_with_more_than_one_name,collapse=",")})'))}
  # num_o_fixed_occurances<-sum(!is.na(source_db$approved_gene_symbol))
  # num_o_fixed_genenames<-length(unique(source_db%>%filter(!is.na(approved_gene_symbol))%>%pull(approved_gene_symbol)))
  # message(glue('Fixed {num_o_fixed_occurances} occurances of {num_o_fixed_genenames} gene names'))
  source_db<-source_db%>%mutate(gene_symbol=ifelse(!is.na(approved_gene_symbol),approved_gene_symbol,gene_symbol))
  return(source_db)
}


# build gencc database
gencc_raw<-readr::read_delim('./data/databases/gencc-submissions_20220327.tsv')
gencc_db<-gencc_raw%>%
  filter(!classification_title%in%c('Refuted Evidence','No Known Disease Relationship','Disputed Evidence'))%>%
  mutate(confidence_level=case_when(
    grepl('Definitive|Strong',classification_title)~3,
    grepl('Supportive|Moderate',classification_title)~2,
    grepl('Limited',classification_title)~1,
  ))%>%
  select(panel_name=disease_title,
         gene_symbol,
         associated_phenotypes=disease_title,
         confidence_level,
         uuid)%>%
  mutate(source='gencc',
         panel_category=NA,
         panel_comment='.',
         panel_url='.',
         panel_joined_name=as.character(glue('{source}|{uuid}|{panel_name}')))%>%
  group_by(panel_name,gene_symbol)%>%
  slice_max(order_by = confidence_level,n=1,with_ties = F)%>%ungroup()
# remove all diseases with only one gene
one_gene_diseases<-gencc_db%>%count(panel_name)%>%filter(n==1)%>%pull(panel_name)
gencc_db<-gencc_db%>%filter(!panel_name%in%one_gene_diseases)
#z<-gencc_db%>%count(panel_name)




