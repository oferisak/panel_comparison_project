hgnc<-readr::read_delim('/media/SSD/Bioinformatics/Databases/hgnc/hgnc_complete_set.txt')

# Get ALL the publications for the given gene list
pubmed_grab_all_gene_publications<-function(gene_list){
  genes_to_test<-data.frame(gene_symbol=gene_list)
  # add entrez_id to gene symbol in genes to test
  genes_to_test<-genes_to_test%>%left_join(hgnc%>%dplyr::select(symbol,entrez_id),by=c('gene_symbol'='symbol'))
  # remove genes without entrez id
  genes_to_test<-genes_to_test%>%filter(!is.na(entrez_id))
  # run entrez search for the given gene ids
  message(glue('running pubmed search for {length(gene_list)} genes..'))
  
  ids<-as.character(genes_to_test$entrez_id)
  nchunks<-round(length(ids)/200)+1
  idChunks<-split(ids,f=c(1:nchunks))
  message(glue('will run {nchunks} search chunks..'))
  entrez_pubmed_res<-NA
  first_chunk<-1
  tmp_chunk_file_name<-glue('./data/tmp_files/tmp_chunk')
  tmp_res_file_name<-glue('./data/tmp_files/tmp_pubmed_res.RData')
  if (file.exists(tmp_chunk_file_name)){
    message('found tmp file. will parse it and continue from that position')
    chunk_df<-readr::read_delim(tmp_chunk_file_name)
    last_done_id<-chunk_df%>%pull(id_num)
    message(glue('the chunk successfully tested was {last_done_id}. will continue from {last_done_id+1}'))
    first_chunk<-last_done_id+1
    load(tmp_res_file_name)
    message(glue('Loaded the pre-ran pubmed search. contains {length(names(entrez_pubmed_res))} genes'))
  }
  for (i in first_chunk:length(names(idChunks))){
    #for (idChunk in idChunks){
    idChunk<-idChunks[[i]]
    message(glue('Analyzing chunk number {i} ({length(idChunk)} genes)'))
    #print(sprintf('Analyzing %s',idChunk))
    chunkLinks  <- entrez_link(db="pubmed", dbfrom="gene",id=idChunk,by_id=TRUE)
    print(class(chunkLinks))
    pubmedIDs<-lapply(chunkLinks, function(x) x$links$gene_pubmed)
    names(pubmedIDs)<-idChunk
    if (is.na(entrez_pubmed_res)){
      entrez_pubmed_res<-pubmedIDs
    }else{
      entrez_pubmed_res<-c(entrez_pubmed_res,pubmedIDs)
    }
    write.table(data.frame(id_num=i,status='DONE'),file=tmp_chunk_file_name,row.names = F,quote = F)
    save(entrez_pubmed_res,file = tmp_res_file_name)
  }
  message('Finished pubmed query. will save the search results for all the genes')
  save(entrez_pubmed_res,file = './data/pre_generated_data/genes_pubmed_search.RData')
  unlink(tmp_chunk_file_name)
  unlink(tmp_res_file_name)
}

# A function that takes in a phenotype name and a list of genes and returns how many publications connect each gene with the phenotype
pubmed_gene_list_vs_phenotype<-function(phenotype_name,gene_list,pre_generated_genes_pubs=NA,with_mesh=F){
  phenotype_pubmed_search_file<-glue('./data/pre_generated_data/{phenotype_name}_pubmed_search.RData')
  if (file.exists(phenotype_pubmed_search_file)){
    message(glue('pubmed search for {phenotype_name} already performed. will load it..'))
    load(phenotype_pubmed_search_file)
    return(gene_phenotype_pubmed_search)
  }
  if (with_mesh){
    phenotype_query<-glue('{phenotype_name}[mesh]')
  }else{
    phenotype_query<-glue('{phenotype_name}')
  }
  message(glue('running pubmed search for {phenotype_name}..'))
  phenotype_search_res<- entrez_search(db="pubmed", term=phenotype_query,retmax=999999)$ids
  message(glue('found {length(phenotype_search_res)} publications for {phenotype_name}'))
  # Get the gene IDs
  genes_to_test<-data.frame(gene_symbol=gene_list)
  # add entrez_id to gene symbol in genes to test
  genes_to_test<-genes_to_test%>%left_join(hgnc%>%dplyr::select(symbol,entrez_id),by=c('gene_symbol'='symbol'))
  # remove genes without entrez id
  genes_to_test<-genes_to_test%>%filter(!is.na(entrez_id))
  
  entrez_pubmed_res<-pre_generated_genes_pubs

  #entrez_pubmed_res<-entrez_link(db="pubmed", dbfrom="gene",id=as.character(genes_to_test$entrez_id),by_id=TRUE)
  # retrieve the pubmed ids from the query result
  #gene_pubmed_ids<-lapply(entrez_pubmed_res, function(x) x$links$gene_pubmed)
  #names(entrez_pubmed_res)<-genes_to_test%>%pull(gene_symbol)
  gene_phenotype_pubmed_search<-NULL
  for (eid in names(entrez_pubmed_res)){
    gene_symbol<-genes_to_test%>%filter(entrez_id==eid)%>%pull(gene_symbol)
    gene_pubmeds<-entrez_pubmed_res[[eid]]
    gene_phenotype_pubmeds<-intersect(phenotype_search_res,gene_pubmeds)
    gene_phenotype_pubmed_search<-gene_phenotype_pubmed_search%>%
      bind_rows(data.frame(gene_symbol,pub_num=length(gene_phenotype_pubmeds),pubmed_ids=paste0(gene_phenotype_pubmeds,collapse=',')))
  }
  save(gene_phenotype_pubmed_search,file=phenotype_pubmed_search_file)
  return(gene_phenotype_pubmed_search)
}

# summarize phenotype vs genes results
summarize_panel_vs_phenotype_pubmed<-function(selected_panels,
                                              panel_genes_pubmed_search,
                                              gtr_with_expert_db,
                                              panelapp_panel=NA,
                                              min_num_of_pub_to_consider_positive=2){
  panel_vs_pubmed<-panel_genes_pubmed_search%>%
    left_join(gtr_with_expert_db%>%filter(panel_joined_name%in%selected_panels)%>%
                dplyr::select(gene_symbol,panel_joined_name))
  if (!is.na(panelapp_panel)){
    panelapp_genes<-gtr_with_expert_db%>%
      filter(panel_joined_name==panelapp_panel)%>%pull(gene_symbol)
    panel_vs_pubmed<-panel_vs_pubmed%>%
      mutate(in_panelapp=gene_symbol%in%panelapp_genes)%>%
      group_by(gene_symbol,pub_num,in_panelapp)
  }else{
    panel_vs_pubmed<-panel_vs_pubmed%>%group_by(gene_symbol,pub_num)
  }
  panel_vs_pubmed<-panel_vs_pubmed%>%
    summarize(n_panels=n())%>%
    mutate(n_panels_rate=n_panels/length(unique(panel_vs_pubmed$panel_joined_name)),
           has_pub=pub_num>=min_num_of_pub_to_consider_positive)
  return(panel_vs_pubmed)
}

# Plotting functions

plot_num_of_publications_per_gene_vs_in_panelapp<-function(panel_vs_pubmed){
  g<-panel_vs_pubmed%>%ggplot(aes(x=in_panelapp,y=pub_num))+
    geom_boxplot(outlier.size = 1)+
    theme_minimal()+
    labs(y='Number of publications',x='Is gene in panelapp')
  print(g)
  return(g)
}

plot_num_of_publications_vs_rate_of_panels<-function(panel_vs_pubmed){
  g<-panel_vs_pubmed%>%ggplot(aes(x=n_panels_rate,y=pub_num,col=in_panelapp))+
    geom_point(alpha=0.5)+
    theme_minimal()+
    ggsci::scale_color_nejm()+
    labs(y='Number of publications',x='Rate of panels')
  print(g)
  return(g)
}

# plot number of genes with publication vs the percent of panels that have the gene
plot_num_of_genes_with_publications_vs_rate_of_panels<-function(panel_vs_pubmed){
  g<-panel_vs_pubmed%>%mutate(has_pub=pub_num>0,
                              n_panels_rate_cat=cut(n_panels_rate,breaks=seq(0,1,0.1),labels=c('0-10%','10-20%','20-30%','30-40%','40-50%','50-60%','60-70%','70-80%','80-90%','90-100%')))%>%
    ggplot(aes(x=n_panels_rate_cat,fill=has_pub))+
    geom_bar()+
    #facet_wrap(in_panelapp~.,nrow = 2,scales='free_y')+
    ggsci::scale_fill_nejm()+
    theme_minimal()+theme(legend.position = 'bottom')+
    labs(x='Percent of panels with gene',y='Number of genes',fill='Phenotype-associated publication')
  print(g)
  g
}

# plot number of genes with publication vs the percent of panels that have the gene and mark genes found in a given panel
plot_num_of_genes_with_publications_vs_rate_of_panels_with_panelapp<-function(panel_vs_pubmed){
  g<-panel_vs_pubmed%>%mutate(has_pub=pub_num>0,
                              gene_group=ifelse(in_panelapp,'PanelApp gene',ifelse(has_pub,'Has publication','No publication')),
                              n_panels_rate_cat=cut(n_panels_rate,breaks=seq(0,1,0.1),labels=c('0-10%','10-20%','20-30%','30-40%','40-50%','50-60%','60-70%','70-80%','80-90%','90-100%')))%>%
    ggplot(aes(x=n_panels_rate_cat,fill=gene_group))+
    geom_bar()+
    #facet_wrap(in_panelapp~.,nrow = 2,scales='free_y')+
    ggsci::scale_fill_nejm()+
    theme_minimal()+theme(legend.position = 'bottom')+
    labs(x='Percent of panels with gene',y='Number of genes',fill='Phenotype-associated publication')
  print(g)
  g
}

# given a gene list generated by some method and a list of genes that match the phenotype by pubmed - 
# calculate sensitivity - what is the rate of published genes that are found in the generated list 
# calculate ppv - rate of generated list genes that are published
calculate_sens_and_ppv_vs_pubmed<-function(generated_gene_list,pubmed_gene_list,precision_recall_naming=T){
  gene_list_vs_pubmed<-list()
  gene_list_vs_pubmed$ngenes=length(generated_gene_list)
  gene_list_vs_pubmed$ppv=sum(generated_gene_list %in% pubmed_gene_list)/length(generated_gene_list)
  gene_list_vs_pubmed$sensitivity=sum(generated_gene_list %in% pubmed_gene_list)/length(pubmed_gene_list)
  gene_list_vs_pubmed$f1=(2*gene_list_vs_pubmed$ppv*gene_list_vs_pubmed$sensitivity)/(gene_list_vs_pubmed$ppv+gene_list_vs_pubmed$sensitivity)
  if (precision_recall_naming){
    gene_list_vs_pubmed<-gene_list_vs_pubmed%>%listr::list_rename(precision=ppv,recall=sensitivity)
  }
  return(gene_list_vs_pubmed)
}
