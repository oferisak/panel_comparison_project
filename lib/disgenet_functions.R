library(disgenet2r)

get_dg_api_key<-function(){
  email<-readline(prompt="Enter email: ")
  password<-readline(prompt="Enter password: ")
  api_key<-disgenet2r::get_disgenet_api_key(email = email,password = password)
  return(api_key)
}

get_disgenet_gene_vs_phenotype_data<-function(api_key,all_genes,database='ALL'){
  
  nchunks<-round(length(all_genes)/500)+1
  geneChunks<-split(all_genes,f=c(1:nchunks))
  message(glue('will run {nchunks} search chunks..'))
  first_chunk<-1
  
  disgenet_chunks<-list()
  for (i in first_chunk:length(names(geneChunks))){
    #for (idChunk in idChunks){
    geneChunk<-geneChunks[[i]]
    message(glue('Analyzing chunk number {i} ({length(geneChunk)} genes)'))
    #print(sprintf('Analyzing %s',idChunk))
    chunk_disgenet <- gene2disease( gene = geneChunk, verbose = TRUE, api_key = api_key,database =database)
    disgenet_chunks[[i]]<-chunk_disgenet
  }
  
  disgenet_df<-NULL
  for (i in 1:length(disgenet_chunks)){
    message(glue('Adding chunk {i}'))
    disgenet_df<-disgenet_df%>%bind_rows(disgenet_chunks[[i]]@qresult)
  }
  
  return(disgenet_df)
}


get_disgenet_gene_vs_evidence_data<-function(api_key,all_genes,database='ALL'){
  
  nchunks<-round(length(all_genes)/500)+1
  geneChunks<-split(all_genes,f=c(1:nchunks))
  message(glue('will run {nchunks} search chunks..'))
  first_chunk<-1
  
  disgenet_chunks<-list()
  for (i in first_chunk:length(names(geneChunks))){
    #for (idChunk in idChunks){
    geneChunk<-geneChunks[[i]]
    message(glue('Analyzing chunk number {i} ({length(geneChunk)} genes)'))
    #print(sprintf('Analyzing %s',idChunk))
    chunk_disgenet <- gene2evidence( gene = geneChunk, verbose = TRUE, api_key = api_key,database =database)
    disgenet_chunks[[i]]<-chunk_disgenet
  }
  
  disgenet_df<-NULL
  for (i in 1:length(disgenet_chunks)){
    message(glue('Adding chunk {i}'))
    disgenet_df<-disgenet_df%>%bind_rows(disgenet_chunks[[i]]@qresult)
  }
  
  return(disgenet_df)
}
