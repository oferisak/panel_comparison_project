# Load all the sources into one table
load_source_dbs <- function() {
  db_files<-list.files('./data/databases/other_sources/')
  source_db<-NULL
  db_versions<-NULL
  for (db_file in db_files){
    db_name<-db_file%>%stringr::str_replace('_db.*','')
    db_version<-db_file%>%stringr::str_match('\\d+-\\d+-\\d+')%>%as.character()
    db_versions<-db_versions%>%bind_rows(data.frame(source=db_name,db_creation_date=db_version))
    message(glue('Parsing {db_name}..'))
    db_table<-readr::read_delim(glue('./data/databases/other_sources/{db_file}'))
    db_table$panel_id<-as.character(db_table$panel_id)
    source_db<-source_db%>%bind_rows(db_table)
  }
  source_db<-source_db%>%mutate(panel_joined_name=glue('{source}|{panel_id}|{panel_name}'))
  return(source_db)
}  

source_db<-load_source_dbs()
# fix the source file
cache("fixed_source_db", {
  # remove red genes from panelapp  
  fixed_source_db<-remove_red_genes_from_panelapp(source_db)
  # remove disputed genes from clingen
  fixed_source_db<-remove_disputed_genes_from_clingen(fixed_source_db)
  # fix "bad" genes
  fixed_source_db<-fix_source_db(fixed_source_db)
  fixed_source_db
})

# summarize the source_db data into panels
get_panels_list_from_source_db <- function(source_db) {
  panels_list<-source_db%>%
    group_by(source,panel_id,panel_name)%>%
    summarize(number_of_genes=n(),
              high_conf_genes=sum(confidence_level==3),
              medium_conf_genes=sum(confidence_level==2),
              low_conf_genes=sum(confidence_level==1),
              genes=paste0(gene_symbol,collapse=', '))%>%
    mutate(number_of_genes_cat=cut(number_of_genes,
                                   breaks=c(0,1,10,50,200,500,1000,100000),
                                   labels=c('1','2-10','11-50','51-200','201-500','501-1000','>1000')),
           panel_joined_name=glue('{source}|{panel_id}|{panel_name}'))
  return(panels_list)
}

gtr_db_file<-'./data/databases/gtr/gtr_panels_db_2022-02-13.csv.gz'
gtr_db<-readr::read_delim(gtr_db_file)
gtr_db<-gtr_db%>%mutate(panel_joined_name=glue('{source}|{panel_id}|{panel_name}'))
gtr_panels_list_file<-'./data/databases/gtr/gtr_panels_list_2022-02-13.csv.gz'
gtr_panels_list<-readr::read_delim(gtr_panels_list_file)

panels_list<-get_panels_list_from_source_db(source_db)