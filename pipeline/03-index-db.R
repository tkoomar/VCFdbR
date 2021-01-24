#! Rscript
args = commandArgs(trailingOnly=TRUE)

while(length(args > 0) ){
  if(args[1] == '--prefix'){
    prefix <- args[2]
    args <- args[-1:-2]
    message("File prefix: ", prefix)
  } else if(args[1] == "--db"){
    db_name <- args[2]
    args <- args[-1:-2]
    message("Input DB: ", db_name)
  } else {
    stop("Unknown argument: ", args[1])
  }
}

if(!exists("prefix")){
  stop("A prefix (name) for the database must be defined with the `--prefix` argument")
}
if(!exists("db_name")){
  stop("A database must be passed with the `--db` argument")
}

suppressPackageStartupMessages(require(VariantAnnotation))
suppressPackageStartupMessages(require(dbplyr))
suppressPackageStartupMessages(require(tidyverse))
suppressPackageStartupMessages(require(magrittr))
suppressPackageStartupMessages(require(DBI))
suppressPackageStartupMessages(require(RSQLite))


## build indicies
message("######\nBUILDING INDICIES\n######")
con <- dbConnect(SQLite(), paste0(db_name))

dbExecute(con, "CREATE INDEX idx_info_variant_id ON variant_info (variant_id)")
message("### done with variant ids on variant_info")

if ('variant_geno' %in% DBI::dbListTables(con)){
  dbExecute(con, "CREATE INDEX idx_geno_variant_id ON variant_geno (variant_id)")
  message("### done with variant ids on variant_geno")
}

if ('variant_impact' %in% DBI::dbListTables(con)){
  dbExecute(con, "CREATE INDEX idx_impact_variant_id ON variant_impact (variant_id)")
  message("### done with variant ids on variant impacts")
  dbExecute(con, "CREATE INDEX idx_impact_symbol ON variant_impact (symbol)")
  message("### done with gene symbols")
  dbExecute(con, "CREATE INDEX idx_impact_gene on variant_impact (gene)")
  message("### done with gene IDs")
  
  message("## Indexing Variant Effects ##")
  dbExecute(con, "CREATE INDEX idx_impact_consequence ON variant_impact (consequence)")
  dbExecute(con, "CREATE INDEX idx_impact_lof ON variant_impact (is_lof)")
  dbExecute(con, "CREATE INDEX idx_impact_exonic ON variant_impact (is_exonic)")
  dbExecute(con, "CREATE INDEX idx_impact_splicing ON variant_impact (is_splicing)")
  dbExecute(con, "CREATE INDEX idx_impact_biotype ON variant_impact (biotype)")
}


message("## Indexing Allele Frequencies ##")
if('af' %in% DBI::dbListFields(con, "variant_info")){dbExecute(con, "CREATE INDEX idx_info_aaf ON variant_info (af)")}
if('an' %in% DBI::dbListFields(con, "variant_info")){dbExecute(con, "CREATE INDEX idx_info_an ON variant_info (an)")}
if('ac' %in% DBI::dbListFields(con, "variant_info")){dbExecute(con, "CREATE INDEX idx_info_ac ON variant_info (ac)")}

dbDisconnect(con)
message("######\nDONE INDEXING\n######")
date()


#### Make Variant Ranges ###
message("######\nMAKING VARIANT RANGES OBJECT\n######")
con <- dbConnect(SQLite(), db_name)
var_info <- tbl(con, 'variant_info')
#var_info <- info.vcf
var_ranges <- tbl(con, 'variant_info') %>%
  select(variant_id, chr, start, end, ref, alt, any_of('geno')) %>%
  collect() %>%
  makeGRangesFromDataFrame( keep.extra.columns = TRUE)

saveRDS(var_ranges, file = paste0(prefix, "-variant-ranges.rds"))


#### Make Gene and Transcript Mapping Table ####
if ('variant_impact' %in% DBI::dbListTables(con)){
  message("######\nMAKING GENE AND TRANSCRIPT MAPPING TABLE\n######")
  var_imp <- tbl(con, 'variant_impact')
  
  gene_dat <- var_imp %>% 
    select(any_of(c('symbol', 'symbol_source', 'gene', 'source', 'feature', 'canonical', 'ensp', 'ccds', 'motif_name', 'feature_type'))) %>% 
    collect() %>% 
    distinct() 
  
  gene_dat[gene_dat == ""] <- NA
  
  DBI::dbWriteTable(
    conn = con,
    name = "gene_map", 
    value = gene_dat)
}
dbDisconnect(con)
message("######\nVCFdb creation ended on\n", date(), "\n######")
