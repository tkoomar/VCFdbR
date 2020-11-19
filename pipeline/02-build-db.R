#! Rscript
args = commandArgs(trailingOnly=TRUE)

## defaults ##
run_parallel <- FALSE

while(length(args > 0) ){
  if(args[1] == '--mode'){
    if(args[2] == 'file'){
      file_mode = TRUE
    } else if(args[2] == 'table'){
      file_mode = FALSE
    } else {
        stop("`--mode` argument requires one of two options: 'file' or 'table'")
      }
    args <- args[-1:-2]
  } else if(args[1] == '--prefix'){
    prefix <- args[2]
    args <- args[-1:-2]
    message("File prefix: ", prefix)
  } else if(args[1] == "--vcf"){
    vcf_name <- args[2]
    args <- args[-1:-2]
    message("Input VCF: ", vcf_name)
  } else if(args[1] == "--chunks-ranges"){
    ranges_name <- args[2]
    args <- args[-1:-2]
    message("Input variant chunk ranges: ", ranges_name)
  } else if(args[1] == "--threads"){
    run_parallel <- TRUE
    threads <- args[2]
    args <- args[-1:-2]
    message("Writing genotypes in parallel (requires furrr package)")
  } else {
    stop("Unknown argument: ", args[1])
  }
}

if(!exists("prefix")){
  stop("A prefix (name) for the database must be defined with the `--prefix` argument")
}
if(!exists("vcf_name")){
  stop("A vcf must be passed with the `--vcf`` argument")
}
if(!exists("ranges_name")){
  stop("A GRanges corresponding to chunks of variants to process must be passed with the `--chunks-ranges` argument")
}
if(exists('threads') ){
  if( is.na(as.integer(threads))){
    stop("If provided, `--threads` must be an integer!")
  }
}

if(!exists('file_mode')){
  stop("You must specify a `--mode` for how genotypes should be stored: either 'file' or 'table'")
} else if(file_mode){
  ## Make the variant genotype directory ##
  geno_dir <- paste0(getwd(), "/", prefix, "-genos")
  dir.create(geno_dir)
  message("Building a database where genotypes are saved to a directory (", geno_dir, "/)")
} else if(!file_mode){
  message("Building a database where genotypes are saved to a table (variant_geno)")
}

require(VariantAnnotation)
require(dbplyr)
require(tidyverse)
require(magrittr)
require(progress)
require(RSQLite)

#### A function to convert genotypes to dosage ####
gt2snp <- function(x){
  case_when(
    str_detect(x, "0\\/\\.|\\.\\/0|0\\/0|0\\|0|^0$")  ~ 0, 
    str_detect(x, c("0\\/1|1\\/0|0\\|1|1\\|0|\\.\\/1|1\\/\\.|\\.\\|1|1\\|\\.|^1$"))  ~ 1, 
    str_detect(x,  c("1\\/1|1\\|1")) ~ 2, 
    TRUE ~ as.numeric(NA)
  )
}

if(!exists('chunk_ranges')){
  chunk_ranges <- read_rds(ranges_name)
}

db_name <- paste0(prefix, ".db")
#### Insert VCF header information ####
message('######\nPARSING VCF HEADER\n######')

vcf_header <- scanVcfHeader(vcf_name)

header <- vcf_header@header
names(header)

lapply(names(header), function(name){
  X <- header[[name]]
  con <- dbConnect(SQLite(), db_name)
  db_insert_into(con = con, 
                 table = name, 
                 as_tibble(X, rownames = 'name'))
  dbDisconnect(con)
})

## add samples ##
con <- dbConnect(SQLite(), db_name)
db_insert_into(con = con, 
               table = 'samples', 
               enframe(vcf_header@samples))
dbDisconnect(con)

#### Check that the geno fields exists ####
message('######\nCHECKING GENOTYPE FIELDS\n')
tmp.vcf <- readVcf(vcf_name, param = chunk_ranges[1])

geno_col_names <- geno(tmp.vcf) %>% names()

## Check to make sure the geno fields declaierd in the header are actually on the data
## Also check that they are of a supported type (SQLite can't handle list columns)
keep_geno_col <- sapply(geno_col_names, function(col_name){
  all_missing <- geno(tmp.vcf)[[col_name]] %>%
    unlist() %>%
    is.na() %>%
    all()
  
  if (all_missing){
    warning("geno column '", col_name, "' appears to be missing. It will be skipped.")
    return(FALSE)
    }
  
  col_class <- type(geno(tmp.vcf)[[col_name]])
  
  if(!(col_class %in% c("character", "integer", "numeric", "double"))){
    warning("geno column '", col_name, "'is of upsupported type '", col_class, ". It will be skipped.")
    return(FALSE)
  }
  
  return(TRUE)
})

geno_col_names <- geno_col_names[keep_geno_col]

message("######\n")


#### Check to see if CSQ exists ####
if(class(tmp.vcf@info$CSQ)!="NULL"){
  csq_exists = TRUE
  csq_cols <- vcf_header@header$INFO$Description[rownames(vcf_header@header$INFO) == 'CSQ'] %>%
    tolower() %>% str_replace_all(" |:|\\.", "_")
} else{
  csq_exists = FALSE
  message("######\nNo VEP annotations detected.\n######")
}
## These are the impacts flagged as 'exonic' from VEP
exonic_impacts = c("stop_gained",
                   "exon_variant",
                   "stop_lost",
                   "frameshift_variant",
                   "initiator_codon_variant",
                   "inframe_deletion",
                   "inframe_insertion",
                   "missense_variant",
                   "protein_altering_variant",
                   "incomplete_terminal_codon_variant",
                   "stop_retained_variant",
                   "5_prime_UTR_premature_start_codon_variant",
                   "synonymous_variant",
                   "coding_sequence_variant",
                   "5_prime_UTR_variant",
                   "3_prime_UTR_variant",
                   "transcript_ablation",
                   "transcript_amplification",
                   "feature_elongation",
                   "feature_truncation")
rm(tmp.vcf)

if(run_parallel){
  require(furrr)
  options(future.globals.maxSize= 10000*1024^2)
  plan(multiprocess)
  options(mc.cores = threads)
}


#### Build Tables ####
message("######\nSTARTING TO BUILD DATABASE\n######")
date()
p <- length(chunk_ranges)
options(progress_enabled = TRUE)
pb <- progress_bar$new(total = p, 
                       clear = FALSE, 
                       force = TRUE,
                       format = ":current/:total chunks completed in :elapsed; eta: :eta"
                       )
index_start <- 0
for(i in 1:p){
  pb$tick()
  
  current_ranges <- chunk_ranges[i]
  
  ## set up looping variables
  .vcf <- readVcf(vcf_name, param = current_ranges)
  
  ## generate indicies for our current variants
  var_ind <- (index_start+1):(index_start + length(.vcf))
  
  index_start <- last(var_ind)
  
  ## VEP CSQ goes in the variant_impact table
  if(csq_exists){
    csq.vcf <- .vcf@info %>%
      as_tibble() %>%
      select(CSQ) %>%
      add_column(variant_id = var_ind) %>%
      unnest(CSQ) %>%
      separate(CSQ, 
               sep = "\\|", 
               into = str_split(csq_cols, pattern = "\\|")[[1]]
      ) %>% 
      separate_rows(consequence,  sep = "&") %>%
      mutate(is_lof = impact == "HIGH" & biotype == "protein_coding", 
             is_splicing = str_detect(consequence, "splice"), 
             is_exonic = biotype %in% exonic_impacts,
             is_intronic = intron != "") %>%
      select(-ends_with("_af"), -any_of(c('clin_sig', 'pheno','somatic', 'pubmed', 'consequence_annotations_from_ensembl_vep__format__allele')))
    
    ## get the clinvar stuff into a better format from the info annotation, rather than VEP's
    clinvar.vcf <- .vcf@info %>%
      as_tibble() %>%
      select(any_of(c("clinvar_sig", "clinvar_disease_name"))) %>%
      as_tibble() %>%
      bind_cols(tibble(variant_id = var_ind))
    
    csq.vcf <- csq.vcf %>%
      {if('clinvar_sig' %in% colnames(clinvar.vcf)){
        left_join(., clinvar.vcf %>% select(variant_id, clinvar_sig) %>% unnest(clinvar_sig)) 
      } else {.}} %>%
      {if('clinvar_disease_name' %in% colnames(clinvar.vcf)){
        left_join(., clinvar.vcf %>% select(variant_id, clinvar_disease_name) %>% unnest(clinvar_disease_name)) 
      } else {.}} %>%
      select(variant_id, everything())
  
  ## other info goes into into the variant_info table
    info.vcf <- tibble(variant_id = var_ind) %>%
      bind_cols(.vcf@rowRanges %>% as_tibble () %>% select(seqnames, start, end) %>% rename('chr' = seqnames)) %>%
      bind_cols(.vcf@fixed %>% as_tibble()) %>%
      bind_cols(.vcf@info %>% as_tibble() %>% select(-any_of(c('CSQ', 'clinvar_sig', 'clinvar_disease_name'))))
    
  } else{
    info.vcf <- tibble(variant_id = var_ind) %>%
      bind_cols(.vcf@rowRanges %>% 
                  as_tibble() %>% 
                  select(seqnames, start, end) %>% 
                  rename('chr' = seqnames)) %>%
      bind_cols(.vcf@fixed %>% as_tibble(), ) %>%
      bind_cols(.vcf@info %>% as_tibble()) 
  }
  
  names(info.vcf) %<>% tolower()
  
  info.vcf <- info.vcf[,!duplicated(colnames(info.vcf))]
  
  ## fix a few weird formatting things on the columns
  ## then add filepaths for variants
  info.vcf <- info.vcf %>%
    mutate(alt = sapply(alt, as.character)) %>%
    mutate_if(function(x){class(x)=="AsIs"}, as.character) 
  
  ## genotypes have to be combined to go into the genotype table
  if (length(geno_col_names) > 0){
    
    .geno_col <- geno_col_names[1]
    geno_col <- enquo(.geno_col)
    geno.vcf <- tibble(group = var_ind, variant_id = var_ind) %>%
      bind_cols(geno(.vcf)[[.geno_col]] %>% as_tibble()) %>%
      gather(sample, !!.geno_col, -variant_id, -group)
    
    for(.geno_col in geno_col_names[-1]){
      message(.geno_col) ## debug
      geno_col <- enquo(.geno_col)
      geno.vcf <- geno.vcf %>%
        bind_cols(
          geno(.vcf)[[.geno_col]] %>% 
            as_tibble() %>%
            gather(sample, !!.geno_col) %>%
            select(-sample) 
        )
    }
    names(geno.vcf) %<>% tolower()
    
    if('GT' %in% geno_col_names){
      geno.vcf <- geno.vcf %>%
        mutate(gt_raw = gt, 
               gt = gt2snp(gt_raw))
    }
    
    
    if(file_mode){
      ## add path to genos to info table
      info.vcf <- info.vcf %>% 
        add_column(geno = str_c(geno_dir, "/", var_ind, ".rds"))
      
      ## Write geno files ##
      geno.vcf <- geno.vcf %>%
        group_by(group) %>%
        nest() 
      
      if(run_parallel){
        tmp <- quietly(future_map2(geno.vcf$data, geno.vcf$group, 
                                   ~write_rds(as.data.frame(.x), str_c(geno_dir, "/", .y, ".rds"))
        )
        )
      } else{
        tmp <- quietly(map2(geno.vcf$data, geno.vcf$group, 
                            ~write_rds(as.data.frame(.x), str_c(geno_dir, "/", .y, ".rds"))
        )
        )
      }
    } else {
      ## write geno to SQLite
      geno.vcf <- geno.vcf %>%
        select(-group) %>%
        arrange(variant_id, sample)
      con <- dbConnect(SQLite(), db_name)
      db_insert_into(con = con, 
                     table = "variant_geno", 
                     geno.vcf) 
      dbDisconnect(con)
    }
    }
  
  ## Write to database
  con <- dbConnect(SQLite(), db_name)
  if(exists('csq.vcf')){
    db_insert_into(con = con, 
                   table = "variant_impact", 
                   csq.vcf) 
  }
  db_insert_into(con = con,
                 table = "variant_info",
                 info.vcf)
  dbDisconnect(con)
}
## end loop
message("######\nDone inserting variants\n#####")
date()
closeAllConnections()

### Make geno dir read-only
if(file_mode){Sys.chmod(geno_dir, mode = "555")}
