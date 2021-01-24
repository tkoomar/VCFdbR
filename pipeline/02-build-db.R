#! Rscript
args = commandArgs(trailingOnly=TRUE)

## defaults ##
run_parallel <- FALSE
multi_gt <- FALSE
end_provided <- FALSE
debug_mode <- FALSE

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
  } else if(args[1] == "--end-chunk"){
    end_provided <- TRUE
    p <- as.integer(args[2])
    args <- args[-1:-2]
    message("Input variant chunk ranges: ", ranges_name)
  } else if(args[1] == "--include-multivalue-gt"){
    multi_gt <- TRUE
    args <- args[-1]
    message("Including GENO Fields with multiple values (may be slow for many samples)")
  } else if(args[1] == "--debug"){
    debug_mode <- TRUE
    args <- args[-1]
    message("Debug mode enabled, producing verbose output")
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

suppressPackageStartupMessages(require(VariantAnnotation))
suppressPackageStartupMessages(require(tidyverse))
suppressPackageStartupMessages(require(dbplyr))
suppressPackageStartupMessages(require(magrittr))
suppressPackageStartupMessages(require(progress))
suppressPackageStartupMessages(require(RSQLite))


if(!exists('chunk_ranges')){
  chunk_ranges <- read_rds(ranges_name)
}

db_name <- paste0(prefix, ".db")

#### If a progress file is detected, resume from it ####
if(!file.exists(paste0(prefix, ".progress.RData"))){
  
  #### A function to convert genotypes to dosage 
  gt2snp <- function(x){
    case_when(
      str_detect(x, "0\\/\\.|\\.\\/0|0\\/0|0\\|0|^0$")  ~ 0, 
      str_detect(x, c("0\\/1|1\\/0|0\\|1|1\\|0|\\.\\/1|1\\/\\.|\\.\\|1|1\\|\\.|^1$"))  ~ 1, 
      str_detect(x,  c("1\\/1|1\\|1")) ~ 2, 
      TRUE ~ as.numeric(NA)
    )
  }
  
  #### Insert VCF header information ####
  message('######\nPARSING VCF HEADER\n######')
  
  vcf_header <- scanVcfHeader(vcf_name)
  
  header <- vcf_header@header
  names(header)
  
  lapply(names(header), function(name){
    X <- header[[name]]
    con <- dbConnect(SQLite(), db_name)
    DBI::dbWriteTable(
      conn = con,
      name = name,
      value = as_tibble(X, rownames = 'name'), 
      append = TRUE)
    DBI::dbDisconnect(con)
  })
  
  ## add samples ##
  con <- dbConnect(SQLite(), db_name)
  DBI::dbWriteTable(
    conn = con,
    name = 'samples',
    value = enframe(vcf_header@samples), 
    append = TRUE)
  DBI::dbDisconnect(con)
  
  #### Check that the geno fields exists ####
  message('######\nCHECKING GENOTYPE FIELDS\n')
  tmp.vcf <- readVcf(vcf_name, param = chunk_ranges[1])
  
  geno_col_names <- geno(tmp.vcf) %>% names()
  
  ## Check to make sure the geno fields declared in the header are actually on the data
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
    
    col_class <- class(geno(tmp.vcf)[[col_name]])
    col_type <- type(geno(tmp.vcf)[[col_name]])
    
    
    if(col_class == "array" & multi_gt){
      message("geno column '", col_name, "' has multiple values per variant (", dim(geno(tmp.vcf)[[col_name]])[3], "). The package `reshape2` is required to parse.")
      suppressPackageStartupMessages(require(reshape2))
    } else if(col_class == "array" & !multi_gt){
      message("geno column '", col_name, "' has multiple values per variant and will be skipped.")
      return(FALSE)
    }
    
    if(!(col_type %in% c("character", "integer", "numeric", "double"))){
      warning("geno column '", col_name, "'is of upsupported type '", col_type, ". It will be skipped.")
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
  
  index_start <- 0
  chunk_start <- 0
} else {
  load(paste0(prefix, ".progress.RData"))
}

if(run_parallel){
  require(furrr)
  options(future.globals.maxSize= 10000*1024^2)
  options(mc.cores = threads)
}


#### Build Tables ####
message("######\nSTARTING TO BUILD DATABASE\n######")
date()

if(!end_provided){
  p <- length(chunk_ranges)
}

pb_total <- p - chunk_start

options(progress_enabled = TRUE)
pb <- progress_bar$new(total = pb_total, 
                       clear = FALSE, 
                       force = TRUE,
                       format = "Starting chunk :current/:total after :elapsed; eta: :eta"
)

chunk_start <- chunk_start + 1
if(p < chunk_start){
  stop("End chunk must be after starting chunk!")
}

for(i in seq(chunk_start,p)){
  pb$tick()
  
  current_ranges <- chunk_ranges[i]
  
  ## set up looping variables
  if(debug_mode){start <- Sys.time()}
  .vcf <- readVcf(vcf_name, param = current_ranges)
  
  ## generate indicies for our current variants
  var_ind <- (index_start+1):(index_start + length(.vcf))
  
  if(debug_mode){
    end <- Sys.time() 
    message("Done Reading Chunk from VCF:") 
    print(end - start)
    start <- end
  }
  
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
        left_join(., clinvar.vcf %>% 
                    select(variant_id, clinvar_sig) %>% 
                    unnest(clinvar_sig), 
                  by = "variant_id") 
      } else {.}} %>%
      {if('clinvar_disease_name' %in% colnames(clinvar.vcf)){
        left_join(., clinvar.vcf %>% 
                    select(variant_id, clinvar_disease_name) %>% 
                    unnest(clinvar_disease_name), 
                  by = "variant_id") 
      } else {.}} %>%
      select(variant_id, everything())
    
    if(debug_mode){
      end <- Sys.time() 
      message("Done Parsing Variant Impact:")
      print(end - start)
      start <- end
    }
    
    ## other info goes into into the variant_info table
    info.vcf <- tibble(variant_id = var_ind) %>%
      bind_cols(.vcf@rowRanges %>% as_tibble () %>% select(seqnames, start, end) %>% rename('chr' = seqnames)) %>%
      bind_cols(.vcf@fixed %>% as_tibble()) %>%
      bind_cols(.vcf@info %>% as_tibble() %>% select(-any_of(c('CSQ', 'clinvar_sig', 'clinvar_disease_name'))))
    
    if(debug_mode){
      end <- Sys.time()
      message("Done Parsing Variant Info:")
      print(end - start)
      start <- end
    }
    
  } else{
    info.vcf <- tibble(variant_id = var_ind) %>%
      bind_cols(.vcf@rowRanges %>% 
                  as_tibble() %>% 
                  select(seqnames, start, end) %>% 
                  rename('chr' = seqnames)) %>%
      bind_cols(.vcf@fixed %>% as_tibble()) %>%
      bind_cols(.vcf@info %>% as_tibble()) 
    
    if(debug_mode){
      end <- Sys.time() 
      message("Done Parsing Variant Info:") 
      print(end - start)  
      start <- end
    }
    
  }
  
  names(info.vcf) %<>% tolower()
  
  info.vcf <- info.vcf[,!duplicated(colnames(info.vcf))]
  
  ## fix a few weird formatting things on the columns
  ## then add filepaths for variants
  
  if(debug_mode){message("fixing 'AsIs' INFO columns")}
  
  info.vcf <- info.vcf %>%
    mutate(alt = sapply(alt, as.character)) %>%
    mutate_if(function(x){class(x)=="AsIs"}, as.character) 
  
  if(debug_mode){
    message("done parsing info")
    end <- Sys.time()
    message("Done Fixing Variant Info:")
    print(end - start)
    start <- end
  }
  
  ## genotypes have to be combined to go into the genotype table
  if (length(geno_col_names) > 0){
    
    .geno_col <- geno_col_names[1]
    geno_col <- enquo(.geno_col)
    geno.vcf <- tibble(group = var_ind, variant_id = var_ind) %>%
      bind_cols(geno(.vcf)[[.geno_col]] %>% as_tibble()) %>%
      gather(sample, !!.geno_col, -variant_id, -group)
    
    if(debug_mode){
      end <- Sys.time()
      message("Done Subsetting Genos:")
      print(end - start)  
      start <- end  
    }
    
    for(.geno_col in geno_col_names[-1]){
      message(.geno_col)
      geno_col <- enquo(.geno_col)
      
      if (class(geno(.vcf)[[.geno_col]]) == "matrix"){
        geno.vcf <- geno.vcf %>%
          bind_cols(
            geno(.vcf)[[.geno_col]] %>% 
              as_tibble() %>%
              gather(sample, !!.geno_col) %>%
              select(-sample) 
          )
      } else if (class(geno(.vcf)[[.geno_col]]) == "array"){
        tmp <- geno(.vcf)[[.geno_col]] 
        rownames(tmp) <- var_ind
        tmp <- tmp %>%
          reshape2::melt() %>%
          select(variant_id = Var1, 
                 sample = Var2, 
                 Var3, 
                 value) %>%
          mutate(Var3 = str_c(.geno_col, "_", Var3)) %>%
          spread(Var3, value)
        
        geno.vcf <- geno.vcf %>%
          left_join(tmp, by = c("variant_id", "sample"))
      }
    }
    names(geno.vcf) %<>% tolower()
    
    if(debug_mode){
      end <- Sys.time()
      message("Done Parsing Geno Fields:")
      print(end - start)
      start <- end
    }
    
    if('GT' %in% geno_col_names){
      geno.vcf <- geno.vcf %>%
        mutate(gt_raw = gt, 
               gt = gt2snp(gt_raw))
    }
    
    if(debug_mode){
      end <- Sys.time()
      message("Done Converting GT to Numeric:")
      print(end - start)
      start <- end
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
        plan(multiprocess)
        tmp <- quietly(future_map2(geno.vcf$data, geno.vcf$group, 
                                   ~write_rds(as.data.frame(.x), str_c(geno_dir, "/", .y, ".rds"))
        )
        )
        plan(sequential)
      } else{
        tmp <- quietly(
          map2(geno.vcf$data, 
               geno.vcf$group, 
               function(.x, .y){
                 write_rds(as.data.frame(.x), str_c(geno_dir, "/", .y, ".rds"))
                 
               }
          )
        )
      }
    } else {
      ## write geno to SQLite
      geno.vcf <- geno.vcf %>%
        select(-group) %>%
        arrange(variant_id, sample)
      con <- dbConnect(SQLite(), db_name)
      DBI::dbWriteTable(
        conn = con, 
        name = "variant_geno", 
        value = geno.vcf, 
        append = TRUE) 
      DBI::dbDisconnect(con)
    }
  }
  
  if(debug_mode){
    end <- Sys.time()
    message("Done Writing Genotypes:")
    print(end - start)
    start <- end
  }
  
  ## Write to database
  con <- dbConnect(SQLite(), db_name)
  if(exists('csq.vcf')){
    DBI::dbWriteTable(
      conn = con, 
      name = "variant_impact",
      value = csq.vcf, 
      append = TRUE) 
  }
  DBI::dbWriteTable(
    conn = con,
    name = "variant_info",
    value = info.vcf, 
    append = TRUE)
  dbDisconnect(con)
  
  if(debug_mode){
  end <- Sys.time() 
  message("Done Inserting Annotations in DB:") 
  print(end - start)
  start <- end
  }
  
  if(debug_mode){message("Done writing/insering variants")}
  
  rm(.vcf, info.vcf, geno.vcf)
  if(exists('csq.vcf')){
    rm(csq.vcf)
  }
  gc(verbose = FALSE)
  
  if(debug_mode){
    end <- Sys.time() 
    message("Done Collecting Garbage:")
    print(end - start)
    start <- end  
  }
  
  chunk_start <- i
  index_start <- last(var_ind)
  save(list = c(
    'index_start', 'chunk_start', 'exonic_impacts', 'geno_col_names', 'gt2snp', 
    'header', 'keep_geno_col', 'ranges_name', 'vcf_header', 'vcf_name',
    if(csq_exists){c('csq_exists', 'csq_cols')}
  ),
       file = paste0(prefix, ".progress.RData"))
}

#### end loop ####
message("######\nDone inserting variants\n#####")
date()
closeAllConnections()

### Make geno dir read-only
# if(file_mode){Sys.chmod(geno_dir, mode = "555")}
