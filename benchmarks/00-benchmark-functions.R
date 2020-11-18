require(tidyverse)
require(dbplyr)
require(DBI)
require(RSQLite)
require(pbapply)
require(parallel)
require(progress)

## Variant  pulling by id function ##
pull_vars_by_id <- function(.db_name, .var_ids, .cores = NA){
  pboptions(type = 'none')
  
  ## check if we have a file database
  if(str_remove(.db_name, ".db") %>% str_c("-genos") %>% dir.exists()){
    file_db <- TRUE
    var_list <- str_c(str_remove(.db_name, ".db") %>% str_c("-genos"), "/", .var_ids, ".rds")
  } else{
    file_db <- FALSE
    var_list = .var_ids
  }
  
  ## split list if multithreaded
  if(!is.na(.cores)){
    var_list <- split(var_list, cut_number(seq_along(var_list), .cores, labels = FALSE))
    cl <- makeForkCluster(.cores)
  } else { cl = NULL  }
  
  if(file_db & is.na(.cores)){
    out <- pblapply(X = var_list, FUN = read_rds) %>%
      bind_rows()
  }  else if(file_db & !is.na(.cores)){
    out <- pblapply(X = var_list, cl = cl, FUN = function(.x){
      map(.x, read_rds) %>%
        bind_rows()
    }) %>% bind_rows()
  } else if(!file_db & !is.na(.cores)) {
    out <- pblapply(X = var_list, cl = cl, FUN = function(ids){
      con <- DBI::dbConnect(RSQLite::SQLite(), .db_name)
      tbl(con, 'variant_geno') %>%
        filter(variant_id %in% ids) %>%
        collect() %>%
        bind_rows()
      DBI::dbDisconnect(con)
    })
  } else {
    con <- DBI::dbConnect(RSQLite::SQLite(), .db_name)
    out <- tbl(con, 'variant_geno') %>%
      filter(variant_id %in% var_list) %>%
      collect() 
    DBI::dbDisconnect(con)
  }
  
  if(!is.na(.cores)){closeAllConnections()}
  return(out)
}

## sampling and timing function
time_var_pull <- function(.db_name, .n,  .cores = NA, .max_vars = 1417043){
  var_ids <- sample(1:.max_vars, .n) %>% sort()
  system.time({
    pull_vars_by_id(.db_name = .db_name, .var_ids = var_ids, .cores = .cores)
  })[3]
}


## get the number of matching variants in each gene ##

filter_test <- function(.symbol, .db_name, .pb = NULL, .af = 0.01, .return_dat = FALSE, .return_con = NA){
  if(!is.null(.pb)){.pb$tick()}
  
  if(!is.na(.return_con)){
    .return_dat = TRUE
    con = .return_con
  } else{
    con <- DBI::dbConnect(RSQLite::SQLite(), .db_name)
    }
  
  filter_time <- system.time({
    tmp_dat <- tbl(con, "variant_impact") %>%
      filter(symbol == .symbol) %>%
      select(variant_id, symbol) %>%
      distinct() %>%
      inner_join(
        tbl(con, "variant_info") %>% 
          select(any_of(c('variant_id', 'geno', 'af'))) %>%
          filter(af < .af), 
        by = "variant_id"
      ) 
    
    if(is.na(.return_con)){
      tmp_dat <- collect(tmp_dat)
      DBI::dbDisconnect(con)
      n_vars = nrow(tmp_dat)
    }
  })  
  
  if(.return_dat){
    return(tmp_dat)
  } else{
    return(tibble(gene = .symbol, filter_time = filter_time[3], n_vars = n_vars))
  }
}

## Function to pull variants ##
pull_geno_test <- function(.symbol, .db_name, .cl = NA, .pb = NULL){
  if(!is.null(.pb)){.pb$tick()}
  
  info_colnames <- DBI::dbConnect(RSQLite::SQLite(), .db_name) %>%
    tbl("variant_info") %>%
    colnames()
  
  tmp_dat <- filter_test(.symbol = .symbol, .db_name = .db_name, .return_dat = TRUE)
  
  pull_time <- system.time({
    if(!is.na(.cl)){
      cl <- makeForkCluster(.cl)
    } else {
      cl <- NULL
    }
    
    if('geno' %in% info_colnames){
        tmp_dat %>%
          mutate(gt = pblapply(X = geno, FUN = read_rds, cl = cl))
      if(!is.na(.cl)){stopCluster(cl)}
      } else if(!is.na(.cl)){
        tmp_dat %>%
          mutate(geno = pblapply(cl = cl, X = tmp_dat$variant_id, FUN = function(x){
              con <- DBI::dbConnect(RSQLite::SQLite(), .db_name)
              tbl(con, 'variant_geno') %>%
                filter(variant_id == x) %>%
                collect()
              DBI::dbDisconnect(con)
            })
          )
        if(!is.na(.cl)){stopCluster(cl)}
      } else {
            con <- DBI::dbConnect(RSQLite::SQLite(), .db_name)
             tbl(con, 'variant_geno') %>%
                filter(variant_id %in% !!tmp_dat$variant_id) %>%
               collect()
          }
  })
  
  
  
  return(pull_time[3])
}


get_genes <- function(.db_name){
  con <- DBI::dbConnect(RSQLite::SQLite(), .db_name)
  
  ## get list genes ####
  genes <- tbl(con, "gene_map") %>%
    filter(symbol_source == "EntrezGene" & feature_type == "Transcript") %>%
    select(symbol) %>%
    distinct() %>%
    select(symbol) %>%
    filter(!is.na(symbol)) %>%
    collect()
  
  ## Done with DB####
  DBI::dbDisconnect(con)
  
  pb <- progress_bar$new(total = nrow(genes),
                         force = TRUE,
                         format = ":percent complete, estimated time remaining: :eta")
  
  
  genes <- genes %>%
    mutate(tmp = map(symbol, filter_test, .pb = pb, .db_name = .db_name)) %>%
    unnest(tmp)
  
  closeAllConnections()
  
  genes <- genes %>%
    select(-gene) %>%
    mutate(bin = cut_width(n_vars,                         
                           width = 500, 
                           center = 250, 
                           dig.lab = 4)
    )
}


pull_variants <- function(db_name, gene_file, sample_size, cores = NA){
  benchmark_dat <- read_csv(gene_file)
  
  ## DB and sample map
  dbs <- tribble(~ db_name, ~n_samples,
                 db_name, sample_size
  )
  
  
  pboptions(type = 'none')
  
  #### run, randomizing the order first ####
  benchmark_dat <- crossing(
    cl = cores, 
    dbs,
    benchmark_dat
  ) %>%
    sample_frac(1L) 
  
  pb <- progress_bar$new(total = nrow(benchmark_dat), 
                         force = TRUE,
                         format = ":percent complete, estimated time remaining: :eta")
  
  out <- benchmark_dat %>%
    mutate(pull_time = 
             pmap_dbl(.l = list(symbol, db_name, cl), 
                      .f = pull_geno_test, 
                      .pb = pb)
    ) %>%
    mutate(cl = replace_na(cl, 1))
  
  return(out)
}