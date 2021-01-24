#! Rscript
args = commandArgs(trailingOnly=TRUE)

file_mode <- NA

while(length(args > 0) ){
  if(args[1] == "--help" | args[1] == "-h" ){
    message(
"This script builds a SQLite representation of a 
normalized VCF in the current working directory. 
** THIS SCRIPT SHOULD NOT BE MOVED **
It relies on other scripts that are located in 
the 'pipeline/' directory. 

Required Arguments:
--prefix [character] 
    Will be the name of the database and other 
    files produced by the pipeline. 
--vcf [character]
    The name of the input VCF. Should not have
    multialleleic sites. If VEP annotatios are
    present, should be named the default 'CSQ'
--mode ['file'|'table']
    In 'table' mode, genotypes are stored in a 
    SQLite table called `variant_geno`. This 
    makes it possible to quickly filter based on 
    genotypes other individual-level metric.
    In 'file' mode, genotypes are saved as 
    separate .rds files in a directory named
    [prefix]-geno. This makes it possible to 
    store data for very large cohorts that would
    violate your filesystem (or SQLite)'s limits
    on the size of individual files. 

Optional Arguments: 
--chunk_size [integer]
    The approximate number of variants to process
    at once. Smaller numbers  use less memory at
    the expense of running slightly slower.
--threads [integer]
    Number of threads to use when writing the
    individual variants to files. Setting this
    high will decrease performace dramatically
    once disk or network speed is saturated.
    Should never be set to more threads than 
    what is available on your machine. 
    Anecdotal testing finds the optimal number
    to be between 2 and 10, depending on many
    factors. Requires the `furrr` package.
--include-multivalue-gt 
    Include GT fields which have multiple values
    for each sample. Can be slow for VCFs with
    many samples, due to the way such fields are
    handeld by VariantAnnotation::readVcf(). If 
    preserving these fields is important, it may
    be worth first separating them in the input 
    VCF.")
    stop()
  } else if(args[1] == '--mode'){
    if(args[2] == 'file'){
      file_mode = TRUE
    } else if(args[2] == 'table'){
      file_mode = FALSE
    } else {
      stop("`--mode` argument requires one of two options: 'file' or 'table'")
    }
    args <- args[-1:-2]
  } else  if(args[1] == '--prefix'){
    prefix <- args[2]
    args <- args[-1:-2]
  } else if(args[1] == "--vcf"){
    vcf_name <- args[2]
    args <- args[-1:-2]
  } else if(args[1] == "--chunk-size"){
    chunk_size <- args[2]
    args <- args[-1:-2]
  } else if(args[1] == "--include-multivalue-gt"){
    multi_gt <- TRUE
    args <- args[-1]
    message("Including GENO Fields with multiple values (may be slow for many samples)")
  } else if(args[1] == "--threads"){
    run_parallel <- TRUE
    threads <- args[2]
    args <- args[-1:-2]
    message("Writing genotypes in parallel (requires furrr package)")
    require(furrr)
  }
  else {
    stop("Unknown argument: ", args[1])
  }
}

if(!exists("prefix")){
  stop("A prefix (name) for the database must be defined with the `--prefix` argument")
}
if(!exists("vcf_name")){
  stop("A vcf must be passed with the `--vcf`` argument")
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


this_file <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) {
    # Rscript
    return(normalizePath(sub(needle, "", cmdArgs[match])))
  } else {
    # 'source'd via R console
    return(normalizePath(sys.frames()[[1]]$ofile))
  }
}

this_dir <- dirname(this_file())

## clear command args before sourcing
commandArgs <- function(...) {}

source(paste0(this_dir, "/pipeline/01-generate-variant-ranges-index.R"))

source(paste0(this_dir, "/pipeline/02-build-db.R"))

source(paste0(this_dir, "/pipeline/03-index-db.R"))
