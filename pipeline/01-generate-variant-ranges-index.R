#! Rscript
args = commandArgs(trailingOnly=TRUE)

chunk_size <- 1000

args = commandArgs(trailingOnly=TRUE)

while(length(args > 0) ){
  if(args[1] == '--prefix'){
    prefix <- args[2]
    args <- args[-1:-2]
  } else if(args[1] == "--vcf"){
    vcf_name <- args[2]
    args <- args[-1:-2]
  } else if(args[1] == "--chunk-size"){
    chunk_size <- args[2]
    args <- args[-1:-2]
  } else {
    stop("Unknown argument: ", args[1])
  }
}

if(!exists("prefix")){
  stop("A prefix (name) for the database must be defined with the `--prefix` argument")
}
if(!exists("vcf_name")){
  stop("A vcf must be passed with the `--vcf` argument")
}

suppressPackageStartupMessages(require(VariantAnnotation))
suppressPackageStartupMessages(require(tidyverse))

if(!file.exists(paste0(prefix, ".progress.RData"))){
#### Generate Chunks of Variants ####
message('######\nLOADING RANGES\n######')

full_ranges <- readVcf(TabixFile(vcf_name),
                       param = ScanVcfParam(fixed = c("ALT"), info = NA, samples = NA)
) %>%
  rowRanges()

## Make sure variants are decomposed ##
n_alts <- sapply(full_ranges$ALT, length)

if(!all(n_alts == 1)){
  multi_sites_file <- paste0(prefix, "-multiallelic-sites.tsv")
  
  names(full_ranges)[n_alts != 1] %>%
    as.data.frame() %>%
    write_tsv(multi_sites_file, col_names = FALSE)
  
  stop("The VCF ", vcf_name, " appears to have multialleleic sites. Sites written out to: ", multi_sites_file )
}

message('######\nMERGING RANGES INTO CHUNKS\n######')
date()

ranges_list <- split(full_ranges, 
                     seqnames(full_ranges) %>% as.numeric()
)

ranges_list <- lapply(ranges_list, function(x){
  chunks <- 1:length(x) %>%
    cut_width(chunk_size, labels = FALSE)
  
  y <- split(x, chunks)
  
  
})

ranges_list <- lapply(ranges_list, function(x){
  lapply(x, function(z){
    GRanges(seqnames = z@seqnames[1],
            ranges = IRanges(start = start(z[1]),
                             end = end(z[length(z)])
            )
    )
  })
})

chunk_ranges <- ranges_list %>% unlist %>% GRangesList %>% unlist
ranges_name <- paste0(prefix, "-chunk-ranges.rds")
saveRDS(chunk_ranges, file = ranges_name)

rm(ranges_list, full_ranges)
gc()

message('######\nDONE IDENTIFYING CHUNKS\nRANGES ARE WRITTEN TO:\n', 
        paste0(prefix, "-chunk-ranges.rds"),
        '\n######')
}