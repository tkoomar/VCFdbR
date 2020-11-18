source("00-benchmark-functions.R")

pull_variants(db_name = "/Dedicated/jmichaelson-genome/1KG/file-1kg.db",
              gene_file = "01-gene-sample-genome.csv",
              sample_size = 2054, 
              cores = c(NA, 2, 4, 8, 16, 32)
) %>%
  write_csv("02-pull-benchmark-genome.csv")

pull_variants(db_name =  "/Dedicated/jmichaelson-genome/1KG/exome/table-1kg.db",
              gene_file = "01-gene-sample-exome.csv",
              sample_size = 1000, 
              cores = NA
) %>%
  write_csv("02-pull-benchmark-exome.csv")
