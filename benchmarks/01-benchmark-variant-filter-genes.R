source("00-benchmark-functions.R")

genes <- get_genes(db_name = "/Dedicated/jmichaelson-genome/1KG/file-1kg.db")
write_csv(genes, "01-filter-benchmark-genome.csv")
## Bin genes by the number of qualifying variants ##
set.seed(129839)
genes_sample <- genes %>%
  filter(n_vars <= 10000) %>%
  group_by(bin) %>%
  sample_n(20) %>%
  write_csv("01-gene-sample-genome.csv")

genes <- get_genes(db_name = "/Dedicated/jmichaelson-genome/1KG/exome/table-1kg.db")
write_csv(genes, "01-filter-benchmark-exome-table.csv")
## Bin genes by the number of qualifying variants ##
set.seed(129839)
genes_sample <- genes %>% 
  filter(n_vars <= 2000) %>%
  group_by(bin) %>%
  sample_n(20, replace = TRUE) %>%
  write_csv("01-gene-sample-exome.csv")

genes <- get_genes(db_name = "/Dedicated/jmichaelson-genome/1KG/exome/table-1kg.db")
write_csv(genes, "01-filter-benchmark-exome-file.csv")

