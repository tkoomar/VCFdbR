source("00-benchmark-functions.R")

## run
time_dat <- tibble(db_name = c("/Dedicated/jmichaelson-genome/1KG/exome/table-1kg.db", 
                               "/Dedicated/jmichaelson-genome/1KG/exome/file-1kg.db", 
                               "/Dedicated/jmichaelson-genome/1KG/file-1kg.db"), 
       db_type = c("table", "file", "file"), 
       data_type = c("exome", "exome", "genome")
       ) %>%
  crossing(
    n_vars = seq(50, 5000, by = 50)
  ) %>%
  crossing(
    cores = c(NA, 2, 4, 8, 16)
  ) %>%
  mutate(bin = cut_width(n_vars,                         
                         width = 500, 
                         center = 250, 
                         dig.lab = 4
                         )
  ) %>%
  sample_frac(1L)

pb <- progress_bar$new(total = nrow(time_dat), 
                       force = TRUE,
                       format = ":percent complete, estimated time remaining: :eta")

time_dat <- time_dat %>%
  mutate(
    time = pmap_dbl(list(.x = db_name, .y = n_vars, .z = cores), function(.x, .y, .z){
      pb$tick()
      time_var_pull(.db_name = .x, 
                    .n = .y, 
                    .cores = .z)
      
    })
  )
  
time_dat <- time_dat %>%
  select(-db_name) %>%
  mutate(cores = if_else(is.na(cores), 1, cores))

write_csv(time_dat, "03-pull-benchmark.csv")



