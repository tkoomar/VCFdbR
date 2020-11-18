library(tidyverse)
library(ggforce)
library(cowplot)
theme_set(
  theme_light(base_size = 15) + 
    theme(strip.background = element_rect(fill = 'grey93'), 
          strip.text = element_text(color = 'grey25'), 
          panel.border = element_blank(), 
          axis.line.x.bottom = element_line(color = 'grey60'),
          axis.line.y.left   = element_line(color = 'grey60'),
          panel.grid  = element_blank()
          )
)

### filttering variants ####
filter_dat <-tibble(
 data_type = c("genome", "exome"), 
 data = list(read_csv("01-filter-benchmark-genome.csv"), read_csv("01-filter-benchmark-exome-table.csv"))
) %>%
  unnest(data)

plot_filter <- function(dat){
  dat %>%
    group_by(bin) %>%
    ungroup() %>%
    ggplot(aes(x = n_vars, y = filter_time)) + 
    geom_point(alpha = 0.05, 
               size = 3,
               color = viridisLite::viridis(9)[5]
    ) + 
    stat_smooth(alpha = 0.5, 
                size = 2, 
                color = viridisLite::viridis(9)[3]
    ) + 
    theme(plot.title = element_text(color = "grey30", hjust = 0.5))
}

filter_a <- filter_dat %>%
  filter(n_vars < 700 & data_type == 'exome') %>%
  plot_filter() + 
  coord_cartesian(xlim = c(0, 500), 
                  ylim = c(0, .3)
                  ) + 
  labs(x = "number of variants", 
       y = "seconds", 
       tag = "A"
  ) 

filter_b <- filter_dat %>%
  filter(n_vars < 10000 & data_type == 'genome') %>%
  plot_filter() + 
  coord_cartesian(xlim = c(0, 10000), 
                  ylim = c(0, 1.5)) + 
  labs(x = "number of variants", 
       y = NULL, 
       tag = "B"
  ) 

plot_grid(plotlist = list(filter_a, filter_b), 
          rel_widths = c(10, 9),
          nrow = 1) +
  ggsave2(filename = "vcfdb-filter-plot.pdf", width = 8, height = 4)

### pulling variants ####
pull_dat <- read_csv("03-pull-benchmark.csv")

plot_pull <- function(dat) {
  dat %>%
    group_by(bin, db_type, data_type) %>%
    ungroup() %>% 
    ggplot(aes(x = n_vars, y = time, color = as.factor(cores))) + 
    geom_point(size = 3, alpha = 0.05) +
    scale_color_viridis_d(option = 'C',
                          begin = 0, 
                          end = 0.85,
                          guide = guide_legend(title = "number of cores", 
                                               title.hjust = 0.5)
    ) + 
    stat_smooth(size = 0.75, se = FALSE) + 
    scale_x_continuous(breaks = scales::pretty_breaks(n = 2)) + 
    theme(legend.position = 'none',
          plot.title = element_text(color = "grey30", hjust = 0.5)
          )
}


pull_a <- pull_dat %>%
  filter(db_type == "table" & data_type == "exome") %>%
  plot_pull() +
  facet_zoom(xlim = c(0, 1000), 
             ylim = c(0, 15),
             zoom.size = 1, 
             show.area = F
             ) + 
  labs(x = NULL,
       y = "seconds", 
       tag = "A") 

pull_b <- pull_dat %>%
  filter(db_type == "file" & data_type == "exome") %>%
  plot_pull() +
  facet_zoom(xlim = c(0, 1000), 
             ylim = c(0, 2),
             zoom.size = 1, 
             show.area = F
  ) + 
  labs(x = NULL, 
       y = "seconds", 
       tag = "B") 

pull_c <- pull_dat %>%
  filter(db_type == "file" & data_type == "genome") %>%
  plot_pull() +
  facet_zoom(xlim = c(0, 1000), 
             ylim = c(0, 20),
             zoom.size = 1, 
             show.area = F
  ) + 
  labs(x = "number of variants", 
       y = "seconds", 
       tag = "C"
       ) 

pull_legend <- get_legend(pull_a +  theme(legend.position = "top"))

plot_grid(plotlist = list(pull_legend, 
                          pull_a , 
                          pull_b , 
                          pull_c
                          ),
          ncol = 1, 
          rel_heights = c(1, 
                          6, 
                          6, 
                          7)
          ) +
  ggsave2("vcfdb-pull-plot.pdf", width = 5, height = 8)


