library(tidyverse)

pa <- read_rds("~/Box/GWAS_Paper/Revisions/Data/prev_abund.rds")

find_core <- function(x){
  tmp <- pa %>% filter(prev >= x) %>% count(ASVID) %>% filter(n==3)
  return(nrow(tmp))
}

find_core_by_site <- function(x){
  tmp <- pa %>% filter(prev >= x) %>% ungroup() %>% count(Site) 
  return(tmp)
}

get_ra <- function(x){
  tmp <- pa %>% filter(prev >= x) %>% add_count(ASVID) %>% filter(n == 3) %>% group_by(Site) %>% summarise(totra = sum(mean_ra)) 
  return(tmp)
}

get_ra_by_site <- function(x){
  tmp <- pa %>% filter(prev >= x) %>% group_by(Site) %>% summarise(totra = sum(mean_ra)) 
  return(tmp)
}

core_ranges <- data.frame(threshold = seq(.2 ,1, by=.02)) %>% 
  mutate(core_size = map_dbl(threshold, find_core))
core_ra <- data.frame(threshold = seq(.2 ,1, by=.02)) %>% 
  mutate(core_size = map(threshold, get_ra)) %>% 
  unnest(core_size)

core_ranges_by_site <- data.frame(threshold = seq(.2 ,1, by=.02)) %>% 
  mutate(core_size = map(threshold, find_core_by_site))
core_ra_by_site <- data.frame(threshold = seq(.2 ,1, by=.02)) %>% 
  mutate(core_size = map(threshold, get_ra_by_site)) %>% 
  unnest(core_size)
core_ranges_by_site %>% 
  unnest(core_size)

ggplot(core_ranges, aes(threshold, core_size)) +
  geom_point()

ggplot(core_ranges_by_site %>% unnest(core_size), aes(threshold, n, color = Site)) +
  geom_line(size = 1) +
  geom_point(data = core_ranges, aes(threshold, core_size), color = "black") +
  scale_color_manual(values = c("darkolivegreen", "gold", "dodgerblue")) +
  labs(x = "Core Prevalence Threshold", y = "Number of ASVS in local core") +
  theme_classic()

ggplot(core_ra_by_site, aes(threshold, totra / 10, color = Site)) +
  geom_line(size = 1) +
  geom_line(data = core_ra, size = 1, linetype = "dashed") +
  scale_color_manual(values = c("darkolivegreen", "gold", "dodgerblue")) +
  labs(x = "Core Prevalence Threshold", y = "% reads attributable to core ASVs") +
  theme_classic()
