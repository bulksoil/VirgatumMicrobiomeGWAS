# Load in the gwas scans from the meta-analysis
## All of this happens on local computer

library(tidyverse)

core_gwas <- read_rds("~/Box/GWAS_Paper/Revisions/Data/core_gwas_metap.rds")

core_leading_snps <- core_gwas %>% 
  filter(pvalue < 5e-08) %>% 
  group_by(ASVID, CHR) %>% 
  mutate(bin = cut_width(POS, width = 50000)) %>% 
  add_count(bin) %>% 
  group_by(bin, ASVID, CHR) %>% 
  top_n(1, -log10(pvalue))

write_rds(core_leading_snps, "~/Box/GWAS_Paper/GWAS/core_leading_snps.rds")

#Perform anova with leading SNPs to calculate variance
## All of this happens on the server because the SNP matrix is huge

library(bigsnpr)
library(tidyverse)
library(switchgrassGWAS)
library(RNOmni)

## Load in SNP information
snp <- snp_attach("/home/edwards/SGMB/Reference/Pvirgatum_V5_GWAS_630g_33M_microbiome.rds")
plants <- snp$fam$sample.ID
markers <- tibble(CHR = snp$map$chromosome, POS = snp$map$physical.pos) %>% mutate(mat_col = 1:nrow(.))

## Load in PCs made from the kinship matrix
## This will be used to control for population structure in the anova
svd <- read_rds("/home/edwards/SGMB/PVirgatum_MAMP/mamp_gwas/SVD_630g_21.5M_SNPs_15PCs.rds")
PCs <- svd$u[,1:5] %>% as_tibble() %>% 
  mutate(PLANT_ID = plants) %>% 
  rename(PC1 = V1, PC2 = V2, PC3 = V3, PC4 = V4, PC5 =  V5)

core_leading_snps <- read_rds("~/SGMB/2019_GWAS/ALL/CORE_GWAS/core_leading_snps.rds") %>% 
  inner_join(markers) %>% 
  mutate(snp_name = paste(CHR, POS, sep = "_"))

## Subsetting the SNP matrix
sub_map <- snp$genotypes[,core_leading_snps$mat_col]
colnames(sub_map) <- paste(core_leading_snps$CHR, core_leading_snps$POS, sep="_")
sub_map <- as_tibble(sub_map) %>% 
  mutate(PLANT_ID = plants)

## Identify the core ASVs
core <- read_rds("../prev_abund.rds") %>% filter(prev >= 0.8) %>% dplyr::count(ASVID) %>% filter(n == 3)

## Load the abundance data and transform it
dat <- read_rds("../microbe_rst.rds") %>%
  mutate(PLANT_ID = factor(PLANT_ID)) %>%
  filter(ra > 0) %>%
  group_by(ASVID, Site) %>%
  mutate(rst2 = RankNorm(ra)) %>%
  inner_join(core, by = c("ASVID")) %>%
  mutate(ASVID2 = paste(Site, ASVID, sep = "_")) %>%
  group_by(PLANT_ID, ASVID, Site) %>%
  summarise(rst = mean(rst2)) 

## Here is the function that actually runs the anova
## The anova itself is specified like this
## aov(rst ~ .*Site, .)
## This means model the abundance (rst) with every column in the dataframe as an independent variable and its interaction with field site
snp_var_fun <- function(ASV){
  tmp <- data.frame(PLANT_ID = plants) %>% as_tibble() %>% 
    left_join(dat %>% filter(ASVID == ASV)) %>% 
    select(ASVID, Site, rst, PLANT_ID) 
  
  snps <- core_leading_snps_0.6 %>%
    filter(ASVID == ASV) %>% 
    pull(snp_name) %>% 
    unique()
  
  snps <- c("PLANT_ID", snps)
  
  sub_map2 <- sub_map %>% select(all_of(snps))
    
  tmp <- tmp %>% 
    inner_join(sub_map2) %>% 
    inner_join(PCs)

  tmp2 <- tmp %>% 
    select(-c(ASVID, PLANT_ID)) %>% 
    aov(rst ~ .*Site, .) %>% tidy() %>% 
    mutate(perv = sumsq / sum(sumsq)) %>% 
    arrange(p.value)
  
  return(tmp2)
}

core_snp_var <- core_leading_snps %>% 
  ungroup() %>% 
  select(ASVID) %>% 
  unique() %>% 
  mutate(aov_res = map(ASVID, ~snp_var_fun(.))) %>% 
  unnest(aov_res)

## At this point it save this data frame and bring it back to my own computer.

## Retrieve the dataframe from the server
core_snp_var <- read_rds("~/Box/GWAS_Paper/GWAS/core_snp_variation.rds")

core_snp_var %>% 
  filter(grepl("Chr", term)) %>% 
  group_by(ASVID) %>% 
  mutate(total_variance = sum(perv)) %>% 
  inner_join(tax %>% dplyr::select(-Seq)) %>% 
  write_tsv("~/Box/GWAS_Paper/Tables/core_snp_variance_0.8.tsv")



