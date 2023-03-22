library(sommer)
library(tidyverse)
library(furrr)
library(RNOmni)

plan(strategy = multiprocess, gc = TRUE)
options(mc.cores = 24)

bc <- read_tsv("/home/edwards/SGMB/2019_GWAS/ALL/barcode_info.tsv", col_names = T) %>% mutate(Fwd = factor(Fwd), Rvs = factor(Rvs), BC = factor(paste0(Fwd, Rvs)))
dat <- read_rds("/home/edwards/SGMB/2019_GWAS/ALL/transformed_counts.rds") %>% mutate(Site = factor(Site), PLANT_ID = factor(PLANT_ID)) %>% 
  mutate(row = factor(str_extract(PLOT_GL, "[CMP][0-9]{2}"))) %>%
  inner_join(bc)
K <- read_rds("/home/edwards/SGMB/2019_GWAS/ALL/Kinship_van_Raden_630_individuals_SNPs_r2_20percent.rds")

gxe_model <- function(x){
  x2 <- x %>% filter(PLANT_ID %in% row.names(K)) %>% mutate(Site = factor(Site)) %>% filter(ra > 0) %>% group_by(Site) %>% mutate(rst = RankNorm(ra))
  mod <- mmer(rst ~ Site + log10(depth), random =~ vs(us(Site), PLANT_ID, Gu = K), rcov = ~vs(ds(Site),units), data = x2 %>% filter(ra > 0), tolparinv = 1e-01, verbose = T)  
  reduced <- mmer(rst ~ Site + log10(depth), random = ~vs(PLANT_ID, Gu=K), rcov = ~vs(ds(Site),units), data = x2 %>% filter(ra > 0), tolparinv = 1e-01, verbose = T)
  reduced2 <- mmer(rst ~ Site + log10(depth), rcov = ~vs(ds(Site),units), data = x2 %>% filter(ra > 0), tolparinv = 1e-01, verbose = T)
  get_p = as.character(anova(reduced, mod)$PrChisq[2] %>% str_split(., pattern = " ") %>% unlist())
  p_gxe= as.numeric(get_p[1])
  get_p <- as.character(anova(reduced2, reduced)$PrChisq[2] %>% str_split(., pattern = " ") %>% unlist())
  p_g= as.numeric(get_p[1])
  type = ifelse(p_gxe < 0.05, "GxE", ifelse(p_g < 0.05, "G", "none"))
  mod_sum <- get_var(mod, reduced, type)
  sites <- paste(unique(x2$Site), collapse="")
  h2 <- extract_h2(mod, sites)
  return(list(modsum = mod_sum, G = p_g, GxE = p_gxe, h2=h2, type = type))
}

gxe_model <- function(x){
  x2 <- x %>% filter(PLANT_ID %in% row.names(K)) %>% mutate(Site = factor(Site)) %>% filter(ra > 0) %>% group_by(Site) %>% mutate(rst = RankNorm(ra))
  mod <- mmer(rst ~ log10(depth) + Site, random =~ vs(us(Site), PLANT_ID, Gu = K), rcov = ~vs(ds(Site),units), data = x2 %>% filter(ra > 0), tolparinv = 1e-01, verbose = T)  
  reduced <- mmer(rst ~ log10(depth) + Site, random = ~vs(PLANT_ID, Gu=K), rcov = ~vs(ds(Site),units), data = x2 %>% filter(ra > 0), tolparinv = 1e-01, verbose = T)
  reduced2 <- mmer(rst ~ log10(depth) + Site, rcov = ~vs(ds(Site),units), data = x2 %>% filter(ra > 0), tolparinv = 1e-01, verbose = T)
  # print(summary(mod))
  get_p = as.character(anova(reduced, mod)$PrChisq[2] %>% str_split(., pattern = " ") %>% unlist())
  p_gxe= as.numeric(get_p[1])
  get_p <- as.character(anova(reduced2, reduced)$PrChisq[2] %>% str_split(., pattern = " ") %>% unlist())
  p_g= as.numeric(get_p[1])
  type = ifelse(p_gxe < 0.05, "GxE", ifelse(p_g < 0.05, "G", "none"))
  mod_sum <- get_var(mod, reduced, type)
  sites <- paste(unique(x2$Site), collapse="")
  h2 <- extract_h2(mod, sites)
  return(list(modsum = mod_sum, G = p_g, GxE = p_gxe, h2=h2, type = type))
}

gxe_model_cs <- function(x){
  E <- diag(length(unique(x$Site)))
  rownames(E) <- colnames(E) <- unique(x$Site)
  EK <- kronecker(E,K, make.dimnames = TRUE)
  
  x2 <- x %>% filter(PLANT_ID %in% row.names(K)) %>% mutate(Site = factor(Site)) %>% filter(tot > 0) %>% group_by(Site) %>% mutate(rst = RankNorm(tot)) %>% ungroup()
  mod <- mmer(rst ~ Site + log10(depth), random =~ vs(PLANT_ID, Gu=K) + vs(Site:PLANT_ID, Gu=EK), rcov = ~units, data = x2, tolparinv = 1e-01, verbose = T)  
  reduced <- mmer(rst ~ Site + log10(depth), random = ~vs(PLANT_ID, Gu=K), rcov = ~units, data = x2, tolparinv = 1e-01, verbose = T)
  reduced2 <- mmer(rst ~ Site + log10(depth), rcov = ~units, data = x2, tolparinv = 1e-01, verbose = T)
  # print(summary(mod))
  get_p = as.character(anova(reduced, mod)$PrChisq[2] %>% str_split(., pattern = " ") %>% unlist())
  p_gxe= as.numeric(get_p[1])
  get_p <- as.character(anova(reduced2, reduced)$PrChisq[2] %>% str_split(., pattern = " ") %>% unlist())
  p_g= as.numeric(get_p[1])
  type = ifelse(p_gxe < 0.001, "GxE", ifelse(p_g < 0.001, "G", "none"))
  #mod_sum <- get_var(mod, reduced, type)
  variances <- get_variances(mod)
  return(list(G = p_g, GxE = p_gxe, variances=variances, type = type))
}

safe_gxe <- possibly(gxe_model_cs, NA_real_)

get_var <- function(mod1, mod2, type){
  sumz <- summary(mod1)$varcomp %>% rownames_to_column("variable") %>% as_tibble()
  return(sumz)
}

get_variances <- function(x){
  variances <- bind_rows(
    vpredict(x, Va ~ (V1) / (V1+V2+V3)),
    vpredict(x, GxE ~ (V2) / (V1+V2+V3)),
    vpredict(x, Genetic ~ (V1+V2) / (V1+V2+V3)),
  )
  return(variances)
}

extract_h2 <- function(mod, sites){
  if(sites == "CM"){
    h2 <- bind_rows(vpredict(mod, h2_m ~ (V1) / ( V1+V4)), 
                    vpredict(mod, h2_c ~ (V3) / (V3+V5))) %>% 
      as.data.frame() %>% 
      mutate(p = pnorm(Estimate / SE, lower.tail = F)*2)
  }
  if(sites == "CMP"){
    h2 <- bind_rows(vpredict(mod, h2_m ~ (V1) / ( V1+V7)), 
                    vpredict(mod, h2_c ~ (V3) / (V3+V8)), 
                    vpredict(mod, h2_p ~ (V6) / (V6+V9))) %>% 
      as.data.frame() %>% 
      mutate(p = pnorm(Estimate / SE, lower.tail = F)*2)
  }
  if(sites == "MP") {
    h2 <- bind_rows(vpredict(mod, h2_m ~ (V1) / ( V1+V4)), 
                    vpredict(mod, h2_p ~ (V3) / (V3+V5))) %>% 
      as.data.frame() %>% 
      mutate(p = pnorm(Estimate / SE, lower.tail = F)*2)
  }
  return(h2)
}

extract_var <- function(x){
  tab <- x$variance %>% rownames_to_column()
  return(tab)
}

get_g <- function(x){
  return(x$G)
}
  
get_gxe <- function(x){
  return(x$GxE)
}

safe_gxe <- possibly(gxe_model, NA_real_)
safe_gxe <- safely(gxe_model)

pull_h2 <- function(mod){
  if(nrow(mod$h2) > 4){
    return(NA)
  }
  return(mod$h2 %>% rownames_to_column("Site"))
}

pull_var <- function(x){
  return(x$modsum)
}

dat2 <- dat %>% group_by(ASVID, Site) %>% filter(sum(ra > 0) / n() >= 0.8) %>% nest() %>% ungroup() %>% add_count(ASVID) %>% filter(n > 1) %>% 
  unnest() %>% group_by(ASVID) %>% nest() %>% ungroup()

sommer_mods <- dat2  %>% 
  mutate(mod = future_map(data, ~safe_gxe(.)))


sommer_mods <- read_rds("~/Box/GWAS_Paper/h2/sommer_mods2.rds")
sommer_mods_cs <- read_rds("~/Box/GWAS_Paper/h2/sommer_mods_cs.rds")
sommer_mods_cs_levels <- read_rds("~/Box/GWAS_Paper/h2/sommer_modls_cs_levels.rds") %>% select(-data)

pvalues_cs <- bind_rows(sommer_mods_cs %>% mutate(Level = "ASV") %>% rename(taxon = ASVID), sommer_mods_cs_levels) %>% 
  mutate(type = map_dbl(mod, ~length(.))) %>% 
  filter(type > 1) %>% 
  mutate(var = map(mod, ~extract_var(.))) %>% 
  mutate(pG = map(mod, ~get_g(.))) %>% 
  mutate(pGxE = map(mod, ~get_gxe(.))) %>% 
  unnest(var) %>% 
  dplyr::select(-c(data, mod, type)) %>% 
  unnest(c(pG, pGxE)) %>% dplyr::select(taxon, Level, pG, pGxE) %>% unique() %>% rename(Va = pG, GxE = pGxE) %>% gather(rowname, pvalue, -c(taxon, Level)) %>% 
  mutate(model = "cs")

pvalues_us <- bind_rows(sommer_mods %>% mutate(Level = "ASV") %>% rename(taxon = ASVID), sommer_mods_cs_levels) %>% 
  mutate(type = map_dbl(mod, ~length(.))) %>% 
  filter(type > 1) %>% 
  mutate(pG = map(mod, ~get_g(.))) %>% 
  mutate(pGxE = map(mod, ~get_gxe(.))) %>% 
  dplyr::select(-c(data, mod, type)) %>% 
  unnest(c(pG, pGxE)) %>% dplyr::select(taxon, Level, pG, pGxE) %>% unique() %>% rename(Va = pG, GxE = pGxE) %>% gather(rowname, pvalue, -c(taxon, Level)) %>% 
  mutate(model = "us")

bind_rows(pvalues_us, pvalues_cs) %>%  
  group_by(taxon, model, Level) %>% mutate(padj = p.adjust(pvalue, "BH")) %>% 
  dplyr::select(taxon, rowname, padj, model, Level) %>%
  filter(Level == "ASV") %>% 
  spread(rowname, padj) %>% 
  mutate(type = ifelse(GxE < 0.1, "GxE", ifelse(Va < 0.1, "G", "none"))) %>% 
  select(taxon, Level, model, type) %>% 
  spread(model, type) %>% ungroup() %>% 
  count(Level, types = paste(cs, us)) %>% 
  separate(types, into = c("cs", "us")) %>% filter(us != "NA") %>% 
  ggplot() +
  geom_segment(aes(x="cs", xend = "us", y = cs, yend = us, size = n)) +
  theme_minimal()

pvalues %>% 
  group_by(Level,rowname) %>% mutate(padj = p.adjust(pvalue, "BH")) %>% 
  dplyr::select(taxon, Level, rowname, padj) %>% 
  spread(rowname, padj) %>% 
  mutate(type = ifelse(GxE < 0.1, "GxE", ifelse(Va < 0.1, "Va", "none"))) %>% 
  filter(Level == "ASV") %>% 
  inner_join(tax, by = c("taxon" = "ASVID")) %>% 
  mutate(Order = ifelse(is.na(Order), paste(as.character(Class), "un"), as.character(Order))) %>% 
  count(type, Order, Phylum, Class) %>% 
  group_by(Order) %>% mutate(nn = sum(n)) %>% filter(nn > 1) %>% 
  mutate(type = fct_relevel(type, "none", "GxE", "G")) %>% 
  ungroup() %>% 
  arrange(Phylum, Class, Order) %>% mutate(id = 1:nrow(.)) %>% 
  mutate(Order = fct_reorder(Order, id, mean)) %>% 
  ggplot(aes(Order, n, fill = type)) + geom_bar(stat="identity", color = "white") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.25)) +
  scale_fill_manual(values = c("grey", "#9EBCDA", "#8856A7"))

bind_rows(sommer_mods_cs %>% mutate(Level = "ASV") %>% rename(taxon = ASVID), sommer_mods_cs_levels) %>% 
  mutate(type = map_dbl(mod, ~length(.))) %>% 
  filter(type > 1) %>% 
  mutate(var = map(mod, ~extract_var(.))) %>% 
  unnest(var) %>% 
  dplyr::select(-c(data, mod, type)) %>% 
  #unnest(c(pG, pGxE)) %>% 
  filter(rowname != "Genetic") %>% group_by(Level, taxon) %>% 
  nest() %>% mutate(sum = map_dbl(data, ~sum(.x$Estimate))) %>% 
  ungroup() %>% 
  mutate(Level = fct_relevel(Level, "Phylum", "Class", "Order", "Family", "Genus", "ASV")) %>% 
  arrange(desc(Level), sum) %>% mutate(id = 1:nrow(.)) %>% 
  unnest() %>% 
  #mutate(unit = fct_reorder(taxon, Estimate, sum)) %>% 
  ggplot(aes(id, Estimate, fill = rowname)) +
  geom_bar(stat = "identity", width = 1, color = "white", size = 0.1) +
  coord_flip() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = c("#9EBCDA", "#8856A7")) +
  theme_classic() +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank())
  
estimates <- bind_rows(sommer_mods_cs %>% mutate(Level = "ASV") %>% rename(taxon = ASVID), sommer_mods_cs_levels %>% filter(n_asv > 10)) %>% 
  mutate(type = map_dbl(mod, ~length(.))) %>% 
  filter(type > 1) %>% 
  mutate(var = map(mod, ~extract_var(.))) %>% 
  unnest(var) %>% 
  dplyr::select(-c(data, mod, type)) %>% 
  #unnest(c(pG, pGxE)) %>% 
  filter(rowname != "Genetic") %>% group_by(Level, taxon)

##########################
## Compound Symmetry model
##########################
## Gxe vs G
estimates %>% 
  mutate(Estimate = Estimate + 0.0001) %>% filter(Estimate > 0) %>% 
  dplyr::select(taxon, Level, rowname, Estimate) %>% spread(rowname, Estimate) %>% 
  mutate(prop = GxE / Va) %>% ungroup() %>% 
  mutate(Level = fct_relevel(Level, "ASV", "Genus", "Family", "Order", "Class", "Phylum")) %>% 
  mutate(Level2 = ifelse(Level == "Phylum", "F", ifelse(Level == "Class", "E", ifelse(Level == "Order", "D", ifelse(Level == "Family", "C", ifelse(Level == "Genus", "B", "A"))))), type = ifelse(GxE > Va, "GxE", "G"), together = paste(type, Level2)) %>% 
  arrange(Level, prop) %>% mutate(id = 1:nrow(.)) %>% 
  ggplot(aes(id, prop, fill = together)) + geom_bar(stat = "identity", width = 1) + 
  scale_y_log10() + coord_flip() + 
  scale_fill_manual(values = c(rev(RColorBrewer::brewer.pal(7, "Purples")[-1]), rev(RColorBrewer::brewer.pal(7, "Blues")[-1]))) + theme_classic()

estimates %>% 
  group_by(Level, rowname) %>% top_frac(0.1, Estimate) %>% ungroup() %>% 
  mutate(Level = factor(Level), Level = fct_relevel(Level, "Phylum", "Class", "Order", "Family", "Genus", "ASV")) %>% 
  ggplot(aes(Level, Estimate, color = rowname)) + 
  geom_boxplot(outlier.size = 0.1) + 
  scale_color_manual(values = c("#9EBCDA", "#8856A7")) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

pvalues %>% 
  group_by(Level,rowname) %>% mutate(padj = p.adjust(pvalue, "BH")) %>% 
  select(taxon, Level, rowname, padj) %>% 
  spread(rowname, padj) %>% 
  mutate(type = ifelse(GxE < 0.1, "GxE", ifelse(Va < 0.1, "Va", "none"))) %>% 
  inner_join(estimates, by = c("taxon", "Level")) %>% 
  filter(type == rowname) %>% 
  group_by(type) %>% top_n(40, Estimate) %>% 
  ungroup() %>%  
  left_join(tax, by = c("taxon" = "ASVID")) %>% 
  mutate(taxon = ifelse(Level == "ASV", paste(Family, taxon), taxon)) %>% 
  mutate(unit = paste0(taxon, " (", Level, ")"), unit = fct_reorder(unit, Estimate, mean)) %>% 
  ggplot(aes(unit, Estimate, fill = type)) + geom_bar(stat = "identity") +
  geom_linerange(aes(ymin = Estimate - SE, ymax = Estimate + SE)) +
  coord_flip() +
  facet_wrap(.~type, scales = "free")

pvalues %>% 
  group_by(Level,rowname) %>% mutate(padj = p.adjust(pvalue, "BH")) %>% 
  select(taxon, Level, rowname, padj) %>% 
  spread(rowname, padj) %>% 
  mutate(type = ifelse(GxE < 0.05, "GxE", ifelse(Va < 0.05, "Va", "none"))) %>% 
  filter(Level == "ASV") %>% 
  inner_join(tax, by = c("taxon" = "ASVID")) %>% 
  mutate(total = 304) %>% 
  mutate(type2 = ifelse(type == "none", "none", "G")) %>% 
  mutate(Order = ifelse(is.na(Order), paste(as.character(Class), "un"), as.character(Order))) %>% 
  count(Order, type2, total) %>% group_by(Order) %>% mutate(Order_total = sum(n)) %>% 
  filter(Order_total > 3) %>% 
  filter(type2 == "G") %>% ungroup() %>% mutate(sig = sum(n)) %>% 
  mutate(p = phyper(n-1, sig, total - sig, Order_total, lower.tail = F)) %>% 
  arrange(p) %>% mutate(padj = p.adjust(p, "BH")) 

pvalues %>% 
  group_by(Level,rowname) %>% mutate(padj = p.adjust(pvalue, "BH")) %>% 
  select(taxon, Level, rowname, padj) %>% 
  spread(rowname, padj) %>% 
  mutate(type = ifelse(GxE < 0.05, "GxE", ifelse(Va < 0.05, "Va", "none"))) %>% 
  filter(Level == "ASV") %>% 
  inner_join(tax, by = c("taxon" = "ASVID")) %>% 
  mutate(total = 304) %>% mutate(type2 = ifelse(type == "none", "none", "G")) %>% 
  count(Family, type2, total) %>% group_by(Family) %>% mutate(Family_total = sum(n)) %>% 
  filter(Family_total > 3) %>% 
  filter(type2 == "G") %>% ungroup() %>% mutate(sig = sum(n)) %>% 
  mutate(p = phyper(n-1, sig, total - sig, Family_total, lower.tail = F)) %>% 
  arrange(p) %>% mutate(padj = p.adjust(p, "BH"))

pvalues %>% 
  group_by(Level,rowname) %>% mutate(padj = p.adjust(pvalue, "BH")) %>% 
  select(taxon, Level, rowname, padj) %>% 
  spread(rowname, padj) %>% 
  mutate(type = ifelse(GxE < 0.1, "GxE", ifelse(Va < 0.1, "Va", "none"))) %>% 
  filter(Level == "ASV") %>% 
  inner_join(tax, by = c("taxon" = "ASVID")) %>% 
  mutate(total = 304) %>% mutate(type2 = ifelse(type == "none", "none", "G")) %>% 
  count(Genus, type2, total) %>% group_by(Genus) %>% mutate(Genus_total = sum(n)) %>% 
  filter(type2 == "G") %>% ungroup() %>% mutate(sig = sum(n)) %>% 
  mutate(p = phyper(n-1, sig, total - Genus_total, Genus_total, lower.tail = F)) %>% 
  arrange(p) %>% mutate(padj = p.adjust(p, "BH"))
  

bind_rows(sommer_mods_cs %>% mutate(Level = "ASV") %>% rename(taxon = ASVID), sommer_mods_cs_levels %>% filter(n_asv > 10)) %>% 
  mutate(type = map_dbl(mod, ~length(.))) %>% 
  filter(type > 1) %>% 
  mutate(var = map(mod, ~extract_var(.))) %>% 
  unnest(var) %>% 
  dplyr::select(-c(data, mod, type)) %>% 
  filter(rowname != "Genetic") %>% group_by(Level, taxon) %>% 
  inner_join(pvalues) %>% 
  group_by(Level, rowname) %>% 
  mutate(padj = p.adjust(pvalue)) %>% 
  select(taxon, Level, rowname, padj) %>% 
  spread(rowname, padj) %>% 
  mutate(type = ifelse(GxE < 0.1, "GxE", ifelse(Va < 0.1, "Va", "none"))) %>% 
  filter(type == "GxE" & Level == "ASV") %>% 
  inner_join(tax, by = c("taxon" = "ASVID")) %>% count(Family) %>% arrange(-n)

bind_rows(sommer_mods_cs %>% mutate(Level = "ASV") %>% rename(taxon = ASVID), sommer_mods_cs_levels) %>% 
  mutate(type = map_dbl(mod, ~length(.))) %>% 
  filter(type > 1) %>% 
  mutate(var = map(mod, ~extract_var(.))) %>% 
  mutate(pG = map(mod, ~get_g(.))) %>% 
  mutate(pGxE = map(mod, ~get_gxe(.))) %>% 
  unnest(var) %>% 
  dplyr::select(-c(data, mod, type)) %>% 
  unnest(c(pG, pGxE)) %>% 
  select(taxon, Level, rowname, Estimate) %>% filter(rowname != "Genetic") %>% 
  spread(rowname, Estimate) %>% 
  mutate(ratio = log2(GxE/Va)) %>% 
  mutate(Level = fct_relevel(Level, "Phylum", "Class", "Order", "Family", "Genus", "ASV")) %>% 
  arrange(Level, ratio) %>% filter(Va != 0 & GxE != 0) %>%  mutate(id = 1:nrow(.)) %>% 
  ggplot(aes(id, ratio)) + geom_bar(stat = "identity", width = 1) +
  facet_grid(.~Level, space = "free_x", scales = "free_x")

######################
## Unstructured models
######################
# Heritability
sommer_mods %>%
  mutate(type = map_dbl(mod, ~length(.))) %>% 
  filter(type > 1) %>% 
  mutate(type = map(mod, ~extract_type(.))) %>% 
  unnest(type) %>% 
  #filter(type == "GxE") %>% 
  mutate(va = map(mod, ~pull_var(.))) %>% 
  dplyr::select(ASVID, va) %>% 
  unnest(va) %>% #add_count(ASVID) %>% filter(n == 3) %>% 
  filter(!grepl("units", variable) & Constraint == "Positive") %>% 
  add_count(ASVID) %>% filter(n == 3) %>% 
  mutate(Site = str_extract(variable, "^.")) %>% 
  mutate(p = pnorm(VarComp / VarCompSE, lower.tail = F)) %>% 
  mutate(padj = p.adjust(p, "BH"), sig = ifelse(padj < 0.1, "sig", "ns")) %>% 
  inner_join(tax) %>% mutate(Phylum2 = ifelse(Phylum == "Proteobacteria", as.character(Class), as.character(Phylum))) %>% 
  arrange(Phylum, Class, Order, Family, Genus, ASVID) %>% 
  mutate(rank = 1:nrow(.)) %>% 
  mutate(Genus = ifelse(grepl("Burkholder", Genus), "Burkholderia", as.character(Genus)),
         Genus = ifelse(grepl("Rhizobium", Genus), "Rhizobium", as.character(Genus))) %>% 
  mutate(unit = ifelse(!is.na(Genus), paste(Genus, ASVID), 
                       ifelse(!is.na(Family), paste(Family, ASVID), 
                              ifelse(!is.na(Order), paste(Order, ASVID), paste(Class, ASVID))))) %>% 
  mutate(unit = fct_reorder(unit, rank, mean), Site = fct_relevel(Site, "P", "C", "M")) %>% 
  mutate(Phylum2 = fct_lump(Phylum2, 7)) %>% 
  ggplot(aes(unit, Site, color = sig, fill = Phylum2, size = VarComp)) +
  geom_point(shape = 21) +
  coord_equal() +
  coord_flip() +
  scale_color_manual(values = c("white", "black")) +
  scale_fill_manual(values = c("#9c003d", "#d43a4a", "#66c2a1", "#f3683f", "#fdae61", "#2e85bd", "#fdde8a", "#5c4aa1", "#abdca1", "grey50", "white", "black")) +
  theme_minimal() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, hjust = 1))

# Va
sommer_mods %>%
  mutate(type = map_dbl(mod, ~length(.))) %>% 
  filter(type > 1) %>% 
  mutate(type = map(mod, ~extract_type(.))) %>% 
  unnest(type) %>% 
  filter(type == "G") %>% 
  mutate(va = map(mod, ~pull_var(.))) %>% 
  dplyr::select(ASVID, va) %>% 
  unnest(va) %>% #add_count(ASVID) %>% filter(n == 3) %>% 
  filter(!grepl("units", variable) & Constraint == "Positive") %>% 
  mutate(Site = str_extract(variable, "^.")) %>% 
  mutate(p = pnorm(VarComp / VarCompSE, lower.tail = F)) %>% 
  mutate(padj = p.adjust(p, "BH"), sig = ifelse(padj < 0.1, "sig", "ns")) %>% 
  inner_join(tax) %>% mutate(Phylum2 = ifelse(Phylum == "Proteobacteria", as.character(Class), as.character(Phylum))) %>% 
  arrange(Phylum, Class, Order, Family, Genus, ASVID) %>% mutate(rank = 1:nrow(.), unit = paste(Genus, ASVID)) %>% 
  mutate(unit = fct_reorder(unit, rank, mean), Site = fct_relevel(Site, "P", "C", "M")) %>% 
  mutate(Phylum2 = fct_lump(Phylum2, 7)) %>% 
  ggplot(aes(unit, Site, color = sig, fill = Phylum2, size = VarComp)) +
  geom_point(shape = 21) +
  coord_equal() +
  scale_color_manual(values = c("white", "black")) +
  scale_fill_manual(values = c("#9c003d", "#d43a4a", "#66c2a1", "#f3683f", "#fdae61", "#2e85bd","#5c4aa1", "#abdca1", "grey50", "white", "black")) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, hjust = 1))

## Covariance from unstructured model
sommer_mods %>%
  mutate(type = map_dbl(mod, ~length(.))) %>% 
  filter(type > 1) %>% 
  mutate(type = map(mod, ~extract_type(.))) %>% 
  unnest(type) %>% 
  filter(type == "G") %>% 
  mutate(va = map(mod, ~pull_var(.))) %>% 
  dplyr::select(ASVID, va) %>% 
  unnest(va) %>% #add_count(ASVID) %>% filter(n == 3) %>% 
  filter(!grepl("units", variable) & Constraint != "Positive") %>% 
  mutate(p = pnorm(VarComp / VarCompSE, lower.tail = F)) %>% 
  ggplot(aes(VarComp, fill = variable)) + 
  geom_density(alpha = 0.25, adjust = 2, color = "white") + 
  scale_fill_manual(values = RColorBrewer::brewer.pal(4, "RdPu")[-1]) + 
  theme_minimal()

sommer_mods %>%
  mutate(type = map_dbl(mod, ~length(.))) %>% 
  filter(type > 1) %>% 
  mutate(type = map(mod, ~extract_type(.))) %>% 
  unnest(type) %>% 
  filter(type == "G") %>% 
  mutate(va = map(mod, ~pull_var(.))) %>% 
  dplyr::select(ASVID, va) %>% 
  unnest(va) %>% #add_count(ASVID) %>% filter(n == 3) %>% 
  filter(!grepl("units", variable) & Constraint != "Positive") %>% 
  mutate(p = pnorm(VarComp / VarCompSE, lower.tail = F)) %>% 
  mutate(padj = p.adjust(p, method = "BH")) %>% 
  filter(padj < 0.1) %>% 
  ggplot(aes(variable, fill = variable)) + geom_bar() + 
  scale_fill_manual(values = RColorBrewer::brewer.pal(4, "RdPu")[-1]) +
  theme_minimal()

sommer_mods %>%
  mutate(type = map_dbl(mod, ~length(.))) %>% 
  filter(type > 1) %>% 
  mutate(type = map(mod, ~extract_type(.))) %>% 
  unnest(type) %>% 
  filter(type == "G") %>% 
  mutate(va = map(mod, ~pull_var(.))) %>% 
  dplyr::select(ASVID, va) %>% 
  unnest(va) %>% #add_count(ASVID) %>% filter(n == 3) %>% 
  filter(!grepl("units", variable) & Constraint != "Positive") %>% 
  mutate(p = pnorm(VarComp / VarCompSE, lower.tail = F)) %>% 
  mutate(padj = p.adjust(p, method = "BH")) %>% 
  filter(padj < 0.1) %>% 
  count(ASVID)

sommer_mods %>%
  mutate(type = map_dbl(mod, ~length(.))) %>% 
  filter(type > 1) %>% 
  mutate(type = map(mod, ~extract_type(.))) %>% 
  unnest(type) %>% 
  filter(type == "G") %>% 
  mutate(va = map(mod, ~pull_var(.))) %>% 
  dplyr::select(ASVID, va) %>% 
  unnest(va) %>% #add_count(ASVID) %>% filter(n == 3) %>% 
  filter(!grepl("units", variable) & Constraint != "Positive") %>% 
  mutate(p = pnorm(VarComp / VarCompSE, lower.tail = F)) %>% 
  mutate(padj = p.adjust(p, method = "BH")) %>% 
  filter(padj < 0.1) %>% 
  count(ASVID) %>% 
  inner_join(tax) %>% 
  count(Order) %>% arrange(-n)

sommer_mods %>%
  mutate(type = map_dbl(mod, ~length(.))) %>% 
  filter(type > 1) %>% 
  mutate(type = map(mod, ~extract_type(.))) %>% 
  unnest(type) %>% 
  filter(type == "G") %>% 
  mutate(va = map(mod, ~pull_var(.))) %>% 
  dplyr::select(ASVID, va) %>% 
  unnest(va) %>% 
  filter(!grepl("units", variable) & Constraint != "Positive") %>% 
  aov(VarComp ~ variable, .) %>% TukeyHSD(.)

mean_h2 <- sommer_mods %>%
  mutate(type = map_dbl(mod, ~length(.))) %>% 
  filter(type > 1) %>% 
  mutate(h2 = map(mod, ~pull_h2(.))) %>% 
  dplyr::select(ASVID, h2) %>% 
  unnest(h2) %>% 
  inner_join(tax) %>% mutate(Phylum2 = ifelse(Phylum == "Proteobacteria", as.character(Class), as.character(Phylum))) %>% 
  arrange(Phylum, Class, Order, Family, Genus, ASVID) %>% mutate(rank = 1:nrow(.), unit = paste(Genus, ASVID)) %>% 
  mutate(unit = fct_reorder(unit, rank, mean)) %>%
  group_by(ASVID, unit, Phylum2) %>% summarise(mean_h2 = mean(Estimate)) %>% 
  #mutate(sig = ifelse(p < 0.05, "sig", "ns")) %>% 
  #filter(sig == "sig") %>% 
  inner_join(gxe_p) %>% mutate(sig = ifelse(gxeP_adj < 0.05, "sig", "ns")) %>% 
  ggplot(aes(unit, mean_h2, fill = sig)) +
  geom_bar(stat = "identity") +
  labs(x = "") +
  theme_minimal() +
  scale_fill_manual(values = c("grey30", "red")) +
  theme(axis.text.x = element_blank()) 

gxe_p <- sommer_mods %>%
  mutate(type = map_dbl(mod, ~length(.))) %>% 
  filter(type > 1) %>% 
  mutate(gxeP = map(mod, ~get_gxe(.))) %>% 
  unnest(gxeP) %>% 
  mutate(gP = map(mod, ~get_g(.))) %>% 
  unnest(gP) %>% 
  separate(gxeP, into = c("gxeP", "nothing"), sep = " ") %>% 
  separate(gP, into = c("gP", "nothing2"), sep = " ") %>% 
  dplyr::select(-c(nothing, nothing2, data, mod, type)) %>% 
  mutate(gxeP_adj = p.adjust(gxeP, "BH"),
         gP_adj = p.adjust(gP, "BH"),
         type = ifelse(gxeP_adj < 0.05, "GxE", ifelse(gP_adj < 0.05, "G", "ns"))) 

core <- sommer_mods %>%
  mutate(type = map_dbl(mod, ~length(.))) %>% 
  filter(type > 1) %>% 
  mutate(h2 = map(mod, ~pull_h2(.))) %>% 
  dplyr::select(ASVID, h2) %>% 
  mutate(l = map_dbl(h2, ~nrow(.))) %>% filter(l <= 3) %>% 
  unnest(h2) %>% 
  count(ASVID) %>% filter(n == 3)


get_top1 <- function(x){
  j <- read_rds(x)
  return(j %>% top_frac(wt = log10p, n = 0.01))
}


mod_h2 <- function(x){
  x2 <- x %>% mutate(rst2 = RNOmni::RankNorm(ra))
  mod <- mmer(rst2 ~ log10(depth), random =~ vs(PLANT_ID, Gu = K), rcov = ~units, data = x2 %>% filter(ra > 0), tolparinv = 1e-01, verbose = T)  
  mod_r <- mmer(rst2 ~ log10(depth), rcov = ~units, data = x2 %>% filter(ra > 0), tolparinv = 1e-01, verbose = T)
  test_stat <- anova(mod, mod_r)
  p = as.character(test_stat$PrChisq[2])
  h2 = 
}

extract_type <- function(x){
  return(x$type)
}

sommer_mods %>% 
  mutate(type = map_dbl(mod, ~length(.))) %>% 
  filter(type > 1) %>% 
  mutate(type = map(mod, ~extract_type(.))) %>% 
  unnest(type) %>% 
  inner_join(tax) %>% 
  count(Phylum, Class, Order, type) %>% 
  arrange(Phylum, Class, Order) %>% mutate(order = 1:nrow(.)) %>% 
  mutate(unit = fct_reorder(Order, order, mean)) %>% 
  mutate(type = fct_relevel(type, "none", "GxE", "G")) %>% 
  filter(!is.na(Order)) %>% 
  ggplot(aes(unit, n, fill = type)) + geom_bar(stat = "identity") + #facet_grid(type ~ .) +
  scale_fill_manual(values = c("#E0ECF4", "#9EBCDA", "#8856A7")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .25))
