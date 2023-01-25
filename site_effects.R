library(tidyverse)
library(emmeans)

## Define functions
get_emm <- function(x){
  emm <- x$emmeans
  return(emm %>% as_tibble())
}

get_contrasts <- function(x){
  cont <- x$contrasts
  return(cont %>% as_tibble())
}

## Load data

dat <- read_rds("~/Box/GWAS_Paper/site/transformed_counts.rds") %>% 
  group_by(ASVID) %>% 
  filter(sum(ra > 0) / n() > 0.5) %>% 
  mutate(rst2 = RNOmni::RankNorm(ra))

## Make models
mods <- dat %>%
  nest() %>% 
  mutate(mod = map(data, ~lm(rst2 ~ log10(depth) + SUBPOP*Site, .)))

## Extract factor effects using anova
site_sub_anova <- mods %>%
  mutate(ano = map(mod, ~broom::tidy(anova(.)))) %>% 
  select(ASVID, ano) %>%
  unnest(ano)

## Analyze site effects
site_emm <- mods %>%
  mutate(emm = map(mod, ~emmeans(., specs = pairwise ~ Site)))


site_means <- site_emm %>%
  mutate(emm2 = map(emm, ~get_emm(.))) %>% 
  select(ASVID, emm2) %>% 
  unnest(emm2)
site_contrasts <- site_emm %>%
  mutate(emm2 = map(emm, ~get_contrasts(.))) %>% 
  select(ASVID, emm2) %>% 
  unnest(emm2)


## Analyze subpop effects
subpop_emm <- mods %>%
  mutate(emm = map(mod, ~emmeans(., specs = pairwise ~ SUBPOP)))

subpop_means <- subpop_emm %>%
  mutate(emm2 = map(emm, ~get_emm(.))) %>% 
  select(ASVID, emm2) %>% 
  unnest(emm2)
subpop_contrasts <- subpop_emm %>%
  mutate(emm2 = map(emm, ~get_contrasts(.))) %>% 
  select(ASVID, emm2) %>% 
  unnest(emm2)

## Analyze interactions
subpop_site_emm <- mods %>%
  mutate(emm = map(mod, ~emmeans(., specs = pairwise ~ SUBPOP:Site)))

subpop_site_means <- subpop_site_emm %>%
  mutate(emm2 = map(emm, ~get_emm(.))) %>% 
  select(ASVID, emm2) %>% 
  unnest(emm2)
subpop_site_contrasts <- subpop_site_emm %>%
  mutate(emm2 = map(emm, ~get_contrasts(.))) %>% 
  select(ASVID, emm2) %>% 
  unnest(emm2)

## Save it
list(site_means = site_means, 
     site_contrasts = site_contrasts, 
     subpop_means = subpop_means, 
     subpop_contrasts = subpop_contrasts, 
     interaction_means = subpop_site_means, 
     interaction_contrasts = subpop_site_contrasts) %>% 
  write_rds("site_subpop_modeling_results.rds")


## Load it into local computer
site_sub_aov <- read_rds("~/Library/CloudStorage/Box-Box/GWAS_Paper/site/site_subpop_mod_anova.rds")
site_sub <- read_rds("~/Library/CloudStorage/Box-Box/GWAS_Paper/site/site_subpop_modeling_results.rds")
tax <- read_rds("~/Library/CloudStorage/Box-Box/GWAS_Paper/total_tax.rds")

## Let's check out the anova stats (this makes Figure 1C)
aov_effects_plot <- site_sub_aov %>% 
  group_by(term) %>% 
  mutate(padj = p.adjust(p.value, "BH")) %>% 
  group_by(ASVID) %>% 
  mutate(prop = sumsq / sum(sumsq)) %>% 
  filter(!grepl("log|Resid", term)) %>% 
  inner_join(tax) %>% ungroup() %>% 
  mutate(Phylum2 = ifelse(Phylum == "Proteobacteria", as.character(Class), as.character(Phylum))) %>% 
  mutate(Phylum2 = fct_lump(Phylum2, 8)) %>% ungroup() %>% 
  mutate(Phylum2 = fct_relevel(Phylum2, "Acidobacteria", "Actinobacteria", "Bacteroidetes", "Chloroflexi", "Alphaproteobacteria", "Deltaproteobacteria", "Gammaproteobacteria", "Verrucomicrobia", "Other")) %>%
  ggplot(aes(Phylum2, prop, color = term)) + 
  geom_hline(yintercept = 0, size = 0.25) + geom_vline(xintercept = -Inf) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, jitter.height = 0), alpha = 0.5, size = 0.5) +
  scale_color_manual(values = RColorBrewer::brewer.pal(4, "Greys")[-1])+
  theme_minimal() +
  labs(x = "") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

## Count and plot ASVs with signficant terms (this makes Figure 1D)
aov_n_asvs <- site_sub_aov %>% 
  group_by(term) %>% 
  mutate(padj = p.adjust(p.value, "BH")) %>% 
  group_by(ASVID) %>% 
  mutate(prop = sumsq / sum(sumsq)) %>% 
  filter(!grepl("log|Resid", term)) %>% 
  inner_join(tax) %>% ungroup() %>% 
  mutate(Phylum2 = ifelse(Phylum == "Proteobacteria", as.character(Class), as.character(Phylum))) %>% 
  mutate(Phylum2 = fct_lump(Phylum2, 8)) %>% ungroup() %>% 
  mutate(Phylum2 = fct_relevel(Phylum2, "Acidobacteria", "Actinobacteria", "Bacteroidetes", "Chloroflexi", "Alphaproteobacteria", "Deltaproteobacteria", "Gammaproteobacteria", "Verrucomicrobia", "Other")) %>% 
  filter(padj < 0.05) %>% 
  count(Phylum2, term) %>% 
  ggplot(aes(Phylum2, -n, fill = term)) +
  geom_hline(yintercept = 0, size = 0.25) + geom_vline(xintercept = -Inf) +
  geom_bar(position = "dodge", stat = "identity", width = 0.5) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(4, "Greys")[-1]) +
  theme_minimal()

gA <- ggplotGrob(aov_effects_plot)
gB <- ggplotGrob(aov_n_asvs)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))

site_sub_aov %>% 
  group_by(term) %>% 
  mutate(padj = p.adjust(p.value, "BH")) %>% 
  group_by(ASVID) %>% 
  mutate(prop = sumsq / sum(sumsq)) %>% 
  filter(!grepl("log", term)) %>% 
  inner_join(tax) %>% ungroup() %>% 
  mutate(Phylum2 = ifelse(Phylum == "Proteobacteria", as.character(Class), as.character(Phylum))) %>% 
  mutate(Phylum2 = fct_lump(Phylum2, 8)) %>% 
  lm(prop ~ Phylum2*term, .) %>% emmeans::emmeans(., specs = pairwise ~ Phylum2:term) %>% 
  .$contrasts %>% as_tibble() %>% 
  filter(grepl("Site.*Site", contrast) & !grepl("SUBPOP", contrast)) %>% 
  arrange(p.value)

## Which microbes are specifically enriched in sites?
site_spec_enr <- site_sub$site_contrasts %>% 
  ungroup() %>% 
  mutate(padj = p.adjust(p.value, "BH")) %>% 
  separate(contrast, into = c("S1", "S2"), sep = " - ") %>% 
  filter(padj < 0.001) %>% 
  mutate(enr = ifelse(estimate < 0, S2, S1)) %>% 
  count(ASVID, enr) %>% filter(n == 2) 


site_spec_depl <- site_sub$site_contrasts %>% 
  ungroup() %>% 
  mutate(padj = p.adjust(p.value, "BH")) %>% 
  separate(contrast, into = c("S1", "S2"), sep = " - ") %>% 
  filter(padj < 0.001) %>% 
  mutate(depl = ifelse(estimate > 0, S2, S1)) %>% 
  count(ASVID, depl) %>% filter(n == 2) 

site_enr_plot <- site_sub$site_means %>% 
  inner_join(site_spec_enr, by = "ASVID") %>% 
  group_by(ASVID) %>% 
  mutate(scaled = scale(emmean)) %>%
  mutate(Site = fct_relevel(Site, "P", "C", "M")) %>% 
  mutate(order = ifelse(enr == "P", 1, ifelse(enr == "C", 2, 3))) %>% 
  ggplot(aes(Site, paste(order, ASVID), fill = scaled)) +
  geom_tile() +
  scale_fill_viridis_c(option = "B") +
  theme_void() 
site_depl_plot <- site_sub$site_means %>% 
  inner_join(site_spec_depl, by = "ASVID") %>% 
  group_by(ASVID) %>% 
  mutate(scaled = scale(emmean)) %>%
  mutate(Site = fct_relevel(Site, "P", "C", "M")) %>% 
  mutate(order = ifelse(depl == "P", 1, ifelse(depl == "C", 2, 3))) %>% 
  ggplot(aes(Site, paste(order, ASVID), fill = scaled)) +
  geom_tile() +
  scale_fill_viridis_c(option = "B") +
  theme_void()

gridExtra::grid.arrange(site_enr_plot, site_depl_plot, nrow = 1)

gA <- ggplotGrob(site_enr_plot)
gB <- ggplotGrob(site_depl_plot)
grid::grid.newpage()
grid::grid.draw(cbind(gA, gB))

## Which microbes by subpop
sub_spec_enr <- site_sub$subpop_contrasts %>% 
  ungroup() %>% 
  mutate(padj = p.adjust(p.value, "BH")) %>% 
  separate(contrast, into = c("S1", "S2"), sep = " - ") %>% 
  filter(padj < 0.01) %>% 
  mutate(enr = ifelse(estimate < 0, S2, S1)) %>% 
  count(ASVID, enr) %>% filter(n == 2) 
sub_spec_depl <- site_sub$subpop_contrasts %>% 
  ungroup() %>% 
  mutate(padj = p.adjust(p.value, "BH")) %>% 
  separate(contrast, into = c("S1", "S2"), sep = " - ") %>% 
  filter(padj < 0.01) %>% 
  mutate(depl = ifelse(estimate > 0, S2, S1)) %>% 
  count(ASVID, depl) %>% filter(n == 2) 

enr_hm <- site_sub$subpop_means %>% 
  inner_join(sub_spec_enr, by = "ASVID") %>% 
  group_by(ASVID) %>% 
  mutate(scaled = scale(emmean)) %>% 
  ggplot(aes(paste(enr, ASVID), SUBPOP, fill = scaled)) +
  geom_tile() +
  scale_fill_viridis_b() +
  coord_equal() +
  labs(x = "") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

depl_hm <- site_sub$subpop_means %>% 
  inner_join(sub_spec_depl, by = "ASVID") %>% 
  group_by(ASVID) %>% 
  mutate(scaled = scale(emmean)) %>% 
  ggplot(aes(paste(depl, ASVID), SUBPOP, fill = scaled)) +
  geom_tile() +
  scale_fill_viridis_b() +
  coord_equal()+
  labs(x = "") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

fam_count_plot <- sub_spec_enr %>% 
  rename(subpop = enr) %>% 
  mutate(type = "Enriched") %>% 
  bind_rows(sub_spec_depl %>% rename(subpop = depl) %>% mutate(type = "Depleted")) %>% 
  inner_join(tax) %>% 
  count(type, subpop, Family) %>% 
  mutate(Family = ifelse(is.na(Family), "Unassigned", Family)) %>% 
  ggplot(aes(paste(n,Family), n)) +
  geom_bar(stat = "identity") + 
  facet_grid(.~type + subpop, scales = "free_x", space = "free_x") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

  

gridExtra::grid.arrange(enr_hm, depl_hm, fam_count_plot)

site_sub$subpop_contrasts %>% 
  ungroup() %>% 
  mutate(padj = p.adjust(p.value, "BH")) %>% 
  separate(contrast, into = c("S1", "S2"), sep = " - ") %>% 
  filter(padj < 0.01) %>% 
  count(ASVID) %>% 
  inner_join(tax) %>% 
  count(Family) %>% arrange(-n)

sub_spec_depl %>% 
  inner_join(tax) %>% 
  count(depl, Family) %>% 
  ggplot(aes(paste(n,Family), n)) +
  geom_bar(stat = "identity") + 
  facet_grid(.~depl, scales = "free_x", space = "free_x") + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

###################
## Figure 2
###################
## What about subpop effects within site?
dat <- read_rds("~/Library/CloudStorage/Box-Box/GWAS_Paper/h2/transformed_counts.rds") %>% 
  group_by(ASVID, Site) %>% 
  filter(sum(ra > 0) / n() > 0.5) %>% mutate(rst2 = RNOmni::RankNorm(ra))

mods <- dat %>%
  nest() %>% 
  mutate(mod = map(data, ~lm(rst2 ~ log10(depth) + SUBPOP, .)))

within_site_emm <- mods %>%
  mutate(emm = map(mod, ~emmeans(., specs = pairwise ~ SUBPOP)))

subpop_means <- within_site_emm %>%
  mutate(emm2 = map(emm, ~get_emm(.))) %>% 
  select(ASVID, emm2) %>% 
  unnest(emm2)
subpop_contrasts <- within_site_emm %>%
  mutate(emm2 = map(emm, ~get_contrasts(.))) %>% 
  select(ASVID, emm2) %>% 
  unnest(emm2)

subpop_contrasts <- read_rds("~/Library/CloudStorage/Box-Box/GWAS_Paper/site/within_site_subpop_contrasts.rds")
subpop_site_means <- read_rds("~/Library/CloudStorage/Box-Box/GWAS_Paper/site/within_site_subpop_means.rds")
prev_abund <- read_rds("~/Library/CloudStorage/Box-Box/GWAS_Paper/GWAS/prev_abund.rds")

## Fig 2B
subpop_contrasts %>% 
  group_by(Site) %>% 
  mutate(padj = p.adjust(p.value, "BH")) %>% 
  filter(padj < 0.1) %>% 
  separate(contrast, into = c("S1", "S2"), sep = " - ") %>% 
  mutate(enr = ifelse(estimate > 0, S1, S2)) %>% 
  count(Site, ASVID, enr) %>% 
  filter(n == 2) %>% 
  count(Site, enr) %>% 
  group_by(Site) %>% 
  mutate(prop = n / sum(n)) %>% 
  ungroup() %>% 
  mutate(Site = fct_relevel(Site, "P", "C", "M")) %>% 
  mutate(enr = fct_relevel(enr, "Gulf", "Atlantic", "Midwest")) %>% 
  ggplot(aes(enr, prop, fill = enr)) + geom_bar(stat = "identity") +
  facet_grid(.~Site) +
  scale_fill_manual(values = c("#FF0000", "#F2AD00", "#00A08A")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

subpop_contrasts %>% 
  group_by(Site) %>% 
  mutate(padj = p.adjust(p.value, "BH")) %>% 
  filter(padj < 0.1) %>% 
  separate(contrast, into = c("S1", "S2"), sep = " - ") %>% 
  mutate(enr = ifelse(estimate > 0, S1, S2)) %>% 
  count(Site, ASVID, enr) %>% 
  filter(n == 2) %>% 
  inner_join(subpop_site_means) %>% 
  mutate(enr = ifelse(enr == "Gulf", "A", ifelse(enr == "Atlantic", "B", "C"))) %>% 
  mutate(unit = paste(enr, ASVID)) %>% 
  mutate(Site = fct_relevel(Site, "P", "C", "M"),
         SUBPOP = fct_relevel(SUBPOP, "Gulf", "Atlantic", "Midwest")) %>% 
  ggplot(aes(SUBPOP, unit, fill = emmean)) +
  geom_tile() +
  facet_wrap(~Site, nrow = 1, scales = "free_y") +
  scale_fill_gradientn(colours = c("white", "grey90", "black")) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

# Fig 2C
subpop_contrasts %>% 
  group_by(Site) %>% 
  mutate(padj = p.adjust(p.value, "BH")) %>% 
  filter(padj < 0.1) %>% 
  separate(contrast, into = c("S1", "S2"), sep = " - ") %>% 
  mutate(enr = ifelse(estimate > 0, S1, S2)) %>% 
  count(Site, ASVID, enr) %>% 
  filter(n == 2) %>% 
  inner_join(prev_abund) %>% 
  group_by(Site) %>% 
  mutate(prop = n / sum(n)) %>% 
  ungroup() %>% 
  mutate(Site = fct_relevel(Site, "P", "C", "M"),
         enr = fct_relevel(enr, "Midwest", "Atlantic", "Gulf")) %>% 
  ggplot(aes(prev, fill = enr)) + 
  #geom_histogram(data = prev_abund %>% filter(prev > 0.5), aes(prev), fill = "grey50", inherit.aes = F, alpha = 0.1) +
  geom_histogram() + 
  facet_grid(enr~Site) +
  scale_fill_manual(values = c("springgreen3", "steelblue", "mediumorchid")) +
  theme_minimal()

subpop_contrasts %>% 
  select(Site, ASVID) %>% distinct() %>% 
  inner_join(prev_abund) %>% 
  ggplot(aes(prev)) +
  geom_histogram() + facet_grid(.~Site)

prev_sample <- function(n, Site_s){
  p <- prev_abund %>% filter(Site == Site_s & prev > 0.5) %>% 
    ungroup() %>% 
    sample_n(size = n, replace = F)
  return(mean(p$prev))
}

g_atx <- data.frame(iter = 1:1000) %>% 
  mutate(perm_prev = map(iter, ~prev_sample(189, "P"))) %>% 
  unnest(perm_prev)
m_atx <- data.frame(iter = 1:1000) %>% 
  mutate(perm_prev = map(iter, ~prev_sample(42, "P"))) %>% 
  unnest(perm_prev)
a_atx <- data.frame(iter = 1:1000) %>% 
  mutate(perm_prev = map(iter, ~prev_sample(37, "P"))) %>% 
  unnest(perm_prev)

g_cmo <- data.frame(iter = 1:1000) %>% 
  mutate(perm_prev = map(iter, ~prev_sample(1, "C"))) %>% 
  unnest(perm_prev)
m_cmo <- data.frame(iter = 1:1000) %>% 
  mutate(perm_prev = map(iter, ~prev_sample(49, "C"))) %>% 
  unnest(perm_prev)
a_cmo <- data.frame(iter = 1:1000) %>% 
  mutate(perm_prev = map(iter, ~prev_sample(10, "C"))) %>% 
  unnest(perm_prev)
g_kmi <- data.frame(iter = 1:1000) %>% 
  mutate(perm_prev = map(iter, ~prev_sample(13, "M"))) %>% 
  unnest(perm_prev)
m_kmi <- data.frame(iter = 1:1000) %>% 
  mutate(perm_prev = map(iter, ~prev_sample(36, "M"))) %>% 
  unnest(perm_prev)
a_kmi <- data.frame(iter = 1:1000) %>% 
  mutate(perm_prev = map(iter, ~prev_sample(11, "M"))) %>% 
  unnest(perm_prev)

bind_rows(g_atx %>% mutate(Site = "P", enr = "Gulf"),
          m_atx %>% mutate(Site = "P", enr = "Midwest"),
          a_atx %>% mutate(Site = "P", enr = "Atlantic"),
          g_cmo %>% mutate(Site = "C", enr = "Gulf"),
          m_cmo %>% mutate(Site = "C", enr = "Midwest"),
          a_cmo %>% mutate(Site = "C", enr = "Atlantic"),
          g_kmi %>% mutate(Site = "M", enr = "Gulf"),
          m_kmi %>% mutate(Site = "M", enr = "Midwest"),
          a_kmi %>% mutate(Site = "M", enr = "Atlantic")) %>% 
  inner_join(mean_prevs) %>% 
  mutate(Site = as.factor(Site), enr = as.factor(enr)) %>% 
  filter(perm_prev > mean_prev) %>% 
  group_by(Site, enr) %>% tally()

mean_prevs <- subpop_contrasts %>% 
  group_by(Site) %>% 
  mutate(padj = p.adjust(p.value, "BH")) %>% 
  filter(padj < 0.1) %>% 
  separate(contrast, into = c("S1", "S2"), sep = " - ") %>% 
  mutate(enr = ifelse(estimate > 0, S1, S2)) %>% 
  count(Site, ASVID, enr) %>% 
  filter(n == 2) %>% 
  inner_join(prev_abund) %>% 
  group_by(enr, Site) %>% 
  summarise(mean_prev = mean(prev))

