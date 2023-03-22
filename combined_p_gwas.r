library(tidyverse)


core_p_values <- read_rds("~/Box/GWAS_Paper/Revisions/Data/core_combined_p_values.rds")
marks <- core_p_values %>% 
  group_by(CHR) %>% 
  summarise(max_pos = max(POS)) %>% 
  ungroup() %>% 
  mutate(POS2 = cumsum(lag(max_pos, default = 0)), subgen = ifelse(grepl("K", CHR), "K", "N"))
cp_anno <- read_tsv("~/Box/GWAS_Paper/GWAS/hard_annos.tsv", col_names = F) %>% 
  mutate(dist = pmin(abs(X2 - X7), abs(X2 - X8)))

tax2 <- tax %>% 
  mutate(Phylum2 = ifelse(Phylum == "Proteobacteria", as.character(Class), as.character(Phylum))) %>% 
  select(ASVID, Phylum2)

sig_bins <- core_p_values %>% dplyr::select(-c(to_add, subgen)) %>% 
  filter(ASVID != "ASV67") %>% 
  group_by(CHR, ASVID) %>% 
  mutate(bin = cut_width(POS, 25000)) %>% 
  group_by(ASVID, CHR, bin) %>% 
  mutate(stat = min(pvalue)) %>% 
  filter(stat < 5e-08) %>% 
  inner_join(marks) %>% 
  inner_join(tax2) 

negs <- core_p_values %>% dplyr::select(-c(to_add, subgen)) %>% 
  filter(ASVID != "ASV67" & pvalue > 5e-08) %>% 
  group_by(CHR) %>% 
  mutate(bin = cut_width(POS, 5000)) %>% 
  group_by(bin, CHR) %>% 
  mutate(mean_p = mean(pvalue), dis = abs(mean_p - pvalue)) %>% 
  arrange(dis) %>% mutate(rank = 1:n()) %>% filter(rank <= 10) %>% 
  inner_join(marks)

##################
#### Figure 4 ####
##################

ggplot(negs %>% filter(dis > 0), aes(POS+POS2, -log10(pvalue))) +
  geom_point(data = . %>% filter(subgen == "K"), aes(POS+POS2, -log10(pvalue)), shape = ".", alpha = 0.2, color = "grey50") +
  geom_point(data = . %>% filter(subgen == "N"), aes(POS+POS2, -log10(pvalue)), shape = ".", alpha = 0.2, color = "grey20") +
  geom_hline(yintercept = -log10(5e-08), linetype = "dashed", color = "grey50") +
  geom_point(data = sig_bins, aes(color = Phylum2), size = 0.25) +
  geom_point(data = sig_bins_soft, aes(color = Phylum2), size = 0.25) +
  geom_hline(yintercept = -Inf) +
  geom_vline(xintercept = -Inf) +
  theme_minimal() +
  labs(x="", y = "-Log10 Pvalue") +
  scale_color_manual(values = c("#9c003d", "#d43a4a", "#66c2a1", "#f3683f", "#2e85bd", "#fdde8aff","#9955ffff", "black","#abdca1")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_blank(), legend.position = "none", text = element_text(size = 15))


## CP go enrichments
genome_go <- read_rds("~/Box/GWAS_Paper/GWAS/genome_go_counts.rds") %>% 
  dplyr::rename(GO_category = GOs)
go_terms <- ontologyIndex::get_ontology("~/Box/GWAS_Paper/GWAS/go-basic.obo")
gos <- go_terms$name %>% as.data.frame() %>% rownames_to_column("GO_category") %>% as_tibble() %>% dplyr::rename(funct = 2)


cp_tail_anno <- read_tsv("~/Box/GWAS_Paper/GWAS/core_asv_anno.bed", col_names = F)
cp_tail_anno %>% dplyr::select(X1, X2, X4, X8, X16) %>% 
  group_by(X1, X2, X4, X8) %>% mutate(GO_category = str_split(X16, ",")) %>% 
  unnest(GO_category) %>% 
  dplyr::select(-X16) %>% 
  distinct() %>% ungroup() %>% na.omit() %>% 
  dplyr::count(GO_category, X4, name = "asv_count") %>% 
  group_by(X4) %>% mutate(asv_total = sum(asv_count)) %>% 
  inner_join(genome_go %>% mutate(genome_total = sum(n))) %>% 
  mutate(p = phyper(asv_count - 1, n, genome_total - n, asv_total, lower.tail = F)) %>% 
  arrange(p) %>% 
  inner_join(gos) %>% 
  dplyr::rename(ASVID = X4) %>% 
  inner_join(tax) %>%
  mutate(tax_unit = ifelse(!is.na(Genus), Genus, ifelse(!is.na(Family), Family, Order))) %>% 
  dplyr::select(ASVID, tax_unit, GO_category, funct, p) %>% 
  mutate(unit = factor(paste(ASVID, GO_category)), unit = fct_reorder(unit, p, max, .desc = T)) %>% 
  ggplot(aes(-log10(p), unit)) + geom_bar(stat = "identity", fill = "grey70") +
  geom_hline(yintercept = -Inf) +
  geom_vline(xintercept = -Inf) +
  theme_minimal()

cp_tail_anno %>% dplyr::select(X1, X2, X4, X8, X16) %>% 
  group_by(X1, X2, X4, X8) %>% mutate(GO_category = str_split(X16, ",")) %>% 
  unnest(GO_category) %>% 
  dplyr::select(-X16) %>% 
  distinct() %>% ungroup() %>% na.omit() %>% 
  dplyr::count(GO_category, X4, name = "asv_count") %>% 
  group_by(X4) %>% mutate(asv_total = sum(asv_count)) %>% 
  inner_join(genome_go %>% mutate(genome_total = sum(n))) %>% 
  mutate(p = phyper(asv_count - 1, n, genome_total - n, asv_total, lower.tail = F)) %>% 
  arrange(p) %>% 
  left_join(gos) %>% 
  dplyr::rename(ASVID = X4) %>% 
  inner_join(tax) %>%
  mutate(tax_unit = ifelse(!is.na(Genus), Genus, ifelse(!is.na(Family), Family, Order))) %>% 
  dplyr::select(ASVID, tax_unit, GO_category, funct, p) %>% 
  mutate(p.adj = p.adjust(p)) %>% ungroup() %>% 
  filter(p.adj < 0.1) %>% dplyr::count(ASVID) %>% arrange(-n)
  write_tsv("~/Library/CloudStorage/Box-Box/GWAS_Paper/Tables/comb_p_go.tsv")
