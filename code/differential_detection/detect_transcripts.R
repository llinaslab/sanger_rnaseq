source("code/utils.R")
load_essentials()
sshhh("limma")
sshhh("org.Pf.plasmo.db")

setwd("~/Dropbox/projects/sanger_rnaseq/")

# Main --------------------------------------------------------------------

# import gene expression
x3d7exp <- readRDS("output/neighboring_genes/gene_reduced_3d7_abund.rds")
xhb3exp <- readRDS("output/neighboring_genes/gene_reduced_hb3_abund.rds")
xitexp  <- readRDS("output/neighboring_genes/gene_reduced_it_abund.rds")

exp <- dplyr::bind_rows(x3d7exp,xhb3exp,xitexp)

rm(x3d7exp,xhb3exp,xitexp)

# import core genes for comparison
core_genes <- readr::read_tsv("data/gene_lists/core_pf3d7_genes.txt",col_names=F)$X1

# unique gene IDs
ug <- unique(exp$gene_id)

# Create data frame of genes that are "on" in each strain
detected <- exp %>%
  dplyr::group_by(gene_id, strain) %>%
  dplyr::summarise(n = sum(TPM >= 5)) %>%
  dplyr::mutate(on = ifelse(n > 0, 1, 0)) %>%
  dplyr::ungroup() %>%
  dplyr::select(-n)

# What genes are "on" in all strains
onall <- detected %>%
  dplyr::filter(on == 1) %>%
  dplyr::group_by(gene_id) %>%
  dplyr::summarise(count = length(unique(strain))) %>%
  dplyr::filter(count == 3) %$%
  gene_id

# What genes are "off" in all strains
offall <- detected %>%
  dplyr::filter(on == 0) %>%
  dplyr::group_by(gene_id) %>%
  dplyr::summarise(count = length(unique(strain))) %>%
  dplyr::filter(count == 3) %$%
  gene_id


# Create venn diagrams ----------------------------------------------------

# create venn diagram for detected and undetected transcripts
on3d7 <- detected %>% filter(strain == "3d7" & on == 1 & gene_id %in% core_genes) %$% gene_id
onhb3 <- detected %>% filter(strain == "hb3" & on == 1& gene_id %in% core_genes) %$% gene_id
onit  <- detected %>% filter(strain == "it" & on == 1 & gene_id %in% core_genes) %$% gene_id

off3d7 <- detected %>% filter(strain == "3d7" & on == 0 & gene_id %in% core_genes) %$% gene_id
offhb3 <- detected %>% filter(strain == "hb3" & on == 0 & gene_id %in% core_genes) %$% gene_id
offit  <- detected %>% filter(strain == "it" & on == 0 & gene_id %in% core_genes) %$% gene_id

X3d7 <- (ug %in% on3d7)
Xhb3 <- (ug %in% onhb3)
Xit  <- (ug %in% onit)
von  <- vennCounts(cbind(X3d7, Xhb3, Xit))

X3d7 <- (ug %in% off3d7)
Xhb3 <- (ug %in% offhb3)
Xit  <- (ug %in% offit)
voff  <- vennCounts(cbind(X3d7, Xhb3, Xit))

# Detect genes specific to strains ----------------------------------------

# Off genes per strain
off_only_3d7 <- detected %>%
  filter(strain == "3d7" & on == 0 & gene_id %nin% offhb3 & gene_id %nin% offit) %$%
  gene_id
off_only_hb3 <- detected %>%
  filter(strain == "hb3" & on == 0 & gene_id %nin% off3d7 & gene_id %nin% offit) %$%
  gene_id
off_only_it <- detected %>%
  filter(strain == "it" & on == 0 & gene_id %nin% offhb3 & gene_id %nin% off3d7) %$%
  gene_id

# On genes per strain
on_only_3d7 <- detected %>%
  filter(strain == "3d7" & on == 1 & gene_id %nin% onhb3 & gene_id %nin% onit) %$%
  gene_id
on_only_hb3 <- detected %>%
  filter(strain == "hb3" & on == 1 & gene_id %nin% on3d7 & gene_id %nin% onit) %$%
  gene_id
on_only_it <- detected %>%
  filter(strain == "it" & on == 1 & gene_id %nin% onhb3 & gene_id %nin% on3d7) %$%
  gene_id


# Annotate lists ----------------------------------------------------------

# import gene names
gene_names <- as.data.frame(org.Pf.plasmoGENENAME)

# annotate with total TPMs
off_3d7hb3_exp <- exp %>%
  filter(gene_id %in% union(off3d7, offhb3) & (strain == "3d7" | strain == "hb3")) %>%
  group_by(gene_id, strain) %>%
  summarise(min = min(TPM), max = max(TPM), mean = mean(TPM), sd = sd(TPM))

off_3d7it_exp <- exp %>%
  filter(gene_id %in% union(off3d7, offit) & (strain == "3d7" | strain == "it")) %>%
  group_by(gene_id, strain) %>%
  summarise(min = min(TPM), max = max(TPM), mean = mean(TPM), sd = sd(TPM))

off_hb3it_exp <- exp %>%
  filter(gene_id %in% union(offhb3, offit) & (strain == "hb3" | strain == "it")) %>%
  group_by(gene_id, strain) %>%
  summarise(min = min(TPM), max = max(TPM), mean = mean(TPM), sd = sd(TPM))

# annotate the off genes
off_3d7hb3_exp <- inner_join(off_3d7hb3_exp, gene_names, by = "gene_id")
off_3d7it_exp  <- inner_join(off_3d7it_exp, gene_names, by = "gene_id")
off_hb3it_exp  <- inner_join(off_hb3it_exp, gene_names, by = "gene_id")


# Write data  -------------------------------------------------------------

cat("\nWriting data to output... ")

write_tsv(x = off_3d7hb3_exp,
          path = "output/differential_detection/off_genes_3d7_hb3.tsv")
write_tsv(x = off_3d7it_exp,
          path = "output/differential_detection/off_genes_3d7_it.tsv")
write_tsv(x = off_hb3it_exp,
          path = "output/differential_detection/off_genes_hb3_it.tsv")
write_tsv(x = inner_join(detected, gene_names, by = "gene_id"),
          path = "output/differential_detection/detected.tsv")

off_hb3_not_3d7 <- off_3d7hb3_exp %>% filter(strain == "hb3" & max < 5) %$% gene_id
off_3d7_not_hb3 <- off_3d7hb3_exp %>% filter(strain == "3d7" & max < 5) %$% gene_id
off_it_not_3d7  <- off_3d7it_exp %>% filter(strain == "it" & max < 5) %$% gene_id
off_3d7_not_it  <- off_3d7it_exp %>% filter(strain == "3d7" & max < 5) %$% gene_id
off_hb3_not_it  <- off_hb3it_exp %>% filter(strain == "hb3" & max < 5) %$% gene_id
off_it_not_hb3  <- off_hb3it_exp %>% filter(strain == "it" & max < 5) %$% gene_id

output <- "output/differential_detection/"

write(off_hb3_not_3d7,paste0(output,"off_genes_3d7_hb3/off_hb3_list.txt"),sep="\n")
write(off_3d7_not_hb3,paste0(output,"off_genes_3d7_hb3/off_3d7_list.txt"),sep="\n")
write(off_it_not_3d7,paste0(output,"off_genes_3d7_it/off_it_list.txt"),sep="\n")
write(off_3d7_not_it,paste0(output,"off_genes_3d7_it/off_3d7_list.txt"),sep="\n")
write(off_hb3_not_it,paste0(output,"off_genes_hb3_it/off_hb3_list.txt"),sep="\n")
write(off_it_not_hb3,paste0(output,"off_genes_hb3_it/off_it_list.txt"),sep="\n")

write(on_only_3d7,paste0(output,"on_only_3d7_list.txt"),sep="\n")
write(on_only_hb3,paste0(output,"on_only_hb3_list.txt"),sep="\n")
write(on_only_it,paste0(output,"on_only_it_list.txt"),sep="\n")

write(off_only_3d7,paste0(output,"off_only_3d7_list.txt"),sep="\n")
write(off_only_hb3,paste0(output,"off_only_hb3_list.txt"),sep="\n")
write(off_only_it,paste0(output,"off_only_it_list.txt"),sep="\n")

saveRDS(von,"output/differential_detection/von.rds")
saveRDS(voff,"output/differential_detection/voff.rds")

cat("Done!\n")
