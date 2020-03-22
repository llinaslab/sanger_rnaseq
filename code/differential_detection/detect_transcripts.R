source("../code/utils.R")
load_essentials()
sshhh("limma")
sshhh("org.Pf.plasmo.db")

# Main --------------------------------------------------------------------

# import gene expression
x3d7exp <- readRDS("../output/neighboring_genes/gene_reduced_3d7_abund.rds")
xhb3exp <- readRDS("../output/neighboring_genes/gene_reduced_hb3_abund.rds")
xitexp  <- readRDS("../output/neighboring_genes/gene_reduced_it_abund.rds")

exp <- dplyr::bind_rows(x3d7exp,xhb3exp,xitexp)

rm(x3d7exp,xhb3exp,xitexp)

# import core genes for comparison
core_genes <- readr::read_tsv("../data/gene_lists/core_pf3d7_genes.txt",col_names=F)$X1

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
on3d7 <- detected %>% dplyr::filter(strain == "3d7" & on == 1 & gene_id %in% core_genes) %$% gene_id
onhb3 <- detected %>% dplyr::filter(strain == "hb3" & on == 1 & gene_id %in% core_genes) %$% gene_id
onit  <- detected %>% dplyr::filter(strain == "it" & on == 1 & gene_id %in% core_genes) %$% gene_id

off3d7 <- detected %>% dplyr::filter(strain == "3d7" & on == 0 & gene_id %in% core_genes) %$% gene_id
offhb3 <- detected %>% dplyr::filter(strain == "hb3" & on == 0 & gene_id %in% core_genes) %$% gene_id
offit  <- detected %>% dplyr::filter(strain == "it" & on == 0 & gene_id %in% core_genes) %$% gene_id

X3d7 <- (ug %in% on3d7)
Xhb3 <- (ug %in% onhb3)
Xit  <- (ug %in% onit)
von  <- vennCounts(cbind(X3d7, Xhb3, Xit))

X3d7 <- (ug %in% off3d7)
Xhb3 <- (ug %in% offhb3)
Xit  <- (ug %in% offit)
voff  <- limma::vennCounts(cbind(X3d7, Xhb3, Xit))

# Detect genes specific to strains ----------------------------------------

# Off genes per strain
off_only_3d7 <- detected %>%
  dplyr::filter(strain == "3d7" & on == 0 & gene_id %nin% offhb3 & gene_id %nin% offit) %$%
  gene_id
off_only_hb3 <- detected %>%
  dplyr::filter(strain == "hb3" & on == 0 & gene_id %nin% off3d7 & gene_id %nin% offit) %$%
  gene_id
off_only_it <- detected %>%
  dplyr::filter(strain == "it" & on == 0 & gene_id %nin% offhb3 & gene_id %nin% off3d7) %$%
  gene_id

# On genes per strain
on_only_3d7 <- detected %>%
  dplyr::filter(strain == "3d7" & on == 1 & gene_id %nin% onhb3 & gene_id %nin% onit) %$%
  gene_id
on_only_hb3 <- detected %>%
  dplyr::filter(strain == "hb3" & on == 1 & gene_id %nin% on3d7 & gene_id %nin% onit) %$%
  gene_id
on_only_it <- detected %>%
  dplyr::filter(strain == "it" & on == 1 & gene_id %nin% onhb3 & gene_id %nin% on3d7) %$%
  gene_id


# Annotate lists ----------------------------------------------------------

# import gene names
gene_names <- as.data.frame(org.Pf.plasmoGENENAME)

# annotate with total TPMs
off_3d7hb3_exp <- exp %>%
  dplyr::filter(gene_id %in% union(off3d7, offhb3) & (strain == "3d7" | strain == "hb3")) %>%
  dplyr::group_by(gene_id, strain) %>%
  dplyr::summarise(min = min(TPM), max = max(TPM), mean = mean(TPM), sd = sd(TPM))

off_3d7it_exp <- exp %>%
  dplyr::filter(gene_id %in% union(off3d7, offit) & (strain == "3d7" | strain == "it")) %>%
  dplyr::group_by(gene_id, strain) %>%
  dplyr::summarise(min = min(TPM), max = max(TPM), mean = mean(TPM), sd = sd(TPM))

off_hb3it_exp <- exp %>%
  dplyr::filter(gene_id %in% union(offhb3, offit) & (strain == "hb3" | strain == "it")) %>%
  dplyr::group_by(gene_id, strain) %>%
  dplyr::summarise(min = min(TPM), max = max(TPM), mean = mean(TPM), sd = sd(TPM))

# annotate the off genes
off_3d7hb3_exp <- dplyr::inner_join(off_3d7hb3_exp, gene_names, by = "gene_id")
off_3d7it_exp  <- dplyr::inner_join(off_3d7it_exp, gene_names, by = "gene_id")
off_hb3it_exp  <- dplyr::inner_join(off_hb3it_exp, gene_names, by = "gene_id")


# Write data  -------------------------------------------------------------

cat("\nWriting data to output... ")

output <- "../output/differential_detection/"
drive_output <- "Shared/Pf RNA-seq manuscript 2017/Supplementary tables/Named tables"
drive_type <- "spreadsheet"

readr::write_tsv(x = off_3d7hb3_exp,
          path = paste0(output,"off_genes_3d7_hb3.tsv"))
readr::write_tsv(x = off_3d7it_exp,
          path = paste0(output,"off_genes_3d7_it.tsv"))
readr::write_tsv(x = off_hb3it_exp,
          path = paste0(output,"off_genes_hb3_it.tsv"))
readr::write_tsv(x = dplyr::inner_join(detected, gene_names, by = "gene_id"),
          path = paste0(output,"detected.tsv"))

googledrive::drive_upload(media=paste0(output,"detected.tsv"),
                          path=drive_output,
                          name="detected",
                          type=drive_type)

off_3d7hb3_exp %>%
  dplyr::filter((strain == "hb3" & max < 5) | (strain == "3d7" & max > 5)) %>%
  dplyr::group_by(gene_id) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::inner_join(off_3d7hb3_exp) %>%
  dplyr::filter(n>1) %>%
  dplyr::select(-n) %>%
  dplyr::arrange(gene_id) %>%
  readr::write_tsv(path=paste0(output,"off_hb3_not_3d7.tsv"))

googledrive::drive_upload(media=paste0(output,"off_hb3_not_3d7.tsv"),
                          path=drive_output,
                          name="off_hb3_not_3d7",
                          type=drive_type)

off_3d7hb3_exp %>%
  dplyr::filter((strain == "3d7" & max < 5) | (strain == "hb3" & max > 5)) %>%
  dplyr::group_by(gene_id) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::inner_join(off_3d7hb3_exp) %>%
  dplyr::filter(n>1) %>%
  dplyr::select(-n) %>%
  dplyr::arrange(gene_id) %>%
  readr::write_tsv(path=paste0(output,"off_3d7_not_hb3.tsv"))

googledrive::drive_upload(media=paste0(output,"off_3d7_not_hb3.tsv"),
                          path=drive_output,
                          name="off_3d7_not_hb3",
                          type=drive_type)

off_3d7it_exp %>%
  dplyr::filter((strain == "it" & max < 5) | (strain == "3d7" & max > 5)) %>%
  dplyr::group_by(gene_id) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::inner_join(off_3d7hb3_exp) %>%
  dplyr::filter(n>1) %>%
  dplyr::select(-n) %>%
  dplyr::arrange(gene_id) %>%
  readr::write_tsv(path=paste0(output,"off_it_not_3d7.tsv"))

googledrive::drive_upload(media=paste0(output,"off_it_not_3d7.tsv"),
                          path=drive_output,
                          name="off_it_not_3d7",
                          type=drive_type)

off_3d7it_exp %>%
  dplyr::filter((strain == "3d7" & max < 5) | (strain == "it" & max > 5)) %>%
  dplyr::group_by(gene_id) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::inner_join(off_3d7hb3_exp) %>%
  dplyr::filter(n>1) %>%
  dplyr::select(-n) %>%
  dplyr::arrange(gene_id) %>%
  readr::write_tsv(path=paste0(output,"off_3d7_not_it.tsv"))

googledrive::drive_upload(media=paste0(output,"off_3d7_not_it.tsv"),
                          path=drive_output,
                          name="off_3d7_not_it",
                          type=drive_type)

off_hb3it_exp %>%
  dplyr::filter((strain == "hb3" & max < 5) | (strain == "it" & max > 5)) %>%
  dplyr::group_by(gene_id) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::inner_join(off_3d7hb3_exp) %>%
  dplyr::filter(n>1) %>%
  dplyr::select(-n) %>%
  dplyr::arrange(gene_id) %>%
  readr::write_tsv(path=paste0(output,"off_hb3_not_it.tsv"))

googledrive::drive_upload(media=paste0(output,"off_hb3_not_it.tsv"),
                          path=drive_output,
                          name="off_hb3_not_it",
                          type=drive_type)

off_hb3it_exp %>%
  dplyr::filter((strain == "it" & max < 5) | (strain == "hb3" & max > 5)) %>%
  dplyr::group_by(gene_id) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::inner_join(off_3d7hb3_exp) %>%
  dplyr::filter(n>1) %>%
  dplyr::select(-n) %>%
  dplyr::arrange(gene_id) %>%
  readr::write_tsv(path=paste0(output,"off_it_not_hb3.tsv"))

googledrive::drive_upload(media=paste0(output,"off_it_not_hb3.tsv"),
                          path=drive_output,
                          name="off_it_not_hb3",
                          type=drive_type)

off_hb3_not_3d7_list <- off_3d7hb3_exp %>% filter(strain == "hb3" & max < 5) %$% gene_id
off_3d7_not_hb3_list <- off_3d7hb3_exp %>% filter(strain == "3d7" & max < 5) %$% gene_id
off_it_not_3d7_list  <- off_3d7it_exp %>% filter(strain == "it" & max < 5) %$% gene_id
off_3d7_not_it_list  <- off_3d7it_exp %>% filter(strain == "3d7" & max < 5) %$% gene_id
off_hb3_not_it_list  <- off_hb3it_exp %>% filter(strain == "hb3" & max < 5) %$% gene_id
off_it_not_hb3_list  <- off_hb3it_exp %>% filter(strain == "it" & max < 5) %$% gene_id

write(off_hb3_not_3d7_list,paste0(output,"off_genes_3d7_hb3/off_hb3_list.txt"),sep="\n")
write(off_3d7_not_hb3_list,paste0(output,"off_genes_3d7_hb3/off_3d7_list.txt"),sep="\n")
write(off_it_not_3d7_list,paste0(output,"off_genes_3d7_it/off_it_list.txt"),sep="\n")
write(off_3d7_not_it_list,paste0(output,"off_genes_3d7_it/off_3d7_list.txt"),sep="\n")
write(off_hb3_not_it_list,paste0(output,"off_genes_hb3_it/off_hb3_list.txt"),sep="\n")
write(off_it_not_hb3_list,paste0(output,"off_genes_hb3_it/off_it_list.txt"),sep="\n")

write(on_only_3d7,paste0(output,"on_only_3d7_list.txt"),sep="\n")
write(on_only_hb3,paste0(output,"on_only_hb3_list.txt"),sep="\n")
write(on_only_it,paste0(output,"on_only_it_list.txt"),sep="\n")

write(off_only_3d7,paste0(output,"off_only_3d7_list.txt"),sep="\n")
write(off_only_hb3,paste0(output,"off_only_hb3_list.txt"),sep="\n")
write(off_only_it,paste0(output,"off_only_it_list.txt"),sep="\n")

saveRDS(von,"../output/differential_detection/von.rds")
saveRDS(voff,"../output/differential_detection/voff.rds")

cat("Done!\n")
