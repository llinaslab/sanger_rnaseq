old_to_new <- readr::read_tsv("../data/gene_lists/clean_pf_new_old.txt",col_names = c("new","old")) %>%
  dplyr::mutate(old=toupper(old))

# read in DeRisi data
read_caro_utrs <- function(file) {

  df <- tibble::as_tibble(rtracklayer::import.bed(file)) %>%
    dplyr::mutate(name=ifelse(stringr::str_count(name,"_") > 1,sub('(^[^_]+_[^_]+)_(.*)$', '\\1 \\2', name),stringr::str_replace(name,"_"," "))) %>%
    tidyr::separate(name, into=c("gene","type"), sep = " ") %>%
    dplyr::mutate(gene=toupper(gene)) %>%
    dplyr::mutate(type=ifelse(type=="5p","5UTR","3UTR")) %>%
    dplyr::inner_join(old_to_new,by=c("gene"="old"))

  return(df)
}

caro_utrs1 <- read_caro_utrs("../data/compare_utrs/GSM1410291_UTRs_1.bed")
caro_utrs2 <- read_caro_utrs("../data/compare_utrs/GSM1410292_UTRs_2.bed")
caro_utrs3 <- read_caro_utrs("../data/compare_utrs/GSM1410293_UTRs_3.bed")
caro_utrs4 <- read_caro_utrs("../data/compare_utrs/GSM1410294_UTRs_4.bed")
caro_utrs5 <- read_caro_utrs("../data/compare_utrs/GSM1410295_UTRs_5.bed")

# read in Sophie's data
sophie_tss <- tibble::as_tibble(rtracklayer::import.gff3("../data/compare_utrs/sorted_Adjalley_Chabbert_TSSs.gff"))

# read in our data
our_utrs <- tibble::as_tibble(rtracklayer::import.gff3("../output/final_utrs/longest_utrs_3d7_plasmodb_compatible.gff"))
our_utrs$Parent <- stringr::str_replace(stringr::str_replace(unlist(our_utrs$Parent), "rna_", ""), "-1", "")

our_tss <- tibble::as_tibble(rtracklayer::import.gff3("../output/ctss_clustering/modified/tag_clusters_annotated_intergenic.gff"))

# combine and compare
compare_utrs1 <- dplyr::inner_join(our_utrs, caro_utrs1, by = c("Parent"="new","type"="type"))
compare_utrs2 <- dplyr::inner_join(our_utrs, caro_utrs2, by = c("Parent"="new","type"="type"))
compare_utrs3 <- dplyr::inner_join(our_utrs, caro_utrs3, by = c("Parent"="new","type"="type"))
compare_utrs4 <- dplyr::inner_join(our_utrs, caro_utrs4, by = c("Parent"="new","type"="type"))
compare_utrs5 <- dplyr::inner_join(our_utrs, caro_utrs5, by = c("Parent"="new","type"="type"))
