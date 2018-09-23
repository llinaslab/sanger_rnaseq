
# import libraries needed
require(readr, quietly = T)
require(dplyr, quietly = T)
require(cowplot, quietly = T)
require(magrittr, quietly = T)
require(stringr, quietly = T)

# FUNCTIONS ###################################################################

# get filtered rpkms
#get_filtered_rpkms <- function(rpkms, norm_rpkms, min_rpkm) {
#
#  filtered_ids <- rpkms %>%
#    gather(tp, exp, -gene_id) %>%
#    group_by(gene_id) %>%
#    summarise(total = sum(exp)) %>%
#    filter(total > min_rpkm) %$%
#    gene_id
#
#  df <- norm_rpkms %>%
#    filter(gene_id %in% filtered_ids) %>%
#    gather(tp, exp, -gene_id) %>%
#    mutate(exp = log2(exp + 0.0005)) %>%
#    spread(tp, exp)
#
#  return(df)
#}

# get filtered tpms
get_filtered_tpms <- function(tpms, min_tpm) {

  filtered_ids <- tpms %>%
    dplyr::group_by(gene_id) %>%
    dplyr::summarise(total = sum(TPM)) %>%
    dplyr::filter(total > min_tpm) %$%
    gene_id

  norm_tpms <- tpms %>%
    dplyr::filter(gene_id %in% filtered_ids) %>%
    dplyr::group_by(gene_id) %>%
    dplyr::summarise(m=mean(TPM)) %>%
    dplyr::ungroup() %>%
    dplyr::inner_join(tpms) %>%
    dplyr::mutate(exp=log2((TPM/m)+0.005)) %>%
    dplyr::select(gene_id,tp,exp) %>%
    tidyr::spread(tp,exp)

  return(norm_tpms)
}

# read in counts file
read_counts <- function(motif_file) {

  df <- tibble::as_tibble(rtracklayer::import.gff3(motif_file)) %>% dplyr::select(seqnames, score, Name) %>%
    dplyr::rename(gene_id=seqnames,motif=Name) %>% dplyr::filter(score >= 0)
  return(df)

}

# get max scores
get_max_scores <- function(counts_df) {

  df <- counts_df %>%
    group_by(motif) %>%
    summarise(max_score = max(score))

  return(df)
}

# get hits
get_hits <- function(counts_df) {

  df <- counts_df %>%
    group_by(gene_id, motif) %>%
    summarise(hits = n())

  return(df)
}

# calculate scores per promoter per motif
calc_scores <- function(counts_df, max_df) {

  df <- inner_join(counts_df, max_df, by = "motif") %>%
    group_by(gene_id, motif) %>%
    summarise(score = sum(score) / max(max_score))

  return(df)
}

# plot scores distribution
plot_scores <- function(scores_df) {
  g <- ggplot(scores_df, aes(x = log2(score + 1))) +
    geom_histogram(color="grey90") +
    ggtitle("Distribution of Scores")
  return(g)
}

# get scores
get_scores <- function(file, count_filter) {

  # calculate scores
  counts <- read_counts(file)
  max    <- get_max_scores(counts)
  hits   <- get_hits(counts)
  scores <- calc_scores(counts, max)

  # remove motifs with less than count_filter hits
  remove <- hits %>%
    group_by(motif) %>%
    summarise(total = sum(hits)) %>%
    filter(total < count_filter) %$% motif

  if(length(remove > 0)) {
    cat("\n")
    message("Some motifs were below the count filter and will be removed from further analysis")
    cat("\n")
    for (m in remove) {
      cat(paste0(m, "\n"))
    }
    # filter for motifs that have less than the threshold
    scores <- scores %>%
      filter(!motif %in% remove)
  }

  # plot distribution of scores
  g <- plot_scores(scores)
  ggsave(plot = g, filename = paste0(str_replace(file, "fimo.gff", ""), "scores_hist.png"))

  # output to file
  scores <- spread(scores, motif, score, fill = 0)
  return(scores)

}

# MAIN FUNCTION ###############################################################

main <- function(file_list, count_filter, tpm_filter) {

  # filter rpkms
  td7 <- get_filtered_tpms(td7, tpm_filter)
  hb3 <- get_filtered_tpms(hb3, tpm_filter)
  it  <- get_filtered_tpms(it, tpm_filter)

  # go through each count file
  for (file in file_list) {

    scores <- get_scores(file, count_filter)

    if(str_count(file, "3d7")) {

      keep <- intersect(scores$gene_id, td7$gene_id)

      exp <- td7 %>%
        filter(gene_id %in% keep)

      scores <- scores %>%
        filter(gene_id %in% keep)

      write.table(x = exp,
                  file = "../output/tfbs_analysis/3d7_tpms.exp", quote = F, sep = "\t", row.names = F)
      write.table(x = scores,
                  file = "../output/tfbs_analysis/3d7_scores.mot", quote = F, sep = "\t", row.names = F)
    }
    else if(str_count(file, "hb3")) {

      keep <- intersect(scores$gene_id, hb3$gene_id)

      exp <- hb3 %>%
        filter(gene_id %in% keep)

      scores <- scores %>%
        filter(gene_id %in% keep)

      write.table(x = exp,
                  file = "../output/tfbs_analysis/hb3_tpms.exp", quote = F, sep = "\t", row.names = F)
      write.table(x = scores,
                  file = "../output/tfbs_analysis/hb3_scores.mot", quote = F, sep = "\t", row.names = F)
    }
    else if(str_count(file, "it")) {

      keep <- intersect(scores$gene_id, it$gene_id)

      exp <- it %>%
        filter(gene_id %in% keep)

      scores <- scores %>%
        filter(gene_id %in% keep)

      write.table(x = exp,
                  file = "../output/tfbs_analysis/it_tpms.exp", quote = F, sep = "\t", row.names = F)
      write.table(x = scores,
                  file = "../output/tfbs_analysis/it_scores.mot", quote = F, sep = "\t", row.names = F)
    }

  }

}

# TESTING #####################################################################

# generate toy data

#

# CALL MAIN ###################################################################

# import expression data
#exp_dir  <- "data/rpkms/sense_mapped_to_3d7"
#td7      <- read_tsv(paste0(exp_dir,"/3D7_stranded_sense_rpkm_all_tps.txt"))
#hb3      <- read_tsv(paste0(exp_dir,"/HB3_stranded_sense_rpkm_all_tps.txt"))
#it       <- read_tsv(paste0(exp_dir,"/IT_stranded_sense_rpkms_all_tps.txt"))
#td7_norm <- read_tsv(paste0(exp_dir,"/3D7_globnorm_stranded_sense_rpkm_all_tps.txt"))
#hb3_norm <- read_tsv(paste0(exp_dir,"/HB3_globnorm_stranded_sense_rpkm_all_tps.txt"))
#it_norm  <- read_tsv(paste0(exp_dir,"/IT_globnorm_stranded_sense_rpkms_all_tps.txt"))

exp_dir  <- "../output/neighboring_genes/"
td7      <- readRDS(paste0(exp_dir,"gene_reduced_3d7_abund.rds"))
hb3      <- readRDS(paste0(exp_dir,"gene_reduced_hb3_abund.rds"))
it       <- readRDS(paste0(exp_dir,"gene_reduced_it_abund.rds"))

# create file lists
#utr_files <- list.files(path = "data/promoters/utrs", pattern = "*.counts", full.names = T)
#atg_files <- list.files(path = "data/promoters/atg", pattern = "*.counts", full.names = T)
#all_files <- c(utr_files, atg_files)

files <- list.files(path = "../output/tfbs_analysis", pattern = "fimo.gff", full.names = T, recursive = T)

# call main
main(file_list = files, count_filter = 100, tpm_filter = 5)
