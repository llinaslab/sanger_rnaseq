df1 <- read_tsv("../data/utrs/homopolymer_analysis/minus_100bp_slop_3d7_5utrs_counts.tsv") # UTR plus 100bp upstream
df2 <- read_tsv("../data/utrs/homopolymer_analysis/minus_100bp_flank_3d7_5utrs_counts.tsv") # only 100 bp upstream of UTR
df3 <- read_tsv("../data/utrs/homopolymer_analysis/3d7_5utrs_counts.tsv") # only UTR

utrs <- as.data.frame(import.gff3("../data/utrs/original_utrs/final.5utrs.3d7.3d7_v3_chr.idc.gff"))
utrs$Parent <- as.character(utrs$Parent)
tmp <- inner_join(utrs, df1, by = c("Parent" = "gene_id"))

tmp %>%
  filter(width < 8000) %>%
  group_by(Parent) %>%
  summarise(width=unique(width),longest=max(longest)) %>%
  ggplot(aes(y=width/longest,x=longest)) + geom_point(alpha=0.5) + stat_smooth(method="lm")

tmp %>%
  filter(width < 8000) %>%
  group_by(Parent) %>%
  summarise(width=unique(width),longest=max(longest)) %>%
  do(tidy(cor.test(.$width/.$longest,.$longest))) %>% View

tmp %>%
  filter(width < 8000 & nuc == "T") %>%
  ggplot(aes(y=width/longest,x=longest)) + geom_point(alpha=0.5) + stat_smooth(method="lm")

tmp %>%
  filter(width < 8000 & nuc == "A") %>%
  do(tidy(cor.test(.$width/.$longest,.$longest))) %>% View
