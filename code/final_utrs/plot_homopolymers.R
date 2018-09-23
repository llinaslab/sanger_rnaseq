df1 <- read_tsv("data/homopolymer_analysis/minus_100bp_slop_3d7_5utrs_counts.tsv")
df2 <- read_tsv("data/homopolymer_analysis/minus_100bp_flank_3d7_5utrs_counts.tsv")
df3 <- read_tsv("data/homopolymer_analysis/3d7_5utrs_counts.tsv")

utrs <- as.data.frame(import.gff3("data/original_utrs/final.5utrs.3d7.3d7_v3_chr.idc.gff"))
utrs$Parent <- as.character(utrs$Parent)
tmp <- inner_join(utrs, df1, by = c("Parent" = "gene_id")) 

tmp %>% filter(width < 8000) %>% group_by(Parent) %>% summarise(width=unique(width),longest=max(longest)) %>% ggplot(aes(y=width/longest,x=longest)) + geom_point(alpha=0.5) + stat_smooth(method="lm")

tmp %>% filter(width < 8000) %>% group_by(Parent) %>% summarise(width=unique(width),longest=max(longest)) %>% do(tidy(cor.test(.$width/.$longest,.$longest))) %>% View
