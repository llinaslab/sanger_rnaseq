# What does the correlation between annotated TSSs and downstream target genes look like?

# import tag clusters
tc_intergenic <- rtracklayer::import.gff3("../output/ctss_clustering/modified/tag_clusters_annotated_intergenic.gff")
tc_exonic     <- rtracklayer::import.gff3("../output/ctss_clustering/modified/tag_clusters_annotated_exons.gff")
tc_intronic   <- rtracklayer::import.gff3("../output/ctss_clustering/modified/tag_clusters_annotated_introns.gff")

# import promoter clusters
pc_intergenic <- rtracklayer::import.gff3("../output/ctss_clustering/modified/promoter_clusters_annotated_intergenic.gff")
pc_exonic     <- rtracklayer::import.gff3("../output/ctss_clustering/modified/promoter_clusters_annotated_exons.gff")
pc_intronic   <- rtracklayer::import.gff3("../output/ctss_clustering/modified/promoter_clusters_annotated_introns.gff")

# abundance estimates
x3d7_abund <- readRDS("../output/neighboring_genes/gene_reduced_3d7_abund.rds")
xhb3_abund <- readRDS("../output/neighboring_genes/gene_reduced_hb3_abund.rds")
xit_abund  <- readRDS("../output/neighboring_genes/gene_reduced_it_abund.rds")

pcg <- tibble::as_tibble(rtracklayer::import.gff3("../data/annotations/PF3D7_codinggenes_for_bedtools.gff"))$ID
get_filtered_ids <- function(abund,tpm_threshold) {
  fabund <- abund %>%
    dplyr::group_by(gene_id) %>%
    dplyr::summarise(f=sum(TPM>=tpm_threshold)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(f>0 & gene_id %in% pcg)
  return(fabund$gene_id)
}

fx3d7 <- get_filtered_ids(x3d7_abund,5)
fxhb3 <- get_filtered_ids(xhb3_abund,5)
fxit  <- get_filtered_ids(xit_abund,5)

tmp <- dplyr::inner_join(as_tibble(intertags) %>% dplyr::mutate(tp=as.integer(tp)),x3d7_abund,by=c("name" = "gene_id","tp"))
filter(tmp, tp==4 & name %in% fx3d7 & as.numeric(tpm.dominant_ctss) >= 2 & TPM >= 5 & (as.numeric(anno_start) - as.numeric(dominant_ctss) <= 2500)) %>%
  mutate(tpm.dominant_ctss=as.numeric(tpm.dominant_ctss)) %>%
  ggplot(aes(x=log2(tpm.dominant_ctss+1),y=log2(TPM+1))) +
  geom_point() +
  stat_smooth()

################################################################################

# Figure 3 Panel A

intertags <- tibble::as_tibble(rtracklayer::import.gff3("../output/ctss_clustering/modified/tag_clusters_annotated_intergenic.gff")) %>%
  group_by(seqnames,start,end,name)

# abundance estimates
x3d7_abund <- readRDS("../output/neighboring_genes/gene_reduced_3d7_abund.rds")
xhb3_abund <- readRDS("../output/neighboring_genes/gene_reduced_hb3_abund.rds")
xit_abund  <- readRDS("../output/neighboring_genes/gene_reduced_it_abund.rds")

pcg <- tibble::as_tibble(rtracklayer::import.gff3("../data/annotations/PF3D7_codinggenes_for_bedtools.gff"))$ID
get_filtered_ids <- function(abund,tpm_threshold) {
  fabund <- abund %>%
    dplyr::group_by(gene_id) %>%
    dplyr::summarise(f=sum(TPM>=tpm_threshold)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(f>0 & gene_id %in% pcg)
  return(fabund$gene_id)
}

fx3d7 <- get_filtered_ids(x3d7_abund,5)
fxhb3 <- get_filtered_ids(xhb3_abund,5)
fxit  <- get_filtered_ids(xit_abund,5)

intertags <- dplyr::inner_join(as_tibble(intertags) %>% dplyr::mutate(tp=as.integer(tp)),x3d7_abund,by=c("name" = "gene_id","tp"))

intertags %>%
  filter(as.numeric(tpm.dominant_ctss) >= 2 & TPM >= 5 & abs(as.numeric(anno_start) - as.numeric(dominant_ctss)) <= 2500) %>%
  group_by(tp,name) %>%
  summarise(n=n()) %>%
  mutate(type=ifelse(n>1,"multi","single")) %>%
  ggplot(aes(x=tp,fill=type)) +
  geom_bar(colour="black",size=1.5) +
  scale_fill_brewer(palette="Greys") +
  ylab("Frequency") +
  xlab("Timepoint") +
  labs(fill="") +
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=28,face="bold"),
        axis.line.x=element_line(colour="black",size=1.5),
        axis.ticks.x=element_line(colour="black",size=1.5),
        axis.line.y=element_line(colour="black",size=1.5),
        axis.ticks.y=element_line(colour="black",size=1.5),
        legend.text=element_text(size=24))

wheretags <- tibble::tibble(var="", type=c(rep("Intergenic",21946),rep("Exon",2082),rep("Intron",241)))

wheretags %>%
  ggplot(aes(x=var,fill=type)) +
  geom_bar(colour="black",size=1.5,width=0.2) +
  scale_fill_brewer(palette="Greys") +
  ylab("Frequency") +
  xlab("") +
  labs(fill="") +
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=28,face="bold"),
        axis.line.x=element_line(colour="black",size=1.5),
        axis.ticks.x=element_blank(),
        axis.line.y=element_line(colour="black",size=1.5),
        axis.ticks.y=element_line(colour="black",size=1.5),
        legend.text=element_text(size=24))

# Figure 3 Panel B

motifs_3d7 <- rtracklayer::import.gff3("../output/tfbs_analysis/promoters_3d7/fimo.gff") %>% tibble::as_tibble()

for (name in unique(motifs_3d7$Name)) {
  n <- motifs_3d7 %>%
    dplyr::filter(Name==name)
  g <- n %>% ggplot(aes(x=start)) + geom_histogram(color="grey70") + xlab("Start Site") + ylab("Frequency") + ggtitle(name) +
    theme(axis.text=element_text(size=28),
          axis.title=element_text(size=32,face="bold"),
          axis.line.x=element_line(colour="black",size=1.5),
          axis.ticks.x=element_line(colour="black",size=1.5),
          axis.line.y=element_line(colour="black",size=1.5),
          axis.ticks.y=element_line(colour="black",size=1.5),
          legend.text=element_text(size=24),
          plot.title=element_text(size=34,face="bold"))
  ggsave(plot=g,filename=paste0("/Users/philippross/Dropbox/projects/figure3/promoters/",name,".svg"))
}


# Figure 3 Panel C

bimotifs_3d7 <- rtracklayer::import.gff3("../output/tfbs_analysis/bidirectional_3d7/fimo.gff") %>% tibble::as_tibble()

df <- tibble::as_tibble(biprom_3d7)
newdf <- tibble::tibble(name=character(),norm_start=numeric())

for (m in unique(bimotifs_3d7$seqnames)) {
  w <- dplyr::filter(df, ID == m)$width
  newdf <- dplyr::bind_rows(newdf, dplyr::filter(bimotifs_3d7, seqnames == m) %>% dplyr::transmute(name=seqnames, norm_start = start / w, motif=Name))
}

for (n in unique(newdf$motif)) {
  b <- newdf %>% dplyr::filter(motif==n)
  g <- b %>% ggplot(aes(x=norm_start)) + geom_line(stat="density",size=1.5) + ggtitle(n) + ylab("Density") + xlab("Normalized Start") +
    theme(axis.text=element_text(size=28),
          axis.title=element_text(size=32,face="bold"),
          axis.line.x=element_line(colour="black",size=1.5),
          axis.ticks.x=element_line(colour="black",size=1.5),
          axis.line.y=element_line(colour="black",size=1.5),
          axis.ticks.y=element_line(colour="black",size=1.5),
          legend.text=element_text(size=24),
          plot.title=element_text(size=34,face="bold"))
  ggsave(plot=g,filename=paste0("/Users/philippross/Dropbox/projects/figure3/bidirectional/",n,".svg"))
}

# Figure 3 Panel D

# This isn't programmable

################################################################################

