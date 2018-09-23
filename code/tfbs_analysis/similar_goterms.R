tibble::as_tibble(rtracklayer::import.gff3("../output/tfbs_analysis/bidirectional_3d7.gff")) %>%
  tidyr::separate(ID,into=c("left","right"),sep="-")

go <- readr::read_tsv("../data/gene_ontology/Pf3D7_go_sept2014.txt",col_names=c("gene","go"))

#  mutate(go = strsplit(go, ",")) %>%
#  unnest(go) %>%
#  group_by(gene) %>%
#  mutate(row = row_number()) %>%
#  spread(row, go)

tibble::as_tibble(rtracklayer::import.gff3("../output/ctss_clustering/modified/tag_clusters_annotated_intergenic.gff")) %>%
  group_by(seqnames,start,end,name)
