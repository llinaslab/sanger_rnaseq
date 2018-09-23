# read in data
# import genome
library(BSgenome.Pfalciparum.PlasmoDB.v24)
library(tidyverse)

# import tag clusters
tc_intergenic <- rtracklayer::import.gff3("../output/ctss_clustering/modified/tag_clusters_annotated_intergenic.gff")
tc_exonic     <- rtracklayer::import.gff3("../output/ctss_clustering/modified/tag_clusters_annotated_exons.gff")
tc_intronic   <- rtracklayer::import.gff3("../output/ctss_clustering/modified/tag_clusters_annotated_introns.gff")

# import promoter clusters
pc_intergenic <- rtracklayer::import.gff3("../output/ctss_clustering/modified/promoter_clusters_annotated_intergenic.gff")
pc_exonic     <- rtracklayer::import.gff3("../output/ctss_clustering/modified/promoter_clusters_annotated_exons.gff")
pc_intronic   <- rtracklayer::import.gff3("../output/ctss_clustering/modified/promoter_clusters_annotated_introns.gff")

########################## for tag clusters

# intergenic
for (i in seq(1,7)) {
  assign(x=paste0("tc_intergenic",i),
         value=as_tibble(tc_intergenic) %>% filter(tp==i) %>% mutate(end=full_end,source="CAGEr",type="tag_cluster") %>% select(-tp,-full_end,-anno_start,-anno_end))
  rtracklayer::export.gff3(GenomicRanges::GRanges(eval(parse(text=paste0("tc_intergenic",i)))),con = paste0("~/Dropbox/people/manuel/tss_data/tso/","tc_intergenic",i,".gff"))
  assign(x=paste0("tc_intergenic",i,"table"),
         value=as_tibble(tc_intergenic) %>%
           dplyr::filter(tp==i) %>%
           dplyr::mutate(start=ifelse(strand=="+",as.numeric(dominant_ctss)-50,as.numeric(dominant_ctss)-48),
                  end=ifelse(strand=="+",as.numeric(dominant_ctss)+50,as.numeric(dominant_ctss)+52)) %>%
           dplyr::select(-tp,-full_end) %>%
           GenomicRanges::GRanges())
  seqs <- tibble(seq=BSgenome::getSeq(BSgenome.Pfalciparum.PlasmoDB.v24,eval(parse(text=paste0("tc_intergenic",i,"table"))),as.character=TRUE))
  assign(x=paste0("tc_intergenic",i,"table"),
         value=dplyr::bind_cols(as_tibble(eval(parse(text=paste0("tc_intergenic",i,"table")))),seqs) %>%
           dplyr::mutate(dist=as.integer(abs(as.numeric(anno_start)-as.numeric(dominant_ctss)))) %>%
           dplyr::select(name,dist,interquantile_width,tpm,tpm.dominant_ctss,seq))
  write_tsv(x=eval(parse(text=paste0("tc_intergenic",i,"table"))),path = paste0("~/Dropbox/people/manuel/tss_data/tso/","tc_intergenic",i,".tsv"))
}

# exonic
for (i in seq(1,7)) {
  assign(x=paste0("tc_exonic",i),
         value=as_tibble(tc_exonic) %>% filter(tp==i) %>% mutate(end=full_end,source="CAGEr",type="tag_cluster") %>% select(-tp,-full_end,-anno_start,-anno_end))
  rtracklayer::export.gff3(GenomicRanges::GRanges(eval(parse(text=paste0("tc_exonic",i)))),con = paste0("~/Dropbox/people/manuel/tss_data/tso/","tc_exonic",i,".gff"))
  assign(x=paste0("tc_exonic",i,"table"),
         value=as_tibble(tc_exonic) %>%
           dplyr::filter(tp==i) %>%
           dplyr::mutate(start=ifelse(strand=="+",as.numeric(dominant_ctss)-50,as.numeric(dominant_ctss)-48),
                         end=ifelse(strand=="+",as.numeric(dominant_ctss)+50,as.numeric(dominant_ctss)+52)) %>%
           dplyr::select(-tp,-full_end) %>%
           GenomicRanges::GRanges())
  seqs <- tibble(seq=BSgenome::getSeq(BSgenome.Pfalciparum.PlasmoDB.v24,eval(parse(text=paste0("tc_exonic",i,"table"))),as.character=TRUE))
  assign(x=paste0("tc_exonic",i,"table"),
         value=dplyr::bind_cols(as_tibble(eval(parse(text=paste0("tc_exonic",i,"table")))),seqs) %>%
           dplyr::mutate(dist=as.integer(abs(as.numeric(anno_start)-as.numeric(dominant_ctss)))) %>%
           dplyr::select(name,dominant_ctss,interquantile_width,tpm,tpm.dominant_ctss,seq))
  write_tsv(x=eval(parse(text=paste0("tc_exonic",i,"table"))),path = paste0("~/Dropbox/people/manuel/tss_data/tso/","tc_exonic",i,".tsv"))
}

# intronic
for (i in seq(1,7)) {
  assign(x=paste0("tc_intronic",i),
         value=as_tibble(tc_intronic) %>% filter(tp==i) %>% mutate(end=full_end,source="CAGEr",type="tag_cluster") %>% select(-tp,-full_end,-anno_start,-anno_end))
  rtracklayer::export.gff3(GenomicRanges::GRanges(eval(parse(text=paste0("tc_intronic",i)))),con = paste0("~/Dropbox/people/manuel/tss_data/tso/","tc_intronic",i,".gff"))
  assign(x=paste0("tc_intronic",i,"table"),
         value=as_tibble(tc_intronic) %>%
           dplyr::filter(tp==i) %>%
           dplyr::mutate(start=ifelse(strand=="+",as.numeric(dominant_ctss)-50,as.numeric(dominant_ctss)-48),
                         end=ifelse(strand=="+",as.numeric(dominant_ctss)+50,as.numeric(dominant_ctss)+52)) %>%
           dplyr::select(-tp,-full_end) %>%
           GenomicRanges::GRanges())
  seqs <- tibble(seq=BSgenome::getSeq(BSgenome.Pfalciparum.PlasmoDB.v24,eval(parse(text=paste0("tc_intronic",i,"table"))),as.character=TRUE))
  assign(x=paste0("tc_intronic",i,"table"),
         value=dplyr::bind_cols(as_tibble(eval(parse(text=paste0("tc_intronic",i,"table")))),seqs) %>%
           dplyr::mutate(dist=as.integer(abs(as.numeric(anno_start)-as.numeric(dominant_ctss)))) %>%
           dplyr::select(name,dominant_ctss,interquantile_width,tpm,tpm.dominant_ctss,seq))
  write_tsv(x=eval(parse(text=paste0("tc_intronic",i,"table"))),path = paste0("~/Dropbox/people/manuel/tss_data/tso/","tc_intronic",i,".tsv"))
}

########################## for promoter clusters

# first create list of distinct promoter clusters and their attributes
# intergenic
fpcs_intregenic <- pc_intergenic %>%
  tibble::as_tibble() %>%
  dplyr::group_by(seqnames,start,full_end,strand,name,dominant_ctss,anno_start) %>%
  dplyr::summarise(tpm=sum(as.numeric(tpm)),tpm_dominant_ctss=sum(as.numeric(tpm.dominant_ctss))) %>%
  dplyr::ungroup()

tmp <- fpcs_intregenic %>%
  dplyr::group_by(seqnames,start,strand,name) %>%
  dplyr::summarise(tpm_dominant_ctss=max(tpm_dominant_ctss)) %>%
  dplyr::ungroup()

fpcs_intregenic <- dplyr::inner_join(fpcs_intregenic,tmp,by=c("seqnames","start","strand","name","tpm_dominant_ctss"))

# exonic
fpcs_exonic <- pc_exonic %>%
  tibble::as_tibble() %>%
  dplyr::group_by(seqnames,start,full_end,strand,name,dominant_ctss,anno_start) %>%
  dplyr::summarise(tpm=sum(as.numeric(tpm)),tpm_dominant_ctss=sum(as.numeric(tpm.dominant_ctss))) %>%
  dplyr::ungroup()

tmp <- fpcs_exonic %>%
  dplyr::group_by(seqnames,start,strand,name) %>%
  dplyr::summarise(tpm_dominant_ctss=max(tpm_dominant_ctss)) %>%
  dplyr::ungroup()

fpcs_exonic <- dplyr::inner_join(fpcs_exonic,tmp,by=c("seqnames","start","strand","name","tpm_dominant_ctss"))

# intronic
fpcs_intronic <- pc_intronic %>%
  tibble::as_tibble() %>%
  dplyr::group_by(seqnames,start,full_end,strand,name,dominant_ctss,anno_start) %>%
  dplyr::summarise(tpm=sum(as.numeric(tpm)),tpm_dominant_ctss=sum(as.numeric(tpm.dominant_ctss))) %>%
  dplyr::ungroup()

tmp <- fpcs_intronic %>%
  dplyr::group_by(seqnames,start,strand,name) %>%
  dplyr::summarise(tpm_dominant_ctss=max(tpm_dominant_ctss)) %>%
  dplyr::ungroup()

fpcs_intronic <- dplyr::inner_join(fpcs_intronic,tmp,by=c("seqnames","start","strand","name","tpm_dominant_ctss"))

# Now write those to a file
# intergenic
exp <- fpcs_intregenic %>%
  dplyr::mutate(end=as.numeric(full_end),source="CAGEr",type="promoter_cluster") %>%
  dplyr::select(-full_end)
rtracklayer::export.gff3(GenomicRanges::GRanges(exp),con = paste0("~/Dropbox/people/manuel/tss_data/tso/","pc_intergenic.gff"))

pc_intergenic_table <- fpcs_intregenic %>%
  dplyr::mutate(full_start=as.integer(start),
                full_end=as.integer(full_end),
                start=ifelse(strand=="+",as.numeric(dominant_ctss)-50,as.numeric(dominant_ctss)-48),
                end=ifelse(strand=="+",as.numeric(dominant_ctss)+50,as.numeric(dominant_ctss)+52)) %>%
  GenomicRanges::GRanges()

seqs <- tibble(seq=BSgenome::getSeq(BSgenome.Pfalciparum.PlasmoDB.v24,pc_intergenic_table,as.character=TRUE))
pc_intergenic_table <- dplyr::bind_cols(tibble::as_tibble(pc_intergenic_table),seqs) %>%
  dplyr::mutate(dist=as.integer(abs(as.numeric(anno_start)-as.numeric(dominant_ctss))),width=full_end-full_start) %>%
  dplyr::select(name,dist,width,tpm,tpm_dominant_ctss,seq)

write_tsv(x=pc_intergenic_table,path = paste0("~/Dropbox/people/manuel/tss_data/tso/","pc_intergenic.tsv"))

# exonic
exp <- fpcs_exonic %>%
  dplyr::mutate(end=as.numeric(full_end),source="CAGEr",type="promoter_cluster") %>%
  dplyr::select(-full_end)
rtracklayer::export.gff3(GenomicRanges::GRanges(exp),con = paste0("~/Dropbox/people/manuel/tss_data/tso/","pc_exonic.gff"))

pc_exonic_table <- fpcs_exonic %>%
  dplyr::mutate(full_start=as.integer(start),
                full_end=as.integer(full_end),
                start=ifelse(strand=="+",as.numeric(dominant_ctss)-50,as.numeric(dominant_ctss)-48),
                end=ifelse(strand=="+",as.numeric(dominant_ctss)+50,as.numeric(dominant_ctss)+52)) %>%
  GenomicRanges::GRanges()

seqs <- tibble(seq=BSgenome::getSeq(BSgenome.Pfalciparum.PlasmoDB.v24,pc_exonic_table,as.character=TRUE))
pc_exonic_table <- dplyr::bind_cols(tibble::as_tibble(pc_exonic_table),seqs) %>%
  dplyr::mutate(dist=as.integer(abs(as.numeric(anno_start)-as.numeric(dominant_ctss))),width=full_end-full_start) %>%
  dplyr::select(name,dominant_ctss,width,tpm,tpm_dominant_ctss,seq)

write_tsv(x=pc_exonic_table,path = paste0("~/Dropbox/people/manuel/tss_data/tso/","pc_exonic.tsv"))

# intronic
exp <- fpcs_intronic %>%
  dplyr::mutate(end=as.numeric(full_end),source="CAGEr",type="promoter_cluster") %>%
  dplyr::select(-full_end)
rtracklayer::export.gff3(GenomicRanges::GRanges(exp),con = paste0("~/Dropbox/people/manuel/tss_data/tso/","pc_intronic.gff"))

pc_intronic_table <- fpcs_intronic %>%
  dplyr::mutate(full_start=as.integer(start),
                full_end=as.integer(full_end),
                start=ifelse(strand=="+",as.numeric(dominant_ctss)-50,as.numeric(dominant_ctss)-48),
                end=ifelse(strand=="+",as.numeric(dominant_ctss)+50,as.numeric(dominant_ctss)+52)) %>%
  GenomicRanges::GRanges()

seqs <- tibble(seq=BSgenome::getSeq(BSgenome.Pfalciparum.PlasmoDB.v24,pc_intronic_table,as.character=TRUE))
pc_intronic_table <- dplyr::bind_cols(tibble::as_tibble(pc_intronic_table),seqs) %>%
  dplyr::mutate(dist=as.integer(abs(as.numeric(anno_start)-as.numeric(dominant_ctss))),width=full_end-full_start) %>%
  dplyr::select(name,dominant_ctss,width,tpm,tpm_dominant_ctss,seq)

write_tsv(x=pc_intronic_table,path = paste0("~/Dropbox/people/manuel/tss_data/tso/","pc_intronic.tsv"))
