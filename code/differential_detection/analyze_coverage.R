setwd("/Users/llinaslab/projects/sanger_rnaseq/analysis/strain_transcript_detection")
source("../../src/utils.R")

load_R_essentials()

options(digits=4)

read_cov <- function(file) {

  df <- read_tsv(file,col_names=F) %>%
    dplyr::select(X9,X10,X11,X12,X13) %>%
    dplyr::rename(att=X9,reads=X10,nonzero=X11,len=X12,fraction=X13)

  df$exon <- apply(df, 1, function(x) {
    str_replace(str_split(str_split(x[["att"]],";")[[1]][1],"=")[[1]][2],"exon_","")
  })

  return(dplyr::select(df,-att))
}

files_3d7 <- list.files("data/cov/",pattern="3d7",full.names=T)
files_hb3 <- list.files("data/cov/",pattern="hb3",full.names=T)
files_it <- list.files("data/cov/",pattern="it",full.names=T)

tp1_3d7 <- read_cov(files_3d7[1]) %>% mutate(strain="3d7",tp="T1")
tp2_3d7 <- read_cov(files_3d7[2]) %>% mutate(strain="3d7",tp="T2")
tp3_3d7 <- read_cov(files_3d7[3]) %>% mutate(strain="3d7",tp="T3")
tp4_3d7 <- read_cov(files_3d7[4]) %>% mutate(strain="3d7",tp="T4")
tp5_3d7 <- read_cov(files_3d7[5]) %>% mutate(strain="3d7",tp="T5")
tp6_3d7 <- read_cov(files_3d7[6]) %>% mutate(strain="3d7",tp="T6")
tp7_3d7 <- read_cov(files_3d7[7]) %>% mutate(strain="3d7",tp="T7")
df_3d7 <- rbind(tp1_3d7,tp2_3d7,tp3_3d7,tp4_3d7,tp5_3d7,tp6_3d7,tp7_3d7)
rm(tp1_3d7,tp2_3d7,tp3_3d7,tp4_3d7,tp5_3d7,tp6_3d7,tp7_3d7)

tp1_hb3 <- read_cov(files_hb3[1]) %>% mutate(strain="hb3",tp="T1")
tp2_hb3 <- read_cov(files_hb3[2]) %>% mutate(strain="hb3",tp="T2")
tp3_hb3 <- read_cov(files_hb3[3]) %>% mutate(strain="hb3",tp="T3")
tp4_hb3 <- read_cov(files_hb3[4]) %>% mutate(strain="hb3",tp="T4")
tp5_hb3 <- read_cov(files_hb3[5]) %>% mutate(strain="hb3",tp="T5")
tp6_hb3 <- read_cov(files_hb3[6]) %>% mutate(strain="hb3",tp="T6")
tp7_hb3 <- read_cov(files_hb3[7]) %>% mutate(strain="hb3",tp="T7")
df_hb3 <- rbind(tp1_hb3,tp2_hb3,tp3_hb3,tp4_hb3,tp5_hb3,tp6_hb3,tp7_hb3)
rm(tp1_hb3,tp2_hb3,tp3_hb3,tp4_hb3,tp5_hb3,tp6_hb3,tp7_hb3)

tp1_it <- read_cov(files_it[1]) %>% mutate(strain="it",tp="T1")
tp2_it <- read_cov(files_it[2]) %>% mutate(strain="it",tp="T2")
tp3_it <- read_cov(files_it[3]) %>% mutate(strain="it",tp="T3")
tp4_it <- read_cov(files_it[4]) %>% mutate(strain="it",tp="T4")
tp5_it <- read_cov(files_it[5]) %>% mutate(strain="it",tp="T5")
tp6_it <- read_cov(files_it[6]) %>% mutate(strain="it",tp="T6")
tp7_it <- read_cov(files_it[7]) %>% mutate(strain="it",tp="T7")
df_it <- rbind(tp1_it,tp2_it,tp3_it,tp4_it,tp5_it,tp6_it,tp7_it)
rm(tp1_it,tp2_it,tp3_it,tp4_it,tp5_it,tp6_it,tp7_it)

coverages <- rbind(df_3d7,df_hb3,df_it)
coverages <- coverages[,c(5,6,7,1,2,3,4)] %>%
  separate(exon, into = c("gene_id","exon"), sep = "-") %>%
  group_by(gene_id,strain, tp) %>%
  summarise(reads=sum(reads),nonzero=sum(nonzero),len=sum(len),fraction=sum(nonzero)/sum(len))

coverages$tp <- factor(coverages$tp, levels = c("T1","T2","T3","T4","T5","T6","T7"))

core_genes <- read_tsv("../../data/gene_lists/core_pf3d7_genes.txt",col_names=F)$X1

mean_coverages <- coverages %>%
  group_by(gene_id,strain) %>%
  summarise(fraction=mean(fraction)) %>%
  spread(strain,fraction) %>%
  mutate(diffhb3 = `3D7` - HB3, diffit = `3D7` - IT)

max_coverages <- coverages %>%
  group_by(gene_id,strain) %>%
  summarise(fraction=max(fraction)) %>%
  spread(strain,fraction) %>%
  mutate(diffhb3 = `3D7` - HB3, diffit = `3D7` - IT)

g1 <- ggplot(mean_coverages,aes(x=diffhb3)) +
  geom_histogram(color="grey90") +
  xlab("Difference from HB3")
g2 <- ggplot(mean_coverages,aes(x=diffit)) +
  geom_histogram(color="grey90") +
  xlab("Difference from IT")

p <- plot_grid(g1,g2,nrow=1)

ggsave(plot = p, filename = "results/differences.pdf")
