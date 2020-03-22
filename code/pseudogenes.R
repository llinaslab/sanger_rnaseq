source("../code/utils.R")
load_essentials()

x3d7exp <- readRDS("../output/neighboring_genes/gene_reduced_3d7_abund.rds")
xhb3exp <- readRDS("../output/neighboring_genes/gene_reduced_hb3_abund.rds")
xitexp  <- readRDS("../output/neighboring_genes/gene_reduced_it_abund.rds")

abund <- dplyr::bind_rows(x3d7exp,xhb3exp,xitexp)

rm(x3d7exp,xhb3exp,xitexp)

pseudo <- readr::read_tsv("../data/pseudogenes/pseudogenes.txt")$`Gene ID`

pseudo_detected <- abund %>%
  dplyr::filter(gene_id %in% pseudo) %>%
  dplyr::group_by(gene_id, strain) %>%
  dplyr::summarise(maxTPM = max(TPM)) %>%
  dplyr::ungroup() %>%
  tidyr::spread(strain,maxTPM) %>%
  dplyr::rename(`Gene ID`=gene_id,`3D7`=`3d7`,HB3=hb3,IT=it)

readr::write_tsv(pseudo_detected,"../output/pseudogenes/daftseq_detected_pseudogenes.tsv")

pseudo_detected_long <- abund %>%
  dplyr::filter(gene_id %in% pseudo) %>%
  dplyr::group_by(gene_id, strain) %>%
  dplyr::summarise(n = sum(TPM >= 1)) %>%
  dplyr::mutate(on = ifelse(n > 0, 1, 0)) %>%
  dplyr::ungroup() %>%
  dplyr::select(-n)

pseudo_detected_long %>% dplyr::group_by(gene_id) %>% dplyr::summarise(n=sum(on)) %$% table(n)

pseudo_detected_long %>% dplyr::filter(on == 1) %$% table(strain)

pseudo_detected_long %>% dplyr::filter(on == 1 & strain == "3d7") %>% dplyr::select(gene_id) %>% readr::write_tsv("..output/pseudogenes/pseudoon3d7.txt")
pseudo_detected_long %>% dplyr::filter(on == 1 & strain == "hb3") %>% dplyr::select(gene_id) %>% readr::write_tsv("../output/pseudogenes/pseudoonhb3.txt")
pseudo_detected_long %>% dplyr::filter(on == 1 & strain == "it") %>% dplyr::select(gene_id) %>% readr::write_tsv("../output/pseudogenes/pseudoonit.txt")

g <- plot_strain_abundances(abund,"PF3D7_0402200") +
  scale_fill_hue(l=55) +
  theme_classic() +
  theme(axis.text.x = element_text(size=24),
        axis.text.y = element_text(vjust=0.5, size=24),
        axis.title.y = element_text(size=28),
        axis.title.x = element_blank(),
        plot.title = element_text(size=28,hjust = 0.5),
        axis.line.x = element_line(size=0.5),
        axis.line.y = element_line(size=0.5),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(size = 0.5),
        axis.ticks.length = unit(0.25, "cm"),
        axis.ticks = element_line(size = 0.5),
        legend.position = "none") +
  fill_colors +
  panel_border(colour="black",size=1) +
  ggtitle("SURFIN 4.1, pseudogene\nPF3D7_0402200")

ggsave(plot=g,filename=paste0("../output/pseudogenes/PF3D7_0402200_abundances.png"),height=6,width=10)

g <- plot_strain_abundances(abund,"PF3D7_0402200") + coord_cartesian(ylim=c(0,25)) +
  scale_fill_hue(l=55) +
  theme_classic() +
  theme(axis.text.x = element_text(size=24),
        axis.text.y = element_text(vjust=0.5, size=24),
        axis.title.y = element_text(size=28),
        plot.title = element_text(size=28,hjust = 0.5),
        axis.line.x = element_line(size=0.5),
        axis.line.y = element_line(size=0.5),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(size = 0.5),
        axis.ticks.length = unit(0.25, "cm"),
        axis.ticks = element_line(size = 0.5),
        legend.position = "none") +
  fill_colors +
  panel_border(colour="black",size=1) +
  ggtitle("SURFIN 4.1, pseudogene\nPF3D7_0402200")

ggsave(plot=g,filename=paste0("../output/pseudogenes/PF3D7_0402200_abundances_zoomed.png"),height=6,width=10)

g <- plot_strain_profiles(abund,"PF3D7_0402200") +
  scale_fill_hue(l=55) +
  theme_classic() +
  theme(axis.text.x = element_text(size=24),
        axis.text.y = element_text(vjust=0.5, size=24),
        axis.title.y = element_text(size=28),
        axis.title.x = element_blank(),
        plot.title = element_text(size=28,hjust = 0.5),
        axis.line.x = element_line(size=0.5),
        axis.line.y = element_line(size=0.5),
        axis.ticks.y = element_line(size = 0.5),
        axis.ticks.length = unit(0.25, "cm"),
        axis.ticks = element_line(size = 0.5),
        legend.position = "none") +
  fill_colors +
  panel_border(colour="black",size=1) +
  ggtitle("SURFIN 4.1, pseudogene\nPF3D7_0402200")

ggsave(plot=g,filename=paste0("../output/pseudogenes/PF3D7_0402200_profiles.png"),height=6,width=10)

g <- plot_strain_abundances(abund,"PF3D7_0402300") +
  scale_fill_hue(l=55) +
  theme_classic() +
  theme(axis.text.x = element_text(size=24),
        axis.text.y = element_text(vjust=0.5, size=24),
        axis.title.y = element_text(size=28),
        axis.title.x = element_blank(),
        plot.title = element_text(size=28,hjust = 0.5),
        axis.line.x = element_line(size=0.5),
        axis.line.y = element_line(size=0.5),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(size = 0.5),
        axis.ticks.length = unit(0.25, "cm"),
        axis.ticks = element_line(size = 0.5),
        legend.position = "none") +
  fill_colors +
  panel_border(colour="black",size=1) +
  ggtitle("RH1\nPF3D7_0402300")

ggsave(plot=g,filename=paste0("../output/pseudogenes/PF3D7_0402300_abundances.png"),height=6,width=10)

g <- plot_strain_profiles(abund,"PF3D7_0402300") +
  scale_fill_hue(l=55) +
  theme_classic() +
  theme(axis.text.x = element_text(size=24),
        axis.text.y = element_text(vjust=0.5, size=24),
        axis.title.y = element_text(size=28),
        axis.title.x = element_blank(),
        plot.title = element_text(size=28,hjust = 0.5),
        axis.line.x = element_line(size=0.5),
        axis.line.y = element_line(size=0.5),
        axis.ticks.y = element_line(size = 0.5),
        axis.ticks.length = unit(0.25, "cm"),
        axis.ticks = element_line(size = 0.5),
        legend.position = "none") +
  fill_colors +
  panel_border(colour="black",size=1) +
  ggtitle("RH1\nPF3D7_0402300")

ggsave(plot=g,filename=paste0("../output/pseudogenes/PF3D7_0402300_profiles.png"),height=6,width=10)


