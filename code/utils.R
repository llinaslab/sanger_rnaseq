# Loading functions -------------------------------------------------------

# Load packages quietly
sshhh <- function(a.package){
  suppressWarnings(suppressPackageStartupMessages(
    library(a.package, character.only=TRUE)))
}

# Load essentials
load_essentials <- function() {
  sshhh("tidyverse")
  sshhh("stringr")
  sshhh("magrittr")
  sshhh("cowplot")
  sshhh("scales")
  sshhh("org.Pf.plasmo.db")
  sshhh("BSgenome.Pfalciparum.PlasmoDB.v24")
}

# Plot Colors -------------------------------------------------------------

myColors <- c("#377EB8", "#E41A1C", "#4DAF4A")
names(myColors) <- levels(c("3d7", "hb3", "it"))
fill_colors <- ggplot2::scale_fill_manual(name = "c", values = myColors)

myColors <- c("#377EB8", "#E41A1C", "#4DAF4A")
names(myColors) <- levels(c("3d7", "hb3", "it"))
outline_colors <- ggplot2::scale_color_manual(name = "c", values = myColors)

# Convenience functions ---------------------------------------------------

"%nin%" <- function(x, y) !(x %in% y)

# Create list from two column data frame
df2list <- function(df, key, value) {
  if((!is.character(key)) | (!is.character(value))) {
	  stop("\nKey and value must be input as character vectors.\n")
	}
  if(length(df[[key]]) != length(unique(df[[key]]))) {
    stop("\nKey column isn't unique! Remove duplicates first.\n")
	}
  whatyouwant <- setNames(as.character(df[[key]]), df[[key]])
	return(whatyouwant)
}

# Retrieve gene names
get_gene_names <- function() {
  require(org.Pf.plasmo.db)
  gene_names <- as.data.frame(org.Pf.plasmoGENENAME)
  gene_names <- setNames(c(gene_names$gene_name), gene_names$gene_id)
  return(gene_names)
}



# Plot functions ----------------------------------------------------------

# Strain Abundances
plot_strain_abundances <- function(df, gid) {

  gene_names <- get_gene_names()

  df %>%
    dplyr::filter(gene_id == gid) %>%
    ggplot(aes(x = toupper(strain), y = as.numeric(TPM), group = tp, fill = toupper(strain))) +
    geom_bar(stat="identity", position="dodge", colour="grey90") +
    ylab("TPM") +
    xlab("Strain") +
    ggtitle(paste(gid, "\n", gene_names[[gid]])) +
    scale_fill_hue(l=55) +
    theme(axis.title.x = element_text(size=14),
          axis.title.y = element_text(size=14),
          axis.text.x = element_text(vjust=0.5, size=14),
          axis.text.y = element_text(vjust=0.5, size=14),
          plot.title = element_text(vjust=0.5, size=16),
          axis.line.x=element_line(size=0),
          axis.line.y=element_line(size=0),
          axis.ticks.x = element_line(size = 0.5),
          axis.ticks.y = element_line(size = 0.5),
          axis.ticks.length = unit(0.25, "cm"),
          axis.ticks = element_line(size = 0.5),
          legend.text = element_text(size = 12),
          legend.title = element_blank()) +
    fill_colors +
    panel_border(colour="black")
}

# Strain Profiles
plot_strain_profiles <- function(df, gid) {

  gene_names <- get_gene_names()

  df %>%
    dplyr::filter(gene_id == gid) %>%
    group_by(gene_id, strain) %>%
    mutate(scaled = (log2(TPM+1)-mean(log2(TPM+1)))/sd(log2(TPM+1))) %>%
    ggplot(aes(x = tp, y = scaled, color = toupper(strain))) +
    ylab("") +
    xlab("") +
    ggtitle(paste(gid, " \n ", gene_names[[gid]])) +
    stat_smooth(aes(group = strain), se = F, size = 1.5) +
    geom_point(aes(group = strain), size = 2) +
    theme(legend.position="bottom",
          legend.direction="vertical",
          legend.title = element_blank(),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          axis.line.x=element_line(size=0),
          axis.line.y=element_line(size=0),
          axis.ticks = element_line(size = 0.5)) +
    outline_colors +
    scale_x_continuous(breaks = c(1,2,3,4,5,6,7), labels = c("T1", "T2", "T3" ,"T4", "T5", "T6", "T7")) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    panel_border(colour="black",remove=F,size=1)
}

#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

# Plot a bunch of genes to one output file
plot_gene_list <- function(func, df, list, output) {

  require(gridExtra)
  # match function name
  plotting_func <- match.fun(func)
  pdf(output, onefile = TRUE)
  # initiate loop
  i <- 1
  j <- 1
  g = list()
  while (i < length(list) + 1) {
    while (j <= 2) {
      if (list[i] %in% df$gene_id) {
        g[[j]] <- plotting_func(df, list[i])
        j <- j + 1
      }
      i <- i + 1
      if (i > length(list)) {
        g[[j]] <- grid.rect(gp=gpar(col="white"))
        break
      }
    }
    multiplot(plotlist=g,cols=1)
    g = list()
    j <- 1
  }
  dev.off()

}
