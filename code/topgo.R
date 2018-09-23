# Options -----------------------------------------------------------------

options(warn=F,dplyr.width=Inf)

# Load packages quietly
sshhh <- function(a.package){
  suppressWarnings(suppressPackageStartupMessages(
    library(a.package, character.only=TRUE)))
}

sshhh("topGO")
sshhh("optparse")
sshhh("readr")
sshhh("dplyr")
sshhh("cowplot")

option_list <- list(
  make_option(c("-l", "--gene_list"), type="character",
              help="Gene list", metavar="character"),
  make_option(c("-t", "--go_terms"), type="character",
              help="GO terms", metavar="character"),
  make_option(c("-a", "--anno"), type="character",
              help="Gene names", metavar="character"),
  make_option(c("-m", "--method"), type="character", default = "weight01",
              help="Method [%default]", metavar="character"),
  make_option(c("-p", "--pvalue"), type="numeric", default = 0.05,
              help="P-value cutoff [%default]", metavar="numeric"),
  make_option(c("-o", "--out_prefix"), type="character", default = "topgo_out",
              help="Output prefix [%default]", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$gene_list)) {
  print_help(opt_parser)
  stop("\nNeed to input file containing gene list\n", call.=FALSE)
}

if (is.null(opt$go_terms)) {
  print_help(opt_parser)
  stop("\nNeed to input file of GO terms\n", call.=FALSE)
}

if (is.null(opt$anno)) {
  print_help(opt_parser)
  stop("\nNeed to input file of gene IDs and descriptions\n", call.=FALSE)
}

# Main --------------------------------------------------------------------

################################################################################
## Read and prepare data                                                       #
################################################################################

ref <- read_tsv(opt$go_terms, col_names = c("id", "go"))
ref.vec <- strsplit(ref$go, split=',', fixed=T)
names(ref.vec) <- ref$id
all.ids <- ref$id

#check to make sure list file isn't empty
if (file.exists(opt$gene_list)) {
	info <- file.info(opt$gene_list)
	if (info$size == 0) {
		message("File is empty...exiting.")
		q(save = "no")
	}
} else {
	message("File doesn't exist...exiting.")
	q(save = "no")
}

vec <- read_tsv(opt$gene_list, col_names = F)$X1

scores <- rep(0, nrow(ref)) # list of scores
names(scores) <- ref$id
scores[ref$id %in% vec] <- 1

desc <- read_tsv(opt$anno,col_names=c("id","desc"))

geneSelectionFun <- function(score){
  return(score >= 1)
}

################################################################################
## Create topGO objects                                                        #
################################################################################

GOdataBP <- new("topGOdata",
                ontology = 'BP',
                allGenes = scores,
                annot = annFUN.gene2GO,
                gene2GO = ref.vec,
                geneSelectionFun = geneSelectionFun,
                nodeSize = 5,
                description = '')
GOdataMF <- new("topGOdata",
                ontology = 'MF',
                allGenes = scores,
                annot = annFUN.gene2GO,
                gene2GO = ref.vec,
                geneSelectionFun = geneSelectionFun,
                nodeSize = 5,
                description = '')
GOdataCC <- new("topGOdata",
                ontology = 'CC',
                allGenes = scores,
                annot = annFUN.gene2GO,
                gene2GO = ref.vec,
                geneSelectionFun = geneSelectionFun,
                nodeSize = 5,
                description = '')

################################################################################
## Run enrichment tests                                                        #
################################################################################

resultTopgoBP <- runTest(GOdataBP,algorithm=opt$method,statistic="Fisher")
resultTopgoMF <- runTest(GOdataMF,algorithm=opt$method,statistic="Fisher")
resultTopgoCC <- runTest(GOdataCC,algorithm=opt$method,statistic="Fisher")

resBP<-GenTable(GOdataBP,topGO=resultTopgoBP,orderBy="topGO",ranksOf="fisher",topNodes=50)
resMF<-GenTable(GOdataMF,topGO=resultTopgoMF,orderBy="topGO",ranksOf="fisher",topNodes=50)
resCC<-GenTable(GOdataCC,topGO=resultTopgoCC,orderBy="topGO",ranksOf="fisher",topNodes=50)

################################################################################
## Write results to output files                                               #
################################################################################

# Write genes in term to file
BPterms <- list()
for (id in resBP$GO.ID) {
  g <- genesInTerm(GOdataBP)[[id]]
  t <- subset(resBP, GO.ID == id)$Term
  BPterms[[paste0(id," - ",t)]] <- desc %>% filter(id %in% g & id %in% vec) %>% data.frame()
}
sink(paste0(opt$out_prefix,"_terms_bp.txt"))
print(BPterms)

MFterms <- list()
for (id in resMF$GO.ID) {
  g <- genesInTerm(GOdataMF)[[id]]
  t <- subset(resMF, GO.ID == id)$Term
  MFterms[[paste0(id," - ",t)]] <- desc %>% filter(id %in% g & id %in% vec) %>% data.frame()
}
sink(paste0(opt$out_prefix,"_terms_mf.txt"))
print(MFterms)

CCterms <- list()
for (id in resCC$GO.ID) {
  g <- genesInTerm(GOdataCC)[[id]]
  t <- subset(resCC, GO.ID == id)$Term
  CCterms[[paste0(id," - ",t)]] <- desc %>% filter(id %in% g & id %in% vec) %>% data.frame()
}
sink(paste0(opt$out_prefix,"_terms_cc.txt"))
print(CCterms)

# Write enrichment results to file
sink(paste0(opt$out_prefix,"_bp.tsv"))
write.table(resBP, sep = "\t", quote=F, row.names=F)
sink(paste0(opt$out_prefix,"_mf.tsv"))
write.table(resMF, sep = "\t", quote=F, row.names=F)
sink(paste0(opt$out_prefix,"_cc.tsv"))
write.table(resCC, sep = "\t", quote=F, row.names=F)

################################################################################
## Create plots of top ten terms                                               #
################################################################################

# Function to order terms for plotting
orderTerms <- function(results) {
  r <- as.data.frame(results) %>% mutate(topGO = as.numeric(topGO))
  o <- order(sort(r$topGO,decreasing=FALSE))
  r$Term <- factor(r$Term,levels=rev(r$Term[o]))
  return(r)
}

# Function to plot top terms
plotTerms <- function(results,num=10) {
  results %>%
  head(num) %>%
    ggplot(aes(x=Term,y=-log10(topGO))) +
    geom_bar(stat="identity",color="black",fill="grey80",size=1) +
    coord_flip() +
    ylab("-Log10(p-value)") +
    theme(axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          axis.line.y = element_blank())
}

pBP <- orderTerms(resBP) %>%
  plotTerms(num=10)
pMF <- orderTerms(resMF) %>%
  plotTerms(num=10)
pCC <- orderTerms(resCC) %>%
  plotTerms(num=10)

ggsave(plot=pBP,filename=paste0(opt$out_prefix,"_bp.pdf"))
ggsave(plot=pMF,filename=paste0(opt$out_prefix,"_mf.pdf"))
ggsave(plot=pCC,filename=paste0(opt$out_prefix,"_cc.pdf"))

q(save="no")
