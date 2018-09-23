# Options -----------------------------------------------------------------

options(stringsAsFactors=F,warn=F)

# Load packages quietly
sshhh <- function(a.package){
  suppressWarnings(suppressPackageStartupMessages(
    library(a.package, character.only=TRUE)))
}

sshhh("readr")
sshhh("dplyr")
sshhh("broom")
sshhh("cowplot")
sshhh("stringr")
sshhh("seqinr")
sshhh("optparse")
sshhh("rtracklayer")

option_list <- list(
  make_option(c("-f", "--fasta"), type="character",
              help="FASTA file", metavar="character"),
  make_option(c("-l", "--length"), type="integer",
              help="Homopolymer length", metavar="integer", default = 15),
  make_option(c("-o", "--output"), type="character",
              help="Output file", metavar="character", default = "homopolymer_counts.tsv")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$fasta)) {
  print_help(opt_parser)
  stop("\nNeed FASTA file as input\n", call.=FALSE)
}

# Main --------------------------------------------------------------------

# read in FASTA file
fa <- read.fasta(opt$fasta,forceDNAtolower=F,as.string=T)

# set nucleotides
nucs <- c("A","T","G","C")

# set homopolymer tract length
k <- opt$length

# create list of homopolymers
hlist <- unlist(
  lapply(nucs, function(nuc) {
    lapply(seq(2,k), function(k) {
      paste(rep(nuc, k),collapse="")
    })
  })
)

# retrieve names of all sequences
# should be in fasta header
ids <- names(fa)

# count how often we see homopolymer tracts in sequences
counts <- lapply(ids, function(id){
    lapply(hlist, function(homopolymer){
    sequence <- fa[[id]][[1]]
    num <- str_count(sequence, homopolymer)
    nuc <- unlist(str_split(homopolymer,""))[1]
    len <- length(unlist(str_split(homopolymer,"")))
    c(gene_id=as.character(id), nuc=as.character(nuc), len=as.integer(len), count=as.integer(num))
    })
})

# Reformat the output so we can compute on it...
# R bullshit...
# suggestions on how to make this step more efficient: http://stackoverflow.com/questions/5942760/most-efficient-list-to-data-frame-method
df <- as.data.frame(matrix(unlist(counts,recursive=F)), nrow=length(unlist(counts[1],recursive=F)))
df <- as.data.frame(t(as.data.frame(df$V1)))
rownames(df) <- NULL
df <- as.data.frame(df)

df$gene_id <- as.character(df$gene_id)
df$nuc <- as.character(df$nuc)
df$len <- as.integer(df$len)
df$count <- as.integer(df$count)

# output data frame
out <- df %>% 
  filter(count > 0) %>%
  group_by(gene_id,nuc) %>% 
  summarise(longest=max(len))

write_tsv(out, path = opt$output)

q(save="no")
