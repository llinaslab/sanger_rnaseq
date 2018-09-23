source("../code/utils.R")

load_R_essentials()
library(rtracklayer)
library(GenomicRanges)

args <- commandArgs(TRUE)

# upload all coverage UTR files
tso_repaired_utrs <- as.data.frame(import.gff3(args[1]))
cov_repaired_utrs <- as.data.frame(import.gff3(args[2]))
org_cov_utrs <- as.data.frame(import.gff3(args[3]))
org_tso_utrs <- as.data.frame(import.gff3(args[4]))
three_prime_utrs <- as.data.frame(import.gff3(args[5]))
include_tso <- args[6]
output <- args[7]

# retrieve IDs from each
tso_repaired_ids <- tso_repaired_utrs$Parent
cov_repaired_ids <- cov_repaired_utrs$Parent
org_cov_ids <- org_cov_utrs$Parent
org_tso_ids <- org_tso_utrs$Parent

# determine coverage repaired UTRs that were not repaired by TSO data
add_to_utrs <- unlist(cov_repaired_ids[cov_repaired_ids %nin% tso_repaired_ids])

# create data frame to combine with TSO repaired UTRs to form combined repaired UTRs
df_to_add <- cov_repaired_utrs %>% filter(Parent %in% add_to_utrs)

# create combined repaired UTRs
combined_utrs <- rbind(tso_repaired_utrs, df_to_add)

# change feature column
combined_utrs$type <- "5UTR"

# only include these if 3D7
if(include_tso == "TRUE") {

  # retrieve combined IDs
  combined_ids <- combined_utrs$Parent

  # determine which original TSO UTRs were not repaired and thus will end up in the final list
  add_to_utrs <- unlist(org_tso_ids[org_tso_ids %nin% combined_ids])

  # create data frame to combine with formerly combined repaired UTRs
  df_to_add <- org_tso_utrs %>% filter(Parent %in% add_to_utrs)

  # add column of NAs to be consistent
  df_to_add$length_added <- 0

  # create a new combined list of UTRs
  combined_utrs <- rbind(combined_utrs, df_to_add)
}

# retrieve combined IDs
combined_ids <- combined_utrs$Parent

# determine which original coverage UTRs were not repaired and are not also TSO UTRs
add_to_utrs <- unlist(org_cov_ids[org_cov_ids %nin% combined_ids])

# create data frame to combine with formerly combined repaired UTRs
df_to_add <- org_cov_utrs %>% filter(Parent %in% add_to_utrs)

# add column of NAs to be consistent
df_to_add$length_added <- 0

# create a new combined list of UTRs
final_five_prime_utrs <- rbind(combined_utrs, df_to_add)

three_prime_utrs$length_added <- 0

final_utrs <- rbind(final_five_prime_utrs, three_prime_utrs)

final_utrs$source <- "Combined_UTRs"

export.gff3(GRanges(final_utrs), output)

# sort output
system(paste("sort -k1,1 -k4,4n", output, "> tmp; cat tmp >", output, "; rm -rf tmp"))
