# These chunks can be shared across multiple documents with knitr::read_chunk().
# http://yihui.name/knitr/demo/externalization/

# ---- knitr-opts-chunk ----
knitr::opts_chunk$set(
  comment = NA,
  fig.align = "center",
  tidy = FALSE,
  fig.path = paste0("figure/", knitr::current_input(), "/"),
  warning = F,
  message = F
)

# ---- load-libraries ----
source("../code/utils.R")
load_essentials()

# ---- session-info ----
sessionInfo()
