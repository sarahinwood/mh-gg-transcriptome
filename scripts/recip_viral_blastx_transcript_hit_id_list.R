#!/usr/bin/env Rscript

#######
# LOG #
#######

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

#############
# LIBRARIES #
#############

library(data.table)

###########
# GLOBALS #
###########

blastx_res <- snakemake@input[["blastx_res"]]

########
# MAIN #
########

blastx_res_table <- fread(blastx_res)
fwrite(unique(blastx_res_table[,list(V1)]), snakemake@output[["transcript_hit_ids"]], col.names = FALSE)

#write log
sessionInfo()