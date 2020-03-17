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

abundance_file <- snakemake@input[["abundance"]]

########
# MAIN #
########

isoform.list <- fread(abundance_file)

#sort by expression
isoforms_by_expression <- isoform.list[,.I[which.max(IsoPct)], by=gene_id][,V1]
#write file sorted by expression
fwrite(isoform.list[isoforms_by_expression,list(transcript_id)],
	snakemake@output[["expression"]], col.names = FALSE)

#sort by length
isoforms_by_length <- isoform.list[,.I[which.max(length)], by=gene_id][,V1]
#write file sorted by length
fwrite(isoform.list[isoforms_by_length,list(transcript_id)],
	snakemake@output[["length"]], col.names = FALSE)

#plot transcript lengths
isoform.list[,hist(length, breaks = 100, xlim=c(0, +5000),
	main = "Transcript Lengths in New MH Transcriptome Assembly",
	xlab = "Transcript Length (bp)")]
sum(isoform.list$length>500)

# write log
sessionInfo()