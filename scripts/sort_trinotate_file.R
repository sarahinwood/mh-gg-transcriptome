library(data.table)

##extracting the most sig annotation per gene (each gene has multiple annotations due to multiple isoforms)

trinotate.report <- fread('output/trinotate/trinotate/trinotate_annotation_report.txt', na.strings = ".")
##split to keep first blastx hit only
trinotate.sorted <- copy(trinotate.report)
trinotate.sorted$blastx_evalue <- tstrsplit(trinotate.report$sprot_Top_BLASTX_hit, "`", fixed=TRUE, keep=c(1))
##split to keep only evalue
trinotate.sorted$blastx_evalue <- tstrsplit(trinotate.sorted$blastx_evalue, "^RecName", fixed=TRUE, keep=c(1))
trinotate.sorted$blastx_evalue <- tstrsplit(trinotate.sorted$blastx_evalue, "E:", fixed=TRUE, keep=c(2))
##reorder to get evalue beside blastx annotation
trinotate.sorted <- trinotate.sorted[,c(1,2,3,17,4,5,6,7,8,9,10,11,12,13,14,15,16)]
##convert evalue column to numeric
trinotate.sorted[,4] <- sapply(trinotate.sorted[,4], as.numeric)
##set order so for each gene the transcript with the lowest evalue hit is at top
setorder(trinotate.sorted, `#gene_id`, `blastx_evalue`, na.last=TRUE)
####extract result with lowest evalue for each gene - what if multiple rows with lowest min?
trinotate.min.eval <- trinotate.sorted[,.SD[which.min(blastx_evalue)], by=`#gene_id`]
##write csv with most sig hit for each gene with annotation
fwrite(trinotate.min.eval, "output/trinotate/trinotate/sorted/most_sig_transcript_blastx_hit_for_each_gene.csv")

##filter out gene ids for unannotated genes
genes_no_annot <- trinotate.report[is.na(sprot_Top_BLASTX_hit),]
##list of unique gene ids from table of genes with no blastx annot
list_ids_no_annot <- list(unique(genes_no_annot$`#gene_id`))
fwrite(list_ids_no_annot, "output/trinotate/trinotate/sorted/ids_genes_no_blastx_annot.txt")

##filter out all genes with viral annots
sorted_virus_annots <- dplyr::filter(trinotate.min.eval, grepl('virus', sprot_Top_BLASTX_hit))
viral_not_transposon <- data.table(dplyr::filter(best_trinotate_virus_annots, !grepl('transposon', sprot_Top_BLASTX_hit)))
fwrite(viral_not_transposon, "output/trinotate/trinotate/sorted/viral_annots_not_transposons.csv")
##table of gene id vs annot
gene_id_viral_annot <- viral_not_transposon[,c(1,3)]
fwrite(gene_id_viral_annot, "output/trinotate/trinotate/sorted/gene_id_vs_viral_annot.csv")
