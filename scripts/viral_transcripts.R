library(data.table)
library(ggplot2)

annotation.report <- fread('output/trinotate/trinotate/trinotate_annotation_report.txt', na.strings = ".")

blastx.results <- annotation.report[!is.na(sprot_Top_BLASTX_hit),.(sprot_Top_BLASTX_hit, `#gene_id`)]
first.blastx.hit <- blastx.results[,tstrsplit(sprot_Top_BLASTX_hit, "`", fixed = TRUE, keep=1), by = `#gene_id`]
split.first.blastx <- first.blastx.hit[,tstrsplit(V1, "^", fixed=TRUE), by=`#gene_id`]
genes.per.taxa <- split.first.blastx[,length(unique(`#gene_id`)), by=V7]
virus_taxa <- dplyr::filter(genes.per.taxa, grepl('viruses', V7))
fwrite(virus_taxa, "output/trinotate/viral/genes_per_taxa_viral.csv")

##plot viral taxa annots
plot.viral.taxa <- fread("output/trinotate/viral/genes_per_taxa_viral_edited_for_plot.csv")

ggplot(plot.viral.taxa, aes(x=reorder(V7, -V1), y=V1))+
  theme(axis.text.x = element_text(angle = 65, hjust = 1, face = "italic")) +
  geom_col()+xlab("Viral Genera")+ylab("Number of BlastX Annotations")


trinotate_virus_annots <- dplyr::filter(annotation.report, grepl('virus', sprot_Top_BLASTX_hit))
viral_not_transposon <- data.table(dplyr::filter(trinotate_virus_annots, !grepl('transposon', sprot_Top_BLASTX_hit)))
fwrite(viral_not_transposon, "output/trinotate/viral/viral_annots.csv")
##write list of viral transcript IDs to pull out of fasta file
viral_transcripts <- viral_not_transposon[,transcript_id]
fwrite(list(viral_transcripts), "output/trinotate/viral/viral_transcript_ids.txt")

##baculovirus genes
baculovirus_genes <- dplyr::filter(viral_not_transposon, grepl('Baculoviridae', sprot_Top_BLASTX_hit))
baculovirus_gene_ids <- baculovirus_genes[,"transcript_id"]
fwrite(list(baculovirus_gene_ids), "output/trinotate/viral/baculoviridae_gene_ids.txt")
fwrite(baculovirus_genes, "output/trinotate/viral/baculoviridae_annots.csv")

bro_n <- dplyr::filter(annotation.report, grepl('Bro-N', Pfam))
fwrite(bro_n, "output/trinotate/viral/bro-n.csv")

fwrite(trinotate_virus_annots, "output/trinotate/viral/viral_annots_trinotate.csv")

##sort out viral annots
dedup_virus_annots <- fread("output/trinotate/viral/dedup_viral_annots_trinotate.csv")
viral_blastx_annots <- dedup_virus_annots[,tstrsplit(sprot_Top_BLASTX_hit, "`", fixed=TRUE, keep=1)]
viral_blastx_table <- viral_blastx_annots[,tstrsplit(V1, "^", fixed=TRUE)]
##table contains 200-odd eukaryote transposon annots - filter for virus in taxa column
virus_blastx_table <- dplyr::filter(viral_blastx_table, grepl('Viruses', V7))
fwrite(virus_blastx_table, "output/trinotate/viral/virus_blastx_annots.csv")
