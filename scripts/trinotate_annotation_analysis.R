library(data.table)
library(VennDiagram)
library(stringr)
library(ggplot2)

#Get annotation report in right format for venn diagram
annotation.report <- fread('output/trinotate/trinotate/trinotate_annotation_report.txt', na.strings = ".")
pfam <- annotation.report[!is.na(Pfam), unique(`#gene_id`)]
blastx <- annotation.report[!is.na(sprot_Top_BLASTX_hit), unique(`#gene_id`)]
kegg <- annotation.report[!is.na(Kegg), unique(`#gene_id`)]
number.genes <- annotation.report[!is.na(`#gene_id`),length(unique(`#gene_id`))]

#Draw Venn Diagram
Set1 <- RColorBrewer::brewer.pal(3, "Set1")
vd <- venn.diagram(x = list("Pfam"=pfam, "BlastX"=blastx, "Kegg"=kegg), filename=NULL, fill=Set1, alpha=0.5, cex = 1, cat.cex=1, lwd=1, main=paste("Total Number of Genes = ", number.genes))
grid.newpage()
grid.draw(vd)

#Sum of genes with any annotation
long.annotationreport <- melt(annotation.report,id.vars = "#gene_id", measure.vars = c("sprot_Top_BLASTX_hit", "sprot_Top_BLASTP_hit", "Pfam",  "eggnog", "Kegg"))
any.annotations <- long.annotationreport[,.(any_annotations = any(!is.na(value))),by=`#gene_id`]
any.annotations[,length(unique(`#gene_id`))]
any.annotations[,sum(any_annotations)]

#Sum of genes with any Transdecoder predicted protein
transdecoder.report <- melt(annotation.report, id.vars="#gene_id", measure.vars = c("prot_id"))
any.transdecoder <- transdecoder.report[,.(transdecoder_annotation = any(!is.na(value))),by=`#gene_id`]
any.transdecoder[,length(unique(`#gene_id`))]
any.transdecoder[,sum(transdecoder_annotation)]
sum(any.transdecoder$transdecoder_annotation==FALSE)

#Annotations per taxa
blastx.results <- annotation.report[!is.na(sprot_Top_BLASTX_hit),.(sprot_Top_BLASTX_hit, `#gene_id`)]
first.blastx.hit <- blastx.results[,tstrsplit(sprot_Top_BLASTX_hit, "`", fixed = TRUE, keep=1), by = `#gene_id`]
split.first.blastx <- first.blastx.hit[,tstrsplit(V1, "^", fixed=TRUE), by=`#gene_id`]
genes.per.taxa <- split.first.blastx[,length(unique(`#gene_id`)), by=V7]
setkey(genes.per.taxa, V1)
print(genes.per.taxa)
fwrite(genes.per.taxa, "output/trinotate/genes_per_taxa.csv")

#meanwhile in excel sort for taxa with most annotations, delete those with low no.
#I'm not interested in, and alter taxa name to just genera

#plot annotations per taxa
plot.genes.per.taxa <- fread("output/trinotate/genes_per_taxa_edited_for_plot.csv")

ggplot(plot.genes.per.taxa, aes(x=reorder(V7, -V1), y=V1))+
  theme(axis.text.x = element_text(angle = 65, hjust = 1, face = "italic")) +
  geom_col()+xlab("Genera")+ylab("Number of BlastX Annotations")

##actin for qpcr?
grep_actin_genes <- dplyr::filter(annotation.report, grepl('actin', sprot_Top_BLASTX_hit))
actin <- dplyr::filter(grep_actin_genes, !grepl('interacting', sprot_Top_BLASTX_hit))
fwrite(actin, "output/trinotate/blastp_grep_actin.csv")
