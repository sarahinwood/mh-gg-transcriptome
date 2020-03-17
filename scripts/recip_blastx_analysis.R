library("data.table")
library("dplyr")
library("ggplot2")

##M.hyp results
mhyp_nr_blastx <- fread("output/recip_blast/nr_blastx/nr_blastx.outfmt3")
setnames(mhyp_nr_blastx, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13"), new=c("transcript_id", "nr_db_id", "%_identical_matches", "alignment_length", "no_mismatches", "no_gap_openings", "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score", "annotation"))
##order so that in event of eval min. tie, which.min takes hit with highest bitscore
setorder(mhyp_nr_blastx, transcript_id, evalue, -bit_score)
##extract result with lowest evalue for each peptide - what if multiple rows with lowest min?
Mh_min_evalues <- mhyp_nr_blastx[,.SD[which.min(evalue)], by=transcript_id]
##filter out virus annotations
Mh_virus <- dplyr::filter(Mh_min_evalues, grepl('virus', annotation))
fwrite(Mh_virus, "output/recip_blast/nr_blastx/viral_annots.csv")
transcript_id_annot <- Mh_virus[,c(1,13)]
fwrite(transcript_id_annot, "output/recip_blast/nr_blastx/viral_transcripts_annots.csv")
