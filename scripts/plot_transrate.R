library(data.table)
library(ggplot2)

transrate.contigs <- fread('output/transrate/Trinity/contigs.csv')

plot.data <- melt(transrate.contigs, id.vars = "contig_name", measure.vars = c("sCnuc", "sCcov", "sCord", "sCseg")) 

facet.order = c("sCnuc"="s(C[nuc])", "sCcov"="s(C[cov])", "sCord"="s(C[ord])", "sCseg"="s(C[seg])")
plot.data[,facet.label:=plyr::revalue(variable, facet.order)]

ggplot(plot.data, aes(x=value))+facet_wrap(~facet.label, labeller=label_parsed)+geom_freqpoly(binwidth=0.01)+
  scale_y_log10()+xlab("TransRate Contig Score")+
  ylab(expression(paste("Number of Contigs ")))
