library(data.table)
library(ggplot2)
fread ("output/trinity_stats/xn50.out.txt")
data.raw<-fread ("output/trinity_stats/xn50.out.txt")
data.raw[,`Ex`:=as.numeric(gsub("E", "", `Ex`))]

ggplot(data.raw, aes(x=`Ex`, y=`ExN50`)) + geom_point()+ scale_y_continuous(limits = c(0, 3000))

