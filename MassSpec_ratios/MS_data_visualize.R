####################################
# Bar Plot for MassSpec-NT data
##################################
#Initialize
setwd("D:/human_mito/T cell MS data")
library(ggplot2)
library(tidyr)
library(ggpubr)
library(stringr)
#################################### Plot
data<-read.table('./Nucleotide_Raw_pmolpermil.txt', sep="\t", header = TRUE)
df<-data[1]

for (n in c('A','C','G','U') ) {print(n)
  ribo=paste(n,'TP',sep="")
  dribo=paste('d',n,'TP',sep="")
  if (n=="U") {dribo='dTTP'}
  df[n]=data[,ribo]/data[,dribo]
  }

#df<-read.table('./NTP Ratio.txt', sep="\t", header = TRUE)
colnames(df)<-c("Cell.Type","rA/dA", "rC/dC", "rG/dG", "rU/dT")
data<-gather(df, key="Nucleotide", value="Ratio", -Cell.Type)
obj<-ggbarplot(data, x = "Cell.Type", y = "Ratio", fill="Nucleotide",
               add = c("mean_sd", "point"),
               #error.plot=c("crossbar"),
               add.params=list(width=0.4),
               palette = c("firebrick1","deepskyblue","gold","green3"),
               position = position_dodge(0.75),
               xlab="", ylab="rNTP/dNTP Ratio",
               width = 0.75,
               size = 1) +
  scale_y_continuous(expand = c(0,0),limits = c(0,1200), breaks = seq(0,1200,100)) +
  theme_classic(base_size = 25)+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 12))
ggpar(obj, legend=c("right"), legend.title = "rNTP/dNTP")
