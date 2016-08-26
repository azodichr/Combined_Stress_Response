library(fpc)
library(reshape2)
library(ggplot2)

# Set arguments
k <- 12
file <- "kmeans"
stresses <-c("High Mg2+","High CO2","Dual")
setwd("~/Desktop/Combined_stress/Niu_MgCO2/02_Clusters/Kmeans/")
open <- "../../Responses_all.txt"

### Generate cluster with kmeans
data <- read.table(open, header=T,row.names = 1, sep="\t")

### Filters out genes that aren't DE under at least 1 condition 
data <- subset(data, abs(st1_logFC) >= 1 | abs(st2_logFC) >= 1 | abs(dual_logFC) >= 1)
names(data) <- stresses
### Run kmeans clustering

c.data <- kmeans(data, k, iter.max = 100, nstart = 50)
#plotcluster(data,c.data$cluster,main="5 Clusters")
c.data$size
cluster <- as.matrix(c.data$cluster)

### Write clusters to .txt files
for (i in 1:k) {
  name <- paste(file, "_" , k,"_", i,".txt", sep="")
  d <- subset(data, cluster==i)
  write.table(d, file = name, quote=FALSE, sep="\t")
}


### Make boxplots of clusters 
cbpalette <- c("#0072B2","indianred1", "#009E73")
for (i in 1:k) {
  d <- subset(data, cluster==i)
  num <- nrow(d)
  name <- paste("plots/",file, "_" , k,"_", i,"_n", num, ".png", sep="")
  melt <- melt(d)
  ggplot(data=melt, aes(x=variable, y=value))+
    geom_boxplot(outlier.shape=NA, fill = cbpalette) + 
    coord_cartesian(ylim=c(-5,5)) +
    geom_hline(yintercept=c(-1,1), linetype='dotted') +
    xlab("") +
    ylab("log2(Fold Change)") +
    theme_bw(base_size=40)
  ggsave(name, width = 8, height = 8)
}

