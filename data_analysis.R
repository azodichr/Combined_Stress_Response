library(ggplot2)
library(ggthemes)
library(reshape)

setwd("Desktop/Dual_dh_norm/03_ByGene_20perc/03_DE/")

############Upload Data by Stress#################
myvars <- c("logFC", "adj.P.Val")
rawdrought <- read.table("contrast_drought.csv", header=TRUE, row.names='gene')[myvars]
rawheat <- read.table("contrast_heat.csv", header=TRUE, row.names='gene')[myvars]
rawdual <- read.table("contrast_dual.csv", header=TRUE, row.names='gene')[myvars]
names(rawdrought) <- c("Drought_logFC", "Drought_adj.P.Val")
names(rawheat) <- c("Heat_logFC", "Heat_adj.P.Val")
names(rawdual) <- c("Dual_logFC", "Dual_adj.P.Val")

###From ordered and combined csv/txt files with all stresses as columns###
At_all <- read.csv("contrast_combined.csv",header=TRUE)
Sb_all <- read.table("~/Desktop/Dual_dh_norm/Sorghum/03_RMA_Contrasts/contrast_combined_ordered.txt", header=TRUE)
ALL <- rbind(At_all, Sb_all)

my.sub <- c('gene',"Drought_logFC", "Heat_logFC","Dual_logFC","species")
ALL_onesig <- subset(ALL, Drought_adj.P.Val <= 0.05 | Heat_adj.P.Val <= 0.05 | Dual_adj.P.Val <= 0.05)[my.sub]
names(ALL_onesig) <- c("gene", "Drought", "Heat", "Combined", 'species')
allsig <- subset(ALL, Drought_adj.P.Val <= 0.05 & Heat_adj.P.Val <= 0.05 & Dual_adj.P.Val <= 0.05)[my.sub]

test <- head(onesig_FC)
testone <- test[order(test$Dual_logFC, test$Heat_logFC),]
ordered <-onesig_FC[order(onesig_FC$Dual_logFC),]


############Scatterplots#################

qplot(data=rawdrought,x=row.names, geom="boxplot")
lis <- c("Drought_logFC", "Heat_logFC", "Dual_logFC")
qplot(row.names, Dual_logFC, data=onesignif)


############Regression Models#################

total <- glm(Dual_logFC ~ Drought_logFC + Heat_logFC + Drought_logFC*Heat_logFC, onesig_FC, family = gaussian)
summary(total)
coef(summary(total))
confint(total, level = 0.95)

############Rasmussen Category Clustering#################

forRasCat <- ALL_onesig  #select which data.frame to sort by Ras Cat.

co1 <- subset(forRasCat, Drought > -1 & Drought < 1 & Heat > -1 & Heat < 1 & Combined >= 1)
co1$RasCat <- "co1"
co2 <- subset(forRasCat, Drought > -1 & Drought < 1 & Heat > -1 & Heat < 1 & Combined <= -1)
co2$RasCat <- "co2"
co3 <- subset(forRasCat, Drought >= 1 & Heat >= 1 & Combined > -1 & Combined < 1)
co3$RasCat <- "co3"
co4 <- subset(forRasCat, Drought <= -1 & Heat <= -1 & Combined > -1 & Combined < 1)
co4$RasCat <- "co4"
co5 <- subset(forRasCat, Drought >= 1 & Heat >= 1 & Combined <= -1) #No At or Sb genes fit this category
co6 <- subset(forRasCat, Drought <= -1 & Heat <= -1 & Combined >= 1)
co6$RasCat <- "co6"
co7 <- subset(forRasCat, Drought <= -1 & Heat > -1 & Heat < 1  & Combined >= 1)
co7$RasCat <- "co7"
co8 <- subset(forRasCat, Drought > -1 & Drought < 1 & Heat <= -1 & Combined >= 1)
co8$RasCat <- "co8"
co9 <- subset(forRasCat, Drought >= 1 & Heat > -1 & Heat < 1  & Combined <= -1)
co9$RasCat <- "co9"
co10 <- subset(forRasCat, Drought > -1 & Drought < 1 & Heat >= 1 & Combined <= -1)
co10$RasCat <- "co10"

ca1 <- subset(forRasCat, Drought > -1 & Drought < 1 & Heat >= 1 & Combined > -1 & Combined < 1)
ca1$RasCat <- "ca1"
ca2 <- subset(forRasCat, Drought > -1 & Drought < 1 & Heat <= -1 & Combined > -1 & Combined < 1)
ca2$RasCat <- "ca2"
ca3 <- subset(forRasCat, Drought >= 1 & Heat > -1 & Heat < 1 & Combined > -1 & Combined < 1)
ca3$RasCat <- "ca3"
ca4 <- subset(forRasCat, Drought <= -1 & Heat > -1 & Heat < 1 & Combined > -1 & Combined < 1)
ca4$RasCat <- "ca4"
ca5 <- subset(forRasCat, Drought >= 1 & Heat <= -1 & Combined > -1 & Combined < 1)
ca5$RasCat <- "ca5"
ca6 <- subset(forRasCat, Drought <= -1 & Heat >= 1 & Combined > -1 & Combined < 1)
ca6$RasCat <- "ca6"

s1 <- subset(forRasCat, Drought >= 1 & Heat >= 1 & Combined >= 1)
s1$RasCat <- "s1"
s2 <- subset(forRasCat, Drought <= -1 & Heat <= -1 & Combined <= -1)
s2$RasCat <- "s2"

p1 <- subset(forRasCat, Drought >= 1 & Heat <= -1 & Combined >= 1)
p1$RasCat <- "p1"
p2 <- subset(forRasCat, Drought <= -1 & Heat >= 1 & Combined <= -1)
p2$RasCat <- "p2"
p3 <- subset(forRasCat, Drought >= 1 & Heat <= -1 & Combined <= -1)
p3$RasCat <- "p3"
p4 <- subset(forRasCat, Drought <= -1 & Heat >= 1 & Combined >= 1)
p4$RasCat <- "p4"

i1 <- subset(forRasCat, Drought >= 1 & Heat > -1 & Heat < 1 & Combined >= 1)
i1$RasCat <- "i1"
i2 <- subset(forRasCat, Drought <= -1 & Heat > -1 & Heat < 1 & Combined <= -1)
i2$RasCat <- "i2"
i3 <- subset(forRasCat, Drought > -1 & Drought < 1  & Heat >= 1 & Combined >= 1)
i3$RasCat <- "i3"
i4 <- subset(forRasCat, Drought > -1 & Drought < 1  & Heat <= -1 & Combined <= -1)
i4$RasCat <- "i4"

nr <- subset(forRasCat, Drought > -1 & Drought < 1 & Heat > -1 & Heat < 1 & Combined > -1 & Combined < 1)
nr$RasCat <- "nr"


All_Cats <- rbind(co1,co2,co3,co4,co6,co7,co8,co9,co10,ca1,ca2,ca3,ca4,ca5,ca6,i1,i2,i3,i4,p1,p2,p3,p4,s1,s2,nr)
All_Cats_mm= melt(All_Cats, id=c('gene','species','RasCat'))


############Boxplots#################

results <- "~/Desktop/Dual_dh_norm/Comparitive_At_Sb/"

cat.graph <- function(df, na.rm=TRUE, ...){
  cat_list<- unique(All_Cats_mm$RasCat)
  for (i in seq_along(cat_list)){
    cat <- cat_list[i]
    atnum <- length(which(All_Cats$RasCat == cat & All_Cats$species == "Arabidopsis"))
    sbnum <- length(which(All_Cats$RasCat == cat & All_Cats$species == "Sorghum"))
    plot <- ggplot(subset(df, df$RasCat==cat_list[i])) + 
      geom_boxplot(aes(x=variable, y=value, fill=variable),outlier.colour="NA") + 
      scale_fill_manual('variable',values=c('Drought'='#56B4E9', 'Heat'='#FF9999', 'Combined'='#009E73')) +
      facet_grid(.~species) +
      ylim(-5,5) +
      theme(axis.title.x=element_blank()) +
      theme(legend.position='none') +
      ylab("Log Fold Change") +
      ggtitle(paste("Rasmussen Category: ", cat," (At n = ", atnum,"; Sb n = ",sbnum," )", sep=""))
    ggsave(plot, file=paste(results,'Boxplots/',cat_list[i],'.pdf', sep=''), scale=1)
  }
}
cat.graph(All_Cats_mm)



############Clustering#################
install.packages("fpc")
library(fpc)
install.packages("jpeg")
install.packages("Cairo")
library(DLL)
jpeg(file="myplot.jpeg")



#generate cluster with kmeans
tocluster <- allsig_FC
k <- 26
c.tocluster <- kmeans(tocluster, k, iter.max = 100, nstart = 50)

#write clusters to a .txt file
list <- as.list(c.tocluster$cluster)
write.table(as.matrix(list), "K160_allClusters.txt", sep="\t")

# See how well the cluster separate from one another
plotcluster(tocluster,c.tocluster$cluster,main=k ,xlim=c(-60,50))
c.tocluster$size

# Get the cluster assignments
c.tocluster_assign <- c.tocluster$cluster

for (i in 1:k) {
  has <- which(c.tocluster_assign==i)
  len <-length(has)
  Matrix <- all[has,]
  MAX <- max(Matrix)
  MIN <- min(Matrix)
  CN <- colnames(Matrix)
  CNN <- length(CN)
  stresses <-c("Drought","Heat","Dual")
  name <- paste("K26 AllSignificantFC: Cluster #",i," (n=",len,")",sep="")
  png(file=paste("K26all_",i,".png",sep=""))
  plot(t(Matrix[1,]),type="l",ylim=c(MIN,MAX),ylab="Expression FC",col="gray",xaxt="n",xlab="",main=name)
  axis(1,at=seq(1,CNN,by=1),labels=stresses)
  for (j in seq(2,dim(Matrix)[1])){
    lines(t(Matrix[j,]),col="gray")
  }
  MED <- sapply(Matrix,median)
  lines(MED,col="blue")
  dev.off()
} 