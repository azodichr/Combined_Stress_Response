## arules: Mining Association Rules in R studio

library(arules)
library(arulesViz)

## Default values
#setwd("~/Desktop/Combined_stress/ARuleMining/")
#num <- 12
#measure <- 'lift'
#file <- "NNU_Ras_df.txt_decisiontree_50"
## 

args = commandArgs(TRUE)
setwd(args[1])
file <- args[2]
num <- as.integer(args[3])
measure <- args[4]

db<-read.table(file, header=T, sep='\t', quote=NULL, row.names=1)

pos <- subset(db, db$Class == 1)
pos$Class <- NULL
pos1 <- as.matrix(pos)

rules <- apriori(pos1, parameter = list(supp=0.15, conf=0.5, target = "rules"))
#top <- inspect(head(sort(rules, by = measure), num))

#plot(rules)

sig <- rules[is.significant(rules, pos1, alpha = 1e-50, adjust = "bonferroni")]

unique <- sig[!is.redundant(sig, measure = measure)]
top_sig_uni <- inspect(head(sort(sig, by = measure), num))


name=paste(basename(args[2]), "_rules.csv", sep="")
write.csv(top_sig_uni, file=name, row.names=FALSE)

