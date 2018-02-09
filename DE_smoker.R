#installing packages
source("https://bioconductor.org/biocLite.R")
biocLite("edgeR")
biocLite("limma")
library("limma")
library("edgeR")
library(RColorBrewer)
my_data <- read.table(file.choose(),check.names= FALSE)
dim(my_data)
matrix1<-as.matrix(my_data)
m1<- cpm(matrix1, log=T)
#heatmap(cor(m1))
#remove lowly expressed genes and visualize distribution
expr_cutoff <- 1.5
keep.exprs <- rowSums(matrix1>expr_cutoff)>=20
data_filtered <- my_data[keep.exprs,]
write.csv(data_filtered, file = "never_smoker_data_filtered_1.csv")
lcpm_filtered <- cpm(data_filtered, log=TRUE)

samples <- ncol(my_data)
col <- brewer.pal(samples, "Paired")
par(mfrow=c(1,2))
plot(density(m1[,1]), col=col[1], lwd=2, ylim=c(0,0.30), las=2, 
     main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:samples){
  den <- density(m1[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
plot(density(lcpm_filtered[,1]), col=col[1], lwd=2, ylim=c(0,0.30), las=2, 
     main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:samples){
  den <- density(lcpm_filtered[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}

#normalize data and visualization
group <- substr(colnames(data_filtered),1,1)
y <- DGEList(counts = data_filtered, group = group)
y_norm <- calcNormFactors(y, method = "TMM")
par(mfrow=c(1,2))
lcpm <- cpm(y, log=TRUE)
boxplot(lcpm, las=2, col=col, main="",xaxt="n")
title(main="A. Unnormalized data",ylab="Log-cpm")
lcpm <- cpm(y_norm, log=TRUE)
boxplot(lcpm, las=2, col=col, main="",xaxt="n")
title(main="B. Normalized data",ylab="Log-cpm")

# differential epxression analysis
set.seed(123)
y_norm <- estimateDisp(y_norm)
#sqrt(y_norm$common.dispersion)
#plotBCV(y_norm)
et <- exactTest(y_norm)
results_edgeR <- topTags(et, n = nrow(data_filtered), sort.by = "none")
head(results_edgeR$table)
write.csv(results_edgeR$table, file = "never_smoker_Results_edgeR.csv")

#MA plot for genes with FDR corrected pvalue smaller than 0.05
results_fdr<- sum(results_edgeR$table$FDR < .01)

write.table(results_fdr, file = "never_smoker_Results_fdr.txt")
#par(mfrow=1)
plotSmear(et, de.tags = rownames(results_edgeR)[results_edgeR$table$FDR < .05])
abline(h = c(-1, 1), col = "blue")

# differential expression patterns for gene "TP53" between two groups
boxplot(as.numeric(data_filtered["DSP", ]) ~ group)

