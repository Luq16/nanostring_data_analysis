
#Function used from NanoNormIter
# imagingQC: Imaging quality flag
# bindingDensityQC: Binding density flag
# limitOfDetectionQC: limit detection flag
# positiveLinearityQC: Linear of positive controls flag
# RUV_total: Function for normalization


## dir is the directory tha contains the raw RCC files

annotation_df <- read.csv("case_study_annotations.csv")

# Create a data frame containing all the raw counts from all raw files.
library(NanoStringQCPro)
files_RCC = list.files(pattern = ".RCC")
n_row = nrow(readRcc(files_RCC[1])$Code_Summary)
n_col = length(files_RCC)
raw_expression = as.data.frame(matrix(nrow = n_row,ncol = length(files_RCC)+2)) # empty Datat frame
colnames(raw_expression)[1:2] = c('Gene','Class')
raw_expression[,1:2] = readRcc(files_RCC[1])$Code_Summary[,c(2,1)] # populate Gene and Class columns in raw_expression

# Create a data frame to store QC flags.
QC_DF = as.data.frame(matrix(nrow = length(files_RCC),ncol = 11))
colnames(QC_DF) = c('BCAC_ID','SampleID','Owner','Comments','Date','GeneRLF','SystemAPF','imagingQC',
                   'bindingDensityQC','limitOfDetectionQC','positiveLinearityQC')# empty data frame



library(NanoNormIter)

#
for (i in 1:length(files_RCC)){
  rcc = readRcc(files_RCC[i])
  raw = rcc$Code_Summary
  raw_expression[,i+2] = as.numeric(raw$Count)
  colnames(raw_expression)[i+2] = strsplit(files_RCC[i],'_')[[1]][1] # populate columns raw_expression with raw counts
  
  # populate QC_DF  with QC flags
  QC_DF[i,2:7] = as.vector(rcc$Sample_Attributes)
  QC_DF$imagingQC[i] = imagingQC(rcc)
  QC_DF$bindingDensityQC[i] = bindingDensityQC(rcc,.05,2.25)
  QC_DF$limitOfDetectionQC[i] = limitOfDetectionQC(rcc)
  QC_DF$positiveLinearityQC[i] = positiveLinQC(rcc)
  
}

#write.csv(raw_expression, "raw_expression.csv")

raw_count = raw_expression[,-c(1:2)]
fData = raw_expression[,c(1:2)]
rownames(raw_count) = fData$Gene # make gene names rowname
HK_genes <- fData$Gene[fData$Class == "Housekeeping"]
QC_DF$HK_Gene_Miss = colSums(raw_count[HK_genes,] == 0)
rownames(fData) = fData$Gene # make gene names rowname
rownames(raw_count) = fData$Gene # make gene names rowname
rownames(QC_DF) = colnames(raw_count) # make sampleid rowname
QC_DF$SampleID = colnames(raw_count)
QC_DF$subjectid = as.factor(rep(c('1','2',
                                  '3','4', '5'),each = 2))
QC_DF$Group = annotation_df$visit



#### Housekeeping control check
library(MASS)
hk_raw = raw_count[HK_genes,]
pval = vector(length = nrow(hk_raw))

for (i in 1:nrow(hk_raw)){
  
  reg = glm.nb(as.numeric(hk_raw[i,]) ~ as.factor(QC_DF$Group))
  pval[i] = coef(summary(reg))[2,4]
  
}

sum(pval <= .05)
####

#### Normalization 
k = 1
vsd = RUV_total(raw_count,QC_DF,fData,k = k)$vsd
set = RUV_total(raw_count,QC_DF,fData,k = k)$set

norm_data <- assay(vsd)
b4_norm <- data.matrix(log2(raw_count[,1:10])) # log transform

#### /Normalization 

####Visualization.

##RLE Plot
library(EDASeq)
#png("RLE_plot.png", width =  8*700, height =  8*400, res = 600)
line = 1
cex = 2
side = 3
adj=-0.05
par(mfrow=c(2,1))
par(cex=0.6)
EDASeq::plotRLE(norm_data, col=rep(1:10, each=1))
mtext("A", side=side, line=line, cex=cex, adj=adj)
EDASeq::plotRLE(b4_norm, col=rep(1:10, each=1))
mtext("B", side=side, line=line, cex=cex, adj=adj)

#dev.off()


# par(mfrow=c(1,2))
# par(cex=0.6)
# plotPCA(unnorm_data, col=rep(1:2, each=5),k=2, isLog=TRUE)
# plotPCA(norm_data, col=rep(1:2, each=5), k=2, isLog=TRUE)

#### /Visualization.

#### Solution to task 3.1.2

## Heat map
library(dplyr)
library(gplots)
library(RColorBrewer)

heat_map <- data.frame(raw_expression$Class, norm_data)
heat_map <- filter(heat_map, heat_map$raw_expression.Class == "Positive" |
                     heat_map$raw_expression.Class == "Negative")
mat_data <- data.matrix(heat_map[,2:ncol(heat_map)]) %>%
  t()
my_palette <- colorRampPalette(c("purple", "#FFFFFF", "orange"))(n = 299)

#png("heat_map.png", width =  8*480, height =  8*480, res = 480)
heatmap.2(mat_data, notecol="black", margins =c(6,9), density.info="none", trace="none", 
          col=my_palette, dendrogram="row", Colv="NA",
          sepcolor="white", sepwidth=c(0.000001,0.0000001), colsep=1:ncol(mat_data),
          rowsep=1:nrow(mat_data),keysize=2,
          ylab = "Genes", xlab = "Samples")
#dev.off()

#### /Solution to task 3.1.2


#### Solution to task 3.2.2
##Boxplot
library(ggplot2)
library(tidyverse)

norm_data_filter <- norm_data[row.names(norm_data) %in% c("MCL1","CXCL1"),] %>%
  t()
file<-read.csv("case_study_annotations.csv")
time_point <- file$visit
DF_new <- data.frame(time_point, norm_data_filter)
DF_new <- gather(DF_new, key = "genes", value = "expression", -time_point)

p <- DF_new %>%
  ggplot(aes(x = genes, y = expression, fill = genes)) +
  stat_boxplot(geom = "errorbar", width = 0.5) +  
  geom_boxplot() +
  stat_summary(fun = mean, geom = "errorbar",aes(ymax = ..y.., ymin = ..y..),
               width = 0.75, size = 0.5, linetype = "solid", colour = "blue") +
  facet_wrap(facets = vars(time_point)) +
  ylab("Expression level") +
  labs(title = "Box plot for MCL1 and CXCL1")

#png("boxPlot.png", width =  8*480, height =  8*300, res = 600)
p
#dev.off()
## Boxplot statistc values
boxPlot_value <- data.frame(ggplot_build(p)$data)
boxPlot_stat <- boxPlot_value[, 2:6]
#write.csv(boxPlot_stat, "boxPlotStat.csv")

#### /Solution to task 3.2.2
