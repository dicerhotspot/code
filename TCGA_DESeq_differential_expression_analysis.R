# Load DESeq2 package

source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
library(DESeq2)
require(DESeq2)

### high_conserved_all_final.txt file contains all hotspot mutations, and WT endometrial cancer cases miRNA expression data ##

countdata.all <- read.table("high_conserved_all_final.txt", header=TRUE, row.names = 1)
countdata <- countdata.all[c(1:563)]

# Remove .bam or .sam from filenames
colnames(countdata) <- gsub("\\.[sb]am$", "", colnames(countdata))
colnames(countdata) <- gsub(".13_mirna_gdc_realn", "", colnames(countdata))

countdata <- ceiling(countdata)
countdata <- as.matrix(countdata)

# change the number in the mutant column for biallelic, non-biallelic and random comparison

(condition <- factor(c(rep("mutant",15), rep("control", 548))))
(treatment <- factor(c(rep("mutant",15), rep("control",548))))

# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
(coldata <- data.frame(row.names=colnames(countdata), condition, treatment))
coldata$condition <- relevel(coldata$condition, "control")
coldata$treatment <- relevel(coldata$treatment, "control")
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design= ~treatment)
dds

dds <- estimateSizeFactors(dds)

dds <- DESeq(dds)

# Get fold-chage and P-values

##### Get resdata for endometrial cancer #####

res <- results(dds)
table(res$padj<0.05)

## Order by adjusted p-value
res <- res[order(res$padj), ]

## Merge with normalized count data
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Geneid"
head(resdata)
resdata$logbaseMean <- with(resdata, log(resdata$baseMean,2))
resdata$log <- -log10(resdata$padj)

# Now separate 5p and 3p arms from the Geneid column 

threep <- resdata[grep("3p", resdata$Geneid),]
threep$threename <- threep$Geneid
fivep <- resdata[grep("5p", resdata$Geneid),]
fivep$fivename <- fivep$Geneid

# Let's make volcano plots; comment out 5p or 3p depending on which plot to generate
# Before that, determine FDR from the Benjamini-Hochberg correction done by DESeq analysis

bhc_5p <- subset(fivep, padj < 0.01)

bhc_3p <- subset(threep, padj < 0.01)

# Now look for bhc_5p$log in the subset to determine FDR cutoff for 5p 

# Now look for bhc_3p$log in the subset to determine FDR cutoff for 3p 


# Let's make volcano plot for 5p differential expression 

ggplot()+
  geom_point(data=fivep, aes(x=fivep$log2FoldChange, y=fivep$log), size=1.5, colour="#00FF00", alpha=1) + theme_bw()+
  geom_text(data=fivep, aes(x=fivep$log2FoldChange, y=fivep$log), label=fivep$Geneid, colour="white", alpha=0.0) +
  geom_text(data=threep, aes(x=threep$log2FoldChange, y=threep$log), label=threep$Geneid, colour="white", alpha=0.0) +
  theme(axis.text.x = element_text(size=16))+
  theme(axis.text.y = element_text(size=16))+
  xlim(-6,6)+ ylim(0,50)+
  geom_hline(aes(yintercept=2.029),color="black", linetype="dashed", size=0.5)+
  annotate(geom="text", x=-4.5, y=3.5, label = paste("FDR=0.01"))+
  xlab("Log2FoldChage") + ylab("-log10 P-value") +
  # title can bet commented out, not always necessary
  ggtitle ("High_conserved_227+53(15) hotspot vs control(548)_5p")


# Let's make volcano plot for 3p differential expression

ggplot()+
  geom_point(data=threep, aes(x=threep$log2FoldChange, y=threep$log), size=1.5, colour="#0000FF", alpha=1) + theme_bw()+
  geom_text(data=fivep, aes(x=fivep$log2FoldChange, y=fivep$log), label=fivep$Geneid, colour="white", alpha=0.0) +
  geom_text(data=threep, aes(x=threep$log2FoldChange, y=threep$log), label=threep$Geneid, colour="white", alpha=0.0) +
  theme(axis.text.x = element_text(size=16))+
  theme(axis.text.y = element_text(size=16))+
  xlim(-6,6)+ ylim(0,50)+
  geom_hline(aes(yintercept=2.008),color="black", linetype="dashed", size=0.5)+
  annotate(geom="text", x=-4.5, y=3.5, label = paste("FDR=0.01"))+
  xlab("Log2FoldChage") + ylab("-log10 P-value") +
  # title can bet commented out, not always necessary
  ggtitle ("High_conserved_227+53(15) hotspot vs control(548)_5p")


# Let's make 5p and 3p distribution plots


# First let's load multiplot function to plot multiple panels

############### Multiplot function ######################

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#################################################################################

fivepmedian <- median(fivep$log2FoldChange, na.rm=TRUE)
p1 <- ggplot()+
  geom_histogram(data=fivep, aes(x=fivep$log2FoldChange), binwidth=0.2, size=0.1, colour="black", alpha=1, fill="#00FF00") + theme_bw()+
  geom_vline(data=fivep, aes(xintercept = median(fivep$log2FoldChange, na.rm=TRUE)),col='red',size=1)+
  annotate(geom="text", x=2, y=25, label = paste("median=", sprintf("%.3f", fivepmedian)))+
  xlim(-6,6)+ xlab("log2FoldChange") + 
  theme(axis.text.x = element_text(size=16))+
  theme(axis.text.y = element_text(size=16))


threepmedian <- median(threep$log2FoldChange, na.rm=TRUE)
p2 <- ggplot()+
  geom_histogram(data=threep, aes(x=threep$log2FoldChange), binwidth=0.2, size=0.1, colour="black", alpha=1, fill="#0000FF") + theme_bw()+
  geom_vline(data=threep, aes(xintercept = median(threep$log2FoldChange, na.rm=TRUE)),col='red',size=1)+
  annotate(geom="text", x=2, y=25, label = paste("median=", sprintf("%.3f", threepmedian)))+
  xlim (-6,6) + xlab("log2FoldChange") +
  theme(axis.text.x = element_text(size=16))+
  theme(axis.text.y = element_text(size=16))

multiplot(p1,p2, cols=1)

## Let's do Wilcox rank sum test 

wilcox.test(fivep$log2FoldChange,threep$log2FoldChange)
