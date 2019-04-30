dicer <- read.table("subset DICER_other from Supplementary Table 3 and log2 mi53 column" header=TRUE)
novel <- read.table("Subset RNAseIIIB_S1344L_biallelic from Supplementary Table 3 and log2 mi53 column", header=TRUE)
hot <- read.table("Subset RNAseIIIB_hotspot  from Supplementary Table 3 and log2 mi53 column", header=TRUE)
mut <- read.table("Subset RNAseIIIB_hotspot_biallelic  from Supplementary Table 3 and log2 mi53 column", header=TRUE)
wt <- read.table ("Subset WT from Supplementary Table 3 and log2 mi53 column", header=TRUE)
dat <- read.table("Subset R944Q column from Supplementary Table 3 and log2 mi53 column", header=TRUE)

samples <- c("WT", "DICER_other", "RNAseIIIB_hotspot", "RNAseIIIB_hotspot_biallelic", "RNAseIIIB_S1344L_biallelic","R944Q") 
dat$new <- factor(dat$type, levels = c("WT", "DICER_other", "RNAseIIIB_hotspot", "RNAseIIIB_hotspot_biallelic", "RNAseIIIB_S1344L_biallelic","R944Q"))
hot$new <- factor(hot$type, levels = c("WT", "DICER_other", "RNAseIIIB_hotspot", "RNAseIIIB_hotspot_biallelic", "RNAseIIIB_S1344L_biallelic","R944Q"))
novel$new <- factor(novel$type, levels = c("WT", "DICER_other", "RNAseIIIB_hotspot", "RNAseIIIB_hotspot_biallelic", "RNAseIIIB_S1344L_biallelic","R944Q"))
mut$new <- factor(mut$type, levels = c("WT", "DICER_other", "RNAseIIIB_hotspot", "RNAseIIIB_hotspot_biallelic", "RNAseIIIB_S1344L_biallelic","R944Q"))
dicer$new <- factor(dicer$type, levels = c("WT", "DICER_other", "RNAseIIIB_hotspot", "RNAseIIIB_hotspot_biallelic", "RNAseIIIB_S1344L_biallelic","R944Q"))
wt$new <- factor(wt$type, levels = c("WT", "DICER_other", "RNAseIIIB_hotspot", "RNAseIIIB_hotspot_biallelic", "RNAseIIIB_S1344L_biallelic","R944Q"))


ggplot()+
  
  #dat, aes(x=dat$new, y=dat$median, color=dat$new)) + 
  geom_point(data=dat,aes(x=dat$new, y=dat$median, color=dat$new), position=position_jitter(width=0.25), alpha=0, size=2) +
  geom_point(data=cont, aes(x=cont$new, y=cont$median), position=position_jitter(width=0.2), alpha=1, shape=16, colour="CornflowerBlue", size=2) +
  geom_point(data=hot, aes(x=hot$new, y=hot$median), position=position_jitter(width=0.2), alpha=1, shape=16, colour="#E69F00", size=2) +
  geom_point(data=novel, aes(x=novel$new, y=novel$median), position=position_jitter(width=0.2), alpha=1, shape=16, colour="DeepPink", size=2) + 
  geom_point(data=mut, aes(x=mut$new, y=mut$median), position=position_jitter(width=0.2), alpha=1, shape=16, colour="#8b0000", size=2) + 
  geom_point(data=dicer, aes(x=dicer$new, y=dicer$median), position=position_jitter(width=0.2), alpha=1, shape=16, colour="aquamarine3", size=2) +
  geom_point(data=wt, aes(x=wt$new, y=wt$median), position=position_jitter(width=0.2), alpha=0.2, shape=16, colour="#999999", size=2) +
  geom_boxplot(data=dat,aes(x=dat$new, y=dat$median, color=dat$new),outlier.shape=NA, alpha=0, color=alpha("darkgray"))+ theme_bw() + 
  ylab ("log2(miR-5p/miR-3p)") + theme(legend.position="none") +
  #scale_x_discrete(labels = data$cancer) + 
  xlab ("") + ylim(-1.5,5) +
  #theme(axis.text.x = element_text(size=16))+
  theme(axis.text.y = element_text(size=16)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank())+
  #theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(axis.title=element_text(size=14,face="bold"))+
  #scale_color_manual(values=c("#999999", "aquamarine3","#E69F00","#8b0000", "DeepPink","CornflowerBlue"))+
  geom_hline(aes(yintercept=0),
             color="red", linetype="dashed", size=0.5)























