library(ggplot2)
library(ggpubr)
library(reshape2)
ggplot(melt(table(tcga_clinical[,c("CANCER_TYPE_ACRONYM","dicer_muts")])), 
       aes(x=factor(CANCER_TYPE_ACRONYM,levels=rownames(data.frame(sort(rowSums(table(tcga_clinical[,c("CANCER_TYPE_ACRONYM","dicer_muts")])[,c(1,2,3,5,6)]) /
                                                                          rowSums(table(tcga_clinical[,c("CANCER_TYPE_ACRONYM","dicer_muts")])),decreasing = T))),
                    labels = paste(rownames(data.frame(sort(rowSums(table(tcga_clinical[,c("CANCER_TYPE_ACRONYM","dicer_muts")])[,c(1,2,3,5,6)]) /
                                                              rowSums(table(tcga_clinical[,c("CANCER_TYPE_ACRONYM","dicer_muts")])),decreasing = T))),": (n=",
                                   tcga_numbers[tcga_numbers$Var1 %in% rownames(data.frame(sort(rowSums(table(tcga_clinical[,c("CANCER_TYPE_ACRONYM","dicer_muts")])[,c(1,2,3,5,6)]) /
                                                                                      rowSums(table(tcga_clinical[,c("CANCER_TYPE_ACRONYM","dicer_muts")])),decreasing = T))),][order(tcga_numbers$Var1),2],")",sep = "")),
           y=value,fill=factor(dicer_muts, levels = c("none",'bialelic','multi_with_hotspot','only_hotspot',"multi","one_mut")))) + theme_pubr()+
  geom_bar(stat = 'identity', position = 'fill') + theme(axis.text.x = element_text(angle=90)) + 
  scale_y_continuous(expand = c(0,0), labels = scales::percent,breaks = c(0,0.05,.1,.15,.25,.50,.75)) +
  scale_fill_manual(values=c('white',"#f03b20","#feb24c", "#ffeda0","#a6bddb","#ece7f2")) + ylab("% Dicer Mutatation (per tumor subytpe)") + 
  theme(legend.position = c(0.5,0.5)) + xlab("") 
