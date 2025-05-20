### load 0505.RData
library(estimate)
library(tidyestimate)
library(dplyr)
library(ggplot2)
library(ggsci)
library(ggpubr)
df_full = cbind(annot_pca$SYMBOL,count)
test_result = estimate_score(df_full, is_affymetrix=FALSE)
test_result$purity_inf = cos(0.61+0.00015*test_result$estimate)
hist(test_result$estimate)
idx_pat_sel = (test_result$estimate < 25000)
idx_new = idx & idx_pat_sel
sum(idx)
sum(idx_new)


to_draw_7 = merge(cor_rna_full, cor_patient_anna, by.x = "ENSEMBL",by.y = "...1")


count = data.frame(t(sort(table(to_draw_7$cor_rna),decreasing = TRUE)))

to_draw_8 = to_draw_7 %>%
  filter(!(cor_rna %in% count$Var2[1:22]))




g_sel = ggplot(to_draw_8, aes(x= cor_rna,y = cor)) +
  geom_point(alpha = 0.1,size = 0.5) +
  stat_cor(method="pearson",label.x = -0.4, label.y=0.9)+
  xlim(-1,1)+
  ylim(-1,1)+
  xlab("PDX model")+
  ylab("patient")+
  theme_classic() +
  coord_fixed()+
  scale_color_npg()+
  #theme(text = element_text( family = 'Arial'))+
  ggtitle("Pearson correlation coefficient related to pseudotime 150-220 \n(after purification)")
g_sel
ggsave("11.3_after_purity_selection.pdf",g_sel,height=4)



g = 
  ggplot(test_result,aes(x=estimate))+  geom_histogram(color="black",fill="white")+ theme_classic() +
  scale_color_npg()+ geom_vline(xintercept=25000,color="red")+
  #theme(text = element_text( family = 'Arial'))+
  xlab("Estimate Score")+
  ylab("Number of Patients")+
  ggtitle("Purity Distribution of Clinical Patient Data Range from \nPseudotime 150–220")
g
ggsave("11.7_estimate_select.pdf",g,height=4,width=4)


to_draw_highlight = meta_data
to_draw_highlight$highlight = idx_new
g = 
  ggplot(to_draw_highlight,aes(x=-PC1,y=PC2,fill=`Sample Type`,alpha=highlight))+geom_point(shape=21,color="black",size=2.5)+
  xlab("PC1: 22% variance")+
  ylab("PC2: 8% variance")+
  theme_classic() +
  guides(alpha = FALSE) +
  coord_fixed()+
  scale_fill_manual(values = c("PRIMARY" = "#e20917", "CRPC" = "#088ac8","NEPC" = "#8f1c5b","NORMAL"="#44984f"),
                     breaks = c("NORMAL","PRIMARY","CRPC","NEPC"))+
  ggtitle("Subset of Patients under Purity Selection ")

g
ggsave("selected_patient_in_PCA.pdf",height=4,width = 6)
