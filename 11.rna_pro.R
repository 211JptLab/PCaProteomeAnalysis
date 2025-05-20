### rna and pro
### cor_merge
library(dplyr)
library(ggrepel)
library(ggsci)
library(ggpubr)

highlight_df_1 <- cor_merge %>% 
  filter(SYMBOL %in% c("TOP2A","BIRC5","CBX2","KIT","FGF2","PAK1","WEE1","CBL","CDH2","DOT1L","USP1","AURKA","USP1","SOX2"))

highlight_df_2 <- cor_merge %>% 
  filter(SYMBOL %in% c("IDH1","CDH1","TRADD","HIPK2","PDCD4","EIF2AK4","FAS","HOXB13","TMPRSS2","AR","NKX3-1","PLPP1"))

highlight_df_3 <- cor_merge %>% 
  filter(SYMBOL %in% c("JUN","KMT1A","HOXC6","COT1","CENPN","COP1","KDM1A","MVP","PHF3","POLA2","POLR2J","SMARCC2","YWHAE"))

highlight_df_5 <- highlight_df_3[c(-8,-11,-14),]
highlight_df_4 <- cor_merge %>% 
  filter(SYMBOL %in% c("FAT4","NF1"))


highlight_df_1$group = "a"
highlight_df_2$group = "b"
highlight_df_5$group = "c"
highlight_df_4$group = "d"
highlight_df = rbind(highlight_df_1,highlight_df_2,highlight_df_5,highlight_df_4)
highlight_df[37,1] = "1433E"


cbPalette <- c("#55b7e6", "#193e8f", "#e53528", "#f09739")

g_m_2 <- ggplot(cor_merge, aes(x=cor_rna, y=cor_pro)) + geom_point(alpha=0.1) + 
  stat_cor(method = "pearson", label.x = -1, label.y = 0.9) +xlim(-1, 1)+ylim(-1,1) + 
  xlab("PDX model (mRNA)")+
  ylab("PDX model (Protein)")+
  theme_classic() +
  coord_fixed()+
  scale_color_npg()+
  geom_point(data=highlight_df, 
             aes(x=cor_rna,y=cor_pro,label=SYMBOL,color=group),size= 2.5,show.legend = FALSE)+
  geom_label_repel(data=highlight_df, 
                   aes(x=cor_rna,y=cor_pro,label=SYMBOL,color=group),size = 3.5,
                   max.overlaps = 25,show.legend = FALSE)+
  scale_colour_manual(values=cbPalette)+
  #theme(text = element_text( family = 'Arial'))+
  ggtitle("Pearson correlation coefficient related to pseudotime 150-220")
print(g_m_2)
ggsave("11.5-mRNA_pro.pdf",height = 6,plot = g_m_2)
