library(readr)
sample_info_all <- read_delim("./PROGRESSION_MERGED/annotations_allmodels.txt", 
                              delim = "\t", escape_double = FALSE, 
                              trim_ws = TRUE)
mean_abl = mean(sample_info_all$PT[1:3])
mean_int = mean(sample_info_all$PT[4:6])
mean_145 = mean(sample_info_all$PT[7:9])
mean_147 = mean(sample_info_all$PT[10:12])
mean_23 = mean(sample_info_all$PT[13:15])
mean_35 = mean(sample_info_all$PT[16:18])
mean_78 = mean(sample_info_all$PT[19:20])
mean_pnp = mean(sample_info_all$PT[21:23])
mean_msk = mean(sample_info_all$PT[24:25])
mean_23R = mean(sample_info_all$PT[26:28])
mean_49 = mean(sample_info_all$PT[29:30])
mean_93 = mean(sample_info_all$PT[31:32])
mean_96R = mean(sample_info_all$PT[33:35])
mean_147R = mean(sample_info_all$PT[36:37])
mean_173.1 = mean(sample_info_all$PT[38:39])
mean_173.2 = mean(sample_info_all$PT[40:41])
mean_176 = mean(sample_info_all$PT[42:43])
mean_96 = mean(sample_info_all$PT[44:45])
mean_LTL310 = mean(sample_info_all$PT[46])
mean_LTL467 = mean(sample_info_all$PT[47])
mean_LTL471 = mean(sample_info_all$PT[48])
mean_LTL508 = mean(sample_info_all$PT[49])
mean_22V = mean(sample_info_all$PT[50:51])                          


pt_11 = c(mean_147,mean_int,mean_78,
          mean_msk,mean_35,mean_abl,
          mean_pnp,mean_22V,mean_23,
          mean_145,mean_23,mean_35,
          mean_145,mean_pnp,mean_msk,
          mean_abl,mean_78,mean_22V,
          mean_147,mean_int)


pt_18 = c(mean_93,mean_147R,mean_23R,mean_23,mean_173.1,
          mean_LTL508,mean_96,mean_96R,mean_LTL471,
          mean_LTL310,mean_LTL467,mean_173.2,
          mean_176,mean_49,mean_78,mean_145,mean_147,
          mean_93,mean_147R,mean_23R,mean_23,
          mean_173.1,mean_LTL508,mean_96,mean_96R,mean_LTL471,
          mean_LTL310,mean_LTL467,mean_173.2,mean_176,
          mean_49,mean_78,mean_145,mean_147)




library(readr)
library(ggvenn)
library(dplyr)
library(ggpubr)
library(ggsci)
library(ggrepel)
library(ggvenn)
library(eulerr)
## get pt info from tmp ##

plex_11 <- read_delim("PILOT/wetransfer-c07d95_BCK/panoply_medainMAD_proteome.gct", 
                      delim = "\t", escape_double = FALSE, 
                      col_names = FALSE, trim_ws = TRUE, skip = 11)
plex_18 <- read_delim("18PLEX/wetransfer_pdx-progression-files_2022-08-26_2220_backup/human.protome.gct_n36x10816.gct",  delim = "\t", escape_double = FALSE, 
                      col_names = FALSE, trim_ws = TRUE, skip = 8)

plex_11$cor = cor((t(plex_11[,19:38])),pt_11,use="pairwise.complete.obs")
plex_18$cor = cor((t(plex_18[,c(17:33,35:51)])),pt_18,use="pairwise.complete.obs")


plex_11$id = gsub("\\..*$", "" , plex_11$X1)
plex_18$id = gsub("\\..*$", "" , plex_18$X1)

a <- list(`Plex 10` = plex_11$id ,
          
          `Plex 18` = plex_18$id )



ggvenn(a)+ggtitle("protein overlap")+ theme(
  panel.background = element_rect(fill = "white"),
  plot.background = element_rect(fill = "white")
)


dat<-c("Batch1" = 2351, 
       "Batch2" = 3029,
       "Batch1&Batch2" = 8465)
png("batch1vsbatch2_2.0.png")
v2 = plot(euler(dat),
     fills = list(fill = c("red", "#008CCA", "#009C47"), alpha = 0.9),
     quantities = list(type = c("counts","percent"),cex=1.5),
     #labels = list(font=3,side="above"),
     legend = list(side = "right",cex=1.5)
)
tags <- v2$children[[2]]$children[[1]]$children$tags$children
tags <- do.call(grid::gList, lapply(tags, function(x) {
  x$children[[2]]$label <- sub(" \\(", "\n(", x$children[[2]]$label)
  x$children[[2]]$just <- NULL
  x$children[[2]]$hjust <- 0.5
  x$children[[2]]$vjust <- 1
  x}))

v2$children[[2]]$children[[1]]$children$tags$children <- tags

v2

dev.off()

#ggsave("protein_overlap.png")

plex_merged = merge(plex_11,plex_18,by="id")


highlight_df_1 <- plex_merged %>% 
  filter(X3.x %in% c("TOP2A","BIRC5","CBX2","KIT","FGF2","PAK1","WEE1","CBL","CDH2","DOT1L","USP1","AURKA","USP1","SOX2"))

highlight_df_2 <- plex_merged %>% 
  filter(X3.x %in% c("IDH1","CDH1","TRADD","HIPK2","PDCD4","EIF2AK4","FAS","HOXB13","TMPRSS2","AR","NKX3-1","PLPP1"))

highlight_df_3 <- plex_merged %>% 
  filter(X3.x %in% c("JUN","KMT1A","HOXC6","COT1","CENPN","COP1","KDM1A","MVP","PHF3","POLA2","POLR2J","SMARCC2","YWHAE"))

highlight_df_5 <- highlight_df_3[c(-8,-11,-14),]
highlight_df_4 <- plex_merged %>% 
  filter(X3.x %in% c("FAT4","NF1"))


highlight_df_1$group = "a"
highlight_df_2$group = "b"
highlight_df_5$group = "c"
highlight_df_4$group = "d"
highlight_df = rbind(highlight_df_1,highlight_df_2,highlight_df_5,highlight_df_4)
highlight_df[30,4] = "1433E"



cbPalette <- c("#55b7e6", "#193e8f", "#e53528", "#f09739")



g_m_2 <- ggplot(plex_merged, aes(x=cor.y, y=cor.x)) + geom_point(alpha=0.1) + 
  stat_cor(method = "pearson", label.x = -0.4, label.y = 0.9) +xlim(-1, 1)+ylim(-1,1) + 
  
  theme_classic() +
  coord_fixed()+
  xlab("Set1 (Protein)")+
  ylab("Set2 (Protein)") +
  scale_color_npg()+
  geom_point(data=highlight_df, 
             aes(x=cor.y,y=cor.x,label=X3.x,color=group),size= 2.5,show.legend = FALSE)+
  geom_label_repel(data=highlight_df, 
                   aes(x=cor.y,y=cor.x,label=X3.x,color=group),size = 3.5,
                   max.overlaps = 25,show.legend = FALSE)+
  scale_colour_manual(values=cbPalette)+
  #theme(text = element_text( family = 'Arial'))+
  ggtitle("Pearson correlation coefficient related to pseudotime 150-220")

g_m_2
ggsave("11.8-pro_2_Batch.pdf", height = 6,plot = g_m_2)
