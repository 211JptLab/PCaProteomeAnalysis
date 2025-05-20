library(readr)
library(ggvenn)
library(dplyr)
library(ggpubr)
library(ggsci)
library(ggrepel)
## get pt info from tmp ##

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

phos_18 <-  read_delim("./18PLEX/wetransfer_pdx-progression-files_2022-08-26_2220_backup/panoply_ptm_normalization_human.phospho.gct_n36x31538-proteome-relative-norm.gct",
                       skip=8,col_names = F, delim="\t")
colnames(phos_18) <- read_delim("./18PLEX/wetransfer_pdx-progression-files_2022-08-26_2220_backup/panoply_ptm_normalization_human.phospho.gct_n36x31538-proteome-relative-norm.gct",
                                skip=2,col_names = F, delim="\t",n_max = 1)
colnames(phos_18[34:67])
pt = rep(c(mean_145,mean_147,mean_147R,mean_173.1,mean_173.2,mean_176,
           mean_23,mean_23R,mean_LTL310,mean_LTL467,mean_LTL471,mean_49,
           mean_LTL508,mean_78,mean_93,mean_96,mean_96R),each=2)


cor_phos_18 = data.frame(site = phos_18$id,
                         NP = phos_18$original_accession_number,
                         SYMBOL = phos_18$geneSymbol,
                         Variable_Sites = phos_18$variableSites,
                         cor_phos = cor(t(phos_18[,34:67]),pt, use="pairwise.complete.obs"),
                         p.value = apply(t(phos_18[,34:67]),2,function(x) cor.test(x,pt)$p.value),
                         NAs = rowSums(is.na(phos_18)))





library(ggfortify)
df <- t(na.omit(phos_18[,34:67]))
pca_ltl=data.frame(df,sample=colnames(phos_18[,34:67]))
pca_ltl$group = substr(pca_ltl$sample, 1, nchar(pca_ltl$sample) - 5)




pca_ltl$group  = gsub("pdx","LuCaP-",pca_ltl$group )
pca_ltl$group  = gsub("lncap_intact","LNCaP",pca_ltl$group )
pca_ltl$group  = gsub("lncap_abl_cast","LNCaP_ABL",pca_ltl$group )
pca_ltl$group  = gsub("rv221","22Rv1",pca_ltl$group )
pca_ltl$group  = gsub("pnpca","PNPCa",pca_ltl$group )
pca_ltl$group  = gsub("LuCaP-310F","LTL-310",pca_ltl$group )
pca_ltl$group  = gsub("LuCaP-467B","LTL-467",pca_ltl$group )
pca_ltl$group  = gsub("LuCaP-508","LTL-508",pca_ltl$group )
pca_ltl$group  = gsub("LuCaP-471","LTL-471",pca_ltl$group )
pca_ltl$group  = gsub("pca1","MSK-PCa1",pca_ltl$group )
pca_ltl$pseudotime =pt
pca_res <- prcomp(df)

to_draw = pca_ltl[,c("group","pseudotime")]
to_draw_2 = unique(to_draw)
to_draw_2$Type = c("3","2","2","3","2","2","2",
                   "2","1","1","1","3","1","2",
                   "3","2","2")
color = c("#8f1c5b","#088ac8","#e20917")
  
g1 = ggplot(to_draw_2, aes(x = pseudotime,y=1,label=group,color=Type)) +
 
  geom_jitter(height=0.01,shape=16, size=12)+
  geom_label_repel(max.overlaps = 25,show.legend = FALSE,
                   segment.color = NA,  # 设置连接线颜色为透明
                   label.size = NA, 
                   color="black",
                   label.padding = unit(0.15, "lines"), # 移除边框
                   fill = NA  )+
  scale_x_continuous(limits = c(150, 250)) +
  scale_y_continuous(limits = c(0.9, 1.1)) +
  geom_hline(yintercept = 1, linetype = "solid", color = "black") +
  scale_color_manual(values = c("1" = "#e20917", "2" = "#088ac8","3" = "#8f1c5b")) +
  theme_minimal() +
  theme(
    legend.position = "none",  # 隐藏图例
    axis.title.y = element_text(angle = 0, vjust = 0.5,size=12),
    axis.text.y  = element_blank(),  # 隐藏y轴文本
    axis.ticks.y = element_blank(),
    
    panel.grid.major = element_blank(),  # 移除主要网格线
    panel.grid.minor = element_blank() # 隐藏y轴刻度
  )+ylab("Set1")+
  geom_vline(xintercept = 150, linetype = "solid", color = "black") +
  geom_vline(xintercept = 250, linetype = "solid", color = "black") 
  g1
ggsave("set1.pdf")




plex_11 <- read_delim("PILOT/wetransfer-c07d95_BCK/panoply_medianMAD_phosphoproteome.gct", 
                      delim = "\t", escape_double = FALSE, 
                      col_names = FALSE, trim_ws = TRUE, skip = 11)
colnames(plex_11) <- read_delim("PILOT/wetransfer-c07d95_BCK/panoply_medianMAD_phosphoproteome.gct",
                                skip=2,col_names = F, delim="\t",n_max = 1)
colnames(plex_11)[24:43]
pt = c(mean_147,mean_int,mean_78,mean_msk,
       mean_35,mean_abl,mean_pnp,mean_22V,
       mean_23,mean_145,mean_23,mean_35,
       mean_145,mean_pnp,mean_msk,mean_abl,
       mean_78,mean_22V,mean_147,mean_int)


cor_phos_10 = data.frame(site = plex_11$id,
                         NP = plex_11$id.description,
                         SYMBOL = plex_11$geneSymbol,
                         Variable_Sites = plex_11$variableSites,
                         cor_phos = cor(t(plex_11[,24:43]),pt, use="pairwise.complete.obs"),
                         p.value = apply(t(plex_11[,24:43]),2,function(x) {
                           if(sum(!is.na(x)) >3)
                           {(cor.test(x,pt)$p.value)}
                           else{return(1)}
                         }),
                         
                         NAs = rowSums(is.na(plex_11)))



library(ggfortify)
df <- t(na.omit(plex_11[,24:43]))
pca_ltl=data.frame(df,sample=colnames(plex_11[,24:43]))
pca_ltl$group = substr(pca_ltl$sample, 7, nchar(pca_ltl$sample) )




pca_ltl$group  = gsub("PDX","LuCaP-",pca_ltl$group )
pca_ltl$group  = gsub("lncap_intact","LNCaP",pca_ltl$group )
pca_ltl$group  = gsub("lncap_abl_cast","LNCaP_ABL",pca_ltl$group )
pca_ltl$group  = gsub("rv221","22Rv1",pca_ltl$group )
pca_ltl$group  = gsub("LuCaP-_PNPCa","PNPCa",pca_ltl$group )
pca_ltl$group  = gsub("LuCaP-310F","LTL-310",pca_ltl$group )
pca_ltl$group  = gsub("LuCaP-467B","LTL-467",pca_ltl$group )
pca_ltl$group  = gsub("LuCaP-508","LTL-508",pca_ltl$group )
pca_ltl$group  = gsub("LuCaP-471","LTL-471",pca_ltl$group )
pca_ltl$group  = gsub("PCa1","MSK-PCa1",pca_ltl$group )
pca_ltl$pseudotime =pt
pca_res <- prcomp(df)

to_draw = pca_ltl[,c("group","pseudotime")]
to_draw_2 = unique(to_draw)
to_draw_2$Type = c("2","2","2","3","2","2",
                   "1","2","2","3")
g2 = ggplot(to_draw_2, aes(x = pseudotime,y=1,label=group,color=Type)) +
  
  geom_jitter(height=0.01,shape=16, size=12)+
  geom_label_repel(max.overlaps = 25,show.legend = FALSE,
                   segment.color = NA,  # 设置连接线颜色为透明
                   label.size = NA, 
                   color="black",
                   label.padding = unit(0.15, "lines"), # 移除边框
                   fill = NA  )+
  scale_x_continuous(limits = c(150, 250)) +
  scale_y_continuous(limits = c(0.9, 1.1)) +
  geom_hline(yintercept = 1, linetype = "solid", color = "black") +
  scale_color_manual(values = c("1" = "#e20917", "2" = "#088ac8","3" = "#8f1c5b")) +
  theme_minimal() +
  theme(
    legend.position = "none",  # 隐藏图例
    axis.title.y = element_text(angle = 0, vjust = 0.5,size=12),
    axis.text.y  = element_blank(),  # 隐藏y轴文本
    axis.ticks.y = element_blank(),
    
    panel.grid.major = element_blank(),  # 移除主要网格线
    panel.grid.minor = element_blank() # 隐藏y轴刻度
  )+ylab("Set2")+
  geom_vline(xintercept = 150, linetype = "solid", color = "black") +
  geom_vline(xintercept = 250, linetype = "solid", color = "black") 
g2

ggsave("set2.pdf")


to_draw_2 = data.frame(group = c("0w","0w","8w","12w","32w","32w"),
                       pseudotime=c(167.52, 169.87, 180.34, 187.18, 220.14, 220.61)
)
to_draw_2$Type = c("1","1","2","2","3","3")
g3 = ggplot(to_draw_2, aes(x = pseudotime,y=1,label=group,color=Type)) +
  
  geom_jitter(height=0.01,shape=16, size=12)+
  geom_label_repel(max.overlaps = 25,show.legend = FALSE,
                   segment.color = NA,  # 设置连接线颜色为透明
                   label.size = NA, 
                   color="black",
                   label.padding = unit(0.15, "lines"), # 移除边框
                   fill = NA  )+
  scale_x_continuous(limits = c(150, 250)) +
  scale_y_continuous(limits = c(0.9, 1.1)) +
  geom_hline(yintercept = 1, linetype = "solid", color = "black") +
  scale_color_manual(values = c("1" = "#e20917", "2" = "#088ac8","3" = "#8f1c5b")) +
  theme_minimal() +
  theme(
    legend.position = "none",  # 隐藏图例
    axis.title.y = element_text(angle = 0, vjust = 0.5,size=12),
    axis.text.y  = element_blank(),  # 隐藏y轴文本
    axis.ticks.y = element_blank(),
    
    panel.grid.major = element_blank(),  # 移除主要网格线
    panel.grid.minor = element_blank() # 隐藏y轴刻度
  )+ylab("LTL331")+
  geom_vline(xintercept = 150, linetype = "solid", color = "black") +
  geom_vline(xintercept = 250, linetype = "solid", color = "black") 
g3
ggsave("set3.pdf")
