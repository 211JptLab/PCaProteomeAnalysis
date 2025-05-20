library(readr)
library(ggvenn)
library(dplyr)
library(ggpubr)
library(ggsci)
library(ggrepel)
library(limma)
library(edgeR)
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

protome_all <-  read_delim("PROGRESSION_MERGED/wetransfer_merged-prostate-progression_2022-12-09_1813_BCK/human.protome_batch_corrected.gct_n54x14382.gct",
                           skip=9,col_names = F, delim="\t")
colnames(protome_all) <- read_delim("PROGRESSION_MERGED/wetransfer_merged-prostate-progression_2022-12-09_1813_BCK/human.protome_batch_corrected.gct_n54x14382.gct",
                                    skip=2,col_names = F, delim="\t",n_max = 1)

pt_pro = c(mean_147,mean_int,mean_78,
           mean_msk,mean_35,mean_abl,
           mean_pnp,mean_22V,mean_23,
           
           mean_145,mean_23,mean_35,
           mean_145,mean_pnp,mean_msk,
           mean_abl,mean_78,mean_22V,
           
           mean_147,mean_int,mean_93,
           mean_147R,mean_23R,mean_23,
           mean_173.1,mean_LTL508,mean_96,
           
           mean_96R,mean_LTL471,mean_LTL310,
           mean_LTL467,mean_173.2,mean_176,
           mean_49,mean_78,mean_145,
           
           mean_147,mean_93,mean_147R,
           mean_23R,mean_23,mean_173.1,
           mean_LTL508,mean_96,mean_96R,
           
           mean_LTL471,mean_LTL310,mean_LTL467,
           mean_173.2,mean_176,mean_49,
           mean_78,mean_145,mean_147)

colnames(protome_all)



library(ggfortify)
df <- t(na.omit(protome_all[16:69]))
pca_ltl=data.frame(df,sample=colnames(protome_all[16:69]))
pca_ltl$group = substr(pca_ltl$sample, 1, nchar(pca_ltl$sample) - 8)

pca_ltl$group  = gsub("pdx","LuCaP-",pca_ltl$group )
pca_ltl$group  = gsub("lncap_intact","LNCaP",pca_ltl$group )
pca_ltl$group  = gsub("lncap_abl_cast","LNCaP-abl",pca_ltl$group )
pca_ltl$group  = gsub("rv221","22Rv1",pca_ltl$group )
pca_ltl$group  = gsub("pnpca","PNPCa",pca_ltl$group )
pca_ltl$group  = gsub("LuCaP-310F","LTL-310",pca_ltl$group )
pca_ltl$group  = gsub("LuCaP-467B","LTL-467",pca_ltl$group )
pca_ltl$group  = gsub("LuCaP-508","LTL-508",pca_ltl$group )
pca_ltl$group  = gsub("LuCaP-471","LTL-471",pca_ltl$group )
pca_ltl$group  = gsub("pca1","MSK-PCa1",pca_ltl$group )
pca_ltl$group 
  
pca_ltl$pseudotime =pt_pro
pca_res <- prcomp(df)


cor_pt = cor.test(pca_res$x[, 1],pt_pro)
cor_pt$p.value
# 2.350107e-25
cor_pt$estimate
# 0.9366894 

g = autoplot(pca_res,data = pca_ltl, colour = 'pseudotime',size=5)+theme(text = element_text( family = 'Arial'))+theme_classic() +
  geom_label_repel( 
    aes(label=group),size = 5,
    max.overlaps = 25,show.legend = FALSE)+
  scale_color_continuous(high = "#132B43", low = "#56B1F7")+
  
  ggtitle("PCA for Protein Set1 and Set2")+
  annotate("text",label="cor: 0.93\np-value: 2.3e-25",x=0.15,y=0.35,size=5)+
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        plot.title = element_text(size = 24),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18)) 
#+guides(color = TRUE)
g
ggsave("merged_pca_new.pdf",g,width = 10,height=10)


