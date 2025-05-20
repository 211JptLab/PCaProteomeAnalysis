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
# pca_ltl$group  = gsub("lncap_intact","LNCaP",pca_ltl$group )
# pca_ltl$group  = gsub("lncap_abl_cast","LNCaP_ABL",pca_ltl$group )
# pca_ltl$group  = gsub("rv221","22Rv1",pca_ltl$group )
# pca_ltl$group  = gsub("pnpca","PNPCa",pca_ltl$group )
pca_ltl$group  = gsub("LuCaP-310F","LTL-310",pca_ltl$group )
pca_ltl$group  = gsub("LuCaP-467B","LTL-467",pca_ltl$group )
pca_ltl$group  = gsub("LuCaP-508","LTL-508",pca_ltl$group )
pca_ltl$group  = gsub("LuCaP-471","LTL-471",pca_ltl$group )
# pca_ltl$group  = gsub("pca1","MSK-PCa1",pca_ltl$group )
pca_ltl$group 

pca_ltl$pseudotime =pt
pca_res <- prcomp(df)
pca_res$x[, 1] <- pca_res$x[, 1] * -1

cor_pt = cor.test(pca_res$x[, 1],pt)
cor_pt$p.value
#5.492384e-09
cor_pt$estimate
# 0.812218 

g = autoplot(pca_res,data = pca_ltl, colour = 'pseudotime',size=5)+theme(text = element_text( family = 'Arial'))+theme_classic() +
  geom_label_repel( 
    aes(label=group),size = 5,
    max.overlaps = 25,show.legend = FALSE)+
  scale_color_continuous(high = "#132B43", low = "#56B1F7")+
  ggtitle("PCA for Phosphorylation in Set1")+
  annotate("text",label="cor: 0.81\np-value: 5.5e-9",x=-0.3,y=0.3,size=5)+
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        plot.title = element_text(size = 24),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18)) 
#+guides(color = TRUE)
g
ggsave("./11-18-18plex_phos_pca_new.pdf",g,width = 10,height=10)




