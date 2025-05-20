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




plex_11 <- read_delim("PILOT/wetransfer-c07d95_BCK/panoply_ptm_normalization_phosphoproteome.gct", 
                      delim = "\t", escape_double = FALSE, 
                      col_names = FALSE, trim_ws = TRUE, skip = 11)
colnames(plex_11) <- read_delim("PILOT/wetransfer-c07d95_BCK/panoply_ptm_normalization_phosphoproteome.gct",
                                skip=2,col_names = F, delim="\t",n_max = 1)
colnames(plex_11)[25:44]
pt = rep(c(mean_22V,mean_abl,mean_int,mean_msk,
           mean_pnp,mean_145,mean_147,mean_23,
           mean_35,mean_78),2)

# cor_test = function(x){
#   if(sum(!is.na(x) >3)
#      {(cor.test(x,pt)$p.value)}
#      else{return(1)}
# }
# }

cor_phos_10 = data.frame(site = plex_11$id,
                         NP = plex_11$original_accession_number,
                         SYMBOL = plex_11$geneSymbol,
                         Variable_Sites = plex_11$variableSites,
                         cor_phos = cor(t(plex_11[,25:44]),pt, use="pairwise.complete.obs"),
                         p.value = apply(t(plex_11[,25:44]),2,function(x) {
                           if(sum(!is.na(x)) >3)
                           {(cor.test(x,pt)$p.value)}
                           else{return(1)}
                         }),
                         
                         NAs = rowSums(is.na(plex_11)))



library(ggfortify)
df <- t(na.omit(plex_11[,25:44]))
pca_ltl=data.frame(df,sample=colnames(plex_11[,25:44]))
pca_ltl$group = substr(pca_ltl$sample, 7, nchar(pca_ltl$sample) )


pca_ltl$group  = gsub("LnCaP_Intact","LNCaP",pca_ltl$group )
pca_ltl$group  = gsub("LnCaP_abl_cast","LNCaP-abl",pca_ltl$group )
pca_ltl$group  = gsub("PDX","LuCaP-",pca_ltl$group )
pca_ltl$group  = gsub("rv221","22Rv1",pca_ltl$group )
pca_ltl$group  = gsub("LuCaP-_PNPCa","PNPCa",pca_ltl$group )
pca_ltl$group  = gsub("LuCaP-310F","LTL-310",pca_ltl$group )
pca_ltl$group  = gsub("LuCaP-467B","LTL-467",pca_ltl$group )
pca_ltl$group  = gsub("LuCaP-508","LTL-508",pca_ltl$group )
pca_ltl$group  = gsub("LuCaP-471","LTL-471",pca_ltl$group )
pca_ltl$group  = gsub("PCa1","MSK-PCa1",pca_ltl$group )
pca_ltl$group 

pca_ltl$pseudotime =pt
pca_res <- prcomp(df)
pca_res$x[, 1] <- pca_res$x[, 1] * -1

cor_pt = cor.test(pca_res$x[, 1],pt)
cor_pt$p.value
#0.01567162
cor_pt$estimate
# 0.5323804

g = autoplot(pca_res,data = pca_ltl, colour = 'pseudotime',size=5)+theme(text = element_text( family = 'Arial'))+theme_classic() +
  geom_label_repel( 
    aes(label=group),size = 5,
    max.overlaps = 25,show.legend = FALSE)+
  scale_color_continuous(high = "#132B43", low = "#56B1F7")+
  ggtitle("PCA for Phosphorylation in Set2")+
  annotate("text",label="cor: 0.53\np-value: 0.02",x=0.25,y=0.3,size=5)+
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        plot.title = element_text(size = 24),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18)) 
#+guides(color = TRUE)
g
ggsave("./11-19-10plex_phos_pca_new.pdf",g,width = 10,height=10)



