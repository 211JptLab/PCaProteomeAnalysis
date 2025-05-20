library(readr)
library(ggvenn)
library(dplyr)
library(ggpubr)
library(ggsci)
library(ggrepel)

mean_0 = mean(167.52,169.87)
mean_8 = 180.34
mean_12 = 187.18
mean_32 = mean(220.14,220.61)
pt_ltl = c()
phos_ltl = read_delim("LTL331/wetransfer_ltl-progression-data_2022-03-29_1755_bkp/ltl.human.phospho.gct_n10x33332-proteome-relative-norm.gct", 
                      delim = "\t", escape_double = FALSE, 
                      col_names = FALSE, trim_ws = TRUE, skip = 9)
colnames(phos_ltl) <- read_delim("LTL331/wetransfer_ltl-progression-data_2022-03-29_1755_bkp/ltl.human.phospho.gct_n10x33332-proteome-relative-norm.gct",
                                 skip=2,col_names = F, delim="\t",n_max = 1)

time_phos = c(32,8,8,12,12,
              0,0,0,32,32)

pt = c(mean_32,mean_8,mean_8,mean_12,mean_12,mean_0,mean_0,mean_0,mean_32,mean_32)


library(ggfortify)
df <- t(na.omit(phos_ltl[,36:45]))
pca_ltl=data.frame(df,sample=colnames(phos_ltl[,36:45]))
pca_ltl$group = c("32w","8w","8w","12w","12w",
                  "0w","0w","0w","32w","32w")


pca_ltl$group 

pca_ltl$pseudotime =pt
pca_res <- prcomp(df)
pca_res$x[, 1] <- pca_res$x[, 1] * -1

cor_pt = cor.test(pca_res$x[, 1],pt)
cor_pt$p.value
#0.001178834
cor_pt$estimate
# 0.8665099 

g = autoplot(pca_res,data = pca_ltl, colour = 'pseudotime',size=5)+theme(text = element_text( family = 'Arial'))+theme_classic() +
  geom_label_repel( 
    aes(label=group),size = 5,
    max.overlaps = 25,show.legend = FALSE)+
  scale_color_continuous(high = "#132B43", low = "#56B1F7")+
  ggtitle("PCA for Phosphorylation in LTL331")+
  annotate("text",label="cor: 0.87\np-value: 0.001",x=0,y=0.3,size=5)+
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        plot.title = element_text(size = 24),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18) ) 
#+guides(color = TRUE)
g
ggsave("11-20-ltl_phos_pca_new.pdf",g,width = 10,height=10)




