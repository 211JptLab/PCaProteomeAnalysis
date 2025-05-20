library(ggplot2)
library(gggenes)
library(dplyr)
library(readr)
library(ggrepel)
library(tidyr)
library(ggforce)
library(AnnotationDbi)

### process dataset ####
mean_0 = mean(167.52,169.87)
mean_8 = 180.34
mean_12 = 187.18
mean_32 = mean(220.14,220.61)

pro_ltl = read_delim("LTL331/wetransfer_ltl-progression-data_2022-03-29_1755_bkp/ltl.human.protome.gct_n10x11322.gct", 
                     delim = "\t", escape_double = FALSE, 
                     col_names = FALSE, trim_ws = TRUE, skip = 9)
colnames(pro_ltl) <- read_delim("LTL331/wetransfer_ltl-progression-data_2022-03-29_1755_bkp/ltl.human.protome.gct_n10x11322.gct",
                                skip=2,col_names = F, delim="\t",n_max = 1)
colnames(pro_ltl)[17:26]
pt_pro = c(mean_8,mean_0,mean_8,mean_0,mean_32,
           mean_12,mean_12,mean_32,mean_0,mean_32)
time_pro = c(8,0,8,0,32,12,12,32,0,32)


# 导入 limma 库
library(limma)

# 读取表达矩阵数据
countMatrix <- pro_ltl[,c(17,19,21,22,23,24,26)]

# 创建一个样本信息表格
sampleTable <- data.frame(
  condition = c("Weeks1", "Weeks1", "Weeks32","Weeks1","Weeks1","Weeks32",  "Weeks32"),
  replicate = c(1, 2, 1, 3,4, 2, 3)
)

# 创建设计矩阵
design <- model.matrix(~condition, data = sampleTable)

# 创建 limma 数据对象
fit <- lmFit(countMatrix, design)
fit <- eBayes(fit)
tab <- topTable(fit, n=Inf, coef=2,sort="none")


fc_pro_ltl = data.frame(site = pro_ltl$id,
                        NP = pro_ltl$id.description,
                        SYMBOL = pro_ltl$geneSymbol,
                        tab)
#write.csv(fc_pro_ltl,"fc_pro_32vsmid.csv",row.names = FALSE)


fc_pro_ltl = data.frame(site = pro_ltl$id,
                        NP = pro_ltl$id.description,
                        SYMBOL = pro_ltl$geneSymbol,
                        fc_pro = tab$logFC)




library(ggplot2)
library(gggenes)
library(dplyr)
library(readr)
library(ggrepel)
library(tidyr)
library(ggforce)
library(stringr)

### process dataset ####
mean_0 = mean(167.52,169.87)
mean_8 = 180.34
mean_12 = 187.18
mean_32 = mean(220.14,220.61)

## cor_pro_ltl

filenames_ltl <- list.files("./LTL331/tabs/", pattern="*.tab", full.names=TRUE)
rna_ltl <- data.frame()
for (i in filenames_ltl){
  item <- read_delim(i, skip=4,col_names = F,  delim="\t")
  if (nrow(rna_ltl) == 0){
    rna_ltl <- rbind(rna_ltl, item[,c(1,4)])
  }else{
    rna_ltl<- cbind(rna_ltl, item[,4])}
}

ttt = substr(filenames_ltl, start=16,stop=10000)
t1 = str_replace(ttt, "_step0_ReadsPerGene.out.tab", "")
colnames(rna_ltl) <- c("gene_id",t1)
sample_info_ltl <- read_delim("./LTL331/tabs/annot_ltl331.txt",delim = "\t", escape_double = FALSE, 
                              trim_ws = TRUE)

sample_info_ltl$PT = c(167.52,
                       169.87,
                       180.34,
                       187.18,
                       220.14,
                       220.61)


library(DESeq2)

# 读取表达矩阵数据
countMatrix <- rna_ltl[,c(4,5,6,7)]

# 创建一个样本信息表格
sampleTable <- data.frame(
  condition = c("Mid Weeks", "Mid Weeks","32 Weeks","32 Weeks"),
  replicate = c(1, 2,  1, 2)
)

# 创建 DESeq2 数据对象
dds <- DESeqDataSetFromMatrix(countData = countMatrix, colData = sampleTable, design = ~ condition)

# 执行标准化和差异表达分析
dds <- DESeq(dds)

# 根据条件比较获取差异表达基因列表
results <- results(dds, contrast = c("condition", "32 Weeks", "Mid Weeks"))

cor_rna_ltl = data.frame(id = rna_ltl$gene_id,
                         fc_rna = results$log2FoldChange)



library(org.Hs.eg.db)
annots <- AnnotationDbi::select(org.Hs.eg.db, keys=gsub("\\..*$", "" , cor_rna_ltl$id), 
                                columns="SYMBOL", keytype="ENSEMBL")
annots = na.omit((annots))
cor_rna_ltl$ENSEMBL = gsub("\\..*$", "" , cor_rna_ltl$id)
cor_rna_ltl_annot = merge(cor_rna_ltl, annots, by = "ENSEMBL", all=F)
fc_rna_ltl = na.omit((cor_rna_ltl_annot))        


#write.csv(fc_rna_ltl,"32vsmid_mRNA.csv",row.names = FALSE)
#write.csv(fc_pro_ltl,"32vsmid_pro.csv",row.names = FALSE)





fc_ltl_merge = merge(fc_rna_ltl,fc_pro_ltl,by = "SYMBOL",all=F)

cor_ltl_merge = fc_ltl_merge

library(ggstatsplot)
library(ggsci)
library(ggpubr)

highlight_df_1 <- cor_ltl_merge %>% 
  filter(SYMBOL %in% c("TOP2A","BIRC5","CBX2","KIT","FGF2","PAK1","WEE1","CBL","CDH2","DOT1L","USP1","AURKA","USP1","SOX2"))

highlight_df_2 <- cor_ltl_merge %>% 
  filter(SYMBOL %in% c("IDH1","CDH1","TRADD","HIPK2","PDCD4","EIF2AK4","FAS","HOXB13","AR"))

highlight_df_3 <- cor_ltl_merge %>% 
  filter(SYMBOL %in% c("JUN","KMT1A","HOXC6","COT1","CENPN","COP1","KDM1A","MVP","PHF3","POLA2","POLR2J","SMARCC2","YWHAE"))

highlight_df_5 <- highlight_df_3[c(-10,-12),]
highlight_df_4 <- cor_ltl_merge %>% 
  filter(SYMBOL %in% c("FAT4","NF1"))



highlight_df_6 <- cor_ltl_merge %>% 
  filter(SYMBOL %in% c("KLK3","KLK2","TMPRSS2","FKBP5","NKX3-1","PLPP1","PMEPA1","PART1","ALDH1A3","STEAP4"))

highlight_df_7 <- cor_ltl_merge %>% 
  filter(SYMBOL %in% c("SYP","CHGA","CHGB","ENO2","CHRNB2","SCG3","SCN3A","PCSK1","ELAVL4"))
# NKX2-1 in TFs
highlight_df_8 <- cor_ltl_merge %>% 
  filter(SYMBOL %in% c("SOX21","POU4F3","POUF4F1","POUF3F1","NKX2-1","FOXA2"))



highlight_df_1$group = "Up Set1/2"
highlight_df_2$group = "Down Set1/2"
# highlight_df_5$group = "3.Only Up in Protein"
# highlight_df_4$group = "4.Only Down in Protein"
highlight_df_6$group = "AR Score"
highlight_df_7$group = "NE Score"
# highlight_df_8$group = "7.TF"
# highlight_df = rbind(highlight_df_1,highlight_df_2,highlight_df_5,highlight_df_4,highlight_df_6,highlight_df_7,highlight_df_8)
highlight_df_8$group = "TF"
highlight_df = rbind(highlight_df_1,highlight_df_2,highlight_df_6,highlight_df_7,highlight_df_8)



#cbPalette <- c("red", "#008CCA", "#009C47", "#8D1F61")

#cbPalette <- c("red", "#11659a", "#69a794", "#fed71a","#d23918","#7c739f","#4c8045")
#cbPalette <- c("#55b7e6", "#193e8f", "#e53528", "#f09739", "#fed71a","#7c739f","#4c8045")
cbPalette <- c("#55b7e6", "#193e8f",  "#fed71a","#7c739f","#4c8045")
#cbPalette <- c("#55b7e6", "#193e8f", "#4c8045")

highlight_df_sel = highlight_df %>%
  filter(group=="TF" | (fc_rna >5 & fc_pro >0) | (fc_rna < -1 & fc_pro < 0) )



# Define cbPalette
cbPalette <- c("Up Set1/2" ="#55b7e6",
               "Down Set1/2" = "#193e8f", 
               "AR Score" = "orange",
               "NE Score" ="#7c739f",
               "TF" = "#4c8045")

# Convert group to factor with specified levels
highlight_df$group <- factor(highlight_df$group, levels = c("Up Set1/2", "Down Set1/2", "AR Score", "NE Score", "TF"))
highlight_df_sel$group <- factor(highlight_df_sel$group, levels = c("Up Set1/2", "Down Set1/2", "AR Score", "NE Score", "TF"))

g_ltl <- ggplot(cor_ltl_merge, aes(x=fc_rna, y=fc_pro)) + 
  geom_point(alpha=0.1) + 
  stat_cor(method = "pearson", label.x = -4, label.y = 4,size=6) + 
  theme_classic() +
  xlab("mRNA") +
  ylab("Protein") +
  geom_point(data=highlight_df, 
             aes(x=fc_rna, y=fc_pro, label=SYMBOL, fill=group), color="white", shape=21, size=2.5, show.legend=TRUE) +
  geom_label_repel(data=highlight_df_sel, 
                   aes(x=fc_rna, y=fc_pro, label=SYMBOL, color=group), size=5, 
                   max.overlaps=25, show.legend=FALSE) +
  scale_color_manual(values=cbPalette) +
  scale_fill_manual(values=cbPalette) +
  
  
  
  ggtitle("Fold change between 32w vs 8/12w")+
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        plot.title = element_text(size = 24),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18)) 
  
print(g_ltl)
ggsave("fc_mRNA_pro_ltl_32vsMid_new.pdf",height = 10,width=10,plot = g_ltl)
