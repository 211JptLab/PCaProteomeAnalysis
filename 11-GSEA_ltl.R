#load 0505.RData
#load to_draw_22
dotPlotListGSEA = function(res_list, out_file, plot_dir="./plots/"){
  all_res = NULL
  for ( k in names(res_list)){
    tmp = res_list[[k]]@result[res_list[[k]]@result$p.adjust < 0.05,]
    if ( nrow(tmp) > 0 ){
      dir.create(plot_dir, showWarnings = F)
      tmp= tmp[order(tmp$p.adjust, decreasing = T),]
      tmp$Description = gsub(pattern = ".*\r: " ,replacement = "" , tmp$Description)
      tmp$Description = gsub(pattern = "_" ,replacement = " " , tmp$Description)
      tmp$Description = gsub(pattern = "^(.{10,30}) " ,replacement = "\\1\n" , tmp$Description)
      tmp$Description = gsub(pattern = "\n(.{10,30}) " ,replacement = "\n\\1\n" , tmp$Description)
      tmp$Description = gsub(pattern = "\n(.{10,30}) " ,replacement = "\n\\1\n" , tmp$Description)
      tmp$Name = k
      all_res = rbind(all_res, tmp)
    } 
  }
  if ( !is.null(all_res) ){
    all_res = all_res[order(all_res$p.adjust),]
    all_res$Name = factor(all_res$Name, levels=names(res_list))
    IDS= unique(all_res$ID)
    pdf(paste0(plot_dir, out_file, ".pdf"), width = max(10, length(names(res_list))))
    ydiv = 20
    for ( i in seq(1 , length(IDS), by=ydiv)){
      to_use = IDS[i:min((i+ydiv-1), length(IDS))]
      
      p <- ggplot(all_res[all_res$ID %in% to_use,])+ geom_point(
        aes(x=Name, y=Description , color=enrichmentScore, size=-log10(p.adjust) ))+ scale_y_discrete(limits=rev)+scale_x_discrete(drop=F)+
        scale_color_gradient2(low = "blue", high = "red", limits = c(min(all_res$enrichmentScore), max(all_res$enrichmentScore)))
      print(p)
    }
    dev.off()
    write.csv(all_res, paste0(plot_dir, out_file, ".csv"), row.names = F)  
  }
}



##### set up ####

h_groups = rbind(
  data.frame(group = "Inflammatory" , hallmark=c(
    "HALLMARK_INFLAMMATORY_RESPONSE",
    "HALLMARK_ALLOGRAFT_REJECTION",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "HALLMARK_INTERFERON_ALPHA_RESPONSE",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "HALLMARK_IL6_JAK_STAT3_SIGNALING",
    "HALLMARK_IL2_STAT5_SIGNALING",
    "HALLMARK_COMPLEMENT",
    "HALLMARK_COAGULATION"
  )  ), data.frame(group = "Damage Response" , hallmark=c(
    "HALLMARK_P53_PATHWAY",
    "HALLMARK_APOPTOSIS",
    "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
    "HALLMARK_DNA_REPAIR",
    "HALLMARK_UV_RESPONSE_UP",
    "HALLMARK_UV_RESPONSE_DN"
  )))

h_groups = rbind(h_groups, data.frame(group = "Cell Proliferation" , hallmark=c(
  "HALLMARK_MITOTIC_SPINDLE",
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_MYC_TARGETS_V2"
)))

h_groups = rbind(h_groups, data.frame(group = "Cell Metabolism" , hallmark=c(
  "HALLMARK_GLYCOLYSIS",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_HYPOXIA",
  "HALLMARK_PEROXISOME",
  "HALLMARK_FATTY_ACID_METABOLISM",
  "HALLMARK_CHOLESTEROL_HOMEOSTASIS",
  "HALLMARK_XENOBIOTIC_METABOLISM",
  "HALLMARK_BILE_ACID_METABOLISM",
  "HALLMARK_HEME_METABOLISM"
)))

h_groups = rbind(h_groups, data.frame(group = "Cell Junction" , hallmark=c(
  "HALLMARK_APICAL_JUNCTION",
  "HALLMARK_APICAL_SURFACE",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"
)))

h_groups = rbind(h_groups, data.frame(group = "Differentiation" , hallmark=c(
  "HALLMARK_ANGIOGENESIS",
  "HALLMARK_ADIPOGENESIS",
  "HALLMARK_MYOGENESIS",
  "HALLMARK_SPERMATOGENESIS",
  "HALLMARK_PANCREAS_BETA_CELLS"
)))

h_groups = rbind(h_groups, data.frame(group = "Signaling" , hallmark=c(
  "HALLMARK_ANDROGEN_RESPONSE",
  "HALLMARK_ESTROGEN_RESPONSE_EARLY",
  "HALLMARK_ESTROGEN_RESPONSE_LATE",
  "HALLMARK_NOTCH_SIGNALING",
  "HALLMARK_WNT_BETA_CATENIN_SIGNALING",
  "HALLMARK_TGF_BETA_SIGNALING",
  "HALLMARK_HEDGEHOG_SIGNALING",
  "HALLMARK_MTORC1_SIGNALING",
  "HALLMARK_PI3K_AKT_MTOR_SIGNALING",
  "HALLMARK_KRAS_SIGNALING_UP",
  "HALLMARK_KRAS_SIGNALING_DN",
  "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
  "HALLMARK_PROTEIN_SECRETION"
)))

dotPlotList2 = function(res_list){
  all_res = NULL
  out=list()
  for ( k in names(res_list)){
    tmp = res_list[[k]]@result
    if ( nrow(tmp) > 0 ){
      tmp= tmp[order(tmp$p.adjust, decreasing = T),]
      tmp$Description = gsub(pattern = ".*\r: " ,replacement = "" , tmp$Description)
      tmp$Description = gsub(pattern = "_" ,replacement = " " , tmp$Description)
      tmp$Description = gsub(pattern = "^(.{10,30}) " ,replacement = "\\1\n" , tmp$Description)
      tmp$Description = gsub(pattern = "\n(.{10,30}) " ,replacement = "\n\\1\n" , tmp$Description)
      tmp$Description = gsub(pattern = "\n(.{10,30}) " ,replacement = "\n\\1\n" , tmp$Description)
      tmp$GeneRatio = unlist(lapply(tmp$GeneRatio, function(x){ as.numeric(strsplit(x, "/")[[1]][[1]]) /as.numeric(strsplit(x, "/")[[1]][[2]])  }))
      tmp$Name = k
      all_res = rbind(all_res, tmp)
    } 
  }
  if ( !is.null(all_res) ){
    all_res = all_res[order(all_res$p.adjust),]
    all_res$Name = factor(all_res$Name, levels=names(res_list))
    IDS= unique(all_res$ID)
    for(grp in unique(h_groups$group) ){
      all_res$signif = ifelse(all_res$p.adjust < 0.05, 1, 0.2)
      p <- ggplot(all_res[all_res$ID %in% h_groups$hallmark[h_groups$group==grp],])+ geom_point(
        aes(x=Name, y=Description , color=-log10(p.adjust), size=GeneRatio,alpha=signif))+ scale_y_discrete(limits=rev)+scale_x_discrete(drop=F)+
        scale_color_gradient(low = "blue", high = "red", limits = c(min(-log10(all_res$p.adjust)), max(-log10(all_res$p.adjust))))+ ggtitle(grp)
      out[[grp]]=p
    }
    out[["all_res"]]=all_res
    
  }
  return(out)
  
}


##### set up finish ###


mrna_32v0 = read.csv("32vs0_mRNA.csv")
pro_32v0 = read.csv("32vs0_pro.csv")
mrna_32vmid = read.csv("32vsmid_mRNA.csv")
pro_32vmid = read.csv("32vsmid_pro.csv")
mrna_midv0 = read.csv("midvs0_mRNA.csv")
pro_midv0 = read.csv("midvs0_pro.csv")




library(clusterProfiler)



H <- NULL
for ( k in unique(h_gs$gs_name)){
  H = rbind(H, data.frame(name = k, id = pro_32v0$SYMBOL[pro_32v0$SYMBOL %in% h_gs$gene_symbol[h_gs$gs_name == k]]))
}

lt = pro_32v0$fc_pro
names(lt) = pro_32v0$SYMBOL
res.gsea_pro_32v0 = GSEA(sort(lt, decreasing=T), TERM2GENE = H, minGSSize = 1, maxGSSize = 1000,seed=123, pvalueCutoff=1, verbose = F)

H <- NULL
for ( k in unique(h_gs$gs_name)){
  H = rbind(H, data.frame(name = k, id = pro_32vmid$SYMBOL[pro_32vmid$SYMBOL %in% h_gs$gene_symbol[h_gs$gs_name == k]]))
}

lt = pro_32vmid$fc_pro
names(lt) = pro_32vmid$SYMBOL
res.gsea_pro_32vmid = GSEA(sort(lt, decreasing=T), TERM2GENE = H, minGSSize = 1, maxGSSize = 1000,seed=123, pvalueCutoff=1, verbose = F)


H <- NULL
for ( k in unique(h_gs$gs_name)){
  H = rbind(H, data.frame(name = k, id = pro_midv0$SYMBOL[pro_midv0$SYMBOL %in% h_gs$gene_symbol[h_gs$gs_name == k]]))
}

lt = pro_midv0$fc_pro
names(lt) = pro_midv0$SYMBOL
res.gsea_pro_midv0 = GSEA(sort(lt, decreasing=T), TERM2GENE = H, minGSSize = 1, maxGSSize = 1000,seed=123, pvalueCutoff=1, verbose = F)





H <- NULL
for ( k in unique(h_gs$gs_name)){
  H = rbind(H, data.frame(name = k, id = mrna_32v0$SYMBOL[mrna_32v0$SYMBOL %in% h_gs$gene_symbol[h_gs$gs_name == k]]))
}

lt = mrna_32v0$fc_rna
names(lt) = mrna_32v0$SYMBOL
res.gsea_rna_32v0 = GSEA(sort(lt, decreasing=T), TERM2GENE = H, minGSSize = 1, maxGSSize = 1000,seed=123, pvalueCutoff=1, verbose = F)

H <- NULL
for ( k in unique(h_gs$gs_name)){
  H = rbind(H, data.frame(name = k, id = mrna_32vmid$SYMBOL[mrna_32vmid$SYMBOL %in% h_gs$gene_symbol[h_gs$gs_name == k]]))
}

lt = mrna_32vmid$fc_rna
names(lt) = mrna_32vmid$SYMBOL
res.gsea_rna_32vmid = GSEA(sort(lt, decreasing=T), TERM2GENE = H, minGSSize = 1, maxGSSize = 1000,seed=123, pvalueCutoff=1, verbose = F)


H <- NULL
for ( k in unique(h_gs$gs_name)){
  H = rbind(H, data.frame(name = k, id = mrna_midv0$SYMBOL[mrna_midv0$SYMBOL %in% h_gs$gene_symbol[h_gs$gs_name == k]]))
}

lt = mrna_midv0$fc_rna
names(lt) = mrna_midv0$SYMBOL
res.gsea_rna_midv0 = GSEA(sort(lt, decreasing=T), TERM2GENE = H, minGSSize = 1, maxGSSize = 1000,seed=123, pvalueCutoff=1, verbose = F)





gsea.results = list()
gsea.results[["32vs0 \n mRNA"]]=res.gsea_rna_32v0@result[,c("ID", "setSize", "enrichmentScore", "NES", "pvalue", "p.adjust", "qvalue", "rank")]
gsea.results[["32vs0 \n Protein"]]=res.gsea_pro_32v0@result[,c("ID", "setSize", "enrichmentScore", "NES", "pvalue", "p.adjust", "qvalue", "rank")]
gsea.results[["32vsmid \n mRNA"]]=res.gsea_rna_32vmid@result[,c("ID", "setSize", "enrichmentScore", "NES", "pvalue", "p.adjust", "qvalue", "rank")]
gsea.results[["32vsmid \n Protein"]]=res.gsea_pro_32vmid@result[,c("ID", "setSize", "enrichmentScore", "NES", "pvalue", "p.adjust", "qvalue", "rank")]
gsea.results[["midv 0 \n mRNA"]]=res.gsea_rna_midv0@result[,c("ID", "setSize", "enrichmentScore", "NES", "pvalue", "p.adjust", "qvalue", "rank")]
gsea.results[["midvs0 \n Protein"]]=res.gsea_pro_midv0@result[,c("ID", "setSize", "enrichmentScore", "NES", "pvalue", "p.adjust", "qvalue", "rank")]

library(ggplot2)
to_plot = NULL
for ( name in names(gsea.results)){
  gsea.results[[name]]$comparison = name
  to_plot = rbind(to_plot, gsea.results[[name]])
}
to_plot = merge(to_plot, h_groups, by.x="ID", by.y="hallmark")
to_plot$HALLMARK = gsub("_", " ", gsub("HALLMARK_", "", to_plot$ID))
#to_plot$comparison = factor(to_plot$comparison , levels=c("GLOBAL", "PRIMARY", "ARPC"))
to_plot$comparison = factor(to_plot$comparison)

# ggplot(to_plot, aes(x = comparison, y= HALLMARK, size= -log10(p.adjust), color = enrichmentScore)) + 
#   geom_point(data=to_plot[to_plot$p.adjust <=0.05,], alpha=1) +
#   geom_point(data=to_plot[to_plot$p.adjust >0.05,], alpha=0.5) +
#   scale_color_gradient2(low="blue", high="red")+ facet_wrap(vars(group), scales = "free") + theme_classic() + 
#   theme( text =element_text(size=8), axis.text.x =element_text(size=7), axis.text.y = element_text(size=5)  )
ggplot(to_plot, aes(x = comparison, y= HALLMARK, size= -log10(p.adjust), color = enrichmentScore)) + 
  geom_point(data=to_plot[to_plot$p.adjust <=0.05,], alpha=1) +
  geom_point(data=to_plot[to_plot$p.adjust >0.05,], alpha=0.1) +
  scale_color_gradient2(low="blue", high="red")+ facet_wrap(vars(group), scales = "free") + theme_classic() + 
  theme( text =element_text(size=8), axis.text.x =element_text(size=7), axis.text.y = element_text(size=5)  )
# ggplot(to_plot, aes(x = comparison, y= HALLMARK, size= -log10(p.adjust), color = enrichmentScore)) + 
#   geom_point(data=to_plot[to_plot$p.adjust <=0.05,], alpha=1) +
#   scale_color_gradient2(low="blue", high="red")+ facet_wrap(vars(group), scales = "free") + theme_classic() + 
#   theme( text =element_text(size=8), axis.text.x =element_text(size=7), axis.text.y = element_text(size=5)  )
ggsave("GSEA_all_ltl.svg", width = 15, height = 15)








#### new version #####
x.1 = res.gsea_rna_32v0
x.2 = res.gsea_pro_32v0
x.3 = res.gsea_rna_32vmid
x.4 = res.gsea_pro_32vmid
x.5 = res.gsea_rna_midv0
x.6 = res.gsea_pro_midv0

library(stringr)
library(dplyr)
library(DOSE)
library(ggplot2)
library(dplyr)
library(stringr)
library(forcats)
## count the gene number for both results
gene_count.x1 <- x.1@result %>% group_by(ID) %>% summarise(count = sum(str_count(core_enrichment, "/")) + 1)
gene_count.x2 <- x.2@result %>% group_by(ID) %>% summarise(count = sum(str_count(core_enrichment, "/")) + 1)
gene_count.x3 <- x.3@result %>% group_by(ID) %>% summarise(count = sum(str_count(core_enrichment, "/")) + 1)
gene_count.x4 <- x.4@result %>% group_by(ID) %>% summarise(count = sum(str_count(core_enrichment, "/")) + 1)
gene_count.x5 <- x.5@result %>% group_by(ID) %>% summarise(count = sum(str_count(core_enrichment, "/")) + 1)
gene_count.x6 <- x.6@result %>% group_by(ID) %>% summarise(count = sum(str_count(core_enrichment, "/")) + 1)

## merge with the original dataframes
dot_df.x1<- left_join(x.1@result, gene_count.x1, by = "ID") %>% mutate(GeneRatio = count/setSize)
dot_df.x2<- left_join(x.2@result, gene_count.x2, by = "ID") %>% mutate(GeneRatio = count/setSize)
dot_df.x3<- left_join(x.3@result, gene_count.x3, by = "ID") %>% mutate(GeneRatio = count/setSize)
dot_df.x4<- left_join(x.4@result, gene_count.x4, by = "ID") %>% mutate(GeneRatio = count/setSize)
dot_df.x5<- left_join(x.5@result, gene_count.x5, by = "ID") %>% mutate(GeneRatio = count/setSize)
dot_df.x6<- left_join(x.6@result, gene_count.x6, by = "ID") %>% mutate(GeneRatio = count/setSize)

## merge the two results
library(clusterProfiler)
merged.res <- as.data.frame(merge_result(list(w32_vs_w0_mRNA=dot_df.x1, w32_vs_w0_Protein=dot_df.x2, 
                                              w32_vs_wmid_mRNA=dot_df.x3, w32_vs_wid_Protein=dot_df.x4, 
                                              wmid_vs_w0_mRNA=dot_df.x5, wmid_vs_w0_Protein=dot_df.x6
                                              )))

## merged.res <- rbind(dot_df.x1, dot_df.x2) #This merging works but it does **not** include source of results (i.e. 'fgsea' or 'dose')

## Set up/downregulation
merged.res$type = "upregulated"
merged.res$type[merged.res$NES < 0] = "downregulated"

merged.res$SIG = -log10(merged.res$p.adj)
# p <- ggplot(merged.res, aes(x = enrichmentScore, y = fct_reorder(Description, SIG))) + 
#   geom_point(aes(size = SIG, color = enrichmentScore)) +
#   theme_bw(base_size = 14) +
#   scale_colour_gradient2(limits=c(-0.7, 0.4), low = "darkblue",high = "#b50404", mid = "white") +
#   ylab(NULL) + facet_grid(.~Cluster) +
#   scale_size(range = c(1,10)) 
# p
# ggsave('GSEA.svg', p, width = 12, height =8)
# readtext("plot.svg")-> temp
# write_clip(temp$text, object_type = "character")


to_show = merged.res 
to_show_all = to_show %>% mutate(alpha =SIG > -log10(0.05) )
to_show_clean = to_show %>% filter(SIG > -log10(0.05))

p <- ggplot(to_show_all, aes(x = NES, y = fct_reorder(Description, SIG))) + 
  geom_point(aes(size = SIG, color = NES,alpha=alpha)) +
  theme_bw(base_size = 14) +
  scale_colour_gradient2(limits=c(-4.5, 4.5), low = "#3953A4", high = "#b50404") +
  ylab(NULL) + facet_grid(.~Cluster) +
  scale_size(range = c(3,12)) +
  ggtitle("GSEA for all")
p
ggsave('GSEA_ltl.svg', p, width =18, height = 30)


p <- ggplot(to_show_clean, aes(x = NES, y = fct_reorder(Description, SIG))) + 
  geom_point(aes(size = SIG, color = NES)) +
  theme_bw(base_size = 14) +
  scale_colour_gradient2(limits=c(-4.5, 4.5), low = "#3953A4", high = "#b50404") +
  ylab(NULL) + facet_grid(.~Cluster) +
  scale_size(range = c(3,12)) +
  ggtitle("GSEA for all (only significant)")
p
ggsave('11-15-GSEA_ltl_sigs.pdf', p, width =18, height = 10)




