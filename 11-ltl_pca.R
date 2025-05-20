library(ggfortify)
df <- t(protome_ltl[17:26])

pt0=mean(167.52,169.87)
pt8=180.34
pt12 = 187.18
pt32 = mean(220.14,220.61)


pca_ltl=data.frame(df,group = factor(c("8w","0w","8w","0w","32w",
                                          "12w","12w","32w","0w","32w"),levels =c("0w","8w","12w","32w")),
                   pseudotime=c(pt8,pt0,pt8,pt0,pt32,pt12,pt12,pt32,pt0,pt32))
pca_res <- prcomp(df)

#cbPalette <- c("red", "#008CCA", "#009C47", "#8D1F61")

cor_pt = cor.test(pca_res$x[, 1],pca_ltl$pseudotime)
cor_pt$p.value
# 2.91076e-07
cor_pt$estimate
# 0.9838614 

g = autoplot(pca_res,data = pca_ltl, colour = 'pseudotime',size=5)+theme(text = element_text( family = 'Arial'))+theme_classic() +
  geom_label_repel( 
    aes(label=group),size = 5,
    max.overlaps = 25,show.legend = FALSE)+
  scale_color_continuous(high = "#132B43", low = "#56B1F7")+
  annotate("text",label="cor: 0.98\np-value: 2.9e-7",x=0.25,y=0.45,size=6)+
  ggtitle("PCA plot for longitudinal model 331 ")+
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        plot.title = element_text(size = 24),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18)) 
g
ggsave("ltl_pca_new.pdf",g,height=10,width=10)

