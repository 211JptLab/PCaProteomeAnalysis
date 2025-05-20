#### save to_draw_22
# save(to_draw_22, file = "1_before_purity_selection.Rdata")
library(ggplot2)
library(ggsci)
library(ggpubr)
library(dplyr)



count = data.frame(t(sort(table(to_draw_22$cor_rna),decreasing = TRUE)))

to_draw_23 = to_draw_22 %>%
  filter(!(cor_rna %in% count$Var2[1:22]))

g = ggplot(to_draw_23, aes(x= cor_rna,y = coef)) +
  geom_point(alpha = 0.1,size=0.5) +
  stat_cor(method="pearson",label.x = -0.4, label.y=0.9,size=6.5)+
  xlim(-1,1)+
  ylim(-1,1)+
  xlab("PDX model")+
  ylab("Unselected patient tumors")+
  theme_classic() +
  theme(
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 22),
    legend.title = element_text(size = 24),
    legend.text = element_text(size = 22),
    plot.title = element_text(size = 24)
  ) +
  coord_fixed()+
  scale_color_npg()+
  #theme(text = element_text( family = 'Arial'))+
  ggtitle("Correlation Gene Expression \n Patient - PDX (PT 150-220)")
g
ggsave("bf_selected.pdf",g,
       width = 10, # PDF宽度（英寸）
       height = 7, # PDF高度（英寸）
       dpi = 300
) # 分辨率

