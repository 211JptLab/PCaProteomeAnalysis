library(dplyr)
library(ggplot2)
library(broom)

# get the first dataframe

set1 = data.frame(pt=pt,PC1 = pca_res$x[,1],set="set1",group = pca_ltl$group )
# get the second dataframe
set2 = data.frame(pt=pt,PC1 = pca_res$x[,1],set="set2",group = pca_ltl$group )
#get the third dataframe
set3 = data.frame(pt=pt,PC1 = pca_res$x[,1],set="LTL331",group = pca_ltl$group )

set = rbind(set1,set2,set3)

model_df <- set %>%
  group_by(set) %>%
  summarise(slope = coef(lm(PC1 ~ pt, data = .))[["pt"]],
            intercept = coef(lm(PC1 ~ pt, data = .))[["(Intercept)"]],
            x = mean(pt),  # 用每组的平均pt值作为文本位置的x坐标
            y = intercept + slope * mean(pt))  
labels <- c("cor: 0.87\np-value: 0.001","cor: 0.81\np-value: 5.5e-9","cor: 0.53\np-value: 0.02")
model_df$labels <- labels

model_df$x = c(180,200,200)
model_df$y = c(50,80,-80)


g = ggplot(set, aes(x = pt, y = PC1, color = set)) + 
  geom_point() +                               
  geom_smooth(method = "lm", se = TRUE,aes(fill=set),show.legend = FALSE) +     
  labs(title = "Correlation between pseudotime and PC1 for \n different sets", x = "pseudotime", y = "PC1")+
  geom_text(data = model_df, aes( x=x,y=y,label = labels), 
            vjust = -1, hjust = 1)+
  theme_minimal() +  
  theme(
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    panel.background = element_blank(),  
    axis.line = element_line(colour = "black") 
  )+
  labs(x=NULL, y=NULL,color="Set",fill=NULL) +
  
  scale_color_manual(values=c("#7abb72","#dd6d5c","#a98cba"),
                    labels = c("LTL331","Set1", "Set2")) +
  scale_fill_manual(values=c("#7abb72","#dd6d5c","#a98cba"),
                     labels = c("LTL331","Set1", "Set2"))



g


ggsave("./Figures_Nov_17/11-21-cor_3_sets.pdf",g,height=6,width=6)
