#first get p18_sel in 18plex_phos_normed_0511.R
ID = p18_sel$ID

# from add_heatbar_0512 get g1
g1 + scale_x_discrete(position = "top") 
#from cmap_inves_new get pdx, pdx_up, pdx_down
p18_sel$overlap_with_normed = p18_sel$ID %in% 
  ID

sum( p18_sel$overlap_with_normed)
# find 2 sigs only in normed data, not in unnormed data

## import pdx and pdx_down in cmap_inves)_new.R


library(stringr)
library(ggplot2)
p18_sel$query <- sapply(str_split(p18_sel$ID, "_"), function(x) ifelse(length(x) > 1, x[2], NA))
p18_sel$profiled <- tolower(p18_sel$query) %in% tolower(pdx$pert_iname)
sum(p18_sel$profiled )
p18_sel$sign_cmap <- tolower(p18_sel$query) %in% tolower(pdx_down$pert_iname)
sum(p18_sel$sign_cmap )

#write.csv(p18_sel,"p18_sel.csv",row.names = FALSE)









p10 <- plex10 %>%
  filter(ID %in% common_elements$ID) %>%
  mutate(per_resv_plex18 = n_inter / plex10[plex10$ID==st,]$setSize) 

p10_aug = merge(p10,group_info,by="ID")
p10_aug$group = as.factor(p10_aug$...6)

p10_sel = p10_aug %>% filter(group != "Not mention")
p10_sel = p10_sel %>%
  group_by(group) %>%
  mutate(n=n())
p10_sel <- p10_sel %>%  arrange( n, group,enrichmentScore)     
p10_sel$y=1:dim(p10_sel)[1]
p10_sel$ID <- factor(p10_sel$ID, levels = p10_sel$ID)



library(tidyr)
p18_sel$percent_overlap = p18_sel$n_inter/p18_sel$setSize

df_long <- p18_sel[,c(1,15,21,23,24)] %>%
  pivot_longer(cols = -c(ID, `...4`), names_to = "variable", values_to = "value")

##here is the start point for re-analysis, just load the RData
# 创建热图

library(ggnewscale)


g1 = ggplot(p18_sel,aes(x= enrichmentScore,y=y)) +   
  geom_rect(aes(xmin=-Inf,xmax=Inf,ymin=0.5,ymax=1.5,fill=(as.character(name_list$group[1]))),alpha=0.01) +
  geom_rect(aes(xmin=-Inf,xmax=Inf,ymin=1.5,ymax=4.5,fill=(as.character(name_list$group[2]))),alpha=0.01) +
  geom_rect(aes(xmin=-Inf,xmax=Inf,ymin=4.5,ymax=8.5,fill=(as.character(name_list$group[3]))),alpha=0.01) +
  geom_rect(aes(xmin=-Inf,xmax=Inf,ymin=8.5,ymax=12.5,fill=(as.character(name_list$group[4]))),alpha=0.01) +
  geom_rect(aes(xmin=-Inf,xmax=Inf,ymin=12.5,ymax=19.5,fill=(as.character(name_list$group[5]))),alpha=0.01) +
  geom_rect(aes(xmin=-Inf,xmax=Inf,ymin=19.5,ymax=26.5,fill=(as.character(name_list$group[6]))),alpha=0.01) +
  geom_rect(aes(xmin=-Inf,xmax=Inf,ymin=26.5,ymax=36.5,fill=(as.character(name_list$group[7]))),alpha=0.01) +
  scale_fill_brewer(palette = 'Dark2', name = 'Group',
                    breaks=c(as.character(name_list$group[1]),
                             as.character(name_list$group[2]),
                             as.character(name_list$group[3]),
                             as.character(name_list$group[4]),
                             as.character(name_list$group[5]),
                             as.character(name_list$group[6]),
                             as.character(name_list$group[7])))+
  
  geom_point(aes(x= enrichmentScore,y=y,size=setSize,color=-log10(p.adjust)))+
  scale_colour_gradient(low="#8D1F61",high="red")+
  labs(color ="Significance for Set1")+
  new_scale_color() +
  geom_point(data=p10_sel,aes(x= enrichmentScore,y=y,size=setSize,color=-log10(p.adjust)))+
  labs(color ="Significance for Set2")+
  theme_classic() +
  ggtitle("Phosphosignatures enriched in Set1 & Set2")+
  labs(y= "Signature", x = "Enrichment Score") +
  
  scale_y_continuous(breaks=p18_sel$y, labels=p18_sel$ID,expand = c(0, 0),position="right"
  )+
  guides(fill = guide_legend(reverse= TRUE, override.aes = list(alpha = 0.1)))
p18_sel$geneExpression = p18_sel$profiled + p18_sel$sign_cmap
p18_sel$geneExpression = as.character((p18_sel$geneExpression))
g2 = 
  ggplot(p18_sel) +
  geom_tile(aes( x=0.1,y =ID, fill = geneExpression),color="black") +
  #scale_fill_gradient(low = "white", high = "red")+
  scale_fill_manual(values = c("0" = "transparent", "1" = "darkgrey","2"="darkgrey"))+
  theme_classic() +
  theme(axis.text.x = element_blank(),       # Hide axis text
        axis.ticks.x = element_blank(),      # Hide axis ticks
        axis.line.x = element_blank(),       # Hide axis lines
        
        legend.position = "left")      +
  labs(y="Target",x="Gene\nset",fill="In CMap")+ 
  scale_fill_manual(values = c("white","lightgrey", "darkgrey"),
                                                                 labels = c("No","Yes", "Yes and Significant"))+
  
  scale_y_discrete(labels = p18_sel$...4)+
  coord_fixed(ratio = 0.1, xlim = c(0.5, 0.6))  # Adjust the ratio and limits



g4 = 
  ggplot(p18_sel) +
  geom_tile(aes( x=0.1,y =ID, fill = overlap_with_normed),color="black") +
  #scale_fill_gradient(low = "white", high = "red")+
  scale_fill_manual(values = c("FALSE" = "transparent", "TRUE" = "lightgreen"))+
  
  theme_classic() +
  theme(axis.text.x = element_blank(),       # Hide axis text
        axis.ticks.x = element_blank(),      # Hide axis ticks
        axis.line.x = element_blank(),       # Hide axis lines
        axis.text.y = element_blank(),       # Hide axis text
        axis.ticks.y = element_blank(),      # Hide axis ticks
        axis.line.y = element_blank(), 
        axis.title.y = element_blank(), 
        legend.position = "none")      +
  labs(x="In\nNormed")+
  
  scale_y_discrete(labels = p18_sel$...4)+
  coord_fixed(ratio = 0.1, xlim = c(0.5, 0.6))  # Adjust the ratio and limits
g4



combined_plot <- g2 + scale_x_discrete(position = "top") +
  g4 + scale_x_discrete(position = "top") +
  g1 + scale_x_discrete(position = "top") 
combined_plot
combined_plot <- g2 + 
  g4 +
  g1 
combined_plot
ggsave("sig_updated_8.0_to_work_final.pdf",combined_plot, width = 12)

g4 = 
  ggplot(p18_sel) +
  geom_tile(aes( x=0.1,y =ID, fill = overlap_with_normed),color="black") +
  #scale_fill_gradient(low = "white", high = "red")+
  scale_fill_manual(values = c("FALSE" = "transparent", "TRUE" = "lightgreen"))+
  
  theme_classic() +
  theme(axis.text.x = element_blank(),       # Hide axis text
        axis.ticks.x = element_blank(),      # Hide axis ticks
        axis.line.x = element_blank(),       # Hide axis lines
        axis.text.y = element_blank(),       # Hide axis text
        axis.ticks.y = element_blank(),      # Hide axis ticks
        
        axis.title.y = element_blank())      +
  labs(x="In\nNormed")+
  
  scale_y_discrete(labels = p18_sel$...4)+
  coord_fixed(ratio = 0.1, xlim = c(0.5, 0.6))  # Adjust the ratio and limits
g4
ggsave("sig_updated_8.0_to_work_legend.pdf",g4, width = 12)



##### the next will be to incorporate the ltl models
mid_vs_0 = read.csv("./Figures_Nov_17/ltl_unnormed_phos/mid_vs_0/phospho_sigs_18plex_new.csv")
u_midv0 = mid_vs_0 %>% filter(p.adjust<0.05)
mid_vs_0 = read.csv("./Figures_Nov_17/ltl_unnormed_phos/32_vs_0/phospho_sigs_18plex_new.csv")
u_32v0 = mid_vs_0 %>% filter(p.adjust<0.05)
mid_vs_0 = read.csv("./Figures_Nov_17/ltl_unnormed_phos/32_vs_mid/phospho_sigs_18plex_new.csv")
u_32vmid = mid_vs_0 %>% filter(p.adjust<0.05)

mid_vs_0 = read.csv("./mid_vs_0/phospho_sigs_18plex_new.csv")
n_midv0 = mid_vs_0 %>% filter(p.adjust<0.05)
mid_vs_0 = read.csv("./32_vs_0/phospho_sigs_18plex_new_32_vs_0.csv")
n_32v0 = mid_vs_0 %>% filter(p.adjust<0.05)
mid_vs_0 = read.csv("./32_vs_mid/phospho_sigs_18plex_new_32_vs_mid.csv")
n_32vmid = mid_vs_0 %>% filter(p.adjust<0.05)


custom_order <- c("wmid_vs_w0", "w32_vs_w0", "w32_vs_wmid")



dat_un <- data.frame(id=1:36,
                  wmid_vs_w0 = p18_sel$ID %in% u_midv0$ID,
                  w32_vs_w0 = p18_sel$ID %in% u_32v0$ID,
                  w32_vs_wmid = p18_sel$ID %in% u_32vmid$ID) %>%
  gather(key = "Variable", value = "Value",wmid_vs_w0:w32_vs_wmid)
dat_un$Variable <- factor(dat_un$Variable, levels = custom_order)

dat_n <- data.frame(id=1:36,
                    wmid_vs_w0 = p18_sel$ID %in% n_midv0$ID,
                    w32_vs_w0 = p18_sel$ID %in% n_32v0$ID,
                    w32_vs_wmid = p18_sel$ID %in% n_32vmid$ID) %>%
  gather(key = "Variable", value = "Value",wmid_vs_w0:w32_vs_wmid)
dat_n$Variable <- factor(dat_n$Variable, levels = custom_order)



g5 = 
  ggplot(dat_un) +
  geom_tile(aes( x=Variable,y =id, fill = Value),color="black")+
  scale_fill_manual(values = c("FALSE" = "white", "TRUE" = "black"))+
  theme_classic() +
  theme(#axis.text.x = element_blank(),       # Hide axis text
        axis.ticks.x = element_blank(),      # Hide axis ticks
        axis.line.x = element_blank(),       # Hide axis lines
        axis.text.y = element_blank(),       # Hide axis text
        axis.ticks.y = element_blank(),      # Hide axis ticks
        axis.title.x = element_blank(),
        axis.title.y = element_blank())      +
  scale_y_discrete(labels = p18_sel$...4)+
  coord_fixed(ratio = 1) +
 theme(axis.text.x = element_text(angle = 65, hjust = 1))


g5
ggsave("./Figures_Nov_17/unnormed_overlap.pdf",g5, width = 12)

g6 = 
  ggplot(dat_n) +
  geom_tile(aes( x=Variable,y =id, fill = Value),color="black")+
  scale_fill_manual(values = c("FALSE" = "white", "TRUE" = "black"))+
  theme_classic() +
  theme(#axis.text.x = element_blank(),       # Hide axis text
    axis.ticks.x = element_blank(),      # Hide axis ticks
    axis.line.x = element_blank(),       # Hide axis lines
    axis.text.y = element_blank(),       # Hide axis text
    axis.ticks.y = element_blank(),      # Hide axis ticks
    axis.title.x = element_blank(),
    axis.title.y = element_blank())      +
  scale_y_discrete(labels = p18_sel$...4)+coord_fixed(ratio = 1) +
  theme(axis.text.x = element_text(angle = 65, hjust = 1))

g6


ggsave("./Figures_Nov_17/normed_overlap.pdf",g6, width = 12)






ggplot(dat_un) +
  geom_tile(aes( x=Variable,y =id, fill = Value),color="black")
+
 # scale_fill_manual(values = c("FALSE" = "white", "TRUE" = "black"))+
  
  geom_polygon(data = dat_n, aes(x = x1, y = y1, group = z),)





df <- tibble::tibble(x = c(LETTERS[1:6], LETTERS[1:5]),
                     y = c(paste0("V", 1:6), paste0("V", 1:5)),
                     group = c(rep("group_1", 6), rep("group_2", 5)))

df1    <- df[!duplicated(interaction(df$x, df$y)),]
df2    <- df[duplicated(interaction(df$x, df$y)),]
df2    <- df[rep(seq(nrow(df)), each = 3),]
df2$x1 <- as.numeric(as.factor(df2$x))
df2$y1 <- as.numeric(as.factor(df2$y))
df2$x1 <- df2$x1 + c(-0.5, 0.5, 0.5)
df2$y1 <- df2$y1 + c(-0.5, -0.5, 0.5)
df2$z  <- rep(seq(nrow(df2)/3), each = 3)

ggplot(df1, aes(x = x, y = y, fill = group)) + 
  geom_tile() +
  geom_polygon(data = df2, aes(x = x1, y = y1, group = z))


dat_un$group = "Unnormed"
dat_n$group = "Normed"
df2    <- dat_n[rep(seq(nrow(dat_n)), each = 3),]
df2$x1 <- as.numeric(as.factor(df2$Variable))
df2$y1 <- as.numeric(as.factor(df2$id))
df2$x1 <- df2$x1 + c(-0.5, 0.5, 0.5)
df2$y1 <- df2$y1 + c(-0.5, -0.5, 0.5)
df2$z  <- rep(seq(nrow(df2)/3), each = 3)

g7 = ggplot(dat_un) +
  geom_tile(aes( x=Variable,y =id, fill = Value),color="black")+
 
  geom_polygon(data = df2, aes(x = x1, y = y1, fill=Value,group = z))+
  scale_fill_manual(values = c("FALSE" = "white", "TRUE" = "black"))+
  theme_classic() +
  theme(#axis.text.x = element_blank(),       # Hide axis text
    axis.ticks.x = element_blank(),      # Hide axis ticks
    axis.line.x = element_blank(),       # Hide axis lines
    axis.text.y = element_blank(),       # Hide axis text
    axis.ticks.y = element_blank(),      # Hide axis ticks
    axis.title.x = element_blank(),
    axis.title.y = element_blank())      +
  scale_y_discrete(labels = p18_sel$...4)+coord_fixed(ratio = 1) +
  theme(axis.text.x = element_text(angle = 65, hjust = 1))
g7
ggsave("./Figures_Nov_17/un_and_normed_overlap.pdf",g7, width = 12)

