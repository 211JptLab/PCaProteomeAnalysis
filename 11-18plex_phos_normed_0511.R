library(readr)
library(dplyr)
library(ggplot2)
library(patchwork)
plex10 <- read_csv("normed_phos_sigs_new//phospho_sigs_pilot_new_1.csv")
plex18 <- read_csv("normed_phos_sigs_new//normed_phospho_sigs_18plex_new_1.csv")
plex9 <- read_csv("normed_phos_sigs_new//phospho_sigs_pilot_pca1.csv")
`Plex 10` <- plex10 %>% filter(p.adjust < 0.05) %>% select(ID)
`Plex 18` <- plex18 %>% filter(p.adjust < 0.05) %>% select(ID)
`Plex 9` <- plex9 %>% filter(p.adjust < 0.05) %>% select(ID)
common_elements <- intersect(intersect(`Plex 10`, `Plex 18`), `Plex 9`)

st = "PERT-P100-DIA2_RESVERATROL"
# p10 <- plex10 %>%
#   filter(ID %in% common_elements$ID) %>%
#   mutate(per_resv_plex10 = n_inter / plex10[plex10$ID==st,]$setSize) %>%
#   select(ID,per_resv_plex10)
# p9 <- plex9 %>%
#   filter(ID %in% common_elements$ID) %>%
#   mutate(per_resv_plex9 = n_inter / plex9[plex9$ID==st,]$setSize) %>%
#   select(ID,per_resv_plex9)
p18 <- plex18 %>%
  filter(ID %in% common_elements$ID) %>%
  mutate(per_resv_plex18 = n_inter / plex18[plex18$ID==st,]$setSize) 

#p = merge(merge(p10,p9,by="ID"),p18,by="ID")
group_info = read_csv("shared normalized.csv")

p18_aug = merge(p18,group_info,by="ID",all.x=TRUE)
#p18_aug$


p18_aug$group = as.factor(p18_aug$...6)



p18_sel = p18_aug %>% filter(group != "Not mention")
p18_sel = p18_sel %>%
  group_by(group) %>%
  mutate(n=n())
p18_sel <- p18_sel %>%  arrange( n, group,enrichmentScore)     
p18_sel$y=1:dim(p18_sel)[1]
p18_sel$ID <- factor(p18_sel$ID, levels = p18_sel$ID)
name_list = p18_sel %>%
  group_by(group) %>%
  summarize(n = n()) %>%
  arrange(n) %>%
  mutate(sum = cumsum(n))

g1 = ggplot(p18_sel,aes(x= enrichmentScore,y=y)) +   
  geom_rect(aes(xmin=-Inf,xmax=Inf,ymin=0.5,ymax=1.5,fill=(as.character(name_list$group[1]))),alpha=0.01) +
  geom_rect(aes(xmin=-Inf,xmax=Inf,ymin=1.5,ymax=2.5,fill=(as.character(name_list$group[2]))),alpha=0.01) +
  geom_rect(aes(xmin=-Inf,xmax=Inf,ymin=2.5,ymax=4.5,fill=(as.character(name_list$group[3]))),alpha=0.01) +
  geom_rect(aes(xmin=-Inf,xmax=Inf,ymin=4.5,ymax=6.5,fill=(as.character(name_list$group[4]))),alpha=0.01) +
  geom_rect(aes(xmin=-Inf,xmax=Inf,ymin=6.5,ymax=10.5,fill=(as.character(name_list$group[5]))),alpha=0.01) +
  geom_rect(aes(xmin=-Inf,xmax=Inf,ymin=10.5,ymax=14.5,fill=(as.character(name_list$group[6]))),alpha=0.01) +
  #geom_rect(aes(xmin=-Inf,xmax=Inf,ymin=26.5,ymax=36.5,fill=(as.character(name_list$group[7]))),alpha=0.01) +
  scale_fill_brewer(palette = 'Dark2', name = 'Group',
                    breaks=c(as.character(name_list$group[1]),
                             as.character(name_list$group[2]),
                             as.character(name_list$group[3]),
                             as.character(name_list$group[4]),
                             as.character(name_list$group[5]),
                      #       as.character(name_list$group[6]),
                            as.character(name_list$group[6])))+
  geom_segment( aes(x=0, xend=enrichmentScore, y=y, yend=y), color="black")+
  geom_point(aes(x= enrichmentScore,y=y,size=setSize,color=-log10(p.adjust)))+
  theme_classic() +
  ggtitle("Normed phosphosignatures enriched in Set1")+
  labs(y= "Signature", x = "Enrichment Score") +
  scale_colour_gradient(low="#8D1F61",high="red")+
  scale_y_continuous(breaks=p18_sel$y, labels=p18_sel$ID,expand = c(0, 0),position="right"
  )+
  guides(fill = guide_legend(reverse= TRUE, override.aes = list(alpha = 0.1)))

g1
#ggsave("plex18_sigs_normed_3.0.png")

View(p18_sel)
p18_sel$percent_overlap = p18_sel$n_inter/p18_sel$setSize
g2 = 
  ggplot(p18_sel) +
  geom_tile(aes( x=0.1,y =ID, fill = percent_overlap),color="black") +
  scale_fill_gradient(low = "white", high = "red")+
  theme_classic() +
  theme(axis.text.x = element_blank(),       # Hide axis text
        axis.ticks.x = element_blank(),      # Hide axis ticks
        axis.line.x = element_blank(),       # Hide axis lines
        axis.title.x = element_blank(),
        legend.position = "left")      +
  labs(y="Target")+
  
  scale_y_discrete(labels = p18_sel$...4)+
  coord_fixed(ratio = 0.1, xlim = c(0.5, 0.6))  # Adjust the ratio and limits
g2+g1
ggsave("11-24-plex18_sigs_normed.pdf",width = 12, height = 8)


