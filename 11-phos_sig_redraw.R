library(ggplot2)
library(readr)
library(dplyr)

k_un  <- read_csv("unnormed_phos_sigs_new/unnormed_phospho_sigs_18plex_new_1.csv")
all_res = k_un %>%
  filter(p.adjust < 0.05)
g_5_m <- ggplot(all_res) + geom_segment( aes(x=0, xend=enrichmentScore, y=reorder(ID,enrichmentScore), yend=reorder(ID,enrichmentScore)), color="black")+
  geom_point(aes(x= enrichmentScore,y=reorder(ID,enrichmentScore),size=setSize,color=-log10(p.adjust)))+
  theme_classic() +
  ggtitle("Phos signatures enriched along disease trajectory in Set1 (unnormed)")+
  labs(y= "Signature", x = "Enrichment Score") +
  scale_colour_gradient(low="#8D1F61",high="red")

print(g_5_m)
ggsave("11-22-18plex_unnormed.pdf",plot=g_5_m,height=10)

k_un  <- read_csv("unnormed_phos_sigs_new/phospho_sigs_pilot_new_1.csv")
all_res = k_un %>%
  filter(p.adjust < 0.05)
g_5_m <- ggplot(all_res) + geom_segment( aes(x=0, xend=enrichmentScore, y=reorder(ID,enrichmentScore), yend=reorder(ID,enrichmentScore)), color="black")+
  geom_point(aes(x= enrichmentScore,y=reorder(ID,enrichmentScore),size=setSize,color=-log10(p.adjust)))+
  theme_classic() +
  ggtitle("Phos signatures enriched along disease trajectory in Set2 (unnormed)")+
  labs(y= "Signature", x = "Enrichment Score") +
  scale_colour_gradient(low="#8D1F61",high="red")

print(g_5_m)
ggsave("11-22-10plex_unnormed.pdf",plot=g_5_m)




k_un  <- read_csv("normed_phos_sigs_new/normed_phospho_sigs_18plex_new_1.csv")
all_res = k_un %>%
  filter(p.adjust < 0.05)
g_5_m <- ggplot(all_res) + geom_segment( aes(x=0, xend=enrichmentScore, y=reorder(ID,enrichmentScore), yend=reorder(ID,enrichmentScore)), color="black")+
  geom_point(aes(x= enrichmentScore,y=reorder(ID,enrichmentScore),size=setSize,color=-log10(p.adjust)))+
  theme_classic() +
  ggtitle("Phos signatures enriched along disease trajectory in Set1 (normed)")+
  labs(y= "Signature", x = "Enrichment Score") +
  scale_colour_gradient(low="#8D1F61",high="red")

print(g_5_m)
ggsave("11-23-18plex_normed.pdf",plot=g_5_m,height=10)


k_un  <- read_csv("normed_phos_sigs_new/phospho_sigs_pilot_new_1.csv")
all_res = k_un %>%
  filter(p.adjust < 0.05)
g_5_m <- ggplot(all_res) + geom_segment( aes(x=0, xend=enrichmentScore, y=reorder(ID,enrichmentScore), yend=reorder(ID,enrichmentScore)), color="black")+
  geom_point(aes(x= enrichmentScore,y=reorder(ID,enrichmentScore),size=setSize,color=-log10(p.adjust)))+
  theme_classic() +
  ggtitle("Phos signatures enriched along disease trajectory in Set2 (normed)")+
  labs(y= "Signature", x = "Enrichment Score") +
  scale_colour_gradient(low="#8D1F61",high="red")

print(g_5_m)
ggsave("11-23-10plex_normed.pdf",plot=g_5_m)