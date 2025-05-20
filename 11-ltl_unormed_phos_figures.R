library(dplyr)
library(readr)
library(ggplot2)

dir.create("ltl_unnormed_phos/", recursive = TRUE, showWarnings = FALSE)



mid_vs_0 = read.csv("./ltl_unnormed_phos/mid_vs_0/phospho_sigs_18plex_new.csv")

all_res = mid_vs_0 %>% filter(p.adjust<0.05)
g_8_m <- ggplot(all_res) + geom_segment( aes(x=0, xend=enrichmentScore, y=reorder(ID,enrichmentScore), yend=reorder(ID,enrichmentScore)), color="black")+
  geom_point(aes(x= enrichmentScore,y=reorder(ID,enrichmentScore),size=setSize,color=-log10(pvalue)))+
  theme_classic() +
  ggtitle("Unnormed Phosphosignatures enriched between mid weeks and 0 weeks")+
  labs(y= "Signature", x = "Enrichment Score") +
  scale_colour_gradient(low="#8D1F61",high="red")

print(g_8_m)
ggsave("./ltl_unnormed_phos/mid_vs_0.pdf",g_8_m)


mid_vs_0 = read.csv("./ltl_unnormed_phos/32_vs_0/phospho_sigs_18plex_new.csv")

all_res = mid_vs_0 %>% filter(p.adjust<0.05)
g_8_m <- ggplot(all_res) + geom_segment( aes(x=0, xend=enrichmentScore, y=reorder(ID,enrichmentScore), yend=reorder(ID,enrichmentScore)), color="black")+
  geom_point(aes(x= enrichmentScore,y=reorder(ID,enrichmentScore),size=setSize,color=-log10(pvalue)))+
  theme_classic() +
  ggtitle("Unnormed Phosphosignatures enriched between 32 weeks and 0 weeks")+
  labs(y= "Signature", x = "Enrichment Score") +
  scale_colour_gradient(low="#8D1F61",high="red")

print(g_8_m)
ggsave("./ltl_unnormed_phos/32_vs_0.pdf",g_8_m)



mid_vs_0 = read.csv("./ltl_unnormed_phos/32_vs_mid/phospho_sigs_18plex_new.csv")

all_res = mid_vs_0 %>% filter(p.adjust<0.05)
g_8_m <- ggplot(all_res) + geom_segment( aes(x=0, xend=enrichmentScore, y=reorder(ID,enrichmentScore), yend=reorder(ID,enrichmentScore)), color="black")+
  geom_point(aes(x= enrichmentScore,y=reorder(ID,enrichmentScore),size=setSize,color=-log10(pvalue)))+
  theme_classic() +
  ggtitle("Unnormed Phosphosignatures enriched between 32 weeks and mid weeks")+
  labs(y= "Signature", x = "Enrichment Score") +
  scale_colour_gradient(low="#8D1F61",high="red")

print(g_8_m)
ggsave("./ltl_unnormed_phos/32_vs_mid.pdf",g_8_m,height=12)
