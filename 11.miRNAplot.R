library(igraph)
library(ggraph)
library(readr)
library(tidyr)
library(dplyr)

up_mirna = read_csv("gene2mir (1).csv")
up_mirna = up_mirna[!duplicated(up_mirna), ]

down_mirna = read_csv("gene2mir (2).csv")
down_mirna = down_mirna[!duplicated(down_mirna), ]

up_enrich = read_csv("up_enrich.csv")
down_enrich = read_csv("down_enrich.csv")

up_list <- up_enrich %>%
  filter(adj_p<0.05)
up_draw <- up_mirna %>%
  filter(ID %in% up_list$ID)


down_list <- down_enrich %>%
  filter(adj_p<0.05)
down_draw <- down_mirna %>%
  filter(ID %in% down_list$ID)



database1 <- read_csv("/Users/jichangzhang/Downloads/hsa/miRNet-mir-gene-hsa-mirecords.csv")
database2 <- read_csv("/Users/jichangzhang/Downloads/hsa/miRNet-mir-gene-hsa-mirtarbase.csv")
database3 <- read_csv("/Users/jichangzhang/Downloads/hsa/miRNet-mir-gene-hsa-tarbase.csv")

datafull <- (rbind(database1, database2, database3))[,c(2,4)]
datafull = datafull[!duplicated(datafull), ]



library(igraph)
library(ggraph)
library(tidyverse)

# 创建一个示例数据框
edges_df <- data.frame(
  source = down_draw$ID,
  target = down_draw$Target
)

# 将数据框转换为 igraph 对象
graph <- graph_from_data_frame(edges_df)
node_df <- data.frame(
  name = V(graph)$name,
  type = ifelse(V(graph)$name %in% edges_df$source, "microRNA", "gene"),
  degree = degree(graph)
)
# 绘制网络图
q=ggraph(graph, layout = "fr") +
  geom_edge_link(alpha=0.3) +
  geom_node_point(aes(color =node_df$type)) +
  scale_color_manual(values = c("red", "blue")) 
#ggsave("test_down_network.png",q)







edges_df <- data.frame(
  source = up_draw$ID,
  target = up_draw$Target
)

# 将数据框转换为 igraph 对象
graph <- graph_from_data_frame(edges_df)

node_df <- data.frame(
  name = V(graph)$name,
  type = ifelse(V(graph)$name %in% edges_df$source, "microRNA", "gene"),
  degree = degree(graph)
) 


# 绘制二分图
q = ggraph(graph) +
  geom_edge_link(alpha=0.05) +
  geom_node_point(aes(color =node_df$type,size=node_df$degree)) +
  scale_color_manual(values = c("red", "blue")) 

#ggsave("test_up_network.png",q)


mi_list = down_enrich$ID[1:6]
FAT4 = datafull %>% 
  filter(symbol == "FAT4") %>%
  filter(mir_id %in% mi_list)

NF1 = datafull %>% 
  filter(symbol == "NF1") %>%
  filter(mir_id %in% mi_list)


FAT4_list = datafull %>% 
  filter(symbol == "FAT4")

NF1_list = datafull %>% 
  filter(symbol == "NF1") 
#write.csv(FAT4_list,"FAT4_list.csv",row.names = FALSE)

#write.csv(NF1_list,"NF1_list.csv",row.names = FALSE)

up_gene_degree <- node_df %>%
  filter(type == "gene")

#write.csv(up_gene_degree,"up_gene_degree.csv",row.names=FALSE)

up_gene_list = c("JUN","KMT1A","HOXC6","COT1","CENPN","COP1","KDM1A","MVP","PHF3","POLA2","POLR2J","SMARCC2","YWHAE")


up_gene_interaction = up_mirna %>%
  filter(Target %in% up_gene_list) %>%
  arrange(Target) %>%
  dplyr::select(c(ID,Target))
#write.csv(up_gene_interaction,"up_gene_interaction.csv",row.names=FALSE)


mirna_list = c("hsa-mir-1-3p","hsa-mir-16-5p","hsa-mir-124-3p","hsa-mir-155-5p","hsa-let-7b-5p","hsa-mir-107")
short_gene_mirna  <- up_gene_interaction %>%
  filter(ID %in% mirna_list)


edges_df <- data.frame(
  source = short_gene_mirna$ID,
  target = short_gene_mirna$Target
)

# 将数据框转换为 igraph 对象
graph <- graph_from_data_frame(edges_df)

idx <- which(names(V(graph)) %in% short_gene_mirna$ID)


V(graph)$shape = "circle"
V(graph)$color = "red"
V(graph)[idx]$shape <- "square"
V(graph)$size = 5
V(graph)[idx]$color = "lightblue"
V(graph)$label = names(V(graph))
q = ggraph(graph,layout = "fr")+geom_edge_link() +
  geom_node_point(aes(shape = shape, size = 3, color =color), stroke = 2)+
  geom_node_label(aes(label = label), repel = TRUE, arrow = arrow(length = unit(0.15, "inches"))) + theme_graph()+
  theme(legend.position = "none", plot.margin = margin(1,1,1,1))


#ggsave("short_list_network.png",q)

q = ggraph(graph)+geom_edge_link() +
  geom_node_point(aes(shape = shape, size = 3, color =color), stroke = 2)+
  geom_node_label(aes(label = label), repel = TRUE, arrow = arrow(length = unit(0.15, "inches"))) + theme_graph()+
  theme(legend.position = "none", plot.margin = margin(1,1,1,1))


#ggsave("short_list_bipart_network.png",q)


MyLO = matrix(0, nrow=vcount(graph), ncol=2)

layer =  V(graph)$color == "red"
## Horizontal position is determined by layer
MyLO[,1] = V(graph)$color == "red"
V(graph)$degree = degree(graph)

## Vertical position is determined by sum of sorted vertex_degree
for(i in 1:2) {
  L  = which(layer == i-1)
  OL = order(V(graph)$degree[L], decreasing=TRUE)
  MyLO[L[OL],2] = cumsum(V(graph)$degree[L][OL])
}



ggraph(graph,layout=MyLO)+geom_edge_link() +
  geom_node_point(aes(shape = shape, size = 3, color =color), stroke = 2)+
  geom_node_label(aes(label = label), repel = TRUE, arrow = arrow(length = unit(0.15, "inches"))) + theme_graph()+
  theme(legend.position = "none", plot.margin = margin(1,1,1,1))
V(graph)$label = gsub("hsa-", "", V(graph)$label)
V(graph)$label[16] = "1433E"
MyLO2 = cbind(MyLO[,2],MyLO[,1])
q = ggraph(graph,layout=MyLO2)+geom_edge_link() +
  geom_node_point(aes(shape = shape, size = 3, color =color), stroke = 2)+
  geom_node_label(aes(label = label), nudge_y=0.1,size=3.5) + theme_graph()+
  theme(legend.position = "none", plot.margin = margin(1,1,1,1))
#V(graph)$label

q
ggsave("11-16-short_list_bipart_network.pdf",height = 3.5,width=12,q)
