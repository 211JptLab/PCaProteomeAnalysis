library(igraph)
library(ggraph)
g <- read.graph("/Users/jichangzhang/Downloads/mirnet1.graphml", format = "graphml")
plot(g)
#V(g)
label_list = c("NF1","FAT4")
idx <- which(names(V(g)) %in% label_list)
idx
V(g)$shape = "circle"
V(g)$color = "red"
V(g)[idx]$shape <- "square"


V(g)$size = 1

deg <- degree(g)
V(g)$label <- ifelse(deg > 1, gsub("hsa-", "", V(g)$name), NA)
V(g)$color <- ifelse(deg > 1, "#009C47", "red")
V(g)$size[c(28,136,137)] = 1.00005
V(g)[idx]$color = "lightblue"
layout <- layout_with_fr(g)


plot(g,layout=layout, vertex.label.dist=1.5,vertex.label.cex=0.3)

graph <- ggraph(g)

graph <- graph + 
  geom_edge_link() +
  geom_node_point(aes(shape = shape, size = size, color =color), stroke = 2)+
  geom_node_label(aes(label = label), size=3.5,repel = TRUE)

graph <- graph + theme_graph()+
  theme(legend.position = "none", plot.margin = margin(1,1,1,1))
   #scale_layout_gem(gamma = 1.2) +
   #geom_node_text_repel(segment.alpha = 0.3, 
  #box.padding = unit(0.4, "lines"),
#force = 2)
graph
ggsave("11-17-2gene_intersect.pdf",graph,height=5,width=7.5)
