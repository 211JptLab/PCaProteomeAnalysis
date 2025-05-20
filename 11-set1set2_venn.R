library(readr)
library(dplyr)
library(patchwork)
library(cowplot)

plex10 <- read_csv("unnormed_phos_sigs_new//phospho_sigs_pilot_new_1.csv")
plex18 <- read_csv("unnormed_phos_sigs_new//unnormed_phospho_sigs_18plex_new_1.csv")

`Plex 10 Unnormed` <- plex10 %>% filter(p.adjust < 0.05) %>% select(ID)
`Plex 18 Unnormed` <- plex18 %>% filter(p.adjust < 0.05) %>% select(ID)


plex10 <- read_csv("normed_phos_sigs_new//phospho_sigs_pilot_new_1.csv")
plex18 <- read_csv("normed_phos_sigs_new//normed_phospho_sigs_18plex_new_1.csv")

`Plex 10 Normed` <- plex10 %>% filter(p.adjust < 0.05) %>% select(ID)
`Plex 18 Normed` <- plex18 %>% filter(p.adjust < 0.05) %>% select(ID)

length(intersect(`Plex 10 Unnormed`,`Plex 10 Normed`))

library(eulerr)

dat<-c("Set2 Unnormed" = 61-21, 
       "Set2 Normed" = 22-21,
       "Set2 Unnormed&Set2 Normed" = 21)
#png("Set2_intersect.png")
v2 = plot(euler(dat),
          fills = list(fill = c("white", "darkgrey"), alpha = 0.9),
          quantities = list(type = c("counts","percent"),cex=1.5),
          #labels = list(font=3,side="above"),
          legend = list(side = "right",cex=1.5)
)
tags <- v2$children[[2]]$children[[1]]$children$tags$children
tags <- do.call(grid::gList, lapply(tags, function(x) {
  x$children[[2]]$label <- sub(" \\(", "\n(", x$children[[2]]$label)
  x$children[[2]]$just <- NULL
  x$children[[2]]$hjust <- 0.5
  x$children[[2]]$vjust <- 1
  x}))

v2$children[[2]]$children[[1]]$children$tags$children <- tags

ggsave("11-25-set2.pdf",v2)

#dev.off()


intersect(`Plex 18 Normed`,`Plex 18 Unnormed`)
dat<-c("Set1 Unnormed" = 72-52, 
       "Set1 Normed" = 66-52,
       "Set1 Unnormed&Set1 Normed" = 52)
#png("Set1_intersect.png")
v2 = plot(euler(dat),
          fills = list(fill = c("white", "darkgrey"), alpha = 0.9),
          quantities = list(type = c("counts","percent"),cex=1.5),
          #labels = list(font=3,side="above"),
          legend = list(side = "right",cex=1.5)
)
tags <- v2$children[[2]]$children[[1]]$children$tags$children
tags <- do.call(grid::gList, lapply(tags, function(x) {
  x$children[[2]]$label <- sub(" \\(", "\n(", x$children[[2]]$label)
  x$children[[2]]$just <- NULL
  x$children[[2]]$hjust <- 0.5
  x$children[[2]]$vjust <- 1
  x}))

v2$children[[2]]$children[[1]]$children$tags$children <- tags

ggsave("11-25-set1.pdf",v2)







dat<-c("Joint Unnormed" = 36-12, 
       "Joint Normed" = 14-12,
       "Joint Unnormed&Joint Normed" = 12)
#png("Joint_intersect.png")
v2 = plot(euler(dat),
          fills = list(fill = c("white", "darkgrey"), alpha = 0.9),
          quantities = list(type = c("counts","percent"),cex=1.5),
          #labels = list(font=3,side="above"),
          legend = list(side = "right",cex=1.5)
)
tags <- v2$children[[2]]$children[[1]]$children$tags$children
tags <- do.call(grid::gList, lapply(tags, function(x) {
  x$children[[2]]$label <- sub(" \\(", "\n(", x$children[[2]]$label)
  x$children[[2]]$just <- NULL
  x$children[[2]]$hjust <- 0.5
  x$children[[2]]$vjust <- 1
  x}))

v2$children[[2]]$children[[1]]$children$tags$children <- tags

ggsave("11-25-joint.pdf",v2)
