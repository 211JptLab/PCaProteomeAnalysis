library(eulerr)
library(readr)
library(dplyr)

library(tidyr)


## unnormed


phos_18 <-  read_delim("./18PLEX/wetransfer_pdx-progression-files_2022-08-26_2220_backup/human.phospho.gct_n36x31538.gct",
                       skip=8,col_names = F, delim="\t")
phos_11 <-  read_delim("./PILOT/wetransfer-c07d95_BCK/panoply_medianMAD_phosphoproteome.gct",
                       skip=11,col_names = F, delim="\t")
phos_ltl = read_delim("LTL331/wetransfer_ltl-progression-data_2022-03-29_1755_bkp/ltl.human.phospho.gct_n10x35268.gct", 
                      delim = "\t", escape_double = FALSE, 
                      col_names = FALSE, trim_ws = TRUE, skip = 9)
v1=phos_18$X1
v2=phos_11$X1
v3=phos_ltl$X1
length(intersect(v1, v2))
length(intersect(v1, v3))
length(intersect(v2, v3))
length(Reduce(intersect, list(v1, v2, v3)))
library(eulerr)

dat<-c("Set1 Unnormed" = length(setdiff(v1, union(v2, v3))), 
       "Set2 Unnormed" = length(setdiff(v2, union(v1, v3))),
       "LTL331 Unnormed"= length(setdiff(v3, union(v1, v2))), 
       "Set1 Unnormed&Set2 Unnormed" = length(intersect(v1, v2))-length(Reduce(intersect, list(v1, v2, v3))),
       "Set1 Unnormed&LTL331 Unnormed" = length(intersect(v1, v3))-length(Reduce(intersect, list(v1, v2, v3))),
       "Set2 Unnormed&LTL331 Unnormed" = length(intersect(v2, v3))-length(Reduce(intersect, list(v1, v2, v3))),
       "Set1 Unnormed&Set2 Unnormed&LTL331 Unnormed"=length(Reduce(intersect, list(v1, v2, v3))))
v2 = plot(euler(dat),
          fills = list(fill = c("white", "darkgrey", "lightgrey"), alpha = 0.9),
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



ggsave("11-26-unnormed_intersect.pdf",v2)

### normed


phos_18 <-  read_delim("./18PLEX/wetransfer_pdx-progression-files_2022-08-26_2220_backup/panoply_ptm_normalization_human.phospho.gct_n36x31538-proteome-relative-norm.gct",
                       skip=8,col_names = F, delim="\t")
phos_11 <-  read_delim("./PILOT/wetransfer-c07d95_BCK/panoply_ptm_normalization_phosphoproteome.gct",
                       skip=11,col_names = F, delim="\t")
phos_ltl = read_delim("LTL331/wetransfer_ltl-progression-data_2022-03-29_1755_bkp/ltl.human.phospho.gct_n10x33332-proteome-relative-norm.gct", 
                      delim = "\t", escape_double = FALSE, 
                      col_names = FALSE, trim_ws = TRUE, skip = 9)


v1=phos_18$X1
v2=phos_11$X1
v3=phos_ltl$X1
length(intersect(v1, v2))
length(intersect(v1, v3))
length(intersect(v2, v3))
length(Reduce(intersect, list(v1, v2, v3)))
library(eulerr)

dat<-c("Set1 Normed" = length(setdiff(v1, union(v2, v3))), 
       "Set2 Normed" = length(setdiff(v2, union(v1, v3))),
       "LTL331 Normed"= length(setdiff(v3, union(v1, v2))), 
       "Set1 Normed&Set2 Normed" = length(intersect(v1, v2))-length(Reduce(intersect, list(v1, v2, v3))),
       "Set1 Normed&LTL331 Normed" = length(intersect(v1, v3))-length(Reduce(intersect, list(v1, v2, v3))),
       "Set2 Normed&LTL331 Normed" = length(intersect(v2, v3))-length(Reduce(intersect, list(v1, v2, v3))),
       "Set1 Normed&Set2 Normed&LTL331 Normed"=length(Reduce(intersect, list(v1, v2, v3))))
v2 = plot(euler(dat),
          fills = list(fill = c("white", "darkgrey", "lightgrey"), alpha = 0.9),
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
ggsave("11-26-normed_intersect.pdf",v2)


