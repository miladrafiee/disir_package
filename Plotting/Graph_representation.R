install.packages("visNetwork")
install.packages("igraph")
library(visNetwork)
library(RColorBrewer)

#################################################################################
############################# Path to input data ################################
#################################################################################
nodes <- read.csv('/Users/i0461476/Desktop/disir_package/results_RA/IL6_IL6R_and_IL6_IL6ST/graph_data/Nodes.csv', header=T, as.is=T)
links <- read.csv('/Users/i0461476/Desktop/disir_package/results_RA/IL6_IL6R_and_IL6_IL6ST/graph_data/Links.csv', header=T, as.is=T)
interactions_raw <- read.csv('/Users/i0461476/Desktop/disir_package/data_test/RA/Subunit_interactions_IL6_IL11.csv', header = FALSE)
interactions_raw <- data.frame(interactions_raw[1,])

##################################################################################
############################# Path to output data ################################
##################################################################################
output_path <- "/Users/i0461476/Desktop/disir_package/results_RA/IL6_IL6R_and_IL6_IL6ST/network_plots/interactions_IL6.html"

##################################################################################
interactions <-as.data.frame(matrix(nrow = nrow(interactions_raw), ncol = 1))
for (i in 1:nrow(interactions_raw)){
  interactions[i, 1] <- paste(interactions_raw[i,1], '+', interactions_raw[i,2], sep = "")
}
colnames(interactions) <- 'V1'
links <- links[links$name %in% interactions$V1,]
vertex_labels <-unique(nodes$name)
noClasses <- length(vertex_labels)
# Edge_color <- rainbow(n = nrow(interactions))
Edge_color <- 'blue'
Edge_color_vec <- vector(length = length(links$name))
for (i in 1:length(links$name)){
  Edge_color_vec[i] <- Edge_color[which(interactions$V1 %in% links$name[i])]
}
nodes$id[nodes$id > noClasses] <- nodes$id[nodes$id > noClasses] - noClasses
links$to[links$to > noClasses] <- links$to[links$to > noClasses] - noClasses
links$from[links$from > noClasses] <- links$from[links$from > noClasses] - noClasses

graph_nodes <- data.frame(id = nodes$id[1:noClasses], label = vertex_labels, value = rep(10, times = noClasses),  font = list(size = 28, bold = TRUE))
noEdges <- nrow(links)
graph_edges <- data.frame(links[c('from', 'to')], value = links$weight, scaling = list(max = 5, min = 0.1), arrows = list(to = list(enabled = TRUE, scaleFactor = 2)), color = Edge_color_vec, smooth = c(TRUE))
ledges <- data.frame(color = Edge_color, label = interactions$V1, arrows = "to", font = list(size = 8, bold = FALSE))
visNetwork(graph_nodes, graph_edges, width = "100%") %>% 
  visLegend(addEdges = ledges, position = "right", width = 0.2, main = list(text = "LR Interactions",
                                                                            style = "font-family:Arial;color:#ff0000;font-size:24px;text-align:center;"), ncol = 2) %>% 
  visIgraphLayout(layout = "layout_in_circle") %>%
  visSave(file = output_path)