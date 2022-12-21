install.packages("RColorBrewer")
library(RColorBrewer)
library(circlize)
library(graphics)
library(grid)

######################################## One threshold #######################################
cord_diagram_mat <- read.csv('/Users/i0461476/Desktop/Papers/NAR_revision/Cord_diagram/RA_IL6/Heatmap.csv', header = TRUE, row.names = 1)
cord_diagram_mat_matrix <- as.matrix(cord_diagram_mat)

grid_colors_all = c(brewer.pal(12,"Paired"), brewer.pal(3,"Set2"))
grid_colors_names_all = colnames(cord_diagram_mat)
names(grid_colors_all) <- grid_colors_names_all

a <- chordDiagram(cord_diagram_mat_matrix, annotationTrack = c("grid"), grid.col = grid_colors_all)
grid_colors_names_col <- c(unique(a$rn),unique(a$cn))
grid_colors_col <- grid_colors_all[names(grid_colors_all) %in% grid_colors_names_col]
legend("right", inset=c(-0.26, 0),
       legend = names(grid_colors_col),
       fill = grid_colors_col,       # Color of the squares
       border = NA,
       title = "Cell type")