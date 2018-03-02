#data.tree
test2<-structure(list(genus = structure(c(4L, 2L, 7L, 8L, 6L, 1L, 3L, 
5L, 5L), .Label = c("Aminobacter", "Bradyrhizobium", "Hoeflea", 
"Hyphomonas", "Mesorhizobium", "Methylosinus", "Ochrobactrum", 
"uncultured"), class = "factor"), family = structure(c(4L, 1L, 
2L, 3L, 5L, 6L, 6L, 6L, 6L), .Label = c("Bradyrhizobiaceae", 
"Brucellaceae", "Hyphomicrobiaceae", "Hyphomonadaceae", "Methylocystaceae", 
"Phyllobacteriaceae"), class = "factor"), order = structure(c(1L, 
2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L), .Label = c("Caulobacterales", 
"Rhizobiales"), class = "factor"), class = structure(c(1L, 1L, 
1L, 1L, 1L, 1L, 1L, 1L, 1L), .Label = "Alphaproteobacteria", class = "factor"), 
    phylum = structure(c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L), .Label = "Proteobacteria", class = "factor")), .Names = c("genus", 
"family", "order", "class", "phylum"), class = "data.frame", row.names = c(NA, 
9L))
test2<-read.csv(file.choose(),sep=";",row.names=1)
library(data.tree)
test2$pathString <- with(test2, 
                               paste(kingdom,
									 phylum,
                                     class,
                                     order,
                                     family,
                                     genus, 
									 species,sep = "/"))

tree_test2 = as.Node(test2)
print(tree_test2)

plot(as.dendrogram(tree_test2),center=TRUE)
SetGraphStyle(tree_test2,rankdir="TB")
SetNodeStyle(tree_test2, shape = "box", fillcolor = "GreenYellow")
#SetNodeStyle(tree_test2$IT, fillcolor = "LightBlue", penwidth = "5px")
plot(tree_test2)

ggdendrogram(as.dendrogram(tree_test2))
test3<-as.dendrogram(tree_test2,center=TRUE)
ggdendrogram(test3)
ddata<-dendro_data(test3, type="rectangle")
p <- ggplot(segment(ddata),center=TRUE) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
  coord_flip() + 
  scale_y_reverse(expand = c(0.2, 0))
p







library(igraph)
library(ggraph)
graph = as.igraph(tree_test2, directed = TRUE, direction = "climb")
graph2 = as.igraph(tree_test2, directed = FALSE)

ggraph(graph, layout = 'kk') + 
  geom_node_text(aes(label = name))+
  geom_edge_link(arrow = arrow(type = "closed", ends = "first",
                               length = unit(0.20, "inches"),
                               angle = 15)) +
  geom_node_point() +
  theme_graph()+
  coord_cartesian(xlim = c(-3,3), expand = TRUE)
  
ggraph(dendrogram,'dendrogram') +
	geom_edge_elbow()
	
	
dendrogram <- as.dendrogram(graph)