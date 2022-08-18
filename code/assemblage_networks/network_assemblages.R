require(tidyverse)
require(tidygraph)
require(igraph)
require(ggplot2)
require(ggraph)
require(visNetwork)
require(lsa)
require(WGCNA)
library(readr)

dat <- read_csv("Catalog_A2.csv")


# Contingency tables of artifacts and locations

prov_assemb <- dat %>%
  group_by(LEVEL_ID, CODE) %>%
  summarise(n = n(), .groups = "drop") %>%
  spread(CODE, n) %>%
  mutate_all(~ replace_na(., 0))

assemb_prov <- dat %>%
  group_by(CODE, LEVEL_ID) %>%
  summarise(n = n(), .groups = "drop") %>%
  spread(LEVEL_ID, n) %>%
  mutate_if(is.numeric, ~ replace_na(., 0))

# Create graph network of diagnostic artifacts by find locations

## Create co-presence matrices from contingency tables 

### Assemblages

assemb_adj <- as.matrix(select(prov_assemb, -LEVEL_ID))
assemb_adj[assemb_adj > 0] <-1 # Binarize matrix
assemb_cos <- cosine(assemb_adj) # Calculate cosine similarity matrix
assemb_thresh <-
  quantile(c(assemb_cos), probs = 0.90) # Quantile threshold for similarity
assemb_cos[assemb_cos < assemb_thresh] <- 0

### Provenience

prov_adj <- as.matrix(select(assemb_prov, -CODE))
prov_adj[prov_adj > 0] <-1 # Binarize matrix
prov_cos <- cosine(prov_adj) # Calculate cosine similarity matrix
prov_thresh <-
  quantile(c(prov_cos), probs = 0.90) # Quantile threshold for similarity
prov_cos[prov_cos < prov_thresh] <-0

## Create assemblage graph

G_assemb <-
  graph_from_adjacency_matrix(
    WGCNA::TOMsimilarity(assemb_cos),
    mode = "undirected",
    diag = FALSE,
    weighted = TRUE
  )

plot(
  G_assemb,
  vertex.size = 4,
  vertex.label = NA,
  layout = layout_nicely(G_assemb)
)

### Remove isolates
G_assemb_isolates <- which(degree(G_assemb)==0)
G_assemb_conn <- delete.vertices(G_assemb, G_assemb_isolates)

### Detect community structure by Louvain method
assemb_partition <- cluster_louvain(G_assemb_conn)

sizes(assemb_partition)
modularity(assemb_partition)

plot(
  assemb_partition,
  G_assemb_conn,
  vertex.size = 4,
  vertex.label = NA,
  layout = layout_with_fr(G_assemb_conn)
)

heatmap(WGCNA::TOMsimilarity(assemb_cos))

## Create provenience graph

G_prov <-
  graph_from_adjacency_matrix(WGCNA::TOMsimilarity(prov_cos),
                              mode = "undirected",
                              diag = FALSE,
                              weighted = TRUE)

plot(
  G_prov,
  vertex.size = 4,
  vertex.label = NA,
  layout = layout_nicely(G_prov)
)

### Remove isolates
G_prov_isolates <- which(degree(G_prov)==0)
G_prov_conn <- delete.vertices(G_prov, G_prov_isolates)

plot(
  G_prov_conn,
  vertex.size = 4,
  vertex.label = NA,
  layout = layout_nicely(G_prov)
)

### Detect community structure by Louvain method
prov_partition <- cluster_louvain(G_prov_conn)

sizes(prov_partition)
modularity(prov_partition)

plot(
  prov_partition,
  G_prov_conn,
  layout = layout_nicely(G_prov_conn),
  vertex.label = NA,
  vertex.size = 4
)

heatmap(WGCNA::TOMsimilarity(prov_cos))
