require(tidyverse)
require(tidygraph)
require(igraph)
require(ggplot2)
require(ggraph)
require(visNetwork)
require(lsa)
require(WGCNA)
library(readr)
require(bipartite)
require(ade4)

dat <- read_csv("Catalog_AC.csv")

# Generate bipartite graph from unique provenience & artifact code pairs
assemblages_bpg <-
  graph.data.frame(unique.data.frame(select(dat, LEVEL_ID, CODE)),
                   directed = FALSE)

V(assemblages_bpg)$type <- bipartite.mapping(assemblages_bpg)$type

# Create incidence matrix from bipartite graph
assemblages_inc <- as_incidence_matrix(assemblages_bpg)

# Extract each mode's adjacency matrix from incidence matrix

prov_adj <-
  assemblages_inc %*% t(assemblages_inc) # Overlap by counts

prov_adj_jacc <- # by Jaccard index
  1 - as.matrix(dist.binary(
    assemblages_inc,
    method = 1,
    upper = TRUE,
    diag = FALSE
  ))

artifact_adj <- t(assemblages_inc) %*% assemblages_inc





