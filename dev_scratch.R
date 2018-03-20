library(GenomicRanges)
#library(gintools)
library(ggplot2)
library(stringr)
library(reshape2)

testDataGr <- gintools:::generate_test_granges(
  positions = c(150, 250, 350, 450), n_reads_p_site = 50, stdev = 2.5)

initialPlot <- ggplot(as.data.frame(testDataGr), aes(x = start)) + geom_bar()

initialRanges <- range(width(reduce(flank(testDataGr, -1), min.gapwidth = 5L)))

stdSitesOutput <- reduce(flank(standardize_intsites(testDataGr), -1), min.gapwidth = 0L)

# identify redundant / bidirectional edges within graphs

get_bidirectional_edge_ids <- function(graph){
  require(dplyr)
  edLogic <- get.data.frame(graph) %>%
    mutate(n = 1:n()) %>%
    group_by(n) %>%
    mutate(eid = paste0(min(from, to), ":", max(from, to))) %>%
    ungroup() %>%
    group_by(eid) %>%
    mutate(count = n()) %>%
    ungroup() %>%
    filter(count == 2)

  if(nrow(edLogic) > 0){
    edges <- unlist(mapply(c, edLogic$from, edLogic$to, SIMPLIFY = FALSE))
    return(sort(get.edge.ids(graph, edges)))
  }else{
    return(c())
  }
}

set.seed(1)
testData <- round(runif(25, 1, 20))
bias <- "upstream"

p0 <- ggplot(data.frame(num = testData), aes(x = num)) + geom_bar()


# find local maxima
## build linear connections to form graph
ir <- IRanges(start = testData, width = 1)
base.ir <- reduce(ir, min.gapwidth = 0L, with.revmap = TRUE)
df <- data.frame(
  pos = start(base.ir),
  count = width(base.ir@elementMetadata$revmap@partitioning)
)

g <- make_empty_graph(n = length(base.ir), directed = TRUE)
pg0 <- plot.igraph(g)

### connect within 1 unit away, bais upstream
hits <- findOverlaps(base.ir, maxgap = 1L, drop.self = TRUE)
el <- data.frame(
  src = queryHits(hits),
  snk = subjectHits(hits),
  srcPos = df$pos[queryHits(hits)],
  snkPos = df$pos[subjectHits(hits)],
  srcCount = df$count[queryHits(hits)],
  snkCount = df$count[subjectHits(hits)])

el <- el[el$srcCount >= el$snkCount,]
#el <- el[ifelse(el$srcCount == el$snkCount, el$srcPos < el$snkPos, TRUE),]
ed <- unlist(with(el, mapply(c, src, snk, SIMPLIFY = FALSE)))

g.1 <- add_edges(g, ed)
pg1 <- plot(
  g.1, edge.arrow.size = 2/3, edge.arrow.width = 1, vertex.size = 5,
  frame = TRUE, asp = 0, edge.curved = FALSE,
  layout = matrix(c(df$pos, df$count), ncol = 2))

### break paths at sink nodes
adj.vrt <- adjacent_vertices(g.1, sinks(g.1), mode = "in")
conflicts <- sapply(adj.vrt, length) == 2
adj.mat <- matrix(
  data = unlist(adj.vrt[conflicts]),
  ncol = 2,
  byrow = TRUE)
br.df <- data.frame(
  snk = sinks(g.1)[conflicts],
  adj.node1 = adj.mat[,1],
  adj.node2 = adj.mat[,2],
  node1.pos = df$pos[adj.mat[,1]],
  node2.pos = df$pos[adj.mat[,2]],
  node1.count = df$count[adj.mat[,1]],
  node2.count = df$count[adj.mat[,2]])

# identify edges to break, reverse logic
break.node <- with(br.df, ifelse(node1.count < node2.count, 1, 2))
if(bias == "upstream"){
  break.node <- with(br.df, ifelse(
    node1.count == node2.count,
    ifelse(node1.pos < node2.pos, 2, 1),
    break.node))
}else if(bias == "downstream"){
  break.node <- with(br.df, ifelse(
    node1.count == node2.count,
    ifelse(node1.pos < node2.pos, 1, 2),
    break.node))
}else{
  stop("Bias variable should be either 'upstream' or 'downstream'.")
}

break.ed <- unlist(mapply(
  c,
  br.df$snk,
  mapply(
    function(i,j) br.df[i,j],
    i = 1:nrow(br.df),
    j = 1+break.node),
  SIMPLIFY = FALSE))
edges.to.break <- get.edge.ids(g.1, break.ed, directed = FALSE)
g.2 <- delete_edges(g.1, edges.to.break)
pg2 <- plot(
  g.2, edge.arrow.size = 2/3, edge.arrow.width = 1, vertex.size = 5,
  frame = TRUE, asp = 0, edge.curved = FALSE,
  layout = matrix(c(df$pos, df$count), ncol = 2))

# resolve bidirectional edges
clus.max.counts <- sapply(split(df$count, clusters(g.2)$membership), max)
bi.nodes <- as.numeric(names(table(as.matrix(
  get.data.frame(g.2)[get_bidirectional_edge_ids(g.2),]))))
adj.bi.nodes <- adjacent_vertices(g.2, bi.nodes, mode = "all")
names(adj.bi.nodes) <- paste0(bi.nodes, ":")
adj.bi.nodes <- unlist(adj.bi.nodes)
bi.node.df <- unique(data.frame(
  node = as.numeric(str_extract(names(adj.bi.nodes), "[0-9]+")),
  adj.node = adj.bi.nodes,
  row.names = NULL))
bi.node.df$adj.node.num <- ifelse(
  duplicated(bi.node.df$node), "adj.node.2", "adj.node.1")
bi.node.df <- dcast(bi.node.df, node ~ adj.node.num, value.var = "adj.node")
bi.node.df$node.count <- df$count[bi.node.df$node]
bi.node.df$adj.1.count <- df$count[bi.node.df$adj.node.1]
bi.node.df$adj.2.count <- df$count[bi.node.df$adj.node.2]
bi.node.df$clus <- clusters(g.2)$membership[bi.node.df$node]
bi.node.df$clus.max.count <- clus.max.counts[bi.node.df$clus]
bi.node.df$term.node <- is.na(bi.node.df$adj.node.2)
bi.node.df$peak.node <- with(bi.node.df, node.count == clus.max.count)

# peak node edges will be determined by bias
bi.node.df$rm.edge <- with(bi.node.df, ifelse(
  peak.node, ifelse(node < adj.node.1,

  )
))


# Connect maxima within reasonable distance (0L for now, skip)

df$mem <- clusters(g.2)$membership
aggregate(pos ~ mem, df, function(x) x[x$count == max(x),])
