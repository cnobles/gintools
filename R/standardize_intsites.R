#' Standardize genomic integration site positions within a dataset
#'
#' \code{standardize_intsites} returns a GRanges object where the integration
#' site positions have been standardized with sites within the gap distance.
#'
#' @description Given a GRanges object ...
#'
#' @usage
#' standardize_intsites(sites)
#'
#' standardize_intsites(sites, min.gap = 1L, sata.gap = 5L)
#'
#' @param sites GRanges Object integration sites to standardize.
#' @param min.gap integer minimum gap to consider combine integration sites.
#' @param sata.gap integer maximum distance to consider combining integration
#' sites.
#'
#' @examples
#' gr <- gintools:::.generate_test_granges()
#'
#' standardize_intsites(gr)
#'
#' @author Christopher Nobles, Ph.D.
#' @export

standardize_intsites <- function(sites, min.gap = 1L, sata.gap = 5L){
  require(gintools)

  # Start by reducing the sites object down to only the "intSite" positions
  # and storing the revmap for going back to the original sites.
  red.sites <- reduce(
    flank(sites, -1, start = TRUE),
    min.gapwidth = 0L,
    with.revmap = TRUE)
  red.sites$siteID <- seq(1:length(red.sites))
  revmap <- as.list(red.sites$revmap)
  red.sites$fragLengths <- sapply(revmap, length) + runif(length(revmap), max = 0.01)
  # Not true if original sites are not dereplicated or unique
  red.hits <- GenomicRanges::as.data.frame(
    findOverlaps(red.sites, maxgap = min.gap, ignoreSelf = TRUE))

  red.hits$q_fragLengths <- red.sites[red.hits$queryHits]$fragLengths
  red.hits$s_fragLengths <- red.sites[red.hits$subjectHits]$fragLengths
  red.hits <- red.hits[red.hits$q_fragLengths >= red.hits$s_fragLengths,]

  # Construct initial graph, the graph will be modified during other parts of
  # the function and can be used to identify clusters as well as 'sources'
  # within those clusters as it's a directed graph. Sources are identified and
  # used to identify the original position given a distribution.
  g <- make_empty_graph(n = length(red.sites), directed = TRUE) %>%
    add_edges(unlist(mapply(
      c, red.hits$queryHits, red.hits$subjectHits, SIMPLIFY = FALSE)))

  sources <- which(Matrix::colSums(get.adjacency(g, sparse = TRUE)) == 0)
  red.sites$clusID <- clusters(g)$membership

  # Identify satalite positions that should be included in clusters up to 5nt
  # away. This portion of the function tries to reach out to up to 5nt from the
  # boundry of a cluster to see if there are any further "satalite" positions
  # that have been annotated. It does this by iterively increasing the size
  # from 2nt to 5nt by 1nt increments.
  lapply(2:sata.gap, function(i){
#    clus.ranges <- unlist(range(split(red.sites, clusters(g)$membership))) #Why can't 'range' work right?
    clus.ranges <- unlist(reduce(GRangesList(split(red.sites, clusters(g)$membership)), min.gapwidth = (i-1)))
#    clus.ranges <- unlist(reduce(split(red.sites, clusters(g)$membership), min.gapwidth = (i-1)))
    sata.hits <- as.data.frame(
      findOverlaps(clus.ranges, maxgap = i, ignoreSelf = TRUE)
    )
    names(sata.hits) <- c("source_clus", "sata_clus")

    red.df <- GenomicRanges::as.data.frame(red.sites)

    if(nrow(sata.hits) > 0){
      clus.data <-
        group_by(red.df, clusID) %>%
        summarize(
          clus_pos_mean = as.integer(mean(start)),
          min_fragLengths = min(fragLengths),
          sum_fragLengths = sum(fragLengths))

      sata.hits <- sata.hits %>%
        mutate(source_pos = clus.data[source_clus,]$clus_pos_mean) %>%
        mutate(sata_pos = clus.data[sata_clus,]$clus_pos_mean) %>%
        mutate(min_q_fragLengths = clus.data[.$source_clus,]$min_fragLengths) %>%
        mutate(min_s_fragLengths = clus.data[.$sata_clus,]$min_fragLengths) %>%
        mutate(q_sumFragLengths = clus.data[.$source_clus,]$sum_fragLengths) %>%
        mutate(s_sumFragLengths = clus.data[.$sata_clus,]$sum_fragLengths) %>%
        mutate(is_upstream = source_pos < sata_pos) %>%
        filter(q_sumFragLengths > s_sumFragLengths) %>%
        filter(as.integer(.$min_q_fragLengths) >= as.integer(.$min_s_fragLengths))

      if(nrow(sata.hits) > 0){
        clus.map <- findOverlaps(clus.ranges, red.sites)
        clus.list <- split(subjectHits(clus.map), queryHits(clus.map))

        sata.hits <- sata.hits %>%
          mutate(source_node = ifelse(
            sata.hits$is_upstream,
            sapply(clus.list[sata.hits$source_clus], last),
            sapply(clus.list[sata.hits$source_clus], first)
          )) %>%
          mutate(sata_node = ifelse(
            is_upstream,
            sapply(clus.list[sata_clus], first),
            sapply(clus.list[sata_clus], last)
          ))

        sata.edges <- unlist(with(
          sata.hits,
          mapply(c, source_node, sata_node, SIMPLIFY = FALSE)
        ))
      }else{
        sata.edges <- c()
      }
    }else{
      sata.edges <- c()
    }
    g <<- add_edges(g, sata.edges)
    sources <<- which(Matrix::colSums(get.adjacency(g, sparse = TRUE)) == 0)
    red.sites$clusID <<- clusters(g)$membership
  })

  # During these steps of expanding the graph, it's possible that clusters grew
  # larger than they needed to be. This could be identified by seeing multiple
  # sources within a single cluster. This may arrize by reaching out and
  # connecting the edge of one smaller cluster to a larger one. If this occures,
  # the path between the sources will be 'snipped' between the nodes that share
  # the largest distance appart from one another.

  sources.p.clus <- split(sources, clusters(g)$membership[sources])

  clus.w.multi.sources <- sources.p.clus[sapply(sources.p.clus, length) > 1]

  if(length(clus.w.multi.sources) > 0){
    adj.pairs <- do.call(c, lapply(clus.w.multi.sources, function(x){
      lapply(1:(length(x)-1), function(i) c(x[i], x[i+1]))
    }))

    sink.nodes <- which(Matrix::rowSums(get.adjacency(g, sparse = TRUE)) == 0)

    edges.to.edit <- data.frame(
      "src_node_i" = sapply(adj.pairs, "[[", 1),
      "src_node_j" = sapply(adj.pairs, "[[", 2)
    ) %>%
      mutate("src_node_i_abund" = as.numeric(red.sites[src_node_i]$fragLengths)) %>%
      mutate("src_node_j_abund" = as.numeric(red.sites[src_node_j]$fragLengths))

    source.paths <- mapply(function(x,y){
      all_simple_paths(as.undirected(g), x, y)},
      edges.to.edit$src_node_i,
      edges.to.edit$src_node_j)

    edges.to.edit <- mutate(edges.to.edit, sink_node = sink.nodes[
      sapply(source.paths, function(x){
        which(sink.nodes %in% x)
      })])

    nodes.adj.to.sink <- bind_rows(lapply(1:nrow(edges.to.edit), function(i){
      sink <- edges.to.edit[i, "sink_node"]
      path <- as.numeric(source.paths[[i]])
      pos <- grep(sink, path)
      data.frame(
        "sink" = rep(sink, 2),
        "adj_node" = c(path[pos-1], path[pos+1])
      )
    })) %>%
      mutate(sink_pos = start(red.sites[sink])) %>%
      mutate(adj_pos = start(red.sites[adj_node])) %>%
      mutate(adj_abund = red.sites[adj_node]$fragLengths) %>%
      mutate(nt_dist = abs(sink_pos - adj_pos)) %>%
      group_by(sink) %>%
      filter(nt_dist == max(nt_dist)) %>%
      filter(adj_abund == min(adj_abund))

    edges.to.edit <- mutate(edges.to.edit, adj_node = nodes.adj.to.sink$adj_node)

    break.edges <- unlist(with(
      edges.to.edit,
      mapply(c, sink_node, adj_node, SIMPLIFY = FALSE)
    ))

    edge.ids.to.break <- get.edge.ids(g, break.edges, directed = FALSE)

    g <- delete_edges(g, edge.ids.to.break)
    sources <- which(Matrix::colSums(get.adjacency(g, sparse = TRUE)) == 0)
    red.sites$clusID <- clusters(g)$membership
  }

  # In the end, sources that are within the range of satalites (5L default),
  # should be grouped together. These sources should be connected by an edge,
  # pointing toward the larger source. The chosen source will have the highest
  # abundance of widths / sonic breaks / fragment lengths.

  near.sources <- findOverlaps(
    red.sites[sources],
    maxgap = sata.gap,
    ignoreSelf = TRUE
  )

  if(length(near.sources) > 0){
    near.src.df <- data.frame(
      node.i = sources[queryHits(near.sources)],
      node.j = sources[subjectHits(near.sources)]
    ) %>%
      mutate(abund.i = red.sites[node.i]$fragLengths) %>%
      mutate(abund.j = red.sites[node.j]$fragLengths) %>%
      filter(abund.i > abund.j)

    edges.to.connect.near.srcs <- unlist(with(
      near.src.df,
      mapply(c, node.i, node.j, SIMPLIFY = FALSE)
    ))

    g <- add.edges(g, edges.to.connect.near.srcs)
    sources <- which(Matrix::colSums(get.adjacency(g, sparse = TRUE)) == 0)
    red.sites$clusID <- clusters(g)$membership
  }

  # As cluster membership has been determined, the remaining section of the
  # function serves to consolidate the information and correct the original
  # sites object with the standardized positions.

  clus.data <- data.frame(
    "clusID" = 1:clusters(g)$no,
    "chr" = seqnames(red.sites[sources]),
    "strand" = strand(red.sites[sources]),
    "position" = start(red.sites[sources]),
#    "width" = width(unlist(range(split(red.sites, clusters(g)$membership))))
    "width" = width(unlist(reduce(GRangesList(split(red.sites, clusters(g)$membership)), min.gapwidth = sata.gap)))
  )

  sites <- sites[unlist(as.list(red.sites$revmap))]
  sites$clusID <- as.numeric(Rle(
    clusters(g)$membership,
    sapply(red.sites$revmap, length)))
  sites$called.pos <- ifelse(strand(sites) == "+", start(sites), end(sites))
  sites$adj.pos <- clus.data[match(sites$clusID, clus.data$clusID), "position"]

  ranges(sites) <- IRanges(
    start = ifelse(strand(sites) == "+", sites$adj.pos, start(sites)),
    end = ifelse(strand(sites) == "+", end(sites), sites$adj.pos)
  )

  sites
}
