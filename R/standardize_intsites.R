#' gr <- GRanges(
#'   seqnames = rep("chr1", 100),
#'   ranges = IRanges(
#'     start = sample(1:200, 100, replace = TRUE),
#'     width = sample(30:200, 100, replace = TRUE)),
#'   strand = rep("+", 100))
#'


standardize_intsites <- function(unstandardized.sites,
                                 std.gap = 1L,
                                 standardize_breakpoints = FALSE,
                                 get.analysis.data = FALSE,
                                 calc.cluster.stats = FALSE){

  #Build data.frame with called info and cluster information
  raw.sites <- unstandardized.sites
  raw.positions <- flank(granges(raw.sites), -1, start = TRUE)
  raw.breakpoints <- flank(granges(raw.sites), -1, start = FALSE)
  raw.positions <- serial_cluster(raw.positions, gaps = c(std.gap, 5L))
  raw.breakpoints <- serial_cluster(raw.breakpoints, gaps = c(std.gap))
  clus.dfr <- data.frame(
    "order" = seq(1:length(raw.sites)),
    "seqnames" = seqnames(raw.sites),
    "strand" = strand(raw.sites),
    "called.pos" = start(raw.positions),
    "called.bp" = end(raw.breakpoints),
    "window.clus" = as.character(as.data.frame(mcols(raw.positions))[,2]),
    "pos.clus" = as.character(as.data.frame(mcols(raw.positions))[,1]),
    "bp.clus" = as.character(as.data.frame(mcols(raw.breakpoints))[,1]),
    stringsAsFactors = FALSE)

  #Using the cluster membership for both ends of the reads, determine the
  #standardized membership for both ends
  clus.list <- lapply(split(clus.dfr, clus.dfr$window.clus), function(clus.dfr){
    bp.pos.clus <- lapply(unique(clus.dfr$bp.clus), function(bp){
      unique(clus.dfr[clus.dfr$bp.clus == bp,]$pos.clus)
    })
    names(bp.pos.clus) <- unique(clus.dfr$bp.clus)
    pos.bp.clus <- lapply(unique(clus.dfr$pos.clus), function(pos){
      unique(clus.dfr[clus.dfr$pos.clus == pos,]$bp.clus)
    })
    names(pos.bp.clus) <- unique(clus.dfr$pos.clus)
    pos.el <- as.matrix(bind_rows(lapply(unique(clus.dfr$pos.clus), function(pos){
      bps <- pos.bp.clus[[pos]]
      subject <- unique(unlist(bp.pos.clus[bps]))
      dfr <- data.frame("query.pos" = pos,
                        "subject.pos" = subject,
                        stringsAsFactors = FALSE
      )
      dfr
    })), ncol = 2)
    pos.graph <- graph.edgelist(pos.el)
    adj.pos.clus <- clusters(pos.graph)$membership
    clus.dfr$adj.pos.clus <- adj.pos.clus[clus.dfr$pos.clus]

    #Determine the standardized position and breakpoints independently
    #Positions
    pos.dfr <- clus.dfr[,c("seqnames", "strand","called.pos",
                           "adj.pos.clus", "bp.clus")]
    pos.list <- split(pos.dfr, pos.dfr$adj.pos.clus)
    std.positions <- as.data.frame(bind_rows(lapply(pos.list, function(dfr){
      dfr <- unique(dfr)
      pos.freq <- as.data.frame(table(dfr$called.pos))
      top.freq <- as.numeric(as.character(
        pos.freq[pos.freq$Freq == max(pos.freq$Freq), "Var1"]
      ))
      std.pos <- ifelse(unique(dfr$strand) == "+", min(top.freq), max(top.freq))
      adj.pos.clus.dfr <- data.frame("adj.pos.clus" = unique(dfr$adj.pos.clus),
                                     "std.pos" = as.integer(std.pos),
                                     stringsAsFactors = FALSE)
      adj.pos.clus.dfr
    })))
    row.names(std.positions) <- std.positions$adj.pos.clus
    clus.dfr$std.pos <- std.positions[clus.dfr$adj.pos.clus, "std.pos"]

    #Breakpoints (optional currently, only recommended on raw ranges or within a
    #replicate, not across replicates or samples)
    if(standardize_breakpoints){
      bp.dfr <- clus.dfr[,c("seqnames", "strand","called.bp",
                            "adj.pos.clus", "bp.clus")]
      bp.list <- split(bp.dfr, bp.dfr$bp.clus)
      std.breakpoints <- as.data.frame(bind_rows(lapply(bp.list, function(dfr){
        dfr <- distinct(dfr)
        bp.freq <- as.data.frame(table(dfr$called.bp))
        top.freq <- as.numeric(as.character(
          bp.freq[bp.freq$Freq == max(bp.freq$Freq), "Var1"]
        ))
        std.bp <- ifelse(unique(dfr$strand) == "+", max(top.freq), min(top.freq))
        adj.bp.clus.dfr <- data.frame("adj.bp.clus" = unique(dfr$bp.clus),
                                      "std.bp" = as.integer(std.bp),
                                      stringsAsFactors = FALSE)
        adj.bp.clus.dfr
      })))
      row.names(std.breakpoints) <- std.breakpoints$adj.bp.clus
      clus.dfr$std.bp <- std.breakpoints[clus.dfr$bp.clus, "std.bp"]
    }else{
      clus.dfr$std.bp <- clus.dfr$called.bp
    }
    clus.dfr
  })
  #Rebuild Granges of standardized sites
  clus.dfr <- do.call(rbind,
                      lapply(1:length(clus.list), function(i) clus.list[[i]]))
  clus.dfr <- arrange(clus.dfr, order)
  std.ranges <- IRanges(start = ifelse(clus.dfr$strand == "+",
                                       clus.dfr$std.pos, clus.dfr$std.bp),
                        end = ifelse(clus.dfr$strand == "+",
                                     clus.dfr$std.bp, clus.dfr$std.pos))
  std.sites <- GRanges(seqnames = clus.dfr$seqnames,
                       ranges = std.ranges,
                       strand = clus.dfr$strand,
                       seqinfo = seqinfo(raw.sites))
  mcols(std.sites) <- mcols(raw.sites)

  if(get.analysis.data){
    std.sites$called.pos <- clus.dfr$called.pos
    std.sites$called.bp <- clus.dfr$called.bp
    std.sites$pos.clus.ori <- clus.dfr$pos.clus
    std.sites$pos.clus.adj <- clus.dfr$adj.pos.clus
    std.sites$bp.clus <- clus.dfr$bp.clus
  }

  if(calc.cluster.stats){
    clus.list <- split(clus.dfr, clus.dfr$adj.pos.clus)
    clus.stats <- as.data.frame(bind_rows(lapply(clus.list, function(clus){
      bp.range.diff <- diff(range(clus$called.bp))
      bp.ori.count <- length(unique(clus$called.bp))
      bp.clus.count <- length(unique(clus$bp.clus))
      bp.mean.diff <- ifelse(nrow(clus) > 1, mean(diff(clus$called.bp)), 0)
      stats.dfr <- data.frame(bp.range.diff, bp.ori.count,
                              bp.clus.count, bp.mean.diff)
      stats.dfr
    })))
    row.names(clus.stats) <- as.numeric(names(clus.list))

    #append the data to std.sites
    std.sites$bp.range.diff <- clus.stats[clus.dfr$adj.pos.clus,"bp.range.diff"]
    std.sites$bp.ori.count <- clus.stats[clus.dfr$adj.pos.clus, "bp.ori.count"]
    std.sites$bp.clus.count <- clus.stats[clus.dfr$adj.pos.clus, "bp.clus.count"]
    std.sites$bp.mean.diff <- clus.stats[clus.dfr$adj.pos.clus, "bp.mean.diff"]
  }

  return(std.sites)
}
