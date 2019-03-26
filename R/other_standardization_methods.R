#' Wrapper around .clusterSites() for standardizing integration site positions
#' given integration site distributions.
#'
#' @description This function coerces GRanges objects of integration site data
#' into a format to work with the .clusterSites() function for standardizing
#' integration sites to a single nucleotide position in the reference geneome
#' given a set or distribution of integration sites.
#'
#' @author Bushman Lab

.standardizeSites <- function(unstandardizedSites){
  if( ! length(unstandardizedSites) > 0 ){
    return(unstandardizedSites)
  }
  #Get called start values for clustering
  unstandardizedSites$Position <- ifelse(strand(unstandardizedSites) == "+", start(unstandardizedSites), end(unstandardizedSites))
  unstandardizedSites$Break <- ifelse(strand(unstandardizedSites) == "+", end(unstandardizedSites), start(unstandardizedSites))
  unstandardizedSites$Score <- 95
  unstandardizedSites$qEnd <- width(unstandardizedSites)

  #Positions clustered by 5L window and best position is chosen for cluster
  standardized <- gintools:::.clusterSites(
    psl.rd = unstandardizedSites,
    weight = rep(1, length(unstandardizedSites))
  )

  start(standardized) <- ifelse(strand(standardized) == "+",
                                standardized$clusteredPosition, standardized$Break)
  end(standardized) <- ifelse(strand(standardized) == "-",
                              standardized$clusteredPosition, standardized$Break)

  standardized$Position <- NULL
  standardized$Break <- NULL
  standardized$Score <- NULL
  standardized$qEnd <- NULL
  standardized$clusteredPosition <- NULL
  standardized$clonecount <- NULL
  standardized$clusterTopHit <- NULL

  sort(standardized)
}

#' Cluster/Correct values within a window based on their frequency given
#' discrete factors
#'
#' Given a group of discrete factors (i.e. position ids) and integer values,
#' the function tries to correct/cluster the integer values based on their
#' frequency in a defined windowsize.
#'
#' @param posID a vector of groupings for the value parameter (i.e. Chr,strand). Required if psl.rd parameter is not defined.
#' @param value a vector of integer with values that needs to corrected or
#' clustered (i.e. Positions). Required if psl.rd parameter is not defined.
#' @param grouping additional vector of grouping of length posID or psl.rd by
#' which to pool the rows (i.e. samplenames). Default is NULL.
#' @param psl.rd a GRanges object returned from \code{\link{getIntegrationSites}}.
#' Default is NULL.
#' @param weight a numeric vector of weights to use when calculating frequency
#' of value by posID and grouping if specified. Default is NULL.
#' @param windowSize size of window within which values should be corrected or
#' clustered. Default is 5.
#' @param byQuartile flag denoting whether quartile based technique should be
#' employed. See notes for details. Default is TRUE.
#' @param quartile if byQuartile=TRUE, then the quartile which serves as the
#' threshold. Default is 0.70.
#' @param parallel use parallel backend to perform calculation with
#' \code{\link{BiocParallel}}. Defaults to TRUE. If no parallel backend is
#' registered, then a serial version is ran using \code{\link{SerialParam}}.
#' Process is split by the grouping the column.
#' @param sonicAbund calculate breakpoint abundance using
#' \code{\link{getSonicAbund}}. Default is FALSE.
#'
#' @note The algorithm for clustering when byQuartile=TRUE is as follows: for
#' all values in each grouping, get a distribution and test if their frequency
#' is >= quartile threshold. For values below the quartile threshold, test if
#' any values overlap with the ones that passed the threshold and is within the
#' defined windowSize. If there is a match, then merge with higher value, else
#' leave it as is. This is only useful if the distribution is wide and polynodal.
#' When byQuartile=FALSE, for each group the values within the defined window
#' are merged with the next highest frequently occuring value, if freuquencies
#' are tied then lowest value is used to represent the cluster. When psl.rd is
#' passed, then multihits are ignored and only unique sites are clustered. All
#' multihits will be tagged as a good 'clusterTopHit'.
#'
#' @return a data frame with clusteredValues and frequency shown alongside with
#' the original input. If psl.rd parameter is defined then a GRanges object is
#' returned with three new columns appended at the end: clusteredPosition,
#' clonecount, and clusterTopHit (a representative for a given cluster chosen
#' by best scoring hit!).
#'
#' @seealso \code{\link{findIntegrations}}, \code{\link{getIntegrationSites}},
#' \code{\link{otuSites}}, \code{\link{isuSites}}, \code{\link{crossOverCheck}},
#' \code{\link{pslToRangedObject}}, \code{\link{getSonicAbund}}
#'
#' @author Nirav Malani
#' @examples
#' \donttest{
#' .clusterSites(posID=c('chr1-','chr1-','chr1-','chr2+','chr15-',
#' 'chr16-','chr11-'), value=c(rep(1000,2),5832,1000,12324,65738,928042),
#' grouping=c('a','a','a','b','b','b','c'))
#' data(psl)
#' psl <- psl[sample(nrow(psl),100),]
#' psl.rd <- getIntegrationSites(pslToRangedObject(psl))
#' psl.rd$grouping <- sub("(.+)-.+","\\1",psl.rd$qName)
#' .clusterSites(grouping=psl.rd$grouping, psl.rd=psl.rd)
#' }
.clusterSites <- function(posID=NULL, value=NULL, grouping=NULL, psl.rd=NULL,
                         weight=NULL, windowSize=5L, byQuartile=FALSE,
                         quartile=0.70, parallel=TRUE, sonicAbund=FALSE) {

  require("BiocParallel")
  # to avoid 'no visible binding for global variable' NOTE during R check #
  posID2 <- freq <- belowQuartile <- isMax <- isClosest <- val <- NULL
  ismaxFreq <- dp <- NULL

  if("pslFile" %in% names(formals()) | "files" %in% names(formals())) {
    if(is.null(files) | length(files)==0) {
      stop("files parameter empty. Please supply a filename to be read.")
    }

    if(any(grepl("\\*|\\$|\\+|\\^",files))) {
      ## vector of filenames
      files <- list.files(path=dirname(files), pattern=basename(files),
                          full.names=TRUE)
    }

    if(length(files)==0) {
      stop("No file(s) found with given paramter in files:", files)
    }
  }

  if("psl.rd" %in% names(formals())) {
    if(!is.null(psl.rd)) {
      if(length(psl.rd)==0 | !is(psl.rd,"GRanges")) {
        stop("psl.rd paramter is empty or not a GRanges object")
      }
    }
  }

  if(parallel & .Platform$OS.type != "windows"){
    dp <- bpparam()
  }else{
    dp <- SerialParam()
  }

  if(is.null(psl.rd)) {
    stopifnot(!is.null(posID))
    stopifnot(!is.null(value))
  } else {

    ## dereplicate & converge sites! ##
    if(!"Position" %in% colnames(mcols(psl.rd))) {
      stop("The object supplied in psl.rd parameter does not have Position ",
           "column in it. Did you run getIntegrationSites() on it?")
    }

    ## make sure column(s) for picking best hits are there! ##
    isthere <- grepl("score", colnames(mcols(psl.rd)), ignore.case=TRUE)
    if(!any(isthere)) {
      message("No 'score' column found in the data. Using 'qEnd' as an ",
              "alternative to pick the best hit!")
      isthere <- grepl("qEnd", colnames(mcols(psl.rd)), ignore.case=TRUE)
      if(!any(isthere)) {
        stop("No 'qEnd' column found in the data either...can't pick the ",
             "cluster best hit without 'qEnd' or 'score' column present in ",
             "the object supplied in psl.rd :(")
      }
    }

    if(length(which(isthere))>1) {
      message("multiple score based columns found: ",
              paste(colnames(mcols(psl.rd))[which(isthere)],
                    collapse=","), " choosing the first one...")
    }
    isthere <- which(isthere)[1]

    good.row <- rep(TRUE, length(psl.rd))
    multi.there <- grepl("isMultiHit", colnames(mcols(psl.rd)),
                         ignore.case=TRUE)
    if(any(multi.there)) { ## see if multihit column exists
      message("Found 'isMultiHit' column in the data. ",
              "These rows will be ignored for the calculation.")
      good.row <- good.row & !mcols(psl.rd)[[which(multi.there)]]
    }

    posIDs <- paste0(as.character(seqnames(psl.rd)),
                     as.character(strand(psl.rd)))
    values <- mcols(psl.rd)$Position
    if(is.null(weight)) {
      ## see if sequences were dereplicated before in the pipeline which
      ## adds counts=x identifier to the deflines
      weight <- suppressWarnings(as.numeric(sub(".+counts=(\\d+)", "\\1",
                                                mcols(psl.rd)$qName)))
      if(all(is.na(weight))) {
        weight <- NULL }
      else {
        weight[is.na(weight)] <- 1
      }
    }

    grouping <- if(is.null(grouping)) {
      rep("group1",length(values))
    } else {
      grouping
    }

    ## for sonic abundance ##
    mcols(psl.rd)$groups <- grouping

    clusters <- .clusterSites(posID=posIDs[good.row],
                             value=values[good.row],
                             grouping=grouping[good.row],
                             weight=weight[good.row],
                             windowSize=windowSize,
                             byQuartile=byQuartile, quartile=quartile)

    message("Adding clustered data back to psl.rd.")
    clusteredValues <- with(clusters,
                            split(clusteredValue,
                                  paste0(posID,value,grouping)))
    groupingVals <- paste0(posIDs, values, grouping)[good.row]
    mcols(psl.rd)$clusteredPosition <- mcols(psl.rd)$Position
    mcols(psl.rd)$clusteredPosition[good.row] <-
      as.numeric(clusteredValues[groupingVals])

    ## add frequency of new clusteredPosition ##
    clusteredValueFreq <- with(clusters,
                               split(clusteredValue.freq,
                                     paste0(posID,value,grouping)))
    mcols(psl.rd)$clonecount <- 0
    mcols(psl.rd)$clonecount[good.row] <-
      as.numeric(clusteredValueFreq[groupingVals])
    rm("clusteredValueFreq","clusteredValues","clusters")
    cleanit <- gc()

    ## pick best scoring hit to represent a cluster ##
    ## make sure to avoid multihit rows! ##
    message("Picking best scoring hit to represent a cluster.")
    groupingVals <- paste0(as.character(seqnames(psl.rd)),
                           as.character(strand(psl.rd)),
                           mcols(psl.rd)$clusteredPosition, grouping)
    bestScore <- tapply(mcols(psl.rd)[[isthere]][good.row],
                        groupingVals[good.row], max)
    isBest <- mcols(psl.rd)[[isthere]][good.row] ==
      bestScore[groupingVals[good.row]]

    ## pick the first match for cases where >1 reads with the same
    ## coordinate had the same best scores ##
    tocheck <- which(isBest)
    res <- tapply(tocheck, names(tocheck), "[[", 1)

    mcols(psl.rd)$clusterTopHit <- FALSE
    mcols(psl.rd)$clusterTopHit[good.row][res] <- TRUE
    mcols(psl.rd)$clusterTopHit[!good.row] <- TRUE

    message("Cleaning up!")
    rm("isBest","bestScore","posIDs","values","groupingVals")
    cleanit <- gc()

    if(sonicAbund) {
      message("Calculating sonic abundance.")
      psl.rd <- getSonicAbund(psl.rd=psl.rd, grouping=mcols(psl.rd)$groups,
                              parallel=parallel)
    }
    mcols(psl.rd)$groups <- NULL

    return(psl.rd)
  }

  # get frequencies of each posID & value combination by grouping #
  groups <- if(is.null(grouping)) { "" } else { grouping }
  weight2 <- if(is.null(weight)) { 1 } else { weight }
  sites <- dplyr::arrange(data.frame(posID, value, grouping=groups,
                              weight=weight2, posID2=paste0(groups, posID),
                              stringsAsFactors=FALSE), posID2, value)
  sites <- plyr::count(sites, c("posID","value","grouping","posID2"), wt_var="weight")
  rm("groups","weight2")

  if(byQuartile) {
    message("Clustering by quartile: ", quartile)
    # obtain the defined quartile of frequency per posID & grouping #
    sites <- arrange(sites, posID2, value, plyr::desc(freq))
    quartiles <- with(sites,
                      tapply(freq, posID2, quantile, probs=quartile, names=FALSE))
    sites$belowQuartile <- with(sites,freq < quartiles[posID2])
    rm(quartiles)

    if(any(sites$belowQuartile)) {
      # for values belowQuartile, see if any within defined windowSize of
      # aboveQuartile
      pos.be <- with(subset(sites,belowQuartile,drop=TRUE),
                     GRanges(IRanges(start=value,width=1),
                             seqnames=posID2, freq=freq))
      pos.ab <- with(subset(sites,!belowQuartile,drop=TRUE),
                     GRanges(IRanges(start=value,width=1),
                             seqnames=posID2, freq=freq))
      pos.overlap <- as.data.frame(as.matrix(findOverlaps(pos.be, pos.ab,
                                                          maxgap=windowSize-1L,
                                                          ignore.strand=TRUE)))

      # for overlapping values, merge them with the biggest respective
      # aboveQuartile site
      pos.overlap$freq <- values(pos.ab[pos.overlap[,"subjectHits"]])$freq
      pos.overlap$isMax <- with(pos.overlap,
                                ave(freq, as.character(queryHits),
                                    FUN=function(x) x==max(x)))
      pos.overlap$isMax <- as.logical(pos.overlap$isMax)
      pos.overlap <- subset(pos.overlap,isMax, drop=TRUE)

      # if there are >1 biggest respective aboveQuartile site, then choose the
      # closest one ... if tied, then use the latter to represent the site
      counts <- xtabs(isMax~queryHits,pos.overlap)
      if(length(table(counts))>1) {
        toFix <- as.numeric(names(which(counts>1)))
        rows <- pos.overlap$queryHits %in% toFix
        pos.overlap$aboveQuartileValue <-
          pos.overlap$belowQuartileValue <- pos.overlap$valueDiff <- 0
        pos.overlap$aboveQuartileValue[rows] <-
          start(pos.ab[pos.overlap[rows,"subjectHits"]])
        pos.overlap$belowQuartileValue[rows] <-
          start(pos.be[pos.overlap[rows,"queryHits"]])
        pos.overlap$valueDiff[rows] <- with(pos.overlap[rows,],
                                            abs(aboveQuartileValue-
                                                  belowQuartileValue))
        mins <- with(pos.overlap[rows,],
                     tapply(valueDiff, as.character(queryHits), min))
        pos.overlap$isClosest <- TRUE
        pos.overlap$isClosest[rows] <- with(pos.overlap[rows,],
                                            valueDiff ==
                                              mins[as.character(queryHits)])
        pos.overlap <- subset(pos.overlap, isMax & isClosest,drop=TRUE)
        rm("counts","mins")
      }

      # trickle the changes back to the original dataframe#
      pos.overlap$clusteredValue <-
        start(pos.ab[pos.overlap[,"subjectHits"]])
      pos.overlap$posID2 <-
        as.character(seqnames(pos.be[pos.overlap[,"queryHits"]]))

      # for cases where no overlap was found, try clustering to themselves #
      rows <- which(!1:length(pos.be) %in% pos.overlap$query)
      loners <- pos.be[rows]
      if(length(loners)>0) {
        times.rep <- values(loners)[["freq"]]
        res <- .clusterSites(rep(as.character(seqnames(loners)),
                                times=times.rep),
                            rep(start(loners),times=times.rep),
                            byQuartile=FALSE)
      }
      pos.overlap <- rbind(pos.overlap[,c("queryHits","clusteredValue")],
                           data.frame(queryHits=rows,
                                      clusteredValue=
                                        as.numeric(res$clusteredValue)))

      sites$clusteredValue <- sites$value
      sites$clusteredValue[sites$belowQuartile][pos.overlap[,"queryHits"]] <-
        pos.overlap$clusteredValue
      stopifnot(any(!is.na(sites$clusteredValue)))
    } else {
      message("No sites found below defined quartile. Try to increase ",
              "the quartile or use standard clustering, byQuartile=FALSE.")
    }
  } else {
    message("Clustering by minimum overlap.")

    sites <- split(sites, sites$grouping)

    sites <- bplapply(sites, function(x) {
      ## find overlapping positions using findOverlaps() using
      ## maxgap adjusted by windowSize!
      sites.gr <- with(x, GRanges(seqnames=posID2, IRanges(start=value, width=1),
                                  strand="*", freq))

      # the key part is drop.self=TRUE,drop.redundant=FALSE..
      # helps overwrite values at later step
      res <- as.data.frame(findOverlaps(sites.gr, drop.self=TRUE,
                                                  drop.redundant=FALSE,
                                                  select="all",
                                                  maxgap=windowSize-1L))
      if(nrow(res)>0) {
        # add accessory columns to dictate decision making!
        # q = query, s = subject, val = value, freq = frequency of query/subject
        res$q.val <- start(sites.gr)[res$queryHits]
        res$s.val <- start(sites.gr)[res$subjectHits]
        res$q.freq <- sites.gr$freq[res$queryHits]
        res$s.freq <- sites.gr$freq[res$subjectHits]
        res$dist <- with(res,abs(q.val-s.val))

        ## do safety checking!
        stopifnot(!any(res$dist>windowSize))

        # favor a lower value where frequence/cloneCount is tied,
        # else use the value of the highest frequency!
        res$val <- with(res, ifelse(q.freq==s.freq,
                                    ifelse(q.val < s.val, q.val, s.val),
                                    ifelse(q.freq >= s.freq, q.val, s.val)))

        # For cases where there are >1 matches between query & subject...
        # find the one with the highest frequency and merge with that.
        # If all frequencies are the same, then use the lowest
        # value to represent the cluster!
        res$maxFreq <- with(res, pmax(q.freq, s.freq))
        res$ismaxFreq <- as.logical(with(res, ave(maxFreq, queryHits,
                                                  FUN=function(x) x==max(x))))
        res$ismaxFreq <- as.logical(res$ismaxFreq)

        ## VIP step...this is what merges high value to low
        ## value for ties in the hash structure below!!!
        res <- dplyr::arrange(res, plyr::desc(queryHits), val)
        clustered <- unique(subset(res,ismaxFreq)[,c("queryHits","val")])
        clustered <- with(clustered, split(val, queryHits))

        ## make sure there is only one entry per hit this is useful in
        ## situations when multiple query & subject are off by 1bp
        ## i.e. queryHits: 1,2,3,4; subjectHits: 1,2,3,4;
        ## vals: 31895692 31895693 31895694 31895695
        clustered <- unlist(sapply(clustered, "[[", 1))

        # trickle results back to sites
        x$clusteredValue <- x$value
        x$clusteredValue[as.numeric(names(clustered))] <- as.numeric(clustered)
        rm("clustered","res")
        cleanit <- gc()
      } else {
        message("No locations found within ", windowSize, "bps for ",
                x$grouping[1], "...no clustering performed!")
        x$clusteredValue <- x$value
      }
      x
    }, BPPARAM=dp)
    sites <- plyr::rbind.fill(sites)
  }

  message("\t - Adding clustered value frequencies.")

  # get frequency of clusteredValue
  counts <- dplyr::select(sites, posID2, clusteredValue, freq) %>%
    plyr::count(., vars = c("posID2","clusteredValue"), wt_var = "freq")

  names(counts)[grep("freq",names(counts),fixed=TRUE)] <- "clusteredValue.freq"
  sites <- merge(sites,counts)

  if(byQuartile) {
    sites <- sites[,c("posID","value","freq","clusteredValue",
                      "clusteredValue.freq","grouping")]
  }

  sites$posID2<-NULL
  if(is.null(grouping)) { sites$grouping<-NULL }
  if(is.null(weight)) { sites$weight<-NULL }

  return(sites)
}
