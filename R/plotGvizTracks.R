#' Generate gene model GRanges object from GTF file
#'
#' @author Charlotte Soneson
#'
#' @param gtfFile GTF file (in Ensembl format)
#'
#' @importFrom rtracklayer import
#' @importFrom S4Vectors mcols
#' @importFrom BiocGenerics subset
#'
#' @export
#'
createGenemodels <- function(gtfFile) {
  genemodels <- rtracklayer::import(gtfFile)
  idx <- match(c("transcript_id", "gene_id", "exon_id"),
               colnames(S4Vectors::mcols(genemodels)))
  colnames(S4Vectors::mcols(genemodels))[idx] <- c("transcript", "gene", "exon")
  S4Vectors::mcols(genemodels)$symbol <- S4Vectors::mcols(genemodels)$transcript
  BiocGenerics::subset(genemodels, type == "exon")
}

#' Plot coverage and gene model tracks
#'
#' @author Charlotte Soneson
#'
#' @param showGene Gene name or symbol for the gene to show in the plot.
#' @param geneModels GRanges object with gene models, typically generated from a
#'   gtf file (e.g. with the \code{createGenemodels} function).
#' @param geneModels2 Second GRanges object with gene models, e.g. those that
#'   were excluded from \code{geneModels}.
#' @param gtfFile GTF file from which \code{geneModels} can be generated if it
#'   is not explicitly provided.
#' @param bwFiles Named vector with paths to bigWig files. These will be used to
#'   generate coverage plots.
#' @param bwCond Named vector corresponding to \code{bwFiles}, giving the group
#'   label for each bigWig file. These will be used to color the coverage plots
#'   by group.
#' @param showChr Name of the chromosome to show in the plot, if \code{showGene}
#'   is not specified.
#' @param minCoord,maxCoord Minimum and maximum genomic coordinate to show, if
#'   \code{showGene} is not specified.
#' @param pdfFilename File name of pdf file to which the plot should be saved.
#'   If set to \code{NULL}, the plot is written to the current graphics device.
#' @param pdfWidth,pdfHeight Width and height of the pdf file to which the plot
#'   should be saved.
#' @param ... Additional arguments (not currently used)
#'
#' @importFrom Gviz GeneRegionTrack GenomeAxisTrack DataTrack plotTracks
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges overlapsAny IRanges
#' @importFrom grDevices dev.off pdf
#' @importFrom GenomeInfoDb seqnames
#' @importFrom S4Vectors mcols
#' @importFrom BiocGenerics start end subset setdiff
#'
#' @export
#'
plotGvizTracks <- function(showGene = NULL, geneModels = NULL, geneModels2 = NULL,
                           gtfFile = NULL, bwFiles = NULL,
                           bwCond = NULL, showChr = NULL,
                           minCoord = NULL, maxCoord = NULL,
                           pdfFilename = NULL, pdfWidth = 7, pdfHeight = 7, ...) {
  options(ucscChromosomeNames = FALSE, envir = .GlobalEnv)

  if (is.null(geneModels) && is.null(gtfFile))
    stop("Either geneModels or gtfFile must be provided")
  if (is.null(geneModels))
    geneModels <- createGenemodels(gtfFile)
  if (!is.null(geneModels))
    stopifnot(all(c("transcript", "gene", "exon") %in%
                    colnames(S4Vectors::mcols(geneModels))))

  ## Create gene region track
  if (!is.null(showGene)) {
    gm <- BiocGenerics::subset(geneModels, tolower(gene) == tolower(showGene) |
                                 tolower(gene_name) == tolower(showGene))
    gm <- BiocGenerics::subset(gm, gene == gene[1])  ## Select only one gene if there are many with the same name
    id <- unique(gm$gene_name)
    idshow <- paste0(id, " (", unique(gm$gene), ")")
    showChr <- unique(GenomeInfoDb::seqnames(gm))[1]
    gm <- BiocGenerics::subset(gm, seqnames == showChr)
    minCoord <- min(BiocGenerics::start(gm)) - 0.15*(max(BiocGenerics::end(gm)) -
                                                       min(BiocGenerics::start(gm)))
    maxCoord <- max(BiocGenerics::end(gm)) + 0.05*(max(BiocGenerics::end(gm)) -
                                                     min(BiocGenerics::start(gm)))

    ## Other features in the considered region
    gmo <- geneModels[IRanges::overlapsAny(
      geneModels,
      GenomicRanges::GRanges(seqnames = showChr,
                             ranges = IRanges::IRanges(start = minCoord,
                                                       end = maxCoord),
                             strand = "*"))]
    # gmo <- gmo[!(gmo %in% gm)]
    # gmo <- gmo[!(IRanges::overlapsAny(gmo, gm, type = "equal"))]
    gmo <- BiocGenerics::subset(gmo, gene != gm$gene[1])

    ## Additional track (e.g., excluded features)
    if (!is.null(geneModels2)) {
      gm2 <- geneModels2[IRanges::overlapsAny(
        geneModels2,
        GenomicRanges::GRanges(seqnames = showChr,
                               ranges = IRanges::IRanges(start = minCoord,
                                                         end = maxCoord),
                               strand = "*")), ]
      ## Excluded features from other genes will be part of the "other genes" track
      gm2other <- BiocGenerics::subset(gm2, gene != gm$gene[1])
      S4Vectors::mcols(gm2other) <-
        S4Vectors::mcols(gm2other)[, match(colnames(S4Vectors::mcols(gmo)),
                                           colnames(S4Vectors::mcols(gm2other)))]
      gmo <- c(gmo, gm2other)
      ## Keep only excluded features from the current gene in this track
      gm2 <- BiocGenerics::subset(gm2, gene == gm$gene[1])
    } else {
      gm2 <- GenomicRanges::GRanges()
    }
  } else {
    if (any(is.null(c(showChr, minCoord, maxCoord)))) {
      stop("Either showGene or genomic region must be provided")
    } else {
      gm <- geneModels[IRanges::overlapsAny(
        geneModels,
        GenomicRanges::GRanges(seqnames = showChr,
                               ranges = IRanges::IRanges(start = minCoord,
                                                         end = maxCoord),
                               strand = "*")), ]
    }
  }
  grtr <- Gviz::GeneRegionTrack(gm, showId = TRUE, col = NULL, fill = "blue",
                                name = ifelse(!is.null(id), id, "Genes"),
                                col.title = "black")
  grtr2 <- Gviz::GeneRegionTrack(gmo, showId = TRUE, col = NULL, fill = "green",
                                 name = "",
                                 col.title = "black")
  grtr3 <- Gviz::GeneRegionTrack(gm2, showId = TRUE, col = NULL, fill = "orange",
                                 name = "",
                                 col.title = "black")

  ## Create genome axis track
  gtr <- Gviz::GenomeAxisTrack()

  if (length(gm) > 0) {
    tracks <- c(gtr, grtr, grtr3, grtr2)
  } else {
    tracks <- c(gtr)
  }

  legend_added <- FALSE
  allcols <- c("#DC050C", "#E8601C", "#7BAFDE", "#1965B0", "#B17BA6",
               "#882E72", "#F1932D", "#F6C141", "#F7EE55", "#4EB265",
               "#90C987", "#CAEDAB", "#777777")[c(1, 3, 5, 10, 13, 8, 4, 6, 12, 2, 7, 9, 11)]

  ## RNAseq coverage tracks
  if (!is.null(bwFiles)) {
    if (is.null(names(bwFiles)))
      stop("RNAseq input file list must be named")
    multiTracks_rnaseq <- lapply(1:length(bwFiles), function(i) {
      assign(paste0("rnaseqtr", i),
             Gviz::DataTrack(range = bwFiles[i],
                             type = "histogram",
                             name = names(bwFiles)[i],
                             chromosome = unique(GenomeInfoDb::seqnames(gm)),
                             col.title = "black",
                             fill = allcols[(as.numeric(as.factor(bwCond)))[i]],
                             col = allcols[(as.numeric(as.factor(bwCond)))[i]],
                             col.histogram = allcols[(as.numeric(as.factor(bwCond)))[i]],
                             fill.histogram = allcols[(as.numeric(as.factor(bwCond)))[i]]
             ))
    })
    tracks <- c(multiTracks_rnaseq, tracks)
  }

  ## Plot tracks
  if (!is.null(pdfFilename))
    pdf(pdfFilename, width = pdfWidth, height = pdfHeight)

  Gviz::plotTracks(tracks, chromosome = showChr,
                   from = minCoord, to = maxCoord,
                   main = ifelse(!is.null(id), idshow, ""),
                   min.width = 0, min.distance = 0, collapse = FALSE, ...)

  if (!is.null(pdfFilename))
    dev.off()
}
