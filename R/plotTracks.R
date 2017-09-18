#' @import Gviz rtracklayer
plotTracks <- function(mygene = NULL, genemodels = NULL, genemodels2 = NULL,
                       gtf_file = NULL, rnaseq_datafiles = NULL,
                       rnaseq_condition = NULL, show_chr = NULL,
                       min_coord = NULL, max_coord = NULL,
                       pdf_filename = NULL, pdf_width = 7, pdf_height = 7, ...) {
  options(ucscChromosomeNames = FALSE, envir = .GlobalEnv)

  if (is.null(genemodels) & is.null(gtf_file))
    stop("Either genemodels or gtf_file must be provided")
  if (is.null(genemodels))
    genemodels <- create_genemodels(gtf_file)
  if (!is.null(genemodels))
    stopifnot(all(c("transcript", "gene", "exon") %in% colnames(mcols(genemodels))))

  ## Create gene region track
  if (!is.null(mygene)) {
    gm <- subset(genemodels, tolower(gene) == tolower(mygene) | tolower(gene_name) == tolower(mygene))
    gm <- subset(gm, gene == gene[1])  ## Select only one gene if there are many with the same name
    id <- unique(gm$gene_name)
    idshow <- paste0(id, " (", unique(gm$gene), ")")
    show_chr <- unique(seqnames(gm))[1]
    gm <- subset(gm, seqnames == show_chr)
    min_coord <- min(start(gm)) - 0.15*(max(end(gm)) - min(start(gm)))
    max_coord <- max(end(gm)) + 0.05*(max(end(gm)) - min(start(gm)))

    ## Other features in the considered region
    gmo <- genemodels[overlapsAny(genemodels,
                                  GRanges(seqnames = show_chr,
                                          ranges = IRanges(start = min_coord,
                                                           end = max_coord),
                                          strand = "*"))]
    gmo <- gmo[!(gmo %in% gm)]

    ## Excluded features
    if (!is.null(genemodels2)) {
      gm2 <- genemodels2[overlapsAny(genemodels2,
                                     GRanges(seqnames = show_chr,
                                             ranges = IRanges(start = min_coord,
                                                              end = max_coord),
                                             strand = "*")), ]
      ## Excluded features from other genes will be part of the "other genes" track
      gm2other <- subset(gm2, gene != gm$gene[1])
      mcols(gm2other) <- mcols(gm2other)[, match(colnames(mcols(gmo)), colnames(mcols(gm2other)))]
      gmo <- c(gmo, gm2other)
      ## Keep only excluded features from the current gene in this track
      gm2 <- subset(gm2, gene == gm$gene[1])
    } else {
      gm2 <- GRanges()
    }
  } else {
    if (any(is.null(c(show_chr, min_coord, max_coord)))) {
      stop("Either mygene or genomic region must be provided")
    } else {
      gm <- genemodels[overlapsAny(genemodels,
                                   GRanges(seqnames = show_chr,
                                           ranges = IRanges(start = min_coord,
                                                            end = max_coord),
                                           strand = "*")), ]
    }
  }
  grtr <- GeneRegionTrack(gm, showId = TRUE, col = NULL, fill = "blue",
                          name = ifelse(!is.null(id), id, "Genes"),
                          col.title = "black")
  grtr2 <- GeneRegionTrack(gmo, showId = TRUE, col = NULL, fill = "green",
                           name = "",
                           col.title = "black")
  grtr3 <- GeneRegionTrack(gm2, showId = TRUE, col = NULL, fill = "orange",
                           name = "",
                           col.title = "black")

  ## Create genome axis track
  gtr <- GenomeAxisTrack()

  if (length(gm) > 0) {
    tracks <- c(gtr, grtr, grtr3, grtr2)
  } else {
    tracks <- c(gtr)
  }

  legend_added <- FALSE
  threecols <- c("#DC050C", "#7BAFDE", "#B17BA6", "#F1932D", "black")
  # twocols <- c(rgb(11, 102, 254, maxColorValue = 255),
  #              rgb(250, 0, 255, maxColorValue = 255))

  ## RNAseq tracks
  if (!is.null(rnaseq_datafiles)) {
    if (is.null(names(rnaseq_datafiles)))
      stop("RNAseq input file list must be named")
    multiTracks_rnaseq <- lapply(1:length(rnaseq_datafiles), function(i) {
      assign(paste0("rnaseqtr", i),
             DataTrack(range = rnaseq_datafiles[i],
                       type = "histogram",
                       name = paste0(rnaseq_condition[names(rnaseq_datafiles)[i]], "_",
                                     gsub("20170530.", "",
                                          gsub("_trimmed_Aligned.sortedByCoord.out", "",
                                               names(rnaseq_datafiles)[i]))),
                       chromosome = unique(seqnames(gm)),
                       col.title = "black",
                       fill = threecols[(as.numeric(as.factor(rnaseq_condition)))[i]],
                       col = threecols[(as.numeric(as.factor(rnaseq_condition)))[i]],
                       col.histogram = threecols[(as.numeric(as.factor(rnaseq_condition)))[i]],
                       fill.histogram = threecols[(as.numeric(as.factor(rnaseq_condition)))[i]]
             ))
    })
    tracks <- c(multiTracks_rnaseq, tracks)
  }

  ## Plot tracks
  if (!is.null(pdf_filename))
    pdf(pdf_filename, width = pdf_width, height = pdf_height)

  plotTracks(tracks, chromosome = show_chr,
             from = min_coord, to = max_coord,
             main = ifelse(!is.null(id), idshow, ""),
             min.width = 0, min.distance = 0, collapse = FALSE, ...)

  if (!is.null(pdf_filename))
    dev.off()
}
