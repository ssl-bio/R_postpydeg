#' Get max and min coverage from a BigWig file for specific region. Modified from http://adomingues.github.io/2016/11/13/max-coverage-in-bigwigs/
#' @name maxCoverageBigWig
#' @param bigwig A bigwig record. BigWig files should be loaded by `rtracklayer::import`.
#' @param chr Char of chromosome name
#' @param x_start Numeric. Coordinate of left end of the plot range.
#' @param x_end Numeric. Coordinate of right end of the plot range.
#' @importFrom GenomicRanges GRanges IRanges subsetByOverlaps
#' @export
maxCoverageBigWig <- function(bigwig, chr, x_start, x_end) {
  # GRanges for the target region
  gr <- GRanges(
    seqnames = chr,
    ranges = IRanges(x_start, x_end),
    strand = c("+", "-")
  )
  # Get BigWig records for the specified region
  ovlp <- subsetByOverlaps(bigwig, gr)
  if (length(ovlp) > 0) {
    max_cov <- max(ovlp$score)
    min_cov <- min(ovlp$score)
  } else {
    max_cov <- 0
    min_cov <- 0
  }
  return(c(min = min_cov, max = max_cov))
}


#' Get max and min coverage from a BigWig file for specific region.
#' @name maxCoverageBigWigList
#' @param bigwigs List of bigwig records. BigWig files should be loaded by `rtracklayer::import`.
#' @param chr Char of chromosome name
#' @param x_start Numeric. Coordinate of left end of the plot range.
#' @param x_end Numeric. Coordinate of right end of the plot range.
#' @export
maxCoverageBigWigList <- function(bigwigs, chr, x_start, x_end) {
    max_min <- lapply(bigwigs, maxCoverageBigWig, chr = chr,
                      x_start = x_start, x_end = x_end)
  max_min <- do.call(rbind, max_min)
  max_cov <- max(max_min[, "max"])
  min_cov <- min(max_min[, "min"])
  return(c(min = min_cov, max = max_cov))
}

#' Invert values and order of elements of a numeric vector while keeping names
#' @name mergeVector
#' @param ivec A string vector with names to construct the output
#' @param sep1 A string to separate the first half of vector elements
#' @param sep2 A string to separate the two halfs before exporting
#' @param parenthesis A pair of symbols to group the half of vector elements
#' @export
mergeVector <- function(ivec,sep1 = " vs ", sep2 = " AND \n",
                        parenthesis=c("(", ")"), single = FALSE) {
    i.sel <-rep(c(TRUE, FALSE),each = length(ivec) / 2)
    part1 <- paste(names(ivec)[i.sel],collapse = sep1)
    if (single) {
        part1 <- gsub(" \\[[0-9]\\]", "", part1)
        imerged <- part1
    } else {
        part2 <- paste(names(ivec)[!i.sel], collapse = sep1)
        imerged <- paste(paste0(parenthesis[1], part1, parenthesis[2]),
                         paste0(parenthesis[1], part2, parenthesis[2]),
                         sep = sep2)
    }
    return(imerged)
}

#' Invert values and order of elements of a numeric vector while keeping names
#' @name invertVector
#' @param x A numeric vector with names
#' @export
invertVector <- function(x) {
  inames <- names(x)
  x <- x[rev(seq_along(x))] * -1
  names(x) <- inames
  return(x)
}

#' Invert the score values of a list of a GenomicRange object
#' @name invertGRscores
#' @param bwl A list of GenomicRange-objects with a score metadata
#' @export
invertGRscores <- function(bwl) {
  for (i in seq_along(bwl)) {
    ibw <- bwl[[i]]
    score(ibw) <- score(ibw) * -1
    bwl[[i]] <- ibw
  }
  return(bwl)
}

#' Returns: Start and end coordinates for the plot
#' @param i_start ranscript start coordinates
#' @param i_end Transcript end coordinates
#' @param i_factor proportion of the transcript range that will be added as offset
#' @param p_width width of the plot in basepairs
#' @name GetPlotCoors
#' @export
GetPlotCoors <- function(i_start, i_end, i_factor = NULL, p_width = NULL) {
  if (i_start < i_end) {
    p_start <- i_start
    p_end <- i_end
  } else {
    p_start <- i_end
    p_end <- i_start
  }
  p_range <- p_end - p_start

  if (!is.null(p_width)) {
    p_offset <- ceiling((p_width - p_range) / 2)
  } else {
    p_offset <- ceiling((p_range * i_factor) / 10) * 10
    p_width <- ceiling(p_range + (p_offset * 2))
  }

  p_start <- i_start-p_offset
  p_end <- p_start+p_width
  results <- list(plot_start = p_start, plot_end = p_end)
  return(results)
}

#' Returns: Coordinates on the GenomeAxisTrack for the marks based on an intial number of segments (ticksn)
#' @param p_start ranscript start coordinates
#' @param p_end Transcript end coordinates
#' @param ticksn Number of segments to divide the transcript length
#' @param i_factor proportion of the transcript range that will be added as offset
#' @param p_width width of the plot in basepairs
#' @name GetPlotMarks
#' @export
GetPlotMarks <- function(i_start, i_end, p_start, p_end, iticks,
                         txPlot =  FALSE, min_ticks = 4, max_ticks = 6) {

    p_width <- abs(p_end - p_start)

    if (p_width < 50) {
        round_f <- 10
        idiv <- 1
    } else {
        round_f <- 100
        idiv <- 10
    }

    t_start <- round_f * ceiling(p_start / round_f)
    t_end <- round_f * floor(p_end / round_f)


    max_ticks <- ceiling(p_width / 10)
    if (max_ticks < iticks) { #Could it cause an infinite loop?
        iticks <- max_ticks
    }
    p_step <- floor((t_end - t_start)/(idiv*iticks))*idiv
    tick_coors <- seq(t_start, t_end, p_step)
    if (txPlot) {
        i.sel <- tick_coors >= i_start & tick_coors <= i_end
        tick_coors <- tick_coors[i.sel]
        if (length(tick_coors) < min_ticks) {
            while (length(tick_coors) < min_ticks) {
                iticks <- iticks + 1
                tick_coors <- GetPlotMarks(i_start, i_end,
                                           t_start, t_end, iticks,
                                           txPlot =  TRUE)
            }
        } else if (length(tick_coors) > max_ticks) {
            while (length(tick_coors) > max_ticks) {
                iticks <- iticks - 1
                tick_coors <- GetPlotMarks(i_start, i_end,
                                           t_start, t_end, iticks,
                                           txPlot =  TRUE)
            }
        }
    }
    ## cat("Final Number of tickmarks  = ",iticks,"\nSpace between marks",p_step,"\n")
    return(tick_coors)
}

#' Draws degradome profile comparing (a pair of) test samples (top panel) against (a pair of) control samples (bottom panel). Depending on the value of the 'plot_peak' parameter (boolean) it will draw a the profile around a narrow region around a peak or it will draw the profile along the transcript range. Peak is defined as the region detected as significantly different by pyDegradome and it is highlighted in yellow. Transcripts mapped on the - strand are reversed to facilitate interpretation
#' @name drawDplot
#' @param p_chr Char of chromosome name
#' @param tx_start Numeric. Coordinate of left end of the transcript.
#' @param tx_end Numeric. Coordinate of right end of the transcript.
#' @param peak_start Numeric vector.Coordinate(s) of left end of the region to be highlighted. If NULL, nothing to be highlinghted (default).
#' @param peak_end Numeric vector.Coordinate(s) of right end of the region to be highlighted. If NULL, nothing to be highlinghted (default).
#' @param tx_ID Gene identification ID, such as ensembl
#' @param sample_list_plot Character vector of sampel names to be plotted. Should be a part of `names(bigwigs)`.
#' @param p_test
#' @param p_ctrl
#' @param ylim Character vector, "all" or "strand". Both will set the scale of all the panels to be the same; the former will consider both strands while the latter will be strand specific.  Null for autoscale (default).
#' @param p_bigwigs List of bigwig records. BigWig files should be loaded by `rtracklayer::import`.
#' @param p_mart a biomart object from which to retrieve the gene information
#' @param p_colors Named vector of color code. This must have names identical to `names(bigwigs)`.
#' @param p_col_hl Charactor (vector) of color code to be used in highlighting.
#' @param ticksn Number of tickmarks drawn on the genomeaxistrack
#' @param p_size Size of plot symbols
#' @param p_alpha Alpha value of plot symbols
#' @import Gviz biomaRt rtracklayer
#' @param ucscnames Logical, if false ucsc enforcing of chromosome names is disabled
#' @param plot_peak Boolean, true: profile is drawn around the peak region (~20nt range); false profile is  drawn in the range of the transcript
#' @param p_font postScript option for the value of the family font
#' @param p_base_col A vector o color values (hex) for the bases. Values have to have a base-name associated and it should take into account whether is RNA or DNA
#' @param ref_genome_seq A DNAStringSet object with all chromosomes for the species under study
#' @param ref_genome A string representing the genome on which the track's ranges are defined (e.g. "TAIR10")
#' @export
#' @importFrom Gviz BiomartGeneRegionTrack DataTrack displayPars GenomeAxisTrack
#' @importFrom Gviz HighlightTrack OverlayTrack plotTracks RNASequenceTrack
drawDplot <- function(p_chr, tx_start, tx_end, peak_start = NULL,
                      peak_end = NULL, tx_ID = NULL, sample_list_plot,
                      p_test = NULL, p_ctrl = NULL, p_strand = NULL,
                      ylim = "strand", p_bigwigs, p_mart = mart,
                      p_colors = "blue", p_col_hl = "#FF000033",
                      ticksn = NULL,  i_factor = NULL, p_width = NULL,
                      p_size = 1, p_alpha = 0.8, ucscnames = TRUE,
                      plot_peak = FALSE, p_font = NULL, p_base_col = NULL,
                      ref_genome_seq = NULL, ref_genome = NULL) {
  # Check if names of color_cond match those of the samples
  if (!all(unlist(lapply(sample_list_plot, is.element, names(p_colors))))) {
    warning("names(p_colors) does not match with sample_list_plot. p_colors will be recycled to cover all samples")
    p_colors <- structure(rep(p_colors, length(sample_list_plot)),
      names = sample_list_plot
    )
  }

  ##Use arbitrary names
    if (!ucscnames) {
        options(ucscChromosomeNames=FALSE)
    }
    if (plot_peak) {
        txPlot <- FALSE
    } else {
        txPlot <- TRUE
    }


    p_coors <- GetPlotCoors(tx_start, tx_end, i_factor, p_width)
    ticks_coors <- GetPlotMarks(i_start = tx_start, i_end = tx_end, p_start = p_coors$plot_start, p_end = p_coors$plot_end, iticks = ticksn, txPlot = txPlot)

  ## Calculate limits; If by p_strand,
  if (ylim == "all") {
    p_ylim <- maxCoverageBigWigList(p_bigwigs, p_chr,
                                  p_coors$plot_start, p_coors$plot_end)
  } else if (ylim == "strand") {
    p_ylim <- maxCoverageBigWigList(p_bigwigs, p_chr,
                                  p_coors$plot_start, p_coors$plot_end)
    if (p_strand == "+") {
        p_ylim[1] <- 0 #Set min value to 0
        test53 <- TRUE
        test35 <- FALSE
        strandD <- FALSE
    } else {
        p_ylim[2] <- 0 #Set max value to 0
        test53 <- FALSE
        test35 <- TRUE
        strandD <- TRUE
        ##Invert values of bigwig data and Y-axis range
        p_bigwigs <- invertGRscores(p_bigwigs)
        p_ylim <- invertVector(p_ylim)
    }
  }
  ## Test if plot around peak is chosen
  if (plot_peak) {
    Sp_fasta_names <- names(ref_genome_seq)
    Sp_chr_names <- sapply(strsplit(Sp_fasta_names, " "),
                           function(x) {
                             x[1]})
    i.sel <- Sp_chr_names %in% p_chr
    ch.seq <- ref_genome_seq[names(ref_genome_seq)[i.sel]]
    names(ch.seq) <- p_chr
    if (p_strand == "+") {
      sTrack <- RNASequenceTrack(#
        ch.seq,
        genome = ref_genome, chromosome = p_chr,
        fontcolor = p_base_col, cex = 0.85, noLetters = FALSE
      )
    } else {
      sTrack <- RNASequenceTrack(#
        ch.seq,
        genome = ref_genome, chromosome = p_chr,
        fontcolor = p_base_col, cex = 0.85, complement = TRUE,
        noLetters = FALSE
      )
    }
    displayPars(sTrack) <- list(fontfamily = p_font, fontfamily.title = p_font)
  }

    tryCatch(
    {
        fm <- Gviz:::.getBMFeatureMap()
        if (!plot_peak) {
            bm <- BiomartGeneRegionTrack(
                chromosome = p_chr, genome = ref_genome,
                start = p_coors$plot_start, end = p_coors$plot_end, biomart = p_mart, size = 1.2,
                name = "Gene model", col.title = "#808080ff", utr5 = "gray",
                utr3 = "gray", protein_coding = "gray30", col.line = NULL,
                cex = 2, filters = list(ensembl_transcript_id = tx_ID),
                cex.title = 0.7, featureMap = fm
            )
        }
    },
    error = function(e) {
      message(paste0(
        "plotBigWig failed to get annotation through biomaRt for: chr = ",
        p_chr, " x_start = ", tx_start, " x_end = ", tx_end, " peak_start = ",
        peak_start, " peak_end = ", peak_end
      ))
      message(e)
      ## bm <- NULL
    }
  )
  AT <- GenomeAxisTrack()
  sample_track <- list()
  for (i.counter in seq_along(sample_list_plot)) {
    i.sample <- sample_list_plot[[i.counter]]
    bw <- p_bigwigs[[i.sample]]
    ## filter bigwig by strand
    sel.strand <- strand(bw) == p_strand
    bw_strand <- bw[sel.strand, ]
    i.lab <- gsub(" \\[[0-9]\\]", "", sample_list_plot[i.counter])
    i.lab <- gsub("At ", "", i.lab)
    sample_track[[i.sample]] <- DataTrack(bw_strand,
                                          chromosome = p_chr,
                                          strand = p_strand,
                                          genome = ref_genome,
                                          ylim = p_ylim,
                                          col = p_colors[i.sample],
                                          col.axis = "white",
                                          background.panel = "#e0e0e07f",
                                          name = i.lab,
                                          fontsize = 6,
                                        # symbols
                                          cex = p_size,
                                          alpha = p_alpha,
                                          alpha.title = 1,
                                          type = c("p", "g"),
                                          legend = TRUE
                                          )
  }

  ot1 <- OverlayTrack(sample_track[sample_list_plot %in% p_test],
                      background.title = "gray30",
                      col.grid = "white",
                      ## fontsize = 14 or 6,
                      cex.axis = 0.7,
                      ## cex.title = 0.8 or 0.7,
                      alpha = 0.75
                      )

  ot2 <- OverlayTrack(sample_track[sample_list_plot %in% p_ctrl],
                      background.title = "gray30",
                      col.grid = "white",
                      ## fontsize = 14 or 6,
                      cex.axis = 0.7,
                      ## cex.title = 0.8 or 0.7,
                      alpha = 0.75
                      )

  if (is.numeric(peak_start) && is.numeric(peak_end)) {
    if (plot_peak) {
      ht <- HighlightTrack(
        trackList = c(ot1, ot2, sTrack), start = peak_start, end = peak_end,
        from = p_coors$plot_start, to = p_coors$plot_end,
        genome = ref_genome, chromosome = p_chr,
        col = "transparent", inBackground = FALSE, fill = p_col_hl, alpha = 0.3
      )
    } else {
      ht <- HighlightTrack(
        trackList = c(
          ot1, ot2, bm, AT
        ), start = peak_start, end = peak_end, chromosome = p_chr,
        col = "transparent", inBackground = FALSE, fill = p_col_hl,
        alpha = 0.68
      )
    }
  } else {
    ht <- NULL
  }
  tryCatch(
    {
      if (!is.null(ht)) {
        if (plot_peak) {
          plotTracks(c(ht, AT),
                     from = p_coors$plot_start, to = p_coors$plot_end,
                     add53 = test53, add35 = test35, labelPos = "below",
                     transcriptAnnotation = "transcript_id", fontsize = 11
                     )
          gc()
        } else {
          plotTracks(c(ht),
                     from = p_coors$plot_start,
                     to = p_coors$plot_end, add53 = test53,
                     add35 = test35, labelPos = "below",
                     transcriptAnnotation = "transcript_id",
                     ticksAt = ticks_coors,
                     fontsize = 14, reverseStrand = strandD)
          gc()

        }
      } else {
        if (plot_peak) {
        plotTracks(c(ht, sTrack, AT),
          from = tx_start, to = tx_end,
          add53 = test53, add35 = test35, labelPos = "below",
          transcriptAnnotation = "transcript_id", fontsize = 11
          )
        gc()
        } else {
          plotTracks(c(ht, bm, AT),
                       from = p_coors$plot_start - 500,
                       to = p_coors$plot_end + 500,
                       add53 = test53, add35 = test35, labelPos = "below",
                       transcriptAnnotation = "transcript_id",
                       ticksAt = ticks_coors, exponent = p_exp,
                       fontsize = 14, reverseStrand = strandD
                     )
          gc()
        }
      }
    },
    error = function(e) {
      message(paste0(
        "plotBigWig failed to plot profile for: chr = ",
        p_chr, " x_start = ", tx_start, " x_end = ", tx_end, " peak_start = ",
        peak_start, " peak_end = ", peak_end
      ))
      message(e)
    }
  )
}
