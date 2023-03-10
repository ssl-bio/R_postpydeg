% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/profile_visualization.R
\name{drawDplot}
\alias{drawDplot}
\title{Draws degradome profile comparing (a pair of) test samples (top panel) against (a pair of) control samples (bottom panel). Depending on the value of the 'plot_peak' parameter (boolean) it will draw a the profile around a narrow region around a peak or it will draw the profile along the transcript range. Peak is defined as the region detected as significantly different by pyDegradome and it is highlighted in yellow. Transcripts mapped on the - strand are reversed to facilitate interpretation}
\usage{
drawDplot(
  p_chr,
  tx_start,
  tx_end,
  peak_start = NULL,
  peak_end = NULL,
  tx_ID = NULL,
  sample_list_plot,
  p_test = NULL,
  p_ctrl = NULL,
  p_strand = NULL,
  ylim = "strand",
  p_bigwigs,
  p_mart = mart,
  p_colors = "blue",
  p_col_hl = "#FF000033",
  ticksn = NULL,
  i_factor = NULL,
  p_width = NULL,
  p_size = 1,
  p_alpha = 0.8,
  ucscnames = TRUE,
  plot_peak = FALSE,
  p_font = NULL,
  p_base_col = NULL,
  ref_genome_seq = NULL,
  ref_genome = NULL
)
}
\arguments{
\item{p_chr}{Char of chromosome name}

\item{tx_start}{Numeric. Coordinate of left end of the transcript.}

\item{tx_end}{Numeric. Coordinate of right end of the transcript.}

\item{peak_start}{Numeric vector.Coordinate(s) of left end of the region to be highlighted. If NULL, nothing to be highlinghted (default).}

\item{peak_end}{Numeric vector.Coordinate(s) of right end of the region to be highlighted. If NULL, nothing to be highlinghted (default).}

\item{tx_ID}{Gene identification ID, such as ensembl}

\item{sample_list_plot}{Character vector of sampel names to be plotted. Should be a part of `names(bigwigs)`.}

\item{p_test}{Plot names of test samples}

\item{p_ctrl}{Plot names of control samples}

\item{ylim}{Character vector, "all" or "strand". Both will set the scale of all the panels to be the same; the former will consider both strands while the latter will be strand specific.  Null for autoscale (default).}

\item{p_bigwigs}{List of bigwig records. BigWig files should be loaded by `rtracklayer::import`.}

\item{p_mart}{a biomart object from which to retrieve the gene information}

\item{p_colors}{Named vector of color code. This must have names identical to `names(bigwigs)`.}

\item{p_col_hl}{Charactor (vector) of color code to be used in highlighting.}

\item{ticksn}{Number of tickmarks drawn on the genomeaxistrack}

\item{p_size}{Size of plot symbols}

\item{p_alpha}{Alpha value of plot symbols}

\item{ucscnames}{Logical, if false ucsc enforcing of chromosome names is disabled}

\item{plot_peak}{Boolean, true: profile is drawn around the peak region (~20nt range); false profile is  drawn in the range of the transcript}

\item{p_font}{postScript option for the value of the family font}

\item{p_base_col}{A vector o color values (hex) for the bases. Values have to have a base-name associated and it should take into account whether is RNA or DNA}

\item{ref_genome_seq}{A DNAStringSet object with all chromosomes for the species under study}

\item{ref_genome}{A string representing the genome on which the track's ranges are defined (e.g. "TAIR10")}
}
\description{
Draws degradome profile comparing (a pair of) test samples (top panel) against (a pair of) control samples (bottom panel). Depending on the value of the 'plot_peak' parameter (boolean) it will draw a the profile around a narrow region around a peak or it will draw the profile along the transcript range. Peak is defined as the region detected as significantly different by pyDegradome and it is highlighted in yellow. Transcripts mapped on the - strand are reversed to facilitate interpretation
}
