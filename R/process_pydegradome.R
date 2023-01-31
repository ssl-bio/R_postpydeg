#'Merges a column (ref_col) from a reference dataframe (i.reads) into an input dataframe (df). An 'indx' is used for merging. The input is divided into matching and non-matching dataframes and returned as a list.
#' @name MatchDFwithRef
#' @param df data.frame to match against a reference
#' @param i.reads reference dataframe
#' @param ref_col Column with values to merge
MatchDFwithRef <- function(df, i.reads, ref_col) {
    ilog <- df$indx %in% i.reads$indx
    df_match <- df[ilog, ]
    df_nonmatch <- df[!ilog, ]
    if (!is.null(df_match) && nrow(df_match) > 0) {
        df_match <- merge(df_match,
                          i.reads[, c("indx", ref_col)],
                          by = "indx",
                          all  = FALSE)[, -1]
    } else {
        df_match <- NULL
    }
    if (is.null(df_nonmatch) || nrow(df_nonmatch) == 0) {
        df_nonmatch <- NULL
    }
    results <- list(df_match = df_match,
                    df_nonmatch = df_nonmatch)
    return(results)
}

#' Merges a column (ref_col) from dataframe B (dfB) into dataframe A (dfA) keeping all entries from the latter. Returns a dataframe
#' @name MergeDFs
#' @param dfA Dataframe to merge into
#' @param dfB Dataframe to merge from
#' @param ref_col Column name (string) to merge
#' @importFrom dplyr dplyr::select
#' @export
MergeDFs <- function(dfA, dfB, ref_col) {
                                        #Create a indx vector to merge
    indx1 <- paste(dfA$tx_name, dfA$chr, dfA$strand, sep = "")
    dfA$indx_merge <- indx1

    indx2 <- paste(dfB$ID, dfB$chr, dfB$strand, sep = "")
    dfB$indx_merge <- indx2

    ref_cols <- c("indx_merge", ref_col)
    dfB1 <- dplyr::select(dfB, all_of(ref_cols))
    dfA_merged <- merge(dfA,  dfB1, by = "indx_merge",
                 all.x = TRUE)[, -1]
    return(dfA_merged)
}

#' Takes the name of a pydegradome comparison ('t_<comparison_1>_c_<comparison_2>') and returns the first element (<comparison_1>). 
#' @name getFirstComparison
#' @param icomp filename based on the comparison between test and control samples
#' @export
getFirstComparison <- function(i.comp) {
    sample_name <- gsubfn::strapplyc(i.comp, "t_(.*)_c_.*")[[1]]
    if (!sample_name %in% names(dg_bigwig_all)) {
        sample_name <- tolower(sample_name)
    }
    return(sample_name)
}

#' Reads a bigwig file and returns a vector with the highest read value over a set of genomic ranges (one or two GRanges objects).
#' @name getMaxRead
#' @param gr_A main GRanges object
#' @param gr_B secondary GRanges object
#' @param f_bigwig reference to a bigwig file
#' @param core Number of threads for the parallel process
#' @importFrom doParallel makeCluster registerDoParallel stopCluster
#' @importFrom foreach foreach
#' @importFrom GenomicRanges subsetByOverlaps
getMaxRead <- function(grA, f_bigwig, grB = NULL, core = 1) {
    cl <- makeCluster(as.numeric(core))
    registerDoParallel(cl)
    on.exit(stopCluster(cl))

    maxread <- foreach(i = seq_along(grA),
                       .packages = c("GenomicRanges"),
                       .combine = "c",
                       .inorder = TRUE) %dopar% {
                           if (is.null(grB)) {
                               scores <- subsetByOverlaps(f_bigwig, grA[i])$score
                           } else {
                               scores <- c(subsetByOverlaps(f_bigwig, grA[i])$score,
                                           subsetByOverlaps(f_bigwig, grB[i])$score)
                           }
                           if (length(scores) > 0) {
                               if (decode(strand(grA[i])) == "+") {
                                   maxread <- max(scores)
                               } else if (decode(strand(grA[i])) == "-") {
                                   maxread <- -min(scores)
                               }
                               maxread
                           } else {
                               ##The selected genomic region has no coverage value
                               ## in the BigWig
                               ## Coverage value is arbitrary set to Zero.
                               maxread <- 0
                           }
                       }
    return(maxread)
}

#' Takes a dataframe with coordinates for a peak region and calculates
#' from a bigwig file the highest read value for the peak.
#' Results are appended to the input dataframe
#' @name getMaxPeak
#' @param df data.frame of pyDegradome result
#' @param i.comp Name of comparison pair (divided into test and ctrl)
#' @param i.cols vector of columns to subset the table to append to the reference file
#' @param idir path to the directory with the reference file
#' @export
getMaxPeak <- function(df, i.comp, i.cols = cols.indx,
                         idir = pydeg_reads_dir) {
    peaks_ref <- file.path(idir,
                           paste0("MaxPeak_reads-",
                                  i.comp))
    cat("Calculating signals in the peak region (max_peak)...\n")
    cat("Processing ", paste0("MaxPeak_reads-", i.comp), "\n")
    if (file.exists(peaks_ref)) {

        cat("\tUsing reference file...\n")
        i.reads <- read.table(peaks_ref, header = TRUE, sep = "\t")
        i.reads <- unique(i.reads)

        ##Subset entries with ID
        df.ID <- df[grepl("^AT", df$ID), ]
        df.ID <- within(df.ID,
                        indx <- paste(ID,
                                      chr,strand,
                                      peak_start,
                                      peak_stop,
                                      sep = ""))
        dfs.ID <- MatchDFwithRef(df.ID,  i.reads, "max_peak")
        df.ID1 <- unique(dfs.ID$df_match)

                                        #Add max_peak to the missing entries
        cat("\t\tAdd max_peak to missing entries...\n")
        if (!is.null(dfs.ID$df_nonmatch) && nrow(dfs.ID$df_nonmatch) > 0) {
            df.ID2 <- unique(dfs.ID$df_nonmatch)
            gr <- GRanges(
                seqnames = df.ID2$chr,
                ranges = IRanges(start = df.ID2$peak_start,
                                 end = df.ID2$peak_stop),
                strand = df.ID2$strand
            )
            sample_name <- getFirstComparison(i.comp)
            bigwig <- dg_bigwig_all[[sample_name]]
            df.ID2$max_peak <- getMaxRead(grA = gr,
                                          f_bigwig = bigwig,
                                          core = env$cores)
            ##---------------------------------------------
                                        #Append the new entries
                                        #to the reference file
            i.reads2 <- dplyr::select(df.ID2,
                                      c(all_of(i.cols), "max_peak"))
            i.reads2 <- within(i.reads2,
                               indx <- paste(ID, chr,
                                             strand,
                                             peak_start,
                                             peak_stop,
                                             sep = ""))
            i.reads <- rbind(i.reads,i.reads2)
            unique(i.reads)
            i.reads <- i.reads[with(i.reads, order(indx)), ]
            
            tryCatch({
                                        #remove indx column
                df.ID2 <- dplyr::select(df.ID2,-indx)
                
                fwrite(i.reads, peaks_ref,sep = "\t", row.names=FALSE)
                cat(paste0("\t\t\tWrote ", basename(peaks_ref), "\n"))
            }, error = function(e) {
                stopCluster(cl)
                message("Writing peaks failed")
                message(e)
            }) #Trycatch
            ##---------------------------------------------
        } else {
            df.ID2 <- NULL
        }
        cat("\t\tFinished adding max_np_gene to missing entries\n")

                                        #Subset entries w/o ID and
                                        #assing max_peak ""
        cat("\t\tWorking on entries without ID...\n")
        df.NonID <- df[!grepl("^AT", df$ID), ]
        if (nrow(df.NonID) > 0) {
            df.NonID$max_peak <- NA
        } else {
            df.NonID <- NULL
        }

                                        #Merge data
        cat("\t\tMerging dataframes...\n")
        df <- rbind(df.ID1, df.ID2, df.NonID)
        df <- unique(df)
    } else {
        cat("\tNo reference was found.\n")
        cat("\t\tCalculating max_peak...\n")
        gr <- GRanges(
            seqnames = df$chr,
            ranges = IRanges(start = df$peak_start, end = df$peak_stop),
            strand = df$strand
        )
        sample_name <- getFirstComparison(i.comp)
        bigwig <- dg_bigwig_all[[sample_name]]
        
        df$max_peak <- getMaxRead(grA = gr,
                                  f_bigwig = bigwig,
                                  core = env$cores)
        ##---------------------------------------------
                                        #Write data to ref file
        i.reads <- dplyr::select(df,
                                 c(all_of(i.cols), "max_peak"))
        i.reads <- unique(i.reads)
        i.reads <- within(i.reads,
                          indx <- paste(ID, chr,
                                        strand, peak_start,
                                        peak_stop,
                                        sep = ""))
        i.reads <- i.reads[with(i.reads, order(indx)), ]
        
        tryCatch({
            fwrite(i.reads, peaks_ref, sep = "\t", row.names = FALSE)
            cat(paste0("\t\t\tWrote ", basename(peaks_ref), "\n"))
        }, error = function(e) {
            stopCluster(cl)
            message("Writing peaks failed")
            message(e)
        }) #Trycatch
        ##---------------------------------------------
    } #Checking if reference file exists
    return(df)
    cat("Finished calculating in peak region (max_peak)\n")
}

#' Returns a list of coordinates for non peak regions.
#' @name getNonPeakCoors
#' @param df A dataframe with peak and gene coordinates
getNonPeakCoors <- function(df) {
    ##Sort by peak_start
    df <- df[order(df$peak_start), ]
    for (i in seq_len(nrow(df))) {
        if (i == 1) {
                                        #Range 1
            tmp_start <- df$gene_region_start[i]
            tmp_stop <- df$peak_start[i] - 1

            non_peak_start <- c(non_peak_start, tmp_start)
            non_peak_stop <- c(non_peak_stop, tmp_stop)
                                        #Range 2
            tmp_start <- df$peak_stop[i] + 1
            tmp_stop <- df$peak_start[i + 1] - 1

            non_peak_start <- c(non_peak_start, tmp_start)
            non_peak_stop <- c(non_peak_stop, tmp_stop)
        } else if (i == nrow(df)) {
                                        #Range N (last)
            tmp_start <- df$peak_stop[i] + 1
            if (tmp_start > df$gene_region_end[i]) {

                tmp_start <- df$gene_region_end[i]
            }
            tmp_stop <- df$gene_region_end[i]

            non_peak_start <- c(non_peak_start, tmp_start)
            non_peak_stop <- c(non_peak_stop, tmp_stop)

        } else {
            tmp_start <- df$peak_stop[i] + 1
            tmp_stop <- df$peak_start[i + 1] - 1
            non_peak_start <- c(non_peak_start, tmp_start)
            non_peak_stop <- c(non_peak_stop, tmp_stop)
        }
    }#nrow df
    results <- list(non_peak_start = non_peak_start,
                    non_peak_stop = non_peak_stop)
    return(results)
}

#' Returns the highest read outside the peak region.
#' Considers regions upstream and downstream of the peak
#' @name getMaxNonPeak
#' @param i.df data frame with PyDegradome results
#' @param i.comp Name of comparison pair (divided into test and ctrl)
#' @param core Number of threads for the parallel process
#' @importFrom GenomicRanges GRanges IRanges subsetByOverlaps
getMaxNonPeak <- function(i.df, i.comp, core = 1) {

                                        # Separate record by with and without gene ID
    df_w_gene <- i.df[!is.na(i.df$ID), ]
    df_wo_gene <- i.df[is.na(i.df$ID), ]

                                        # Peak range can go outside genes.
                                        #In that case, set start =0 end = 0

    ## Define GRanges upstream the peak region
    gr_f <- GRanges(seqnames = df_w_gene$chr,
                    ranges = IRanges(#
                        start = ifelse(#
                            df_w_gene$gene_region_start < df_w_gene$peak_start,
                            df_w_gene$gene_region_start,  0),
                        end = ifelse(#
                            df_w_gene$gene_region_start < df_w_gene$peak_start,
                            df_w_gene$peak_start - 1,  0)),
                    strand = df_w_gene$strand)

    ## Define GRanges downstream the peak region
    gr_r <- GRanges(seqnames = df_w_gene$chr,
                    ranges = IRanges(#
                        start = ifelse(#
                            df_w_gene$peak_stop < df_w_gene$gene_region_end,
                            df_w_gene$peak_stop + 1, 0),
                        end =  ifelse(#
                            df_w_gene$peak_stop < df_w_gene$gene_region_end,
                            df_w_gene$gene_region_end, 0)),
                    strand = df_w_gene$strand)

    sample_name <- getFirstComparison(i.comp)
    n_bigwig <- dg_bigwig_all[[sample_name]]

    max_non_peak <-  getMaxRead(grA = gr_f,  grB = gr_r,
                                f_bigwig = n_bigwig,
                                core = env$cores)

    df_w_gene$max_np_gene <- as.numeric(max_non_peak)
    if (nrow(df_wo_gene) > 0) {
        df_wo_gene$max_np_gene <- ""
    }else {
        df_wo_gene <- NULL
    }
    return(rbind(df_w_gene, df_wo_gene))
}

#' Takes a dataframe with coordinates for a single peak,
#' gene start and end positions. With those coordinates it process a bigwig file
#' and gets highest read value outside the peak region.
#' Results are appended to the input dataframe
#' @name getMaxNonPeak_single
#' @param i.df data.frame of pyDegradome result
#' @param i.comp Name of comparison pair (divided into test and ctrl)
#' @param cols.indx A vector of column names used to construct the index.
#' @importFrom dplyr dplyr::select
#' @importFrom data.table fwrite
#' @export
getMaxNonPeak_single <- function(df, i.comp, cols.indx,
                               idir = pydeg_reads_dir) {
    ## i.comp <- i.comp
    peaks_ref <- file.path(idir,
                           paste0("MaxNonPeak_reads-",
                                  i.comp))


    cat("Calculating signals outside peak region (max_np_gene)...\n")
    cat("Processing ", paste0("MaxNonPeak_reads-", i.comp), "\n")
    if (file.exists(peaks_ref)) {
        cat("\tUsing reference file...\n")
        i.reads <- read.table(peaks_ref, header = TRUE, sep = "\t")
        i.reads <- unique(i.reads)

        ##Subset entries with ID
        df.ID <- df[grepl("^AT", df$ID), ]#!is.na(df$ID)
        if (!is.null(df.ID) && nrow(df.ID) > 0) {
            df.ID <- within(df.ID,
                            indx <- paste(ID, chr, strand,
                                          peak_start, peak_stop,
                                          sep = ""))

            dfs.ID <- MatchDFwithRef(df.ID, i.reads,"max_np_gene")
            df.ID1 <- dfs.ID$df_match
            df.ID2 <- dfs.ID$df_nonmatch
        } else {
            df.ID <- NULL
        }
                                        #Add max_np_gene to the missing entries
        cat("\t\tAdd max_np_gene to missing entries...\n")
        if (!is.null(df.ID2) && nrow(df.ID2) > 0) {
            df.ID2 <- getMaxNonPeak(df.ID2, i.comp, core = env$cores)
            ##---------------------------------------------
                                        #Append the new entries
                                        #to the reference file
            i.reads2 <-  dplyr::select(df.ID2,
                                       c(all_of(cols.indx), "max_np_gene"))
            i.reads2 <- within(i.reads2,
                               indx <- paste(ID, chr,strand,
                                             peak_start, peak_stop,
                                             sep = ""))
            i.reads <- rbind(i.reads, i.reads2)
            i.reads <- unique(i.reads)
            i.reads <- i.reads[with(i.reads, order(indx)), ]
            
            tryCatch({
                                        #remove indx column
                df.ID2 <- dplyr::select(df.ID2,-indx)
                fwrite(i.reads, peaks_ref, sep = "\t", row.names = FALSE)
                cat(paste0("\t\t\tWrote ",  basename(peaks_ref), "\n"))
            }, error = function(e) {
                message("Writing peaks failed")
                message(e)
            }) #Trycatch
            ##---------------------------------------------
        } else {
            df.ID2 <- NULL
        }#df.ID2 missing max_np_gene entries
        cat("\t\tFinished adding max_np_gene to missing entries\n")
                                        #Subset entries w/o ID and
                                        #assing max_np_gene ""
        cat("\t\tWorking on entries without ID...\n")
        df.NonID <- df[!grepl("^AT", df$ID), ]
        if (nrow(df.NonID) > 0) {
            df.NonID$max_np_gene <- NA
        }else {
            df.NonID <- NULL
        }

                                        #Merge data
        cat("\t\tMerging dataframes...\n")

        df <- rbind(df.ID1, df.ID2, df.NonID)
        df <- unique(df)

    } else {
        cat("\tNo reference was found.\n")
        cat("\t\tCalculating max_np_gene...\n")
        df <- getMaxNonPeak(df, i.comp, core = env$cores)
        ##---------------------------------------------
                                        #Write data to ref file
        i.reads <- dplyr::select(df,
                                 c(all_of(cols.indx), "max_np_gene"))
        i.reads <- unique(i.reads)
        i.reads <- within(i.reads,
                          indx <- paste(ID, chr,strand,
                                        peak_start, peak_stop,
                                        sep = ""))
        i.reads <- i.reads[with(i.reads, order(indx)), ]
        
        tryCatch({
            fwrite(i.reads, peaks_ref, sep = "\t", row.names=FALSE)
            cat(paste0("\t\t\tWrote ",  basename(peaks_ref), "\n"))
        }, error = function(e) {
            stopCluster(cl)
            message("Writing peaks failed")
            message(e)
        }) #Trycatch
        ##---------------------------------------------
    }#Checking reference
    return(df)
    cat("Finished calculating signals outside peak region (max_np_gene)\n")
}

#' Takes a data frame with coordinates for two or more peaks.
#' Then, using a dataframe with non-peak coordinates, process a bigwig file
#' to get the highest read value in every segment (region).
#' @name getMaxNonPeakRegion
#' @param df A data.frame with non peak coordinates of entries with more than one peak in the gene
#' @param bigwig Filename of the bigwig file
#' @param core Number of threads for the parallel process
#' @param i.comp Name of comparison pair (divided into test and ctrl)
#' @import doParallel foreach
#' @importFrom GenomicRanges GRanges IRanges
#' @importFrom data.table fwrite
#' @export
getMaxNonPeakRegion <- function(df, bigwig, i.comp, core = env$cores) {
                                        #Columns for index & comparison
    i.cols  <- c("ID",
                 "chr",
                 "strand",
                 "non_peak_start",
                 "non_peak_stop")

    ## bigwig <- bigwig ##<<- Global assignment may be needed?

    peaks_ref <- file.path(pydeg_reads_dir,
                           paste0("MaxNonPeak_reads_dup-",
                                  i.comp))

    cat("\t\tCalculating signals outside peak region (max_non_peak_region Peaks +1)...\n")
    cat("\t\tProcessing ", paste0("MaxNonPeak_reads-", i.comp), "\n")
    if (file.exists(peaks_ref)) {
        cat("\tUsing reference file...\n")
        i.reads <- read.table(peaks_ref, header = TRUE, sep = "\t")
        i.reads <- unique(i.reads)
                                        #Subset entries with ID
        df.ID <- df[grepl("^AT", df$ID), ]
        if (nrow(df.ID) > 0) {
            df.ID <- within(df.ID,
                            indx <- paste(ID, chr,strand,
                                          non_peak_start,
                                          non_peak_stop,
                                          sep = ""))

            dfs.ID <- MatchDFwithRef(df.ID, i.reads,"max_non_peak_region")
            df.ID1 <- dfs.ID$df_match
        } else {
            df.ID <- NULL
        }#Check if df.ID has rows

        ## Add max_non_peak_region to the missing entries
        cat("\t\tAdd max_non_peak_region to missing entries...\n")
        if (!is.null(dfs.ID$df_nonmatch) && nrow(dfs.ID$df_nonmatch) > 0) {
            df.ID2 <- dfs.ID$df_nonmatch


            ##==================================================
            cat("\t\tdefining GRanges...\n")

            ##GRanges for the target region
            gr <- GRanges(
                seqnames = df.ID2$chr,
                ranges = IRanges(start = df.ID2$non_peak_start,
                                 end = df.ID2$non_peak_stop),
                strand = df.ID2$strand
            )

            cat("\t\tCalculating Non Peak Max...\n")

            df.ID2$max_non_peak_region <- getMaxRead(grA = gr,
                                                     f_bigwig = bigwig,
                                                     core = env$cores)

            cat("\t\tFinished Non Peak Max calculation\n")
            ##==================================================
            ## Append the new entries to the reference file
            i.reads2 <- dplyr::select(df.ID2,
                                      c(all_of(i.cols),
                                        "max_non_peak_region"))
            i.reads2 <- within(i.reads2,
                               indx <- paste(ID, chr,strand,
                                             non_peak_start,
                                             non_peak_stop,
                                             sep = ""))
            i.reads <- rbind(i.reads, i.reads2)
            i.reads <- unique(i.reads)
            i.reads <- i.reads[with(i.reads, order(indx)), ]
            
            tryCatch({
                                        #remove indx column
                df.ID2 <- dplyr::select(df.ID2,-indx)
                fwrite(i.reads, peaks_ref, sep = "\t", row.names = FALSE)
                cat(paste0("\t\t\tWrote ",  basename(peaks_ref), "\n"))
            }, error = function(e) {
                message("Writing peaks failed")
                message(e)
            }) #Trycatch
            ##---------------------------------------------
        } else {
            df.ID2 <- NULL
        }#df.ID2 missing max_non_peak_region entries
        cat("\t\tFinished adding max_non_peak_region to missing entries\n")
                                        #Subset entries w/o ID and
                                        #assing max_non_peak_region ""
        cat("\t\tWorking on entries without ID...\n")
        df.NonID <- df[!grepl("^AT", df$ID), ]
        df.NonID <- NULL

                                        #Merge data
        cat("\t\tMerging dataframes...\n")

        df <- rbind(df.ID1, df.ID2, df.NonID)
        df <- unique(df)

    } else {
        cat("\tNo reference was found.\n")
        cat("\t\tCalculating max_non_peak_region...\n")
        ##==================================================
        ## df$max_non_peak_region <- MaxNonPeak(df, i.comp)
        cat("\t\tdefining GRanges...\n")
                                        #GRanges for the target region
        gr <- GRanges(
            seqnames = df$chr,
            ranges = IRanges(start = df$non_peak_start,
                             end = df$non_peak_stop),
            strand = df$strand
        )

        cat("\t\tCalculating Non Peak Max...\n")

        df$max_non_peak_region <- getMaxRead(grA = gr,
                                             f_bigwig = bigwig,
                                             core = env$cores)

        cat("\t\tFinished Non Peak Max calculation\n")
        ##==================================================
        ##---------------------------------------------
                                        #Write data to ref file
        i.reads <- dplyr::select(df,
                                 c(all_of(i.cols),
                                   "max_non_peak_region"))
        i.reads <- unique(i.reads)
        i.reads <- within(i.reads,
                          indx <- paste(ID,
                                        chr,
                                        strand,
                                        non_peak_start,
                                        non_peak_stop,
                                        sep = ""))
        i.reads <- i.reads[with(i.reads, order(indx)), ]
        
        tryCatch({
            fwrite(i.reads, peaks_ref,sep = "\t", row.names = FALSE)
            cat(paste0("\t\t\tWrote ", basename(peaks_ref), "\n"))
        }, error = function(e) {
            message("Writing peaks failed")
            message(e)
        }) #Trycatch
        ##---------------------------------------------
    }#Checking reference
    return(df)
    cat("Finished calculating signals outside peak region (max_non_peak_region)\n")
}

#' Takes the output of getMaxNonPeakRegion and return the highest read
#' along the gene (loops through isoforms)
#' @name getMaxNonPeak_multiple
#' @param pydeg_np_region Ouput from getMaxNonPeakRegion
#' @param dup_ID tx_names from transcripts to select the max read from
#' @param core Number of cores for this function
#' @import doParallel foreach
#' @export
getMaxNonPeak_multiple <- function(pydeg_np_region, dup_ID, core = 1) {
                                        #Define parallel parameters
    cl <- makeCluster(as.numeric(core))
    registerDoParallel(cl)
    on.exit(stopCluster(cl))
                                        #Loop over all the IDS
    pydeg_np_gene <- foreach(j = seq_along(dup_ID),
                             .packages = c("data.table"),
                             .combine = "rbind",
                             .inorder = TRUE) %dopar% {
                                 ## single ID
                                 i.id <- dup_ID[j]
                                 ## subset of regions for i.id
                                 tmp_pydeg_np_gene <- pydeg_np_region[pydeg_np_region$ID==i.id, ]
                                        #calculate absolute max_non_peak
                                 tmp_max_np_gene <- max(tmp_pydeg_np_gene$max_non_peak_region)
                                 ##-------------------------
                                        #create df to merge max_np_gene_gene
                                        #single line tmp_df_np_merge
                                 pydeg_np_gene <- data.table(#
                                     chr = tmp_pydeg_np_gene$chr[1],
                                     strand = tmp_pydeg_np_gene$strand[1],
                                     ID = tmp_pydeg_np_gene$ID[1],
                                     max_np_gene = tmp_max_np_gene)

                                 cat(paste0("Finished ID = ", i.id, "\n"))
                                 pydeg_np_gene
                             }#Loop over ID values
    return(pydeg_np_gene)
}

#' Splits a dataframe based on whether the values of a column are duplicated or unique.
#' Returns a list of dataframes
#' @name splitDFbyCol
#' @param df Dataframe with peak and gene coordinates
#' @param ref_col Column used to judge duplicates
#' @param keep_ref Boolean for keeping or removing 'ref_col'
#' @export
splitDFbyCol <- function(df, ref_col, keep_ref = FALSE) {
                                        #get table of duplicated values
    n_occur_id <- data.frame(table(df[[ref_col]]))
    colnames(n_occur_id) <- c("Var1", "Freq")
                                        #logical vect of duplicated values
    i.log.dup <- df[[ref_col]] %in% n_occur_id$Var1[n_occur_id$Freq > 1]

    results <- list(df_dup = df[i.log.dup, ],
                    df_uniq = df[!i.log.dup, ])
    if (!keep_ref) {
        results <- lapply(results, function(x) {
            dplyr::select(x, -ref_col)
        })
    }
    return(results)
}

#' Takes a dataframe with start and stop coordinates and calculates (from a bigwig file)
#' the highest read value. The function loops over a list of bigwig files
#' corresponding to replicates of a particular condition (test or control). The boolean 'gene' is used to define the columns used to get coordinates and the name of the reference file If TRUE, coordinates will be taken from the gene_region_ column and the reference will be have the root name 'Max_reads_tx-'; if FALSE coordinates
#' @name maxReadPairs
#' @param f_df data frame with pydegradome processed output
#' @param f_bigwigs List of bigwig records. BigWig files should be loaded by `rtracklayer::import
#' @param f_pairs pairs of genotype B concentration to be compared
#' @param f_input_dir directory of the pooled file
#' @param gene boolean used to define which columns and references will be used. TRUE to calculate control max read and FALSE to calculate peak max read (test)
#' @param f_core Number of threads for the nested funtion (getMaxread)
#' @export
maxReadPairs <- function(f_df, f_bigwigs, f_pairs, f_input_dir,
                         gene=FALSE, f_core = 1) {
    if (gene) {
        start_col <- "gene_region_start"
        end_col <- "gene_region_end"
        ref_name <- "Max_reads_tx-"
        ref_col <- "max_read_tx"
        imsg <- "control"
    } else {
        start_col <- "peak_start"
        end_col <- "peak_stop"
        ref_name <- "Peak_reads_test-"
        ref_col <- "max_peak"
        imsg <- "test"
    }
    
    setDF(f_df)
    f_df <- within(f_df,
                   indx <- paste(tx_name,
                                 chr,strand,
                                 eval(parse(text=start_col)),
                                 eval(parse(text=end_col)),
                                 sep = ""))
    isel <- !duplicated(f_df$indx)
    f_df <- f_df[isel, ]
    f_cols <- colnames(f_df)
    out_list <- list()
    
  
    for (i in seq_along(f_pairs)) {
        f_sample <- f_pairs[i]
        f_bigwig <- f_bigwigs[[f_sample]]

        reads_ref <- file.path(f_input_dir,
                               paste0(ref_name,
                                      f_sample))

        if (file.exists(reads_ref)) {
            cat("\tUsing reference file...\n")
            f_reads <- read.table(reads_ref, header = TRUE,  sep = "\t")
            dfs.ID <- MatchDFwithRef(f_df, f_reads, ref_col)
            df1 <- dfs.ID$df_match
            df2 <- dfs.ID$df_nonmatch

            cat("\t\tAdd peak signal to missing entries...\n")
            if (!is.null(df2) && nrow(df2) > 0) {
                f_gr <- GRanges(seqnames = df2$chr,
                                ranges = IRanges(start = df2[[start_col]],
                                               end = df2[[end_col]]),
                                strand = df2$strand)
                df2[[ref_col]] <- getMaxRead(grA = f_gr,
                                           f_bigwig = f_bigwig,
                                           core = f_core)

                cat("\t\tFinished peak signal calculation on", imsg,  "sample\n")
                f_reads2 <- df2
                f_reads <- rbind(f_reads, f_reads2)
                f_reads <- unique(f_reads)
                f_reads<- f_reads[with(f_reads, order(indx)), ]
                fwrite(f_reads,reads_ref,sep = "\t", row.names = FALSE)
                df2 <- dplyr::select(df2, -indx)
                cat(paste0("\t\t\tWrote ", basename(reads_ref), "\n"))
            } else {
                df2 <- NULL
            }
            cat("\t\tMerging dataframes...\n")

            df.out <- rbind(df1, df2)
            df.out$sample <- f_sample
            out_list[[f_sample]] <- unique(df.out)

        } else {
            cat("\tNo reference was found.\n")
            cat("\t\tCalculating peak signal on",imsg,"sample...\n")
            f_gr <- GRanges(seqnames = f_df$chr,
                                ranges = IRanges(start = f_df[[start_col]],
                                               end = f_df[[end_col]]),
                                strand = f_df$strand)
            f_df[[ref_col]] <- getMaxRead(grA = f_gr,
                                        f_bigwig = f_bigwig, core = f_core)

            cat("\t\tFinished peak signal calculation on",imsg,"sample\n")
            f_reads <- dplyr::select(f_df,
                                     all_of(c(f_cols,ref_col, "indx")))
            f_reads <- unique(f_reads)
            f_reads<- f_reads[with(f_reads, order(indx)), ]
            
            fwrite(f_reads,reads_ref,sep = "\t", row.names=FALSE)
            cat(paste0("\t\t\tWrote ", basename(reads_ref), "\n"))

            df.out <- dplyr::select(f_df,-indx)
            df.out$sample <- f_sample
            out_list[[f_sample]] <- unique(df.out)
        }
        cat("Finished calculating peak signals on",imsg,"sample: ",
            f_sample, "\n")
    }
    return(out_list)
}
