#' Save pdf of degradome profiles around detected peaks from a processed pyDegradome result table.
#' @name plotProfilesAllPydegResults
#' @param pydeg data.frame of (processed) pyDegradome output.
#' @param core Number of thread for multi-threading.
#' @param min_plot_width Numeric. Minimum plotting range.
#' @param sample_list Character vector of sampel names to be plotted. Should be a part of `names(bigwigs)`.
#' @param bigwigs List of bigwig records. BigWig files should be loaded by `rtracklayer::import`.
#' @param col_cond Named vector of color code. This must have names identical to `names(bigwigs)`.
#' @param plotrange Character. Define range of the plot.'transcript' for whole transcript or 'peak' for around peaks. `tr_range` must be given for 'transcript'
#' @param tr_range Granges object for transcript ranges. Result of `GenomicFeatures::transcripts(txdb)`
#' @param col_hl Charactor (vector) of color code to be used in highlighting.
#' @param out_dir Path to output directory.
#' @import doParallel
#' @export
plotProfilesAllPydegResults <- function(pydeg,
                                     sample_list,
                                     bigwigs,
                                     bigwigs2 = NULL,
                                     out_dir,
                                     min_plot_width=1000,
                                     col_cond="blue",
                                     col_hl="#00ff0099",
                                     plotrange = "peak",
                                     tr_range = NULL,
                                     core=1) {

    if(plotrange == "transcript" && is.null(tr_range)) {
        stop("tr_range must be given when plotrange = 'transcript'")
    }

                                        #Plot only when the signal has gene annotation
    pydeg <- pydeg[!is.na(pydeg$ID),]

                                        #Generate profile plot(pdf) for each peaks detected
    cl <- makeCluster(as.numeric(core))
    registerDoParallel(cl)
    on.exit(stopCluster(cl))

                                        # foreach(i=sample_list, .export=ls(envir=parent.frame()),.packages = c("ribotools")) %do% {SaveGeneReadDrv(i)}
                                        # foreach(i=1:nrow(pydeg), .export=ls(envir=parent.frame()),.packages = c("cid7degradomeR","Gviz","stringr","biomaRt","GenomicFeatures")) %dopar%{
    foreach(i=1:nrow(pydeg), .export= c("pydeg",
                                        "plotrange",
                                        "min_plot_width",
                                        "tr_range",
                                        "sample_list",
                                        "bigwigs",
                                        "bigwigs2",
                                        "col_cond",
                                        "col_hl"),

                                        # .packages = c("cid7degradomeR","Gviz","stringr","biomaRt","GenomicFeatures")) %dopar%{
            .packages = c("cid7degradomeR",
                          "stringr",
                          "GenomicFeatures"),
            .errorhandling = 'remove') %dopar%{
                                        # for(i in 1:nrow(pydeg)) {
                                        # for(i.record in 1:nrow(pydeg)) {
                cat("\tPlotting ",i,"/",nrow(pydeg),"\r")
                                        #Get gene ID and plot all peaks detected on the gene.
                gene_ID <- pydeg[i,"ID"]
                if(!is.na(gene_ID) && (gene_ID!="NA")) {
                                        #gene_ID could have multiple genes saparated by ";". Get all records on those genes.
                    gene_ID_list <- unique(unlist(str_split(gene_ID,";")))

                                        # Peak positions
                                        # sel <- (pydeg$ID %in% gene_ID_list)
                    sel <-grepl(paste(gene_ID_list,collapse = "|"), pydeg$ID)

                    sel[is.na(sel)] <- FALSE
                    p_start <- pydeg[sel,"peak_start"]
                    p_end <- pydeg[sel,"peak_stop"]

                                        # Plotting range
                    if(plotrange =="peak") {
                                        #To have minimum x plot range
                        max_peak_distace <- max(p_end)-min(p_start)
                        plot_offset <- max(40,round((min_plot_width-max_peak_distace)/2))
                        plot_start <- max(0,min(p_start)-plot_offset)
                        plot_end <- max(p_end)+plot_offset
                    } else if (plotrange == "transcript") {
                        
                        gr <- tr_range[tr_range$tx_name %in% gene_ID_list ,]
                        if(length(gr)>0) {
                            plot_start <- min(start(gr))
                            plot_end <- max(end(gr))
                        }else {
                            plot_start <- max(min(p_start) -  min_plot_width/2,0)
                            plot_end <- max(p_end) + min_plot_width/2
                        }

                    } else {
                        stop("invalid plotrange. It must be \'peak\' or \'transcript\'.")
                    }

                                        #Plot
                    tryCatch(
                    {
                        my.plot <- file.path(out_dir,
                                             paste0(str_replace(gene_ID,
                                                                ";",
                                                                "_"),
                                                    ".pdf"))
                        if(!file.exists(my.plot)) {
                            pdf(my.plot)
                            plotBigWig(chr =pydeg[i,"chr"],
                                       p_start = p_start,
                                       p_end = p_end,
                                       x_start =  plot_start,
                                       x_end = plot_end,
                                       sample_list = sample_list,
                                       ylim="all",
                                       bigwigs = bigwigs,
                                       col_cond = col_cond,
                                       col_hl = col_hl)
                            if(!is.null(bigwigs2)) {
                                plotBigWig(chr =pydeg[i,"chr"],
                                           p_start = p_start,
                                           p_end = p_end,
                                           x_start =  plot_start,
                                           x_end = plot_end,
                                           sample_list = names(bigwigs2),
                                           ylim="all",
                                           bigwigs = bigwigs2,
                                           col_cond = col_cond,
                                           col_hl = col_hl)
                            }#bigwigs2
                            dev.off()
                        }#If file exists
                    },error=function(e) {
                        message(paste0("plotBigWig failed to complete plotting for: chr=",pydeg[i,"chr"]," plot_start=",plot_start," plot_end=",plot_end, " p_start=",p_start," p_end=",p_end))
                        message(e)
                        bm <- NULL
                    }
                    ) #Trycatch
                }#if NA
            }
}


#' Get max and min coverage from a BigWig file for specific region. Modified from http://adomingues.github.io/2016/11/13/max-coverage-in-bigwigs/
#' @name maxCoverageBigWig
#' @param bigwig A bigwig record. BigWig files should be loaded by `rtracklayer::import`.
#' @param chr Char of chromosome name
#' @param x_start Numeric. Coordinate of left end of the plot range.
#' @param x_end Numeric. Coordinate of right end of the plot range.
#' @import GenomicRanges
#' @export
maxCoverageBigWig <- function(bigwig, chr, x_start, x_end, strand=c("+", "-")) {
                                        #GRanges for the target region
    gr <- GRanges(
        seqnames=chr,
        ranges=IRanges(x_start, x_end),
        strand=strand
    )
                                        #Get BigWig records for the specified region
    ovlp <- subsetByOverlaps(bigwig, gr)
    if (length(ovlp) > 0) {
        max_cov <- max(ovlp$score)
        min_cov <- min(ovlp$score)
    } else {
        max_cov <- 0
        min_cov <- 0
    }
    return(c(min=min_cov,max=max_cov))
}


#' Merges a given column of a dataframe with a reference based on an index 'indx' column
#' Returns a list of dataframes one matched with the reference an another with nonmatches
#' @name MatchDFwithRef
#' @param df data.frame to match against a reference
#' @param ref reference dataframe
#' @param ref_col Column with values to merge 
#' @export
MatchDFwithRef <- function(df, i.reads, ref_col) {
    ## i.reads <- read.table(ref,header = TRUE, sep = "\t")
    ## i.reads <- unique(i.reads)
    ## df <- within(df, indx <- paste(ID,chr,strand,non_peak_start,non_peak_stop,
                                   ## sep=""))  
    ilog <- df$indx %in% i.reads$indx
    df_match <- df[ilog,]
    df_nonmatch <- df[!ilog,]
    if (!is.null(df_match) && nrow(df_match) > 0) {
        df_match <- merge(df_match,
                          i.reads[, c("indx",ref_col)],
                          by="indx",
                          all =FALSE)[,-1]
    } else {
        df_match <- NULL
    }
    if (is.null(df_nonmatch) || nrow(df_nonmatch) == 0) {
        df_nonmatch <- NULL
    }
    results <- list(df_match=df_match,
                    df_nonmatch=df_nonmatch)
    return(results) 
}

#' Gets the first element of the comparison used for the pydegradome analysis
#' @name SampleName
#' @param icomp filename based on the comparison between test and control samples
#' @export
SampleName <- function(i.comp) {
  sample_name <- gsubfn::strapplyc(i.comp,"t_(.*)_c_.*")[[1]]
    if (!sample_name %in% names(dg_bigwig_all)) {
      sample_name <- lowerCaseSampleName(sample_name)
    }
  return(sample_name)
}


#' Merges a dataframe with intermediary results into the main dataframe 
#' @name MergeDFs
#' @param dfA Dataframe to merge into
#' @param dfB Dataframe to merge from
#' @param ref_col Column name (string) to merge
#' @export
MergeDFs <- function(dfA, dfB, ref_col) {
                                        #Create a indx vector to merge
    indx1 <- paste(dfA$tx_name,dfA$chr,dfA$strand,sep="")
    dfA$indx_merge <- indx1
    
    indx2 <- paste(dfB$ID,dfB$chr,dfB$strand,sep="")
    dfB$indx_merge <- indx2
    
    ref_cols <- c("indx_merge", ref_col)
    dfB1 <- dplyr::select(dfB,all_of(ref_cols))
    dfA <- merge(dfA,  dfB1, by="indx_merge",
                 all.x = TRUE)[,-1]
    
    dfA <- dplyr::select(dfA,-indx_dup)#indx_merge
    return(dfA)
}

#'Calculates the max read in a series of GRanges using a foreach loop
#' @name GetMaxRead
#' @param gr GRanges object
#' @param bigwig reference to a bigwig file
#' @import doParallel foreach iterators
#' @export
GetMaxRead <- function(grA, grB=NULL, f_bigwig, core=1) {
    cl <- makeCluster(as.numeric(core))
    registerDoParallel(cl)
    on.exit(stopCluster(cl))
    
    maxread <- foreach(i=seq_along(grA),
                       .export=c("grA","grB","f_bigwig"),
                       ## .packages = c("cid7degradomeR","GenomicRanges"),
                       .packages = c("GenomicRanges"),
                       .combine = "c",
                       .inorder = TRUE) %dopar% {
                         if (is.null(grB)) {
                           scores <- subsetByOverlaps(f_bigwig, grA[i])$score
                         } else {
                           scores <- c(subsetByOverlaps(f_bigwig, grA[i])$score,
                                     subsetByOverlaps(f_bigwig, grB[i])$score)
                         }
                         
                         if(decode(strand(grA[i])) == "+") {
                           maxread <- max(scores)
                         }else if(decode(strand(grA[i])) == "-") {
                           maxread <- -min(scores)
                         }
                         maxread
                       }
    return(maxread)
  }

#'Calculate highest peak read (Adapted from the addMaxNonPeakSignal function)
#' @name MaxPeakBigWig
#' @param i.df data.frame of pyDegradome result
#' @param core Number of cores available
#' @param i.comp Name of comparison pair (divided into test and ctrl)
#' @import doParallel foreach iterators
#' @export
MaxPeakBigWig <- function(df, i.comp, core=1) {
                                        #GRanges for the target region
    gr <- GRanges(
        seqnames=df$chr,
        ranges=IRanges(start=df$peak_start, end=df$peak_stop),
        strand=df$strand
    )
                                        # Import bigwig data based on i.comp variable
    sample_name <- gsubfn::strapplyc(i.comp,"t_(.*)_c_.*")[[1]]
    if (!sample_name %in% names(dg_bigwig_all)) {
        sample_name <- lowerCaseSampleName(sample_name)
    }
    bigwig <- dg_bigwig_all[[sample_name]]
    
    peak_max <- GetMaxRead(grA=gr, f_bigwig=bigwig, core=env$core)
    return(peak_max)
}


#' Add maxNonPeakRegion to pyDegradome result data.frame
#' @name addMaxNonPeakSignal
#' @param x data.frame of pyDegradome result
#' @param ann data.frame of GFF annotation
#' @param core Number of cores available
#' @import doParallel foreach iterators
#' @export
addMaxNonPeakSignal <- function(i.df, i.comp, core=1) {

                                        # Separate record by with and without gene ID
    df_w_gene <- i.df[!is.na(i.df$ID),]
    df_wo_gene <- i.df[is.na(i.df$ID),]

                                        # Peak range can go outside genes.
                                        #In that case, set start =0 end = 0
    gr_f <- GRanges(seqnames = df_w_gene$chr,
                    ranges = IRanges(#
                        start = ifelse(#
                            df_w_gene$gene_region_start < df_w_gene$peak_start,
                            df_w_gene$gene_region_start, 0),
                        end = ifelse(#
                            df_w_gene$gene_region_start < df_w_gene$peak_start,
                                     df_w_gene$peak_start -1, 0)),
                    strand = df_w_gene$strand)

    gr_r <- GRanges(seqnames = df_w_gene$chr,
                    ranges = IRanges(#
                        start = ifelse(#
                            df_w_gene$peak_stop < df_w_gene$gene_region_end,
                            df_w_gene$peak_stop + 1,0),
                        end =  ifelse(#
                            df_w_gene$peak_stop < df_w_gene$gene_region_end,
                            df_w_gene$gene_region_end,0)),
                    strand = df_w_gene$strand)

    sample_name <- SampleName(i.comp)
    n_bigwig <- dg_bigwig_all[[sample_name]]

    max_non_peak <-  GetMaxRead(grA=gr_f, grB=gr_r, f_bigwig=n_bigwig, core=env$core)
    
    df_w_gene$max_np_gene <- as.numeric(max_non_peak)
    if(nrow(df_wo_gene)>0) {
        df_wo_gene$max_np_gene <- ""
    }else {
        df_wo_gene <- NULL
    }
    return(rbind(df_w_gene,df_wo_gene))
}


#' Returns a list of coordinates for non peak regions.
#' @name GetNonPeakCoors
#' @param df A dataframe with peak and gene coordinates
#' @export
GetNonPeakCoors <- function(df) {
  ## non_peak_start <- NULL
  ## non_peak_stop <- NULL

  ##Sort by peak_start
  df <- df[order(df$peak_start),]
  for (i in 1:nrow(df)) {
    if ( i == 1 ) {
                                        #Range 1
      tmp_start <- df$gene_region_start[i]
      tmp_stop <- df$peak_start[i]-1
      
      non_peak_start<-c(non_peak_start,tmp_start)
      non_peak_stop<-c(non_peak_stop,tmp_stop)
                                        #Range 2
      tmp_start <- df$peak_stop[i]+1
      tmp_stop <- df$peak_start[i+1]-1
      
      non_peak_start<-c(non_peak_start,tmp_start)
      non_peak_stop <- c(non_peak_stop,tmp_stop)
    } else if(i==nrow(df)) {
                                        #Range N (last)
      tmp_start <- df$peak_stop[i]+1
      if(tmp_start > df$gene_region_end[i]) {
        
        tmp_start <- df$gene_region_end[i]
      }
      tmp_stop <- df$gene_region_end[i]
      
      non_peak_start <- c(non_peak_start,tmp_start)
      non_peak_stop <- c(non_peak_stop,tmp_stop)
      
    } else {
      tmp_start <- df$peak_stop[i]+1
      tmp_stop <- df$peak_start[i+1]-1
      non_peak_start <- c(non_peak_start,tmp_start)
      non_peak_stop <- c(non_peak_stop,tmp_stop)
    }
  }#nrow df
  results <- list(non_peak_start = non_peak_start,
                  non_peak_stop = non_peak_stop)
  return(results)
}


#'Calculate highest peak read (Adapted from the addMaxNonPeakSignal function)
#' @name addNonPeak
#' @param df_dup_ID_sub data.frame (subset) of entries with more than one peak in the gene
#' @param core Number of cores available
#' @param i.comp Name of comparison pair (divided into test and ctrl)
#' @import doParallel foreach
#' @export
addNonPeak <- function(df_dup_ID_sub, df_dup_sub, bigwig, i.comp, core=env$core) {
                                        #Vector of columns
    i.cols  <- c("ID",
                 "chr",
                 "strand",
                 "non_peak_start",
                 "non_peak_stop")
    
                                        #Define parallel parameters

    cl <-  makeCluster(as.numeric(core))
    registerDoParallel(cl)
    on.exit(stopCluster(cl))
                                        #Loop over all the IDS
    ## for (i.id in df_dup_ID) {
                                        #remove previous log file
    if (file.exists(file.path(pydeg_log_dir,"Log-dups.txt"))) {
                                        #Delete file if it exists
        file.remove(file.path(pydeg_log_dir,"Log-dups.txt"))
    }

    df_np_merge <- foreach(j=seq_along(df_dup_ID_sub),
                           .export=c("pydeg_reads_dir","pydeg_log_dir","i.cols"),
                           .packages=c("data.table","GenomicRanges","dplyr"),
                           .combine = "rbind",
                           .inorder = TRUE) %dopar% {
                               i.id <- df_dup_ID_sub[j]
                               
                                        #Temp df for a single id
                                        #to define all non peaks
                               temp_df <- df_dup_sub[df_dup_sub$ID==i.id,]
                              
                               ##==================================================
                               ## sink(file.path(pydeg_log_dir,"Log-dups.txt"),append=FALSE)
                               cat("Working on ID",i.id,"\n")
                               cat("\t\tObtaining NonPeak coordinates\n")
                               np_coors <- GetNonPeakCoors(temp_df)
                               tmp_df_np_coor <- data.table(#
                                   chr=temp_df$chr[1],
                                   strand=temp_df$strand[1],
                                   ID=temp_df$ID[1],
                                   non_peak_start=np_coors$non_peak_start,
                                   non_peak_stop=np_coors$non_peak_stop)
                               cat("\t\tFinished obtaining NonPeak coordinates\n")
                               ##==================================================
                               peaks_ref <- file.path(pydeg_reads_dir,
                                                      paste0("MaxNonPeak_reads_dup-",
                                                             i.comp))
                               
                               
                               cat("Calculating signals outside peak region (max_non_peak_region Peaks +1)...\n")
                               cat("Processing ", paste0("MaxNonPeak_reads-",i.comp), "\n")
                               if(file.exists(peaks_ref)) {
  cat("\tUsing reference file...\n")
  i.reads <- read.table(peaks_ref,header = TRUE, sep = "\t")
  i.reads <- unique(i.reads)
                                        #Subset entries with ID
  df.ID <- tmp_df_np_coor[grepl("^AT",tmp_df_np_coor$ID),]#!is.na(df$ID)
  
  if(nrow(df.ID) > 0) {
      df.ID <- within(df.ID,
                      indx <- paste(ID,chr,strand,non_peak_start,non_peak_stop,
                                    sep=""))
      idfs <- MatchDFWithRef(df.ID,i.reads,"max_non_peak_region")
      df.ID1 <- idfs$df_match
  }else {
    df.ID <- NULL
  }
                                        #Add max_non_peak_region to the missing entries
  cat("\t\tAdd max_non_peak_region to missing entries...\n")
  df.ID2 <- idfs$df_nonmatch
  
  if(nrow(df.ID2) > 0) {
    ##==================================================
    ## df.ID2$max_non_peak_region <- MaxNonPeak(df.ID2,i.comp)
    cat("\t\tdefining GRanges...\n")
                                        #GRanges for the target region
    gr <- GRanges(
      seqnames=df.ID2$chr,
      ranges=IRanges(start=df.ID2$non_peak_start,
                     end=df.ID2$non_peak_stop),
      strand=df.ID2$strand
    )
    
    cat("\t\tCalculating Non Peak Max...\n")
    df.ID2$max_non_peak_region <- GetMaxRead(grA=gr,f_bigwig=bigwig, core=env$core)
    cat("\t\tFinished Non Peak Max calculation\n")
    ##==================================================
    ##---------------------------------------------
                                        #Append the new entries
                                        #to the reference file
    ## i.reads2 <- df.ID2[,c(i.cols,"max_non_peak_region")]
    i.reads2 <- dplyr::select(df.ID2,
                              c(all_of(i.cols),
                                "max_non_peak_region"))
    
    
    i.reads2 <- within(i.reads2, indx <- paste(ID,chr,strand,non_peak_start,non_peak_stop,
                                         sep=""))
    i.reads <- rbind(i.reads,i.reads2)
    
    tryCatch({
                                        #remove indx column
      df.ID2 <- dplyr::select(df.ID2,-indx)
      ##fwrite(i.reads,peaks_ref,sep="\t",row.names=FALSE)
      cat(paste0("\t\t\tWrote ",
                 basename(peaks_ref),"\n"))
    },error=function(e) {
      stopCluster(cl)
      message("Writing peaks failed")
      message(e)
    }) #Trycatch
    ##---------------------------------------------
  }else {
    df.ID2 <- NULL
  }#df.ID2 missing max_non_peak_region entries
  cat("\t\tFinished adding max_non_peak_region to missing entries\n")
                                        #Subset entries w/o ID and
                                        #assing max_non_peak_region ""
  cat("\t\tWorking on entries without ID...\n")
  df.NonID <- tmp_df_np_coor[!grepl("^AT",
                                    tmp_df_np_coor$ID),
                             ]#!is.na(df$ID)
  df.NonID <- NULL
  
                                        #Merge data
  cat("\t\tMerging dataframes...\n")
  
  tmp_df_np_coor <- rbind(df.ID1,df.ID2,df.NonID)
  tmp_df_np_coor <- unique(tmp_df_np_coor)
  
} else {
  cat("\tNo reference was found.\n")
  cat("\t\tCalculating max_non_peak_region...\n")
  ##==================================================
  ## tmp_df_np_coor$max_non_peak_region <- MaxNonPeak(tmp_df_np_coor,i.comp)
  cat("\t\tdefining GRanges...\n")
                                        #GRanges for the target region
  gr <- GRanges(
    seqnames=tmp_df_np_coor$chr,
    ranges=IRanges(start=tmp_df_np_coor$non_peak_start,
                   end=tmp_df_np_coor$non_peak_stop),
    strand=tmp_df_np_coor$strand
  )
  
  cat("\t\tCalculating Non Peak Max...\n")
  tmp_df_np_coor$max_non_peak_region <- GetMaxRead(grA=gr, f_bigwig=bigwig, core=env$core)
  cat("\t\tFinished Non Peak Max calculation\n")
  ##==================================================
  ##---------------------------------------------
                                        #Write data to ref file
  i.reads <- dplyr::select(tmp_df_np_coor,
                           c(all_of(i.cols),
                             "max_non_peak_region"))
  i.reads <- unique(i.reads)
  i.reads <- within(i.reads, indx <- paste(ID,chr,strand,non_peak_start,non_peak_stop,
                                         sep=""))
  tryCatch({
    ##fwrite(i.reads,peaks_ref,sep="\t",row.names=FALSE)
    cat(paste0("\t\t\tWrote ",peaks_ref,"\n"))
  },error=function(e) {
    stopCluster(cl)
    message("Writing peaks failed")
    message(e)
  }) #Trycatch
  ##---------------------------------------------
}#Checking reference
                               cat("Finished calculating signals outside peak region (max_non_peak_region)\n")
                               
                               ##==================================================
                               tmp_df_np_coor <- within(tmp_df_np_coor, indx <- paste(ID,chr,strand,non_peak_start,non_peak_stop,
                                         sep=""))
                               ##-------------------------
                                        #calculate absolute max_non_peak
                               tmp_max_np_gene <- max(tmp_df_np_coor$max_non_peak_region)
                               ##-------------------------
                                        #create df to merge max_np_gene_gene
                                        #single line tmp_df_np_merge
                               df_np_merge <- data.table(#
                                   chr=tmp_df_np_coor$chr[1],
                                   strand=tmp_df_np_coor$strand[1],
                                   ID=tmp_df_np_coor$ID[1],
                                   max_np_gene=tmp_max_np_gene)
                               ##df_np_merge = rbind(df_np_merge,tmp_df_np_merge)
                               
                               
                               ## sink(file.path(pydeg_log_dir,"Log-dups.txt"),
                               ##      append=TRUE)
                               cat(paste0("Finished ID = ",i.id, "\n"))
                               ## sink()
                               df_np_merge
                           }#FOR LOOP [i.id in df_dup_ID_sub]
    return(df_np_merge)
}


#'Calculate highest read outside peak (Adapted from the addMaxNonPeakSignal function)
#' @name addMaxNonPeakRegion
#' @param df data.frame (subset) of entries with more than one peak in the gene
#' @param bigwig Filename of the bigwig file
#' @param core Number of cores available
#' @param i.comp Name of comparison pair (divided into test and ctrl)
#' @import doParallel foreach
#' @export
addMaxNonPeakRegion <- function(df, bigwig, i.comp, core=env$core) {
                                        #Columns for index & comparison
    i.cols  <- c("ID",
                 "chr",
                 "strand",
                 "non_peak_start",
                 "non_peak_stop")

    bigwig <<- bigwig
    
    peaks_ref <- file.path(pydeg_reads_dir,
                           paste0("MaxNonPeak_reads_dup-",
                                  i.comp))
    
    
    cat("\t\tCalculating signals outside peak region (max_non_peak_region Peaks +1)...\n")
    cat("\t\tProcessing ", paste0("MaxNonPeak_reads-",i.comp), "\n")
    if(file.exists(peaks_ref)) {
        cat("\tUsing reference file...\n")
        i.reads <- read.table(peaks_ref,header = TRUE, sep = "\t")
        i.reads <- unique(i.reads)
                                        #Subset entries with ID
        df.ID <- df[grepl("^AT",df$ID),]#is.na(df$ID)
        if(nrow(df.ID) > 0) {
            df.ID <- within(df.ID,
                            indx <- paste(ID,chr,strand,non_peak_start,non_peak_stop,
                                          sep=""))
           
            dfs.ID <- MatchDFwithRef(df.ID,i.reads,"max_non_peak_region")
            df.ID1 <- dfs.ID$df_matchdfs.ID <- MatchDFwithRef(df.ID,i.reads,"max_non_peak_region")
            df.ID1 <- dfs.ID$df_match
        } else {
            df.ID <- NULL
        }#Check if df.ID has rows

                                        #Add max_non_peak_region to the missing entries
        cat("\t\tAdd max_non_peak_region to missing entries...\n")
        if (!is.null(dfs.ID$df_nonmatch) && nrow(dfs.ID$df_nonmatch) > 0) {
        df.ID2 <- dfs.ID$df_nonmatch
        
  
            ##==================================================
            ## df.ID2$max_non_peak_region <- MaxNonPeak(df.ID2,i.comp)
            cat("\t\tdefining GRanges...\n")
                                        #GRanges for the target region
            gr <- GRanges(
                seqnames=df.ID2$chr,
                ranges=IRanges(start=df.ID2$non_peak_start, end=df.ID2$non_peak_stop),
                strand=df.ID2$strand
            )
            
            cat("\t\tCalculating Non Peak Max...\n")
            
            df.ID2$max_non_peak_region <- GetMaxRead(grA=gr, f_bigwig=bigwig, core=env$core)
           
            cat("\t\tFinished Non Peak Max calculation\n")
            ##==================================================
            ##---------------------------------------------
                                        #Append the new entries
                                        #to the reference file
            ## i.reads2 <- df.ID2[,c(i.cols,"max_non_peak_region")]
            i.reads2 <- dplyr::select(df.ID2,
                                      c(all_of(i.cols),
                                        "max_non_peak_region"))
            
            
            ## i.reads2 <- unique(i.reads2)
            i.reads2 <- within(i.reads2,
                               indx <- paste(ID,chr,strand,non_peak_start,non_peak_stop,
                                   sep=""))
            i.reads <- rbind(i.reads,i.reads2)
            i.reads <- unique(i.reads)
            tryCatch({
                                        #remove indx column
                df.ID2 <- dplyr::select(df.ID2,-indx)
                fwrite(i.reads,peaks_ref,sep="\t",row.names=FALSE)
                cat(paste0("\t\t\tWrote ",peaks_ref,"\n"))
            },error=function(e) {
                stopCluster(cl)
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
        df.NonID <- df[!grepl("^AT",df$ID),]#!is.na(df$ID)
        df.NonID <- NULL
        ## if(nrow(df.NonID) > 0) {
        ##     df.NonID$max_non_peak_region <- NA
        
        ## }else {
        ##     df.NonID <- NULL
        ## }
        
                                        #Merge data
        cat("\t\tMerging dataframes...\n")
        
        df <- rbind(df.ID1,df.ID2,df.NonID)
        df <- unique(df)
        
    } else {
        cat("\tNo reference was found.\n")
        cat("\t\tCalculating max_non_peak_region...\n")
        ##==================================================
        ## df$max_non_peak_region <- MaxNonPeak(df,i.comp)
        cat("\t\tdefining GRanges...\n")
                                        #GRanges for the target region
        gr <- GRanges(
            seqnames=df$chr,
            ranges=IRanges(start=df$non_peak_start, end=df$non_peak_stop),
            strand=df$strand
        )
        
        cat("\t\tCalculating Non Peak Max...\n")
       
        df$max_non_peak_region <- GetMaxRead(grA=gr, f_bigwig=bigwig, core=env$core)
        
        cat("\t\tFinished Non Peak Max calculation\n")
        ##==================================================
        ##---------------------------------------------
                                        #Write data to ref file
        i.reads <- dplyr::select(df,
                                 c(all_of(i.cols),
                                   "max_non_peak_region"))
        i.reads <- unique(i.reads)
        i.reads <- within(i.reads,
                          indx <- paste(ID,chr,strand,non_peak_start,non_peak_stop,
                                   sep=""))
        tryCatch({
            fwrite(i.reads,peaks_ref,sep="\t",row.names=FALSE)
            cat(paste0("\t\t\tWrote ",peaks_ref,"\n"))
        },error=function(e) {
            stopCluster(cl)
            message("Writing peaks failed")
            message(e)
        }) #Trycatch
        ##---------------------------------------------
    }#Checking reference
    return(df)
    cat("Finished calculating signals outside peak region (max_non_peak_region)\n")
}


#'Calculate highest read outside peak (Genewise) Depends on the results of addMaxNonPeakRegion
#' @name addMaxNonPeakGene
#' @param pydeg_np_region data.frame with non-peak regions and max reads within them
#' @param dup_ID tx_names from transcripts to select the max read from
#' @param core Number of cores for this function
#' @import doParallel foreach
#' @export
 addMaxNonPeakGene <- function(pydeg_np_region, dup_ID, core=1) {
                                        #Define parallel parameters
     cl <- makeCluster(as.numeric(core))
     registerDoParallel(cl)
     on.exit(stopCluster(cl))
                                        #Loop over all the IDS 
     pydeg_np_gene <- foreach(j=seq_along(dup_ID),
                              .export=c("dup_ID", "pydeg_np_region"),
                              .packages=c("data.table"),
                              .combine = "rbind",
                              .inorder = TRUE) %dopar% {
                                        #single ID
                                  i.id <- dup_ID[j]

                                        #subset of regions for i.id
                                  tmp_pydeg_np_gene <- pydeg_np_region[pydeg_np_region$ID==i.id,]
                                        #calculate absolute max_non_peak
                                  tmp_max_np_gene <- max(tmp_pydeg_np_gene$max_non_peak_region)
                                  ##-------------------------
                                        #create df to merge max_np_gene_gene
                                        #single line tmp_df_np_merge
                                  pydeg_np_gene <- data.table(#
                                      chr=tmp_pydeg_np_gene$chr[1],
                                      strand=tmp_pydeg_np_gene$strand[1],
                                      ID=tmp_pydeg_np_gene$ID[1],
                                      max_np_gene=tmp_max_np_gene)
                                  
                                  cat(paste0("Finished ID = ",i.id, "\n"))
                                  ## sink()
                                  pydeg_np_gene
                              }#Loop over ID values
     return(pydeg_np_gene)
}


#'Calculate highest peak read (Adapted from the addMaxNonPeakSignal function)
#' @name addNonPeakSingle
#' @param df_dup_ID_sub data.frame (subset) of entries with more than one peak in the gene
#' @param core Number of cores available
#' @param i.comp Name of comparison pair (divided into test and ctrl)
#' @import doParallel foreach
#' @export
addNonPeakSingle <- function(df_dup_ID_sub, df_dup_sub, bigwig, i.comp) {
                                        #Vector of columns
    i.cols  <- c("ID",
                 "chr",
                 "strand",
                 "non_peak_start",
                 "non_peak_stop")
                                        #Loop over all the IDS
    ## for (i.id in df_dup_ID) {
                                        #remove previous log file
    if (file.exists(file.path(pydeg_log_dir,"Log-dups.txt"))) {
                                        #Delete file if it exists
        file.remove(file.path(pydeg_log_dir,"Log-dups.txt"))
    }


    df_np_total <- data.table()
    for(j in seq_along(df_dup_ID_sub)) {
        
        i.id <- df_dup_ID_sub[j]

                                        #Temp df for a single id
                                        #to define all non peaks
        temp_df <- df_dup_sub[df_dup_sub$ID==i.id,]
                                        #Order based on peak start
        temp_df <- temp_df[order(temp_df$peak_start),]
        ##==================================================
        ## sink(file.path(pydeg_log_dir,"Log-dups.txt"),
        ##      append=TRUE)
        cat("\t\tObtaining NonPeak coordinates\n")
        ## tmp_df_np_coor <- NonPeakCoor(temp_df)
        non_peak_start <- NULL
        non_peak_stop <- NULL
        for (i.row in 1:nrow(temp_df)) {
            
            if ( i.row == 1 ) {
                                        #Range 1
                tmp_start <- temp_df$gene_region_start[i.row]
                tmp_stop <- temp_df$peak_start[i.row]-1

                non_peak_start <- c(non_peak_start,tmp_start)
                non_peak_stop <- c(non_peak_stop,tmp_stop)
                                        #Range 2
                tmp_start <- temp_df$peak_stop[i.row]+1
                tmp_stop <- temp_df$peak_start[i.row+1]-1

                non_peak_start <- c(non_peak_start,tmp_start)
                non_peak_stop <- c(non_peak_stop,tmp_stop)
            } else if(i.row==nrow(temp_df)) {
                                        #Range N (last)
                tmp_start <- temp_df$peak_stop[i.row]+1
                if(tmp_start > temp_df$gene_region_end[i.row]) {
                    
                    tmp_start <- temp_df$gene_region_end[i.row]
                }
                tmp_stop <- temp_df$gene_region_end[i.row]

                non_peak_start <- c(non_peak_start,tmp_start)
                non_peak_stop <- c(non_peak_stop,tmp_stop)

            } else {
                tmp_start <- temp_df$peak_stop[i.row]+1
                tmp_stop <- temp_df$peak_start[i.row+1]-1
                non_peak_start <- c(non_peak_start,tmp_start)
                non_peak_stop <- c(non_peak_stop,tmp_stop)
            }
        }#nrow temp_df
        tmp_df_np_coor <- data.table(chr=temp_df$chr[1],
                                     strand=temp_df$strand[1],
                                     ID=temp_df$ID[1],
                                     non_peak_start=non_peak_start,
                                     non_peak_stop=non_peak_stop)
        cat("\t\tFinished obtaining NonPeak coordinates\n")
        ##==================================================
        ## tmp_df_np_coor <- addNonPeakRefDups(tmp_df_np_coor,i.comp)
        

        peaks_ref <- file.path(pydeg_reads_dir,
                               paste0("MaxNonPeak_reads_dup-",
                                      i.comp))


        cat("Calculating signals outside peak region (max_non_peak_region Peaks +1)...\n")
        cat("Processing ", paste0("MaxNonPeak_reads-",i.comp), "\n")
        if(file.exists(peaks_ref)) {
            
            cat("\tUsing reference file...\n")
            i.reads <- read.table(peaks_ref,header = TRUE, sep = "\t")
            i.reads <- unique(i.reads)
                                        #Subset entries with ID
            df.ID <- tmp_df_np_coor[grepl("^AT",tmp_df_np_coor$ID),]#!is.na(df$ID)

            if(nrow(df.ID) > 0) {
                df.ID <- within(df.ID,
                                indx <- paste(ID,chr,strand,non_peak_start,non_peak_stop,
                                   sep=""))

                                        #Logical vector of indeces with data in reference file
                i.wread <- df.ID$indx %in% i.reads$indx

                cat("\t\tMerging max_non_peak_region from Ref...\n")
                df.ID1 <- df.ID[i.wread,]
                ## df.ID1 <- unique(df.ID1)
                df.ID1 <- merge(df.ID1,
                                i.reads[,
                                        c("indx",
                                          "max_non_peak_region")],
                                by="indx",
                                all =FALSE)[,-1]#all.x = TRUE
            }else {
                df.ID <- NULL
            }
                                        #Add max_non_peak_region to the missing entries
            cat("\t\tAdd max_non_peak_region to missing entries...\n")
            df.ID2 <- df.ID[!i.wread,]

            if(nrow(df.ID2) > 0) {
                ##==================================================
                ## df.ID2$max_non_peak_region <- MaxNonPeak(df.ID2,i.comp)
                cat("\t\tdefining GRanges...\n")
                                        #GRanges for the target region
                gr <- GRanges(
                    seqnames=df.ID2$chr,
                    ranges=IRanges(start=df.ID2$non_peak_start, end=df.ID2$non_peak_stop),
                    strand=df.ID2$strand
                )
                
                cat("\t\tCalculating Non Peak Max...\n")
                df.ID2$max_non_peak_region <- GetMaxRead(grA=gr,
                                                         f_bigwig=bigwig,
                                                         core=env$core)
                cat("\t\tFinished Non Peak Max calculation\n")
                ##==================================================
                ##---------------------------------------------
                                        #Append the new entries
                                        #to the reference file
                ## i.reads2 <- df.ID2[,c(i.cols,"max_non_peak_region")]
                i.reads2 <- dplyr::select(df.ID2,
                                          c(all_of(i.cols),"max_non_peak_region"))
                ## i.reads2 <- unique(i.reads2)
                i.reads2 <- within(i.reads2,
                             indx <- paste(ID,chr,strand,non_peak_start,non_peak_stop,
                                   sep="")) 
                i.reads <- rbind(i.reads,i.reads2)
                i.reads <- unique(i.reads)
                tryCatch({
                    ## i.reads <- apply(i.reads,2,as.character)
                                        #remove indx column
                    df.ID2 <- dplyr::select(df.ID2,-indx)
                    fwrite(i.reads,peaks_ref,sep="\t",row.names=FALSE)
                    cat(paste0("\t\t\tWrote ",peaks_ref,"\n"))
                },error=function(e) {
                    
                    stopCluster(cl)
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
            df.NonID <- tmp_df_np_coor[!grepl("^AT",tmp_df_np_coor$ID),]#!is.na(df$ID)
            df.NonID <- NULL
            ## if(nrow(df.NonID) > 0) {
            ##     df.NonID$max_non_peak_region <- NA

            ## }else {
            ##     df.NonID <- NULL
            ## }

                                        #Merge data
            cat("\t\tMerging dataframes...\n")

            tmp_df_np_coor <- rbind(df.ID1,df.ID2,df.NonID)
            tmp_df_np_coor <- unique(tmp_df_np_coor)

        } else {
            cat("\tNo reference was found.\n")
            cat("\t\tCalculating max_non_peak_region...\n")
            ##==================================================
            ## tmp_df_np_coor$max_non_peak_region <- MaxNonPeak(tmp_df_np_coor,i.comp)
            cat("\t\tdefining GRanges...\n")
                                        #GRanges for the target region
            gr <- GRanges(
                seqnames <- tmp_df_np_coor$chr,
                ranges <- IRanges(start=tmp_df_np_coor$non_peak_start, end=tmp_df_np_coor$non_peak_stop),
                strand <- tmp_df_np_coor$strand
            )
                                        # Import bigwig data based on i.comp variable
            ## sample_name <- gsubfn::strapplyc(i.comp,"t_(.*)_c_.*")[[1]]
            ## if (!sample_name %in% names(dg_bigwig_all)) {
            ##     sample_name <- lowerCaseSampleName(sample_name)}
            ## bigwig <- dg_bigwig_all[[sample_name]]

            cat("\t\tCalculating Non Peak Max...\n")
            tmp_df_np_coor$max_non_peak_region <- GetMaxRead(#
                grA=gr, f_bigwig=bigwig, core=env$core)
            cat("\t\tFinished Non Peak Max calculation\n")
            ##==================================================
            ##---------------------------------------------
                                        #Write data to ref file
            i.reads <- dplyr::select(tmp_df_np_coor,
                                     c(all_of(i.cols),"max_non_peak_region"))
            i.reads <- unique(i.reads)
            i.reads <- within(i.reads,
                         indx <- paste(ID,chr,strand,non_peak_start,non_peak_stop,
                                   sep=""))
            tryCatch({
                ## i.reads <- apply(i.reads,2,as.character)
                fwrite(i.reads,peaks_ref,sep="\t",row.names=FALSE)
                cat(paste0("\t\t\tWrote ",peaks_ref,"\n"))
            },error=function(e) {
                stopCluster(cl)
                message("Writing peaks failed")
                message(e)
            }) #Trycatch
            ##---------------------------------------------
        }#Checking reference
        cat("Finished calculating signals outside peak region (max_non_peak_region)\n")

        ##==================================================
        tmp_df_np_coor <- within(tmp_df_np_coor,
                     indx <- paste(ID,chr,strand,non_peak_start,non_peak_stop,
                                   sep=""))
        ##-------------------------
                                        #calculate absolute max_non_peak
        tmp_max_np_gene <- max(tmp_df_np_coor$max_non_peak_region)
        ##-------------------------
                                        #create df to merge max_np_gene_gene
                                        #single line tmp_df_np_merge
        df_np_merge <- data.table(chr=tmp_df_np_coor$chr[1],
                                  strand=tmp_df_np_coor$strand[1],
                                  ID=tmp_df_np_coor$ID[1],
                                  max_np_gene=tmp_max_np_gene)
        df_np_total <- rbind(df_np_total,df_np_merge)

        
        ## sink(file.path(pydeg_log_dir,"Log-dups.txt"),
        ##      append=TRUE)
        cat(paste0("Finished ID = ",i.id, "\n"))
        ## sink()
    }#i.id in df_dup_ID_sub
    return(df_np_total)
}


#'Calculate highest read outside peak (peak = 1)
#' @name checkNonPeakRefOne
#' @param i.df data.frame of pyDegradome result
#' @param i.comp Name of comparison pair (divided into test and ctrl)
#' @param cols.indx A vector of column names used to construct the index.
#' @export
checkNonPeakRefOne <- function(df, i.comp, cols.indx) {
    i.comp <<- i.comp
    peaks_ref <- file.path(pydeg_reads_dir,
                           paste0("MaxNonPeak_reads-",
                                  i.comp))
    

    cat("Calculating signals outside peak region (max_np_gene)...\n")
    cat("Processing ", paste0("MaxNonPeak_reads-",i.comp), "\n")
    if(file.exists(peaks_ref)) {
        cat("\tUsing reference file...\n")
        i.reads <- read.table(peaks_ref,header = TRUE, sep = "\t")
        i.reads <- unique(i.reads)
        
      ##Subset entries with ID
      df.ID <- df[grepl("^AT",df$ID),]#!is.na(df$ID)
        if(!is.null(df.ID) && nrow(df.ID) > 0) {
            df.ID <- within(df.ID,
                         indx <- paste(ID,chr,strand,peak_start,peak_stop,
                                       sep=""))
            
            dfs.ID <- MatchDFwithRef(df.ID,i.reads,"max_np_gene")
            df.ID1 <- dfs.ID$df_match
            df.ID2 <- dfs.ID$df_nonmatch
        } else {
            df.ID <- NULL
        }
                                        #Add max_np_gene to the missing entries
        cat("\t\tAdd max_np_gene to missing entries...\n")
        if(!is.null(df.ID2) && nrow(df.ID2) > 0) {
            df.ID2 <- addMaxNonPeakSignal(df.ID2, i.comp, core=env$core)
            ##---------------------------------------------
                                        #Append the new entries
                                        #to the reference file
            i.reads2 <-  dplyr::select(df.ID2,
                                       c(all_of(cols.indx),"max_np_gene"))
            ## i.reads2 <- df.ID2[,c(cols.indx,"max_np_gene")]
            ## i.reads2 <- unique(i.reads2)
            i.reads2 <- within(i.reads2,
                               indx <- paste(ID,chr,strand,peak_start,peak_stop,
                                   sep=""))
            i.reads <- rbind(i.reads,i.reads2)
            i.reads <- unique(i.reads)
            tryCatch({
                ## i.reads <- apply(i.reads,2,as.character)
                                        #remove indx column
                df.ID2 <- dplyr::select(df.ID2,-indx)
                fwrite(i.reads,peaks_ref,sep="\t",row.names=FALSE)
                cat(paste0("\t\t\tWrote ",peaks_ref,"\n"))
            },error=function(e) {
                stopCluster(cl)
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
        df.NonID <- df[!grepl("^AT",df$ID),]#!is.na(df$ID)
        if(nrow(df.NonID) > 0) {
            df.NonID$max_np_gene <- NA
        }else {
            df.NonID <- NULL
        }

                                        #Merge data
        cat("\t\tMerging dataframes...\n")

        df <- rbind(df.ID1,df.ID2,df.NonID)
        df <- unique(df)

    } else {
        cat("\tNo reference was found.\n")
        cat("\t\tCalculating max_np_gene...\n")
        df <- addMaxNonPeakSignal(df,i.comp,core = env$core)
        ##---------------------------------------------
                                        #Write data to ref file
        i.reads <- dplyr::select(df,
                                 c(all_of(cols.indx),"max_np_gene"))
        ## i.reads <- df[,c(cols.indx,"max_np_gene")]
        i.reads <- unique(i.reads)
        i.reads <- within(i.reads,
                               indx <- paste(ID,chr,strand,peak_start,peak_stop,
                                   sep=""))
        tryCatch({
            ## i.reads <- apply(i.reads,2,as.character)
            fwrite(i.reads,peaks_ref,sep="\t",row.names=FALSE)
            cat(paste0("\t\t\tWrote ",peaks_ref,"\n"))
        },error=function(e) {
            stopCluster(cl)
            message("Writing peaks failed")
            message(e)
        }) #Trycatch
        ##---------------------------------------------
    }#Checking reference
    return(df)
    cat("Finished calculating signals outside peak region (max_np_gene)\n")
}



#'Calculate reads for a given region. Checks a reference file to speed up
#' @name checkPeakRef
#' @param df data.frame of pyDegradome result
#' @param i.comp comparison
#' @export
checkPeakRef <- function(df,i.comp,i.cols=cols.indx,
                         idir=pydeg_reads_dir) {
    i.comp <<- i.comp
    df <<- df
    peaks_ref <- file.path(idir,
                           paste0("MaxPeak_reads-",
                                  i.comp))
    cat("Calculating signals in the peak region (max_peak)...\n")
    cat("Processing ", paste0("MaxPeak_reads-",i.comp), "\n")
    if(file.exists(peaks_ref)) {
        
        cat("\tUsing reference file...\n")
        i.reads <- read.table(peaks_ref,header = TRUE, sep = "\t")
        i.reads <- unique(i.reads)
                                        #Subset entries with ID
        df.ID <- df[grepl("^AT",df$ID),]#is.na(df$ID)
        df.ID <- within(df.ID,
                               indx <- paste(ID,chr,strand,peak_start,peak_stop,
                                             sep=""))
        dfs.ID <- MatchDFwithRef(df.ID, i.reads, "max_peak")
        df.ID1 <- unique(dfs.ID$df_match)

                                        #Add max_peak to the missing entries
        cat("\t\tAdd max_peak to missing entries...\n")
        if (!is.null(dfs.ID$df_nonmatch) && nrow(dfs.ID$df_nonmatch) > 0) {
            df.ID2 <- unique(dfs.ID$df_nonmatch)
            df.ID2$max_peak <- MaxPeakBigWig(df.ID2, i.comp,core = env$core)
            ##---------------------------------------------
                                        #Append the new entries
                                        #to the reference file
            ## i.reads2 <- df.ID2[,c(i.cols,"max_peak")]
            i.reads2 <- dplyr::select(df.ID2,
                                      c(all_of(i.cols),"max_peak"))
            ## i.reads2 <- unique(i.reads2)
            i.reads2 <- within(i.reads2,
                               indx <- paste(ID,chr,strand,peak_start,peak_stop,
                                   sep=""))
            i.reads <- rbind(i.reads,i.reads2)
            unique(i.reads)
            tryCatch({
                ## i.reads <- apply(i.reads,2,as.character)
                                        #remove indx column
                df.ID2 <- dplyr::select(df.ID2,-indx)
                fwrite(i.reads,peaks_ref,sep="\t",row.names=FALSE)
                cat(paste0("\t\t\tWrote ",peaks_ref,"\n"))
            },error=function(e) {
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
        df.NonID <- df[!grepl("^AT",df$ID),]#is.na(df$ID)
        if(nrow(df.NonID) > 0) {
            df.NonID$max_peak <- NA
        }else {
            df.NonID <- NULL
        }

                                        #Merge data
        cat("\t\tMerging dataframes...\n")
        df <- rbind(df.ID1,df.ID2,df.NonID)
        df <- unique(df)
    } else {
        cat("\tNo reference was found.\n")
        cat("\t\tCalculating max_peak...\n")
        df$max_peak <- MaxPeakBigWig(df, i.comp, core = env$core)
        ##---------------------------------------------
                                        #Write data to ref file
        i.reads <- dplyr::select(df,
                                 c(all_of(i.cols),"max_peak"))
        ## i.reads <- df[,c(i.cols,"max_peak")]
        i.reads <- unique(i.reads)
        i.reads <- within(i.reads,
                               indx <- paste(ID,chr,strand,peak_start,peak_stop,
                                   sep=""))
        tryCatch({
            ## i.reads <- apply(i.reads,2,as.character)
            fwrite(i.reads,peaks_ref,sep="\t",row.names=FALSE)
            cat(paste0("\t\t\tWrote ",peaks_ref,"\n"))
        },error=function(e) {
            stopCluster(cl)
            message("Writing peaks failed")
            message(e)
        }) #Trycatch
        ##---------------------------------------------
    } #Checking if reference file exists
    return(df)
    cat("Finished calculating in peak region (max_peak)\n")
}


#' Find genes with two peaks based on repetition of the ID value
#' @name SplitDFbyNpeaks
#' @param df Dataframe with peak and gene coordinates
#' @param ref_col Column used to judge duplicates
#' @export
SplitDFbyNpeaks <- function(df,ref_col) {
                                        #get table of duplicated values
  n_occur_id <- data.frame(table(df[[ref_col]]))
  colnames(n_occur_id) <- c("Var1","Freq")
                                        #logical vect of duplicated values
  i.log.dup <- df[[ref_col]] %in% n_occur_id$Var1[n_occur_id$Freq > 1]

  results <- list(df_dup=df[i.log.dup,],
                  df_uniq=df[!i.log.dup,])
  return(results)
}

#'Get the mean value of peak between replicates in the test sample
#' @name maxPtest
#' @param f_df data frame with pydegradome processed output
#' @param f_bigwigs bigwigs List of bigwig records. BigWig files should be loaded by `rtracklayer::import
#' @param f_pairs pairs of genotype B concentration to be compared
#' @param f_input_dir directory of the pooled file
#' @param f_core Number of threads for the nested funtion (GetMaxread)
#' @export
maxPtest <-function(f_df, f_bigwigs, f_pairs, f_input_dir, f_core=1) {
    
    f_cols <- colnames(f_df)
    max.peak.l <- list()
    for (i in seq_along(f_pairs)) {
        
        f_sample <- f_pairs[i]
        bigwig <- f_bigwigs[[f_sample]]
        
        reads_ref <- file.path(f_input_dir,
                               paste0("Peak_reads_test-",
                                      f_sample))
        
        if(file.exists(reads_ref)) {
            cat("\tUsing reference file...\n")
            f_reads <- read.table(reads_ref,header = TRUE, sep = "\t")
            f_df <- within(f_df,
                         indx <- paste(tx_name,chr,strand,peak_start,peak_stop,
                                       sep=""))
            dfs.ID <- MatchDFwithRef(f_df,f_reads,"max_peak")
            df1 <- dfs.ID$df_match
            df2<- dfs.ID$df_nonmatch

            cat("\t\tAdd peak signal to missing entries...\n")
            if(!is.null(df2) && nrow(df2) > 0) {
                
                f_gr <- GRanges(seqnames=df2$chr,
                              ranges=IRanges(start=df2$peak_start,
                                             end=df2$peak_stop),
                              strand=df2$strand)
                df2$max_peak <- GetMaxRead(grA = f_gr, f_bigwig = bigwig, core=f_core)
               
                cat("\t\tFinished peak signal calculation on test sample\n")
                f_reads2 <- df2
                f_reads <- rbind(f_reads,f_reads2)
                f_reads <- unique(f_reads)
                df2 <- dplyr::select(df2,-indx)
                fwrite(f_reads,reads_ref,sep="\t",row.names=FALSE)
                cat(paste0("\t\t\tWrote ",reads_ref,"\n"))
            } else {
                df2 <- NULL
            }
            cat("\t\tMerging dataframes...\n")
            
            df.out <- rbind(df1,df2)
            df.out$sample <- f_sample
            max.peak.l[[f_sample]] <- unique(df.out)
            
        } else {
            cat("\tNo reference was found.\n")
            cat("\t\tCalculating peak signal on test sample...\n")
            f_df <- within(f_df,
                         indx <- paste(tx_name,chr,strand,peak_start,peak_stop,
                                       sep=""))
            f_gr <- GRanges(seqnames=f_df$chr,
                          ranges=IRanges(start=f_df$peak_start,
                                         end=f_df$peak_stop),
                          strand=f_df$strand)
            f_df$max_peak <- GetMaxRead(grA = f_gr, f_bigwig = bigwig, core=f_core)
            
            cat("\t\tFinished peak signal calculation on test sample\n")
            f_reads <- dplyr::select(f_df, all_of(c(f_cols,"max_peak","indx")))
            f_reads <- unique(f_reads)
            
            fwrite(f_reads,reads_ref,sep="\t",row.names=FALSE)
            cat(paste0("\t\t\tWrote ",reads_ref,"\n"))
            
            df.out <- dplyr::select(f_df,-indx)
            df.out$sample <- f_sample
            max.peak.l[[f_sample]] <- unique(df.out)
        }
        cat("Finished calculating peak signals on test sample: ",
            f_sample, "\n")
    }
    return(max.peak.l)
}




#'Get the mean value of the highest signal found in each of the replicates in control samples
#' @name maxTxctrl
#' @param tx_id Vector with transcripts from which the reads should be obtained
#' @param f_bigwigs bigwigs List of bigwig records. BigWig files should be loaded by `rtracklayer::import
#' @param f_pairs pairs of genotype B concentration to be compared
#' @param f_gr_tr subset of the GRanges object tx_range using pydeg_tx
#' @param f_input_dir directory of the pooled file
#' @param f_core Number of threads for the nested function (GetMaxread)
#' @export
maxTxctrl <- function(tx_id, f_bigwigs, f_pairs, f_gr_tr, f_input_dir, f_core=1) {
    i.list <- list()
    for (i in seq_along(f_pairs)) {
        f_sample <- f_pairs[i]
        bigwig <- f_bigwigs[[f_sample]]
        
                                        #create max read df
        pydeg_df <- data.frame(tx_name=tx_id, #1
                               indx=paste0(tx_id,f_sample))
        
        reads_ref <- file.path(f_input_dir,
                               paste0("Max_reads_tx-",
                                      f_sample))
        
        
        if(file.exists(reads_ref)) {

            cat("\tUsing reference file...\n")
            f_reads <- read.table(reads_ref,header = TRUE, sep = "\t")

            ##logical vector of indices in ref file
            i.wread <- pydeg_df$indx %in% f_reads$indx

            dfs.ID <- MatchDFwithRef(pydeg_df, f_reads, "max_read_tx")
            pydeg_df1 <- dfs.ID$df_match
            pydeg_df2 <- dfs.ID$df_nonmatch
            
            cat("\t\tAdd max_read_tx to missing entries...\n")
            if(!is.null(pydeg_df2) && nrow(pydeg_df2) > 0) {
                f_gr  <-  f_gr_tr[!i.wread,]
                pydeg_df2$max_read_tx <- GetMaxRead(grA = f_gr, f_bigwig = bigwig, core=f_core)

                cat("\t\tFinished max_read_tx calculation on control sample\n")
                f_reads2 <- pydeg_df2
                f_reads <- rbind(f_reads,f_reads2)
                f_reads <- unique(f_reads)
                pydeg_df2 <- dplyr::select(pydeg_df2,-indx)
                fwrite(f_reads,reads_ref,sep="\t",row.names=FALSE)
                cat(paste0("\t\t\tWrote ",reads_ref,"\n"))
            } else {
                pydeg_df2 <- NULL
            }
            cat("\t\tMerging dataframes...\n")
            
            pydeg_df <- rbind(pydeg_df1,pydeg_df2)
            pydeg_df$sample <- f_sample
            i.list[[f_sample]] <- unique(pydeg_df)
            
        } else {
            cat("\tNo reference was found.\n")
            cat("\t\tCalculating max_read_tx on control sample...\n")
            f_gr <- f_gr_tr
            pydeg_df$max_read_tx <- GetMaxRead(grA = f_gr, f_bigwig = bigwig, core=env$core)
            
            cat("\t\tFinished max_read_tx calculation on control sample\n")
            f_reads <- pydeg_df
            f_reads <- unique(f_reads)
            
            fwrite(f_reads,reads_ref,sep="\t",row.names=FALSE)
            cat(paste0("\t\t\tWrote ",reads_ref,"\n"))
            
            pydeg_df <- dplyr::select(pydeg_df,-indx)
            pydeg_df$sample <- f_sample
            i.list[[f_sample]] <- unique(pydeg_df)
        }#Test if there is reference file
        cat("Finished calculating signals along the transcript for control sample: ",
            f_sample, "\n")
        
    }#for i (loop over replicates)
    return(i.list)
}
