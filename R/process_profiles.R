#' Omit last record from a Granges object considering the strandness.
#' @name omitLastExon
#' @param x `Granges` object of exons in a transcript
#' @return `Granges` object with lenght of `length(x) - 1`, after omiting the last (3'end in transcript) record in the transcript
omitLastExon <- function(x){
  ret <- x[-which.max(x$exon_rank)]
  return(ret)
}


#' Calculate mean of relative distribution from a matrix or data.frame of counts(position x count).
#' @name getRelativeMeanDist
#' @param x Numeric data.frame or matrix.
#' @param thr Numeric. Threshold for total count. Records whose total counts are less than `thr` will be excluded from the calculation.
#' @return Numeric vector of normalized mean count. The last value is the number of exon counted.
#' @export
getRelativeMeanDist <- function(x,thr = 1) {
  total_count <- apply(x,1,sum)
  # Omit exon without signal
  sel <- total_count >= thr
  x <- x[sel,]
                                        #Test if matrix has more than one row
                                        #if not return a null value
  i.test <- dim(x)
	     if(!is.null(i.test))
	     {
                 x <- sweep(x,MARGIN = 1,FUN = "/",apply(x,1,sum))

                                        # Add number of exon counted to the last
                 ret <- c(apply(x,2,mean),sum(sel))
                 return(ret)
             }else{
                 return(NULL)}
}

#' Variation of getRelativeMeanDist. In addition to calculate the mean distribution it also outputs the raw counts as a data frame. Information about the transcript is collected. This was written to get the transcript id for those peaks of interest
#' @name getRelativeMeanDistList 
#' @param x Numeric data.frame or matrix.
#' @param thr Numeric. Threshold for total count. Records whose total counts are less than `thr` will be excluded from the calculation.
#' @return List: First element, numeric vector of normalized mean count. The last value is the number of exon counted. Second element, data frame or raw counts that have at least one read over the range (-50:50)
#' @export
getRelativeMeanDistList <- function(x, l=-50, thr = 1 ) {
    total_count <- apply(x, 1, sum)
    sel <- total_count >= thr
    x <- x[sel, ]
    i.test <- dim(x)
    if(!is.null(i.test))
    {
        y <- as.data.frame(x)
        colnames(y) <- seq(l,dim(x)[2]-(abs(l)+1))
        y$txidcoor <- rownames(y)
        y <- reshape2::melt(y, id.var="txidcoor", variable.name="pos")
        y$pos <- as.numeric(as.character(y$pos))
        y <- y[order(y$txidcoor),]
        z <- sweep(x, MARGIN = 1, FUN = "/", apply(x, 1, sum))
        means <- c(apply(z, 2, mean), sum(sel))
        i.list <- list(means=means,df=y)
        return(i.list)
    }else{return(NULL)}
}


#' Variation of getRelativeMeanDist. It calculates the mean distribution while checking that the input matrix is not empty
#' @name getRelativeMeanDist2
#' @param x Numeric data.frame or matrix.
#' @param thr Numeric. Threshold for total count. Records whose total counts are less than `thr` will be excluded from the calculation.
#' @return List of numeric vectors of normalized mean counts. The last value is the number of exon counted.
#' @export
getRelativeMeanDist2 <- function (x, thr = 1 ) {
    total_count <- apply(x, 1, sum)
    sel <- total_count >= thr
    x <- x[sel, ]
    i.test <- dim(x)
    if(!is.null(i.test))
    {
        z <- sweep(x, MARGIN = 1, FUN = "/", apply(x, 1, sum))
        means <- c(apply(z, 2, mean), sum(sel))
        return(means)
    }else{return(NULL)}
}

#' Variation of getRelativeMeanDist. It also outputs the raw counts as a data frame. Information about the transcript is collected. This was written to get the transcript id for those peaks of interest
#' @name getAbsCountsDistList 
#' @param x Numeric data.frame or matrix.
#' @param thr Numeric. Threshold for total count. Records whose total counts are less than `thr` will be excluded from the calculation.
#' @return List: data frames of raw counts that have at least one read over the range (-50:50)
#' @export
 getAbsCountsDistList <-
function (x, l = -50, thr = 1) 
{
    total_count <- apply(x, 1, sum)
    sel <- total_count >= thr
    x <- x[sel, ]
    i.test <- dim(x)
    if (!is.null(i.test)) {
        y <- as.data.frame(x)
        colnames(y) <- seq(l, dim(x)[2] - (abs(l) + 1))
        y$txidcoor <- rownames(y)
        y <- reshape2::melt(y, id.var = "txidcoor", variable.name = "pos")
        y$pos <- as.numeric(as.character(y$pos))
        y <- y[order(y$txidcoor), ]
        return(y)
    }
    else {
        return(NULL)
    }
}

#' Extract profiles of specific regions from a list of matrixes of profiles.
#' @param psite_prof  `psite_prof` object (list[by transcript] of matrix[read count of position x sample]) from Ribo-seq analysis.
#' @param gr `Granges` to define the area of interest.
#' @name getProfileFromGrange
#' @return List[by region] of matrix[read count of position x sample].
#' @importFrom pbapply pbmapply
#' @import GenomicRanges
#' @export
getProfileFromGrange <- function(psite_prof, gr){

  ret <- pbmapply(FUN = function(seq,start,end){psite_prof[[seq]][start:end,]},
                  as.character(seqnames(gr)),
                  start(gr),
                  end(gr),
                  SIMPLIFY = F

  )
  names(ret) <- paste(as.character(seqnames(gr)),start(gr),end(gr), sep = "_")
  ret <- ret[!sapply(ret,is.null)]

  return(ret)
}
