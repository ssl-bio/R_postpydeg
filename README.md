# R\_postpydeg


## Description

Different functions to quantify the number of reads in bigwig files. All of them are based on the function `getMaxRead` which returns the highest read value in a given interval (defined as one or two GRanges objects, grA and grB).

```R
## Based on http://adomingues.github.io/2016/11/13/max-coverage-in-bigwigs/

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
```

The different functions arise depending on the number or intervals and the meaning of the defined interval. Regarding the latter, it could represent a *peak* identified by the script `PyDegradome` (<a href="#citeproc_bib_item_1">Gaglia, Rycroft, and Glaunsinger 2015</a>) or a region outside it. For a given gene, if a single peak is present there will be two outside regions (one upstream and one downstream) but when two or more peaks are present then the number of outside regions will increase. This is relevant when trying to obtain the highest read value outside the peak region(s) and this results in the need of additional functions.

Because querying a bigwig file is a relatively slow process, every time an interval is read it creates an entry in a reference file which will be read on a next round in order to avoid re-reading the bigwig file. This is particularly useful if different PyDegradome settings are being tested as this results in a number of identical regions from each test.

Ultimately, values obtained from the different functions are used to calculate ratios between the read found within peaks and those outside it; these in turn are used to classify the peaks into different categories. For further details are described in a blog [post](https://ssl-blog.netlify.app/posts/degradome-analysis/degradome-code/#peak-annotation-and-computing-of-ratios)


## Installation

Open an R session in the cloned directory and run the following commands.

```R
library(devtools)
devtools::document()
devtools::build()
devtools::install()
```

Note that it requires the library `devtools` to be previously installed. Other dependencies to run the analysis are described [here](https://ssl-blog.netlify.app/posts/degradome-analysis/degradome-code/#r-package-installation)


## References
  <div class="csl-entry"><a id="citeproc_bib_item_1"></a>Gaglia, Marta Maria, Chris H. Rycroft, and Britt A. Glaunsinger. 2015. “Transcriptome-Wide Cleavage Site Mapping on Cellular mRNAs Reveals Features Underlying Sequence-Specific Cleavage by the Viral Ribonuclease SOX.” Edited by Pinghui Feng. <i>PLOS Pathogens</i> 11 (12): e1005305. doi:<a href="https://doi.org/10.1371/journal.ppat.1005305">10.1371/journal.ppat.1005305</a>.</div>
</div>
