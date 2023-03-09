

# seqArchR
<!-- badges: start -->
[![DOI](https://zenodo.org/badge/188449833.svg)](https://zenodo.org/badge/latestdoi/188449833)
[![codecov](https://codecov.io/gh/snikumbh/seqArchR/branch/main/graph/badge.svg?token=NEjCGuOUlW)](https://codecov.io/gh/snikumbh/seqArchR)
<!-- bioc badges: start -->
  [![Bioc release status](http://www.bioconductor.org/shields/build/release/bioc/seqArchR.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/seqArchR)
  [![Bioc downloads rank](https://bioconductor.org/shields/downloads/release/seqArchR.svg)](http://bioconductor.org/packages/stats/bioc/seqArchR/)
  [![Bioc support](https://bioconductor.org/shields/posts/seqArchR.svg)](https://support.bioconductor.org/tag/seqArchR)
  [![Bioc history](https://bioconductor.org/shields/years-in-bioc/seqArchR.svg)](https://bioconductor.org/packages/release/bioc/html/seqArchR.html#since)
  [![Bioc dependencies](https://bioconductor.org/shields/dependencies/release/seqArchR.svg)](https://bioconductor.org/packages/release/bioc/html/seqArchR.html#since)
  <!-- bioc badges: end -->
  
<!-- [![R build status](https://github.com/snikumbh/seqArchR/workflows/R-CMD-check/badge.svg)](https://github.com/snikumbh/seqArchR/actions) -->
<!-- badges: end -->


seqArchR is an unsupervised, non-negative matrix factorization (NMF)-based algorithm for discovery of sequence architectures de novo.
Below is a schematic of seqArchR's algorithm.

<img src="https://github.com/snikumbh/seqArchR/blob/main/vignettes/seqArchR_algorithm_1080p_cropped.gif" width="550" align="center">


## Installation

### Python scikit-learn dependency
This package requires the Python module scikit-learn. Please see installation instructions [here](https://scikit-learn.org/stable/install.html).


### To install this package, use 

```r
if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")   
}

remotes::install_github("snikumbh/seqArchR", build_vignettes = FALSE)
``` 



### Usage
```r
# load package
library(seqArchR)
library(Biostrings)


# Creation of one-hot encoded data matrix from FASTA file
# You can use your own FASTA file instead
inputFastaFilename <- system.file("extdata", "example_data.fa", 
                                  package = "seqArchR", 
                                  mustWork = TRUE)

# Specifying dinuc generates dinucleotide features
inputSeqsMat <- seqArchR::prepare_data_from_FASTA(inputFastaFilename,
                                                  sinuc_or_dinuc = "dinuc")

inputSeqsRaw <- seqArchR::prepare_data_from_FASTA(inputFastaFilename, 
                                               raw_seq = TRUE)

nSeqs <- length(inputSeqsRaw)
positions <- seq(1, Biostrings::width(inputSeqsRaw[1]))

# Set seqArchR configuration
# Most arguments have default values
seqArchRconfig <- seqArchR::set_config(
        parallelize = TRUE,
        n_cores = 2,
        n_runs = 100,
        k_min = 1,
        k_max = 20,
        mod_sel_type = "stability",
        bound = 10^-6,
        chunk_size = 100,
	result_aggl = "ward.D",
	result_dist = "euclid",
        flags = list(debug = FALSE, time = TRUE, verbose = TRUE,
                     plot = FALSE)
        )

#
### Call/Run seqArchR
seqArchRresult <- seqArchR::seqArchR(config = seqArchRconfig,
                               seqs_ohe_mat = inputSeqsMat,
                               seqs_raw = inputSeqsRaw,
                               seqs_pos = positions,
                               total_itr = 2,
			       set_ocollation = c(TRUE, FALSE))

```


# Contact
Comments, suggestions, enquiries/requests are welcome! Feel free to email sarvesh.nikumbh@gmail.com or [create an new issue](https://github.com/snikumbh/seqArchR/issues/new)
