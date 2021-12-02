

# seqArchR
<!-- badges: start -->
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![DOI](https://zenodo.org/badge/188449833.svg)](https://zenodo.org/badge/latestdoi/188449833)
[![Build status](https://travis-ci.org/snikumbh/seqArchR.svg?branch=master)](https://travis-ci.org/snikumbh/seqArchR)
[![Codecov test coverage](https://codecov.io/gh/snikumbh/seqArchR/branch/master/graph/badge.svg)](https://codecov.io/gh/snikumbh/seqArchR?branch=master)
[![R build status](https://github.com/snikumbh/seqArchR/workflows/R-CMD-check/badge.svg)](https://github.com/snikumbh/seqArchR/actions)
<!-- badges: end -->

Note: _This package is currently under development. So, please bear with me while I put the final blocks together. Thanks for your understanding!_ 


seqArchR is an unsupervised, non-negative matrix factorization (NMF)-based algorithm for discovery of sequence architectures de novo.
Below is a schematic of seqArchR's algorithm.

<img src="https://github.com/snikumbh/seqArchR/blob/master/vignettes/archR_algorithm_1080p_cropped.gif" width="550" align="center">


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
archRresult <- seqArchR::seqArchR(config = archRconfig,
                               seqs_ohe_mat = inputSeqsMat,
                               seqs_raw = inputSeqsRaw,
                               seqs_pos = positions,
                               total_itr = 2,
			       set_ocollation = c(TRUE, FALSE))

```


# Contact
Comments, suggestions, enquiries/requests are welcome! Feel free to email sarvesh.nikumbh@gmail.com or [create an new issue](https://github.com/snikumbh/seqArchR/issues/new)
