

# seqArchR
<!-- badges: start -->
[![R build status](https://github.com/snikumbh/seqArchR/workflows/R-CMD-check/badge.svg)](https://github.com/snikumbh/seqArchR/actions)
[![codecov](https://codecov.io/gh/snikumbh/seqArchR/branch/main/graph/badge.svg?token=NEjCGuOUlW)](https://codecov.io/gh/snikumbh/seqArchR)
<!-- badges: end -->

Note: _This package is currently under development. So, please bear with me while I put the final blocks together. Thanks for your understanding!_ 


seqArchR is an unsupervised, non-negative matrix factorization (NMF)-based algorithm for discovery of sequence architectures de novo.
Below is a schematic of seqArchR's algorithm.


<img src="reference/figures/seqArchR_algorithm_1080p_cropped.gif" width="550" align="center">


## Installation

### Python scikit-learn dependency
This package requires the Python module scikit-learn. See installation instructions [here](https://scikit-learn.org/stable/install.html).


### To install this package, use 

```r
if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")   
}

remotes::install_github("snikumbh/seqArchR", build_vignettes = FALSE)
``` 

If the above command produces errors, see this troubleshooting [note](https://github.com/snikumbh/seqArchR#troubleshooting-installation) below.


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
        n_cores = 4,
        n_iterations = 100,
        k_min = 1,
        k_max = 20,
        mod_sel_type = "stability",
        bound = 10^-8,
        inner_chunk_size = 500,
        flags = list(debug = FALSE, time = TRUE, verbose = TRUE,
                     plot = FALSE)

#
### Call/Run seqArchR
seqArchRresult <- seqArchR::seqArchR(config = seqArchRconfig,
                               seqs_ohe_mat = inputSeqsMat,
                               seqs_raw = inputSeqsRaw,
                               seqs_pos = positions,
                               threshold_itr = 2)

```

## Troubleshooting Installation

1. Error during installation (with `remotes::install_github` or `devtools::install_github`. If the error occurs during `** testing if installed package can be loaded`, try adding the `--no-test-load` option as `remotes::install_github("snikumbh/seqArchR", INSTALL_opts = c("--no-test-load"))`


# Contact
Comments, suggestions, enquiries/requests are welcome! Feel free to email sarvesh.nikumbh@gmail.com or [create an new issue](https://github.com/snikumbh/seqArchR/issues/new)
