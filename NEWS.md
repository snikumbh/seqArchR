# seqArchR 1.1.1

* (User-facing) Fixed: Too many warnings from Python's scikit-learn were 
  being printed to console (warning about an expected change in future 
  scikit-learn version)

# seqArchR 0.99.339

* seqArchR available on Bioconductor

# seqArchR 0.99.0

## New features
* `viz_seqs_acgt_mat()` function can now add a legend via new arguments 
   add_legend and use_legend.
* Package name change from `archR` to `seqArchR`


## Breaking changes
* (User-facing) Function name `archR_set_config()` changed to `set_config()`
* (User-facing) Function name `viz_seqs_acgt_mat_from_seqs()` changed to 
  `viz_seqs_acgt_mat()`
* (User-facing) Function `runArchRUI()` that launched a Shiny app is moved to 
  a different package coming up in the future
