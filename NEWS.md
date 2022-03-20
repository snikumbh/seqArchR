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

# archR 0.1.8

## New features
* Custom ordering of individual clusters within collated clusters. Earlier this 
followed label ordering as per `cutree` output
* Improved messages for model selection using cross-validation
* (User-facing) The agglomeration and distance method used for collation of 
clusters can be set as part of configuration settings

## Breaking changes
* (User-facing) Argument name `inner_chunk_size` changed to `chunk_size` in 
  `archR_set_config()`


# archR 0.1.7

## New features
* Accounting and correcting for over-fitting improves clusters/architectures' 
quality
* Improved messages
* (User-facing) New function `collate_clusters()`

## Breaking changes
* (User-facing) tolerance no longer an argument in archR_set_configuration


# archR 0.1.6

## New features
* Can choose to fix aspect ratio for visualization functions
  - `plot_ggheatmap()`
  - `plot_ggseqlogo()`
  - `viz_bas_vec_heatmap_seqlogo()`
  - `viz_bas_vec_heatmap()`
  - `viz_bas_vec_seqlogo()`
  - `plot_arch_for_clusters()`
  - `plot_ggseqlogo_of_seqs()`
* Ability to detect occurrence of just for sake clustering when fetching from 
  HAC

# archR 0.1.5

## New features
* (User-facing) Getter functions to fetch archR result object elements.
  - `seqs_str()` to fetch DNA sequences for a particular cluster, iteration etc.
  - `get_clBasVec()` to fetch basis vectors element at given iteration.
  - `get_clBasVec_k()` to fetch number of basis vectors/clusters.
  - `get_clBasVec_m()` to fetch basis vectors.
  - `get_seqClLab()` to fetch cluster Ids for each sequence.
  
## Breaking changes
* (User-facing) Arguments of all functions consistently use underscores.


# archR 0.1.2 

## Bug fixes
* Fix staged installation error. This enables seamless installation on macOS.

## New features
* CI with GHA.
* (User-facing) archR result list has new elements:
  - `timeInfo` if timeFlag is TRUE;
  - `clustSol` storing combined clusters from the last iteration of archR.
* (User-facing) Architectures/clusters sequence logos can be plotted with 
`auto` y-axis limits for information content (0-max instead of 0-2).
* (User-facing) TFBSTools moved from dependency to suggests. This eases 
installation of archR by reducing possibility of failed installation.


# archR 0.99.3
## New features
* Accepts FASTA sequences as DNAStringSet object(s)

