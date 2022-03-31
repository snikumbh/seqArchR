# @title One-hot encode
#
# @description One-hot encode a given DNA sequence.
#
# @param givenSeq A single DNA sequence as character/string.
#
# @return The one-hot encoded sequence.
#
# This func is a combined func instead of one for each sinuc, dinuc and trinuc
.one_hot_encode <- function(givenSeq, k = 1){
    if(k > 3) stop("Only mono-, di- and trinucleotides are supported")
    dna_alph <- Biostrings::DNA_BASES
    use_alph <- switch(k,
                   dna_alph,
                   get_dimers_from_alphabet(dna_alph),
                   get_trimers_from_alphabet(dna_alph)
                )
    seqlen <- length(givenSeq)
    if(seqlen < 1) stop("Empty or NULL sequences")
    use_colnames <- .get_feat_names(alph = dna_alph, k = k, seqlen = seqlen)
    ohe <- .form_ohe(k = k, givenSeq, use_alph, seqlen, use_colnames)
}
## =============================================================================

.form_ohe <- function(k = 1, gseq, alph, seqlen, set_colnames){
    ##
    ohe <- matrix(rep(0, length(alph) * seqlen), nrow = 1, byrow = TRUE)
    ##
    modseq <- unlist(lapply(seq_len(seqlen - (k-1)), function(x) {
        retVal <- gseq[x]
        if(k == 2) retVal <- paste0(gseq[x], gseq[x + 1])
        if(k == 3) retVal <- paste0(gseq[x], gseq[x + 1], gseq[x + 2])
        retVal
    }))
    ##
    for (i in seq_along(alph)) {
        ohe[, (i - 1) * seqlen + which(modseq == alph[i])] <- 1
    }
    colnames(ohe) <- set_colnames
    ohe
}
## =============================================================================

.get_feat_names <- function(alph = c('A', 'C', 'G', 'T'), k=1, seqlen){
    stopifnot(k<=3)
    if(k==1) kmers <- alph
    if(k==2) kmers <- get_dimers_from_alphabet(alph)
    if(k==3) kmers <- get_trimers_from_alphabet(alph)
    use_feat_names <- paste(rep(kmers, each = seqlen), seq_len(seqlen),
                            sep=".")
    use_feat_names
}
## =============================================================================



# @title One-hot decode
#
# @description One-hot decode a given one-hot encoded DNA sequence.
#
# @param oneHotEncodedSeqV A single one-hot encoded sequence vector.
#
# @return The one-hot decoded sequence of ACGTs.
#
.one_hot_decode <- function(oneHotEncodedSeqV) {

    dna_alphabet <- c("A", "C", "G", "T")
    seqlen <- length(oneHotEncodedSeqV)/length(dna_alphabet)
    decodedSeq <- rep("Q", seqlen)
    for (alpha_char in seq_along(dna_alphabet)) {
        cutp <- seqlen
        startp <- (alpha_char * cutp) - cutp + 1
        endp <- (alpha_char * cutp)
        decodedSeq[which(oneHotEncodedSeqV[startp:endp] == 1)] <-
            dna_alphabet[alpha_char]
    }
    return(paste0(decodedSeq, collapse = ""))
}
## =============================================================================

#' @title Get one-hot encoded sequences
#'
#' @description Get the one-hot encoding representation of the given sequences.
#'
#' @param seqs A \code{\link[Biostrings]{DNAStringSet}} object holding the
#' given DNA sequences
#' @param sinuc_or_dinuc character string, 'sinuc' or 'dinuc' to select for
#' mono- or dinucleotide profiles.
#'
#' @return A sparse matrix of sequences represented with one-hot-encoding
#' @family input functions
#' @seealso \code{\link{prepare_data_from_FASTA}} for generating one-hot
#' encoding of sequences from a FASTA file
#' @importFrom Matrix Matrix
#'
#' @examples
#'
#' fname <- system.file("extdata", "example_data.fa.gz",
#'                         package = "seqArchR", mustWork = TRUE)
#'
#'
#' rawSeqs <- prepare_data_from_FASTA(fasta_fname = fname,
#'                         raw_seq = TRUE)
#'
#' seqs_dinuc <- get_one_hot_encoded_seqs(seqs = rawSeqs,
#'                                        sinuc_or_dinuc = "dinuc")
#'
#' @export
get_one_hot_encoded_seqs <- function(seqs, sinuc_or_dinuc = "sinuc") {
    #
    if(is.null(seqs) || length(seqs) == 0){
        stop("Empty or NULL sequences")
    }
    seqs_split_as_list <-
        base::strsplit(as.character(seqs), split = NULL)
    if (length(seqs_split_as_list) > 0) {
        if (sinuc_or_dinuc == "sinuc") {
            .msg_pstr("Generating dinucleotide profiles", flg=TRUE)
            encoded_seqs <- lapply(seqs_split_as_list,
                                    .one_hot_encode, k = 1)
        } else if (sinuc_or_dinuc == "dinuc") {
            .msg_pstr("Generating dinucleotide profiles", flg=TRUE)
            encoded_seqs <-
                lapply(seqs_split_as_list, .one_hot_encode, k = 2)
        }  else if (sinuc_or_dinuc == "trinuc") {
            .msg_pstr("Generating trinucleotide profiles", flg=TRUE)
            encoded_seqs <-
                lapply(seqs_split_as_list, .one_hot_encode, k = 3)
        }

        encoded_seqs <- do.call(rbind, encoded_seqs)
        encoded_seqs <- t(encoded_seqs)
        ## Use Matrix package, store as sparse matrix
        ## Helps reduce object size and computation time
        encoded_seqs <- Matrix::Matrix(encoded_seqs, sparse = TRUE)
        ##
        return(encoded_seqs)
    }
}
## =============================================================================

# @title Assert attributes of sequences
#
# @description Assert the attributes of the sequences provided. This includes
# checking for (1) the length of the sequences, (2) characters in the
# sequences.
#
# @param givenSeqs DNA sequences as a list.
#
#
# @return nothing. Only prints a warning to the screen.
# @importFrom Biostrings width
#
.assert_seq_attributes <- function(givenSeqs) {
    # Check that all sequences are of same length
    seqs_split_as_list <-
        base::strsplit(as.character(givenSeqs), split = NULL)
    length_vals <- levels(as.factor(Biostrings::width(givenSeqs)))
    char_levels <- levels(as.factor(unlist(seqs_split_as_list)))
    dna_alphabet <- c("A", "C", "G", "T")
    if (0 %in% length_vals) {
        # Checking sequences of length 0
        stop("Found ", which(0 == length_vals),
                " sequence(s) of length zero\n")
        #
    } else if (length(levels(as.factor(length_vals))) > 1) {
        # Checking all sequences are of same length
        stop("Sequences are of different length\n")
        #
    } else if (any(!(char_levels %in% dna_alphabet))) {
        # Check for non-alphabet characters
        # Raise either an error or just warn!
        warning(c("Non DNA-alphabet character in the sequences: ",
                char_levels), immediate. = TRUE)

    } else {
        # All OK!
        .msg_pstr("Sequences OK, ", levels(length_vals)[1], flg=TRUE)
    }
}
## =============================================================================

#' @title
#' Generate one-hot encoding of sequences given as FASTA file
#'
#' @description
#' Given a set of sequences in a FASTA file this function returns a sparse
#' matrix with one-hot encoded sequences.
#' In this matrix, the sequence features are along rows, and sequences along
#' columns. Currently, mono- and dinucleotide features for DNA sequences are
#' supported. Therefore, the length of the feature vector is 4 and 16 times
#' the length of the sequences (since the DNA alphabet is four characters)
#' for mono- and dinucleotide features respectively.
#'
#' @param fasta_fname Provide the name (with complete path) of the input
#' FASTA file.
#'
#' @param raw_seq TRUE or FALSE, set this to TRUE if you want the raw sequences.
#'
#' @param sinuc_or_dinuc character string, 'sinuc' or 'dinuc' to select for
#' mono- or dinucleotide profiles.
#'
#' @return A sparse matrix of sequences represented with one-hot-encoding.
#' @family input functions
#' @seealso \code{\link{get_one_hot_encoded_seqs}} for directly using a
#' DNAStringSet object
#' @importFrom Biostrings DNAStringSet
#'
#' @examples
#'
#' fname <- system.file("extdata", "example_data.fa.gz",
#'                         package = "seqArchR", mustWork = TRUE)
#'
#' # mononucleotides feature matrix
#' rawSeqs <- prepare_data_from_FASTA(fasta_fname = fname,
#'                         sinuc_or_dinuc = "sinuc")
#'
#' # dinucleotides feature matrix
#' rawSeqs <- prepare_data_from_FASTA(fasta_fname = fname,
#'                         sinuc_or_dinuc = "dinuc")
#'
#' # FASTA sequences as a Biostrings::DNAStringSet object
#' rawSeqs <- prepare_data_from_FASTA(fasta_fname = fname,
#'                         raw_seq = TRUE)
#'
#' @export
prepare_data_from_FASTA <- function(fasta_fname, raw_seq = FALSE,
                                    sinuc_or_dinuc = "sinuc") {
    if (!file.exists(fasta_fname)) {
        stop("File not found, please check if it exists")
    }

    if (raw_seq) {
        givenSeqs <- Biostrings::readDNAStringSet(
            filepath = fasta_fname, format = "fasta", use.names = TRUE)
        givenSeqs <- Biostrings::DNAStringSet(toupper(givenSeqs))
        return(givenSeqs)
    } else {
        #
        givenSeqs <- Biostrings::readDNAStringSet(
            filepath = fasta_fname, format = "fasta", use.names = FALSE)
        givenSeqs <- Biostrings::DNAStringSet(toupper(givenSeqs))
        .assert_seq_attributes(givenSeqs)
        .msg_pstr("Read ", length(givenSeqs), " sequences", flg=TRUE)
        #
        oheSeqs <- get_one_hot_encoded_seqs(givenSeqs,
                                            sinuc_or_dinuc = sinuc_or_dinuc)
        return(oheSeqs)
    }
}
## =============================================================================

