context("Auxiliary Functions I")


test_that("fetch cutree by hcorder ordering", {

    hcorder <- c(6, 7, 1, 10, 4, 2, 5, 8, 9, 3)
    clust_list <- list(sort(hcorder[1:3]), sort(hcorder[4:8]), sort(hcorder[9]),
                       sort(hcorder[10]))
    expect_identical(fetch_cutree_by_hc_order(clust_list = clust_list,
                                              hcorder = hcorder),
                     list(c(6,7,1), c(10, 4, 2, 5, 8), c(9), c(3)))
})

test_that("collate_clusters handles empty globClustAssignments", {
    Obj <- list(clustering = list(k = 5, sizes = rep(5,5), order = 1:25))
    tempList <- vector("list", 5)
    globClustAssignments <- lapply(seq_along(tempList),
                                   function(x){
                                       if ( x != 3) {
                                            tempList[[x]] <- round(10*runif(5))
                                       }
                                   })
    expect_error(collate_clusters(Obj, globClustAssignments),
                 "Cluster assignments variable has a 0-length entry")
})

test_that("unfurl_nodeList handles example nodeList well (vec version)", {
    nodeList <- list(1, c(2, 5), 6, 7)
    unfurled_nodeList <- list(1, 2, 5, 6, 7)
    ans <- .unfurl_nodeList(nodeList = nodeList, vec_ver=TRUE)
    expect_identical(unfurled_nodeList, ans)

})

test_that("unfurl_nodeList handles example nodeList well", {
    nodeList <- list(c(1,3,4), list(2, 5), c(6,7))
    unfurled_nodeList <- list(c(1,3,4), 2, 5, c(6,7))
    ans <- .unfurl_nodeList(nodeList = nodeList, vec_ver=FALSE)
    expect_identical(unfurled_nodeList, ans)

})


test_that("unfurl_nodeList handles empty/null nodeList", {
    nodeList <- NULL
    expect_error(.unfurl_nodeList(nodeList = nodeList, vec_ver=FALSE))
    nodeList <- list()
    expect_error(.unfurl_nodeList(nodeList = nodeList, vec_ver=FALSE))
})




test_that("update_cluster_labels handles empty collatedClustAssignments", {

    toyClustLabels <- c("9", "22", "20", "23", "12", "22", "1", "16", "20",
                        "9", "21", "10", "18", "8", "13", "1", "17", "18",
                        "11", "14", "11", "9", "25", "3", "12"
                        # "16", "1",
                        # "19", "15", "2", "21", "10", "6", "22", "8", "6",
                        # "25", "10", "12", "5", "20", "17", "7", "14", "8",
                        # "21", "7", "18", "7", "4", "2", "14", "4", "23", "4",
                        # "24", "13", "3", "15", "5", "13", "19", "5", "11", "23",
                        # "24", "24", "15", "16", "17", "2", "6", "19", "25", "3"
                        )

    collection <- as.numeric(c("3", "15", "24", "14", "19", "13", "1", "5", "12", "9",
                    "11", "8", "4", "17", "20", "16", "25", "10", "2", "7",
                    "6", "18", "21", "22", "23"))
    tempList <- vector("list", 5)
    collatedClustAssignmentsErr <- lapply(seq_along(tempList),
                                   function(x){
                                       if ( x != 3) {
                                           tempList[[x]] <- round(10*runif(5))
                                       }
                                   })

    idx <- rep(1:5,5)
    collatedClustAssignments <- lapply(seq_along(tempList),
                                   function(x){
                                       tempList[[x]] <- collection[idx == x]
                                   })
    ## this tests an assertion
    expect_error(.update_cluster_labels(toyClustLabels,
                                        collatedClustAssignmentsErr),
                 "Cluster assignments variable has a 0-length entry")
    ####
    foo <- c("2", "4", "1", "3", "3", "1", "5", "2", "5", "3", "1", "4",
                 "1", "4", "2", "1", "4", "2", "5", "5", "3", "4", "5", "3", "2")
    ## this tests the updated cluster labels
    ## -- the updated labels should have 5 unique cluster labels since
    ## globClustAssignments has 5 clusters
    ## --

    foo_ans <- .update_cluster_labels(toyClustLabels, collatedClustAssignments)

    expect_equal(foo_ans, foo)

})

test_that("get_seqs_clust_list works", {
    seqsClustLabels <- c("2", "4", "1", "3", "3", "1", "5", "2",
                       "5", "3", "1", "4", "1", "4", "2", "1",
                       "4", "2", "5", "5", "3", "4", "5", "3", "2")
    ansList <- seqArchR::get_seqs_clust_list(seqsClustLabels)
    expect_length(ansList, 5)
    })



# test_that("update_cluster_labels handles inconsistent #sequences", {
#     seqsClustLabels <- rep("0-1-2-3", 20)
#     tempList <- vector("list", 5)
#     globClustAssignments <- lapply(seq_along(tempList),
#                                    function(x){
#                                     tempList[[x]] <- round(10*runif(5))
#                                    })
#     expect_error(.update_cluster_labels(seqsClustLabels, globClustAssignments),
#     "Number of sequences in seqsClustLabels and clustAssignments not equal")
# })



test_that("get seq clusters as list works fine", {
    set.seed(1234)
    toyClustLabels <- sample(x = rep(as.character(1:14),3), size = 42, replace = FALSE)
    seqsClustLabels <- toyClustLabels
    expAns <- 14
    expect_equal(length(get_seqs_clust_list(seqsClustLabels)), 14)
    expect_error(get_seqs_clust_list(NULL))
})




test_that("prepare_chunks handles negative chunkSize", {
    seqsClustLabels <- rep("0-1-2-3", 20)
    tempList <- vector("list", 5)
    globClustAssignments <- lapply(seq_along(tempList),
                                   function(x){
                                       tempList[[x]] <- round(10*runif(5))
                                   })
    expect_error(.prepare_chunks(seqsClustLabels, -25),
                 "'chunk_size' should be > 0")
})





test_that("clustering of factors handles all-zero featuresMatrix", {

    fMat <- matrix(rep(0,1000), ncol = 5)
    expect_error(.handle_clustering_of_factors(fMat),
                   "All zeroes as factors")
})


test_that("handle_clustering_of_factors handles NA in featuresMatrix", {

    fMat <- matrix(rep(0,1000), ncol = 5)
    fMat[2,2] <- NA
    expect_error(.handle_clustering_of_factors(fMat),
                    "Factors have NA")
})

test_that("handle_clustering_of_factors handles improper flags", {

    fMat <- matrix(rep(0,1000), ncol = 5)
    fMat[2,2] <- NA
    expect_error(.handle_clustering_of_factors(fMat, flags = list()),
                    "")
})




#
# test_that("Q^2 computation: handling original non-matrix", {
#     matA <- c(rnorm(20)) # matrix(rnorm(20), nrow = 4)#rnorm(20), nrow = 4)
#     reMatA <- matrix(rnorm(20), nrow = 4)
#     expect_error(.compute_q2(matA, reMatA), "not of type matrix")
# })
#
# test_that("Q^2 computation: handling reconstructed non-matrix", {
#     matA <- matrix(rnorm(20), nrow = 4)
#     reMatA <- c(rnorm(20)) # matrix()
#     expect_error(.compute_q2(matA, reMatA), "not of type matrix")
# })

# test_that("NMF solved correctly", {
#
#
#
# })
