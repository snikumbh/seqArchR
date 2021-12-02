# Author: snikumbh

# Perform NMF in Python
# A barebones function
# 


from sklearn.decomposition import NMF


def perform_nmf_func(givenMat, nPatterns, nIter=200, givenAlpha = 0, 
                                    givenL1_ratio = 1, seed_val = 3456):
    #
    model = NMF(n_components=nPatterns,
                solver='cd',
                init='nndsvd', 
                alpha=givenAlpha, 
                max_iter=nIter,
                tol=1e-03,
                l1_ratio=givenL1_ratio, 
                random_state=seed_val)
    #
    W = model.fit_transform(givenMat)
    H = model.components_
    return W, H
