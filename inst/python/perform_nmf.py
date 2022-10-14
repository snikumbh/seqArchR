# Author: snikumbh

# Perform NMF in Python
# A barebones function
# 

import sklearn
from sklearn.decomposition import NMF
from packaging import version

def perform_nmf_func(givenMat, nPatterns, nIter=200, givenAlpha = 0, 
                                    givenL1_ratio = 1, seed_val = 3456):
    #
    if(version.parse(sklearn.__version__) < version.parse('1.0')):
        model = NMF(n_components=nPatterns,
                solver='cd',
                init='nndsvd', 
                alpha=givenAlpha,
                max_iter=nIter,
                tol=1e-03,
                l1_ratio=givenL1_ratio, 
                random_state=seed_val)
    else:
        model = NMF(n_components=nPatterns,
                solver='cd',
                init='nndsvd', 
                alpha_H=givenAlpha,
                alpha_W=givenAlpha,
                max_iter=nIter,
                tol=1e-03,
                l1_ratio=givenL1_ratio, 
                random_state=seed_val)
    #
    W = model.fit_transform(givenMat)
    H = model.components_
    return W, H

## Future versions of scikit-learn
## will require alpha_H and alpha_W
## See full warning below.
## 
## In the future, we could ask for the minimum version 
## of scikit-learn to be >= 1.2
##
# alpha_H=givenAlpha,
# alpha_W=givenAlpha,
## FutureWarning: `alpha` was deprecated in version 1.0 and will be removed 
## in 1.2. Use `alpha_W` and `alpha_H` instead.
