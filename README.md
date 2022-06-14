## README

Code for faster approximation algorithms and more practical lower bounds for correlation clustering and cluster deletion, by rounding lower bounds to strong triadic closure edge labeling problems. To accompany ICML 2022 paper:

"Correlation Clustering via Strong Triadic Closure Labeling: Fast Approximation Algorithms and Practical Lower Bounds"

This repo contains some, but not all of the datasets used in the experiments. Other graphs may be accessed via the SNAP repository or the suitesparse matrix collection.

Details for how to download and standardize all snap graphs for experiments can be found in the folder data/snap-graphs, in file "standardize-snap-graphs.jl".

## CE vs CD

Experiments and implementations for cluter deletion algorithms are typically marked with "CD". Experiments and implementations for cluster editing are typically marked with "CE", or in some cases "CC" as this is equivalent to complete unweighted Correlation Clustering.


## Reproducing results

In the Experiments folder, there is a subfolder for experiments in each section of the paper. 

You can also find demo files to give an example of running algorithms.