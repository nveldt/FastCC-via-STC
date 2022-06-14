
# Datasets
----------


fb-graph-names.txt 
snap-graph-names.txt

These are the names of the graph files for Facebook100 and the SNAP collection respectively. In order to reproduce these experiments, these datasets need to be obtained separately, standardized, and stored in a folder. Standardization just means that they are converted into simple graphs (e.g., with directions and weights removed). Instructions for converting these to simple graphs are provided in the "standardize-snap-graphs.jl" file in the data folder.

# Running experiments
---------------------

All results for the Facebook100 graphs can be obtained by running:

"run-louvain-fb.jl"

This runs both our MFP algorithm (applying the randomized rounding procedure once), and the LambdaLouvain algorithm (with the fastest settings, only visiting each node once). This runs within a few minutes on a laptop.

For snap graphs, we have two separate .jl files for running MFP and LambdaLouvain, as these take longer.

"run-snap-just-louvain.jl"

Runs just the Louvain method, with the fastest settings.

"run-snap-just-mfp.jl"

Runs the MFP algorithm with only one pivot for rounding.

We use the fastest settings for obtaining a clustering from each method, in order to provide as far of a comparison as possible. We found that running algorithms for more settings tends to only improve a posteriori approximation ratios slightly.

The file "run-louvain-fb-longer-round.jl" runs experiments again on Facebook100 graphs, with algorithmic settings that take longer but lead to improved approximation guarantees.

# Plotting Results
------------------

See "plot-louvain-fb-snap.jl" for plots of results.
