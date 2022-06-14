README

Instructions for reproducing experiments on PACE Challenge Graphs.

## Datasets

The PACE graphs can be downloaded at

https://fpt.akt.tu-berlin.de/pace2021/exact.tar.gz
https://fpt.akt.tu-berlin.de/pace2021/heur.tar.gz

The first link is to the set of benchmark graphs for the competition between exact algorithms for cluster editing, and the second link is for the benchmark graphs for the competition between heuristic methods. 

The files “read-pace-data.jl” and “read-gr.jl” in the data folder can be used to store these in a .mat format which we use for our experiments. 

## Running KaPoCE

In order to run KaPoCE, this software needs to be separately downloaded and installed which may take a significant amount of configuration. A copy of the code with its original license is contained in the folder “include/cluster_editing”. This includes a few changes to the make files that were needed to run it on our platform. These changes and this version of the code may not work for you.

The PACE challenge graphs should be downloaded and installed in the “include/cluster_editing/instances” folder.

The files “run_kapoce_exact.jl” and “run_kapoce_heuristics.jl” in the folder 

“include/cluster_editing/build” 

illustrate how this software was run using a Julia wrapper. These will not work unless you download and install KaPoCE. 


## Running Louvain and MFP

The files “run-kapoce-exact.jl” and “run-kapoce-heuristics.jl” run LambdaLouvain (with the fastest settings) and MFP-CE on the exact and heuristic graphs. File “run-pace-stcpluslp.jl” runs our LP-STC+ algorithm, which comes with tighter lower bounds but is slower.

Results are contained in the folders: “pace_exact” and “pace_heuristic”.

## Plotting results

First, unzip the heur_out.zip folder.

Then, run “plot-pace.jl”.