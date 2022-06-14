using MAT
using DataStructures
using StatsBase
using SparseArrays

include("../../src/cc_lp_relaxations.jl")
include("run-cc-functions.jl")
include("../../src/faster_many_pivot.jl")

tinytestgraphs = [
    "KarateA";
    "dolphinsA";
    "lesmisA";
    "polbooksA";
    "adjnounA";
    "footballA";
]

graphs = [
    "Harvard500A";
    "Erdos991A"; 
    "celegansneuralA";
    "Netscience";
    "celegansmetabolicA";
    "RogetA";
    "SmaGriA";
    "standard_emailA";
    "standard_polblogsA";
    "standard_ca-GrQc";
    "standard_caHepThA";
    "standard_EmailEnronA";
    "standard_condmat2005A";
    "standard_ca-AstroPhA";
    "standard_loc-Brightkite";
    # "Caltech36";
    # "Reed98";
    # "Simmons81";
    # "Haverford76";
    # "Swarthmore42";
    # "Amherst41";
    # "Bowdoin47";
    # "Rice31";
    # "Lehigh96";
]


graphs = tinytestgraphs
graphset = 1:length(graphs)

## MFP cc
for i = graphset
    pivtimes = 50
    graph = graphs[i]
    F = matread("../../data/smallgraphs/$graph.mat")
    A = F["A"]
    n = size(A,1)
    A = sparse(A)
    println("$graph")
    runtimes, results = run_matchflippivot_cc(A,pivtimes)
    matwrite("cc-experiment-results/$(graph)_mfp_cc_det.mat", Dict("runtimes"=>runtimes,"results"=>results,"pivtimes"=>pivtimes))
end

# combined: this computes both the CC LP relaxation, as well as the STC+ LP relaxation along the way
for i = graphset
    tl = 1800; pivtimes = 50;  times206 = 10; 
    graph = graphs[i]
    F = matread("../../data/smallgraphs/$graph.mat")
    A = F["A"]; n = size(A,1); A = sparse(A)
    println("$graph")
 
    m = sum(A)/2
    if m < 50000
        checkcc = true
    else
        checkcc = false
    end
    runtimes_cclp, results_cclp, runtimes_stclp, results_stclp = run_lazycc_stcplus_combined(A,pivtimes,times206, tl,checkcc)
    matwrite("cc-experiment-results/$(graph)_stcplus_cc_combined.mat", Dict("runtimes_cclp"=>runtimes_cclp,
        "results_cclp"=>results_cclp,"pivtimes"=>pivtimes,"times206"=>times206,
        "runtimes_stclp"=>runtimes_stclp,"results_stclp"=>results_stclp))
end


# stc lp
# for i = graphset
#     tl = 1800; pivtimes = 50; graph = graphs[i]
#     F = matread("../../data/smallgraphs/$graph.mat")
#     A = F["A"]; n = size(A,1); A = sparse(A)
#     println("$graph")
#     times206 = 10
#     runtimes, results = run_stcplus_lp(A,pivtimes,times206,tl)
#     matwrite("cc-experiment-results/$(graph)_stcplus_lp_faster.mat", Dict("runtimes"=>runtimes,"results"=>results,"pivtimes"=>pivtimes,"times206"=>times206))
# end

# cc lp lazy
# for i = graphset
#     tl = 1800; graph = graphs[i]
#     F = matread("../../data/smallgraphs/$graph.mat")
#     A = F["A"]; n = size(A,1); A = sparse(A)

#     times206 = 10
#     runtimes, results = run_lazycc_lp(A,times206,tl)
#     matwrite("cc-experiment-results/$(graph)_lazycc_lp.mat", Dict("runtimes"=>runtimes,"results"=>results,"times206"=>times206))
# end






