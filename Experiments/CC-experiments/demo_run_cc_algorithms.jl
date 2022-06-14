using MAT
using SparseArrays

include("../../src/cc_lp_relaxations.jl")
include("../../src/faster_many_pivot.jl")
include("../../src/mfp_cc_fast.jl")
include("run-cc-functions.jl")


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


# Pick a graph from the list
i = 3
graph = tinytestgraphs[i]

tl = 1800           # time limit
pivtimes = 100      # number of times to run pivot in rounding step
   
F = matread("../../data/smallgraphs/$graph.mat")
A = F["A"]
n = size(A,1)
A = sparse(A)


## MFP-CC results
runtimes_mfp, results_mfp = run_matchflippivot_cc(A,pivtimes)
lb_mfp = results_mfp[1]
ub_mfp = results_mfp[2]
ub_mfp_det = results_mfp[4] # deterministic rounding
approx_mfp = ub_mfp/lb_mfp

## LP-STC+ results and LP-CC results
checkcc = true
tl = 1800
times206 = 10
pivtimes = 50
runtimes_cclp, results_cclp, runtimes_stclp, results_stclp = run_lazycc_stcplus_combined(A,pivtimes,times206, tl,checkcc)

