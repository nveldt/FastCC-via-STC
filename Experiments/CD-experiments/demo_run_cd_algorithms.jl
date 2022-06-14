using MAT
using SparseArrays

include("../../src/cd_lp_relaxations.jl")
include("../../src/cd_rounding.jl")
include("../../src/faster_many_pivot.jl")
include("run-cd-functions.jl")


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
Is, Js, Vs = findnz(triu(A))
Elist = [Is Js]

## MFP-CD results
clus_mfp, runtimes_mfp, results_mfp = run_matchflippivot_cd(A,Elist,pivtimes)
lb_mfp = results_mfp[1]
ub_mfp = results_mfp[2]
approx_mfp = ub_mfp/lb_mfp


## LP-STC results
clus_stclp, runtimes_stclp, results_stclp = run_stclp(A,Elist,pivtimes,tl)
lb_stclp = results_stclp[1]
ub_stclp = results_stclp[2]
approx_stclp = ub_stclp/lb_stclp

 
## LP-CD results
clus_cdlp, runtimes_cdlp, results_cdlp = run_stclp(A,Elist,pivtimes,tl)
lb_cdlp = results_cdlp[1]
ub_cdlp = results_cdlp[2]
approx_cdlp = ub_cdlp/lb_cdlp
