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



for i = graphset
    tl = 1800
    pivtimes = 100
    graph = graphs[i]

    # Load graph
    F = matread("../../data/smallgraphs/$graph.mat")
    A = F["A"]
    n = size(A,1)
    A = sparse(A)
    Is, Js, Vs = findnz(triu(A))
    Elist = [Is Js]
    println("$graph mfp")
    clus, runtimes, results = run_matchflippivot_cd(A,Elist,pivtimes)
    matwrite("cd-experiment-results/$(graph)_mfp.mat", Dict("clus"=>clus,
        "runtimes"=>runtimes,"results"=>results))
end

## STC LP
for i = graphset
    tl = 1800; pivtimes = 100; graph = graphs[i]
    F = matread("../../data/smallgraphs/$graph.mat")
    A = F["A"]; A = sparse(A); Is, Js, Vs = findnz(triu(A)); Elist = [Is Js]
    println("$graph stc")

    clus, runtimes, results = run_stclp(A,Elist,pivtimes,tl)
    matwrite("cd-experiment-results/$(graph)_stclp.mat", Dict("clus"=>clus,
        "runtimes"=>runtimes,"results"=>results))
end

## CD LP
for i = graphset
    tl = 1800; pivtimes = 100; graph = graphs[i]
    F = matread("../../data/smallgraphs/$graph.mat")
    A = F["A"]; A = sparse(A); Is, Js, Vs = findnz(triu(A)); Elist = [Is Js]
    println("$graph cdlp")

    clus, runtimes, results = run_cdlp(A,Elist,pivtimes,tl)
    matwrite("cd-experiment-results/$(graph)_cdlp.mat", Dict("clus"=>clus,
        "runtimes"=>runtimes,"results"=>results))
end

## STC ILP
# for i = graphset
#     tl = 1800; pivtimes = 100; graph = graphs[i]
#     F = matread("../../data/smallgraphs/$graph.mat")
#     A = F["A"]; A = sparse(A); Is, Js, Vs = findnz(triu(A)); Elist = [Is Js]
    
#     clus, runtime, results = run_stc_opt(A,Elist,tl)
#     matwrite("cd-experiment-results/$(graph)_stc_opt.mat", Dict("clus"=>clus,
#         "runtime"=>runtime,"results"=>results))
# end


## CD ILP
# for i = graphset
#     tl = 1800; graph = graphs[i]
#     F = matread("../../data/smallgraphs/$graph.mat")
#     A = F["A"]; A = sparse(A); Is, Js, Vs = findnz(triu(A)); Elist = [Is Js]
    
#     clus, runtime, results = run_cd_opt(A,Elist,tl)
#     matwrite("cd-experiment-results/$(graph)_cd_opt.mat", Dict("clus"=>clus,
#         "runtime"=>runtime,"results"=>results))
# end