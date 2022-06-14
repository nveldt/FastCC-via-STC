using MAT
include("../Exp-CC/run_alg_functions.jl")
include("../../src/stcplus_lp_relaxation.jl")
for j = 1:2:57

    M = matread("../../data/heurmat/heur$j.mat")
    A = Float64.(M["A"])
    n = size(A,1)
    m = sum(A)/2
    sparsity = round(m/binomial(n,2),digits = 2)

    num_tri, num_ow = count_wedges_triangles(A)
    println("$i $graph \t $n \t $m \t $num_ow \t $num_tri")
    graph_stats = [n;m;num_ow;num_tri]

    tl = 1800; pivtimes = 1; 
    times206 = 1
    runtimes, results = run_stcplus_lp_fasternotfaster(A,pivtimes,times206,tl)
    # runtimes = [lb_time; time_round; time_round_206]
    # results = [stclp_bd, ub_stclp, ub_stclp_206, status, iscclp]
    matwrite("pace_heuristic/heur$(j)_stcplus_lp.mat", Dict("runtimes"=>runtimes,"results"=>results,"pivtimes"=>pivtimes,"times206"=>times206))

   
end