using MAT
include("../Large-graph-experiments/run_mfp_functions.jl")
include("../../src/run_louvain_cc_cd.jl")

for j = 1:2:199

    M = matread("../../data/exactmat/ex$j.mat")
    A = Float64.(M["A"])
    n = size(A,1)
    m = sum(A)/2
    sparsity = round(m/binomial(n,2),digits = 2)

    num_tri, num_ow = count_wedges_triangles(A)
    println("$i $graph \t $n \t $m \t $num_ow \t $num_tri")
    graph_stats = [n;m;num_ow;num_tri]


    ## MFP methods
    # runtimes = [lb_time; time_rand; time_rand_avg; time_piv]
    # results = [bd; ub_rand_best; ub_rand_avg; ub_piv]
    pivtimes = 1
    runtimes_cc, results_cc = run_matchflippivot_cc(A,pivtimes)

    # runtimes = [timer; timer_r]
    # results = [bd, r_match, comp]
    runtimes_cd, results_cd = run_matchflippivot_cd(A,pivtimes)

    cd_t = sum(runtimes_cd)
    cc_t = runtimes_cc[1]+runtimes_cc[4]
    cd_r = runtimes_cd[2]
    cc_r = runtimes_cc[4]
    cc_rat = results_cc[4]/results_cc[1]
    cd_rat = results_cd[2]/results_cd[1]

    # now try lambda-louvain
    numtimes = 1
    maxits = 1
    ll_cc, llcc_time = run_ll_cc(A,numtimes,maxits)

    llcc_rat = ll_cc/results_cc[1]

    ll_cd, llcd_time = run_ll_cd(A,numtimes,maxits)

    llcd_rat = ll_cd/results_cd[1]

    louv_results = [ll_cc, ll_cd]
    louv_time = [llcc_time, llcd_time]

    println("\n $graph")
    println("Total time: $cd_t \t $cc_t \t")
    println("Round time: $cd_r \t $cc_r \t")
    println("Louvain time: $llcc_time $llcd_time")
    println("MFP Ratios:  $cd_rat \t $cc_rat \t")
    println("Louvain Ratio: $llcc_rat $llcd_rat")

    matwrite("pace_exact/exact$(j)_mfp_louvain.mat", Dict("graph_stats"=>graph_stats,
        "runtimes_cc"=>runtimes_cc,"results_cc"=>results_cc,"runtimes_cd"=>runtimes_cd,
        "results_cd"=>results_cd,"louv_results"=>louv_results, "louv_time"=>louv_time))

end