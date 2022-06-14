include("../../src/cc_lp_relaxations.jl")
include("../../src/cd_lp_relaxations.jl")
include("../../src/cd_rounding.jl")
include("../../src/faster_many_pivot.jl")
include("../../src/mfp_cc_fast.jl")

"""
Run MFP (cluster editing version) for one graph A.
"""
function run_matchflippivot_cc(A,pivtimes)
    
    # Get lower bound
    tic = time()
    Edel, Eadd, bd = maximal_disjoint_openwedge_cc_dict(A)
    lb_time = time()-tic

    # Randomized rounding
    tic = time()
    ub_rand_best, c_rand, ub_rand_avg, time_rand_avg = STCplus_to_CC_round_rand_faster(A,Eadd, Edel,pivtimes)
    time_rand = time()-tic

    # Get pivot results, to compare against
    II,JJ = findnz(triu(A))
    Elist = [II JJ]
    tic = time()
    ub_piv, c_piv = many_pivot_faster(A,Elist,pivtimes)
    time_piv = time()-tic

    runtimes = [lb_time; time_rand; time_rand_avg; time_piv]
    results = [bd; ub_rand_best; ub_rand_avg; ub_piv]

    return runtimes, results
end



"""
Run MFP (cluster deletion version) for one graph A.
"""
function run_matchflippivot_cd(A,pivtimes)
    # Get lower bounds
    tic = time()
    bd, Elist, soln_match = maximal_disjoint_openwedge_rounding(A)
    timer = time()-tic

    # Round the matching
    tic = time()
    r_match, clus, comp = CD_LP_round(A,Elist,soln_match,pivtimes)
    timer_r = time()-tic

    runtimes = [timer; timer_r]
    results = [bd, r_match, comp]

    return runtimes, results
end