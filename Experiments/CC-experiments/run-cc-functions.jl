# These are functions that take in a graph, compute different types of lower bounds, and rounds them in different ways.

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


    # Deterministic rounding
    tic = time()
    ub_det, c_det = STCplus_to_CC_round_det(A,Eadd,Edel)
    time_det = time()-tic
    
    # Get pivot results, to compare against
    II,JJ = findnz(triu(A))
    Elist = [II JJ]
    tic = time()
    ub_piv, c_piv = many_pivot_faster(A,Elist,pivtimes)
    time_piv = time()-tic

    @assert(ub_piv == check_cc_obj(A,c_piv))
    runtimes = [lb_time; time_rand; time_rand_avg; time_det; time_piv]
    results = [bd; ub_rand_best; ub_rand_avg; ub_det; ub_piv]

    return runtimes, results
end

function run_lazycc_lp(A,times206, tl)
    outflag = false
    LP = true
    verb = true
    tic = time()
    status, bd, Elist, cclp_soln = LazyExactCCGurobi(A;outputflag = outflag,verbose = verb,LP = LP, timelimit = tl)
    lb_time = time()-tic
    
    tic = time()
    cclp_round, c = many_206_faster(A,cclp_soln,times206)
    time_rd = time()-tic

    runtimes = [lb_time; time_rd]
    results = [bd, cclp_round, status]

    return runtimes, results
end

function run_lazycc_stcplus_combined(A,pivtimes,times206, tl,checkcc)
    outflag = false
    LP = true
    verb = true
    tic = time()
    status_stclp, lb_stclp, stclp_soln, lb_time_stclp, status_cclp, lb_cclp, Elist, cclp_soln = LazyCCandSTCplus(A;outputflag = outflag, verbose = verb,LP = LP,timelimit = tl)
    lb_time_cclp = time()-tic

    # Does it satisfy CC constraints?
    if checkcc
        iscclp = is_CC_feasible(A,Elist,stclp_soln)
    else
        iscclp = -1
    end
    
    # round LP-CC
    tic = time()
    cclp_round, c = many_206_faster(A,cclp_soln,times206)
    time_rd = time()-tic

    runtimes_cclp = [lb_time_cclp; time_rd]
    results_cclp = [lb_cclp, cclp_round, status_cclp]

    # round STCplus LP standard
    tic = time()
    ub_stclp, c, avg_stclp, avgtime_stclp = STCplus_to_CC_round_LP(A,Elist,stclp_soln,pivtimes)
    time_round = time()-tic

    # Round using 2.06 approximation for CC (Chawla et al. 2015). No guarantees, but works well
    ub_stclp_206, c = many_206_faster(A,stclp_soln,times206)
    time_round_206 = time()-tic

    runtimes_stclp = [lb_time_stclp; time_round; time_round_206]
    results_stclp = [lb_stclp, ub_stclp, ub_stclp_206, status_stclp, iscclp]

    return runtimes_cclp, results_cclp, runtimes_stclp, results_stclp
end

function run_lazycc_opt(A,tl)
    outflag = false
    LP = false
    verb = true
    tic = time()
    status, cc_opt, Elist, cc_soln = LazyExactCCGurobi(A;outputflag = outflag,verbose = verb,LP = LP, timelimit = tl)
    timer = time()-tic

    return timer, cc_opt, status
end


