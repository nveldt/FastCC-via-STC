# A few packaged up functions for running a single experiment


"""
Run MFP (cluster deletion version) for one graph A.
"""
function run_matchflippivot_cd(A,Elist,pivtimes)
    # Get lower bounds
    tic = time()
    bd, Elist_match, soln_match = maximal_disjoint_openwedge_rounding(A)
    timer = time()-tic

    # Round the matching
    tic = time()
    r_match, clus, comp = CD_LP_round(A,Elist,soln_match,pivtimes)
    timer_r = time()-tic
    # println("Upper bound done")
    @assert(isCDfeasible(A,clus))

    # println("Done checking feasibility for CD ")
    runtimes = [timer; timer_r]
    results = [bd, r_match, comp]

    return clus, runtimes, results
end

function run_stclp(A,Elist,pivtimes,tl)
    tic = time()
    stat_stclp, bd, Elist_stc, soln_stclp = CDSTC_Gurobi(A; clusterdeletion = false, LP = true,timelimit = tl,outputflag = false)
    timer = time()-tic

    # Same as CD-LP?
    iscdlp = check_cd_constraints(A,Elist,soln_stclp)

    # Round STCLP
    tic = time()
    r_stclp, clus, comp = CD_LP_round(A,Elist,soln_stclp,pivtimes)
    timer_r = time()-tic

    if stat_stclp
        # If here, then the LP solver converged
        @assert(isCDfeasible(A,clus))
    end

    runtimes = [timer; timer_r]
    results = [bd, r_stclp, comp,stat_stclp,iscdlp]

    return clus, runtimes, results
end


# 2-approximation, full cluster deletion LP relaxation
function run_cdlp(A,Elist,pivtimes,tl)
    tic = time()
    stat_cdlp, bd, Elist_cd, soln_cdlp = CDSTC_Gurobi(A; clusterdeletion = true, LP = true,timelimit = tl,outputflag = false)
    timer = time()-tic
    @assert(Elist_cd == Elist)

    tic = time()
    r_cdlp, clus, comp = CD_LP_round(A,Elist,soln_cdlp,pivtimes)
    timer_r = time()-tic
    if stat_cdlp
        @assert(isCDfeasible(A,clus))
    end

    runtimes = [timer, timer_r]
    results = [bd, r_cdlp, comp, stat_cdlp]

    return clus, runtimes, results
end

# Find optimal solution to the minSTC problem
function run_stc_opt(A,Elist,tl)
    tic = time()
    stat_stc, stc_opt, Elist_stc, soln_stc_opt = CDSTC_Gurobi(A; clusterdeletion = false, LP = false,timelimit = tl)
    runtime = time()-tic
    @assert(Elist_stc == Elist)

    # Check if answer also is same as CD-ILP
    iscd = check_cd_constraints(A,Elist,soln_stc_opt)
    
    # Can be used as a 2-approx for CD.
    r_stc_opt, clus, comp = CD_LP_round(A,Elist,soln_stc_opt,1)

    results = [stc_opt, r_stc_opt, comp, stat_stc, iscd]

    return clus, runtime, results
end

# Find optimal solution to cluster deletion 
function run_cd_opt(A,Elist,tl)
    tic = time()
    stat_cd, cd_opt, Elist_cd, soln_cd_opt = CDSTC_Gurobi(A; clusterdeletion = true, LP = false,timelimit = tl)
    runtime = time()-tic
    @assert(Elist_cd == Elist)

    r_cd_opt, clus, comp = CD_LP_round(A,Elist,soln_cd_opt,1)
    
    if stat_cd
        @assert(round(Int64,r_cd_opt) == round(Int64,cd_opt))
    end

    results = [cd_opt, comp, stat_cd]

    return clus, runtime, results
end
