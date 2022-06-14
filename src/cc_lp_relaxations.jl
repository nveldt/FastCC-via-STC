# Approximation algorithms for correlatioin clustering
# By rounding STC+ labeling lower bounds

using Gurobi
using SparseArrays
using LinearAlgebra

gurobi_env = Gurobi.Env()

include("helpers.jl")

"""
CC_STCplus_Gurobi
(complete unweighted correlation clustering & strong triadic closure with edge additions)

Use Gurobi WITHOUT JuMP.jl to solve correlation clustering or
minimum weakness strong triadic closure labeling with edge additions,
either optimally, or just the canonical LP relaxation.

This is quite a bit faster than using the JuMP package.

Input:
    A = adjacency matrix (unsigned, unweighted)

Optional parameters:
    correlationclustering = true:
        solve correlation clustering objective and not just stc+

    LP = true:
        just solve the LP relaxation, not the exact ILP

    Timelimit: upper limit on how long to let Gurobi run

    Outputflag: whether or not to show Gurobi solver status during solve

Output:
    obj = output objective
    Elist = ordered list of edges in graph
    Evals = variable value for each edge
"""
function CC_STCplus_Gurobi(A; correlationclustering = true, outputflag = false, LP = false, timelimit = Inf)
    n = size(A,1)
    W = ones(n,n)
    for i = 1:n
        for j = 1:n
            if A[i,j] == 0
                W[i,j] = -1
            end
        end
    end
    # build varmap and obj
    vmap = -ones(Int,n,n)
    obj = Float64[]
    nvars = 0
    Elist = zeros(Int64,round(Int64,n*(n-1)/2),2)
    for j=1:n
        for i=1:j-1
            nvars += 1
            Elist[nvars,1] = i
            Elist[nvars,2] = j
            vmap[i,j] = nvars-1
            vmap[j,i] = nvars-1
            push!(obj, W[i,j])
        end
    end
    aptr = Ref{Ptr{Cvoid}}()
    if LP
        vtypes = repeat(GRB_CONTINUOUS, nvars)
        ub = ones(nvars)
        err = GRBnewmodel(gurobi_env, aptr, "CCLP", nvars, obj, C_NULL, ub, vtypes, C_NULL)

    else
        vtypes = repeat(GRB_BINARY, nvars)
        err = GRBnewmodel(gurobi_env, aptr, "ExactCC", nvars, obj, C_NULL, C_NULL, vtypes, C_NULL)
    end

    m = aptr[]
    GRBsetdblparam(GRBgetenv(m), GRB_DBL_PAR_TIMELIMIT, timelimit)
    GRBsetintparam(GRBgetenv(m), "OutputFlag",outputflag)



    try
        if correlationclustering

            # all triangle inequality constraints
            cind = Int32[0,0,0]
            cval = Float64[0,0,0]
            for i = 1:n
                for j = i+1:n
                    for k = j+1:n
                        cind[1] = vmap[min(j,k),max(j,k)] # xij
                        cind[2] = vmap[min(i,k),max(i,k)] # xik
                        cind[3] = vmap[min(i,j),max(i,j)] # xjk

                        # xij - xik - xjk <= 0.0
                        cval[1] = 1
                        cval[2] = -1
                        cval[3] = -1
                        error = GRBaddconstr(m, 3, cind, cval, GRB_LESS_EQUAL, 0.0, C_NULL)

                        # -xij + xik - xjk <= 0.0
                        cval[1] = -1
                        cval[2] = 1
                        cval[3] = -1
                        error = GRBaddconstr(m, 3, cind, cval, GRB_LESS_EQUAL, 0.0, C_NULL)

                        # -xij - xik + xjk <= 0.0
                        cval[1] = -1
                        cval[2] = -1
                        cval[3] = 1
                        error = GRBaddconstr(m, 3, cind, cval, GRB_LESS_EQUAL, 0.0, C_NULL)
                    end
                end
            end
        else
            T,W = get_wedges_triangles(A)
            # not correlation clustering, just stc+
            # so the only triangle inequality contraints we need
            # are x_ij + x_ik >= x_jk if (i,j,k) is an open wedge centered at i
            cind = Int32[0,0,0]
            cval = Float64[0,0,0]
            for t in W
                i = t[1]    # This is always the center
                j = t[2]
                k = t[3]
                cind[1] = vmap[min(i,j),max(i,j)] # xij
                cind[2] = vmap[min(i,k),max(i,k)] # xik
                cind[3] = vmap[min(j,k),max(j,k)] # xik
                # x[i,j] + x[i,k] >= x[j,k] so xjk - xij - xik <= 0
                cval[1] = -1
                cval[2] = -1
                cval[3] = 1
                error = GRBaddconstr(m, 3, cind, cval, GRB_LESS_EQUAL, 0.0, C_NULL)
            end

        end

    GRBoptimize(m)
    stat = Ref{Int32}(0)
    GRBgetintattr(m, GRB_INT_ATTR_STATUS, stat)
    #Status codes: https://www.gurobi.com/documentation/9.5/refman/optimization_status_codes.html#sec:StatusCodes
    status = stat[]
    # println("Status = $status")
    if status == 2
        optimal = true
    else
        optimal = false
    end
    robj = Ref{Float64}(0.0)
    GRBgetdblattr(m, GRB_DBL_ATTR_OBJVAL, robj)

    obj = robj[] + (nvars - sum(A)/2)   # adjust by number of non-edges to get the actual objective
    soln = zeros(nvars)
    GRBgetdblattrarray(m, GRB_DBL_ATTR_X, 0, nvars, soln)

    # Elist is a linearization of the node pairs
    # soln[k] is the x variable value for node-par Elist[k]
    # Sanity check if you ever need it for computing obj
    # obj = 0
    # for k = 1:size(Elist,1)
    #     i = Elist[k,1]
    #     j = Elist[k,2]
    #     if A[i,j] == 0
    #         obj += 1-soln[k]
    #     else
    #         obj += soln[k]
    #     end
    # end
    return optimal, obj, Elist, soln
 finally
     #@show "freemodel here!"
     GRBfreemodel(m)
 end
end


"""
Find maximal node-pair disjoint set of open wedges in A.

This lower bounds cluster editing.
"""
function maximal_disjoint_openwedge_cc(A)
    n = size(A,1)
    Neighbs = ConstructAdj(A,n)
    Is = Vector{Int64}()
    Js = Vector{Int64}()
    B = copy(A)
    for i = 1:n
        N = Neighbs[i]
        for jj = 1:length(N)
            j = N[jj]
            ij1 = max(i,j)
            ij2 = min(i,j)
            if B[ij1,ij2] == 2
                continue
            end
            for kk = jj+1:length(N)
                k = N[kk]
                ik1 = max(i,k)
                ik2 = min(i,k)

                jk1 = max(j,k)
                jk2 = min(j,k)

                if B[ik1,ik2] == 2 || A[k,j] > 0 || B[jk1,jk2] == 2
                    continue

                else
                    # if we reach here: (i,j) and (i,k) are edges, but (j,k) is not
                    # So (i,j,k) is an open wedge centered at i.
                    # Also, neither (i,j) nor (i,k) nor (j,k) already appears in our open wedge set.
                    push!(Is,ij2)
                    push!(Js,ij1)

                    push!(Is,ik2)
                    push!(Js,ik1)

                    push!(Is,jk2)
                    push!(Js,jk1)

                    B[ij1,ij2] = 2
                    B[ik1,ik2] = 2
                    B[jk1,jk2] = 2

                    # Now exit this and start with a new j, as we have already
                    # used edge (i,j) and we don't need to check for more
                    # edges involving (i,j)
                    break
                end
            end
        end
    end
    Welist = [Is Js]    # list of the node-pairs in the wedgelist
    bd = length(Is)/3   # lower bound on the cluster editing objective
    return Welist, round(Int64,bd)
end

# confirm this is a set of pair disjoint node pairs
function is_pair_disjoint(Welist)

    seenpairs = Dict()
    for k = 1:size(Welist,1)
        i = Welist[k,1]
        j = Welist[k,2]
        if haskey(seenpairs,(min(i,j),max(i,j)))
            return false
        else
            seenpairs[(min(i,j),max(i,j))] = true
        end
    end
    return true

end

"""
Given a set of node-pairs to flip,
(the output of maximal_disjoint_openwedge_cc(A))
this returns a vector soln that is feasible for
the ILP formulation of the STC+ problem.

The idea is to make the output of our match-flip-pivot
technique fit the output of our functions that use Gurobi
to solve an LP or ILP, so that they can be rounded in the same way.

A is the adjacency matrix for a graph, and
Elist is a linearization of its node pairs.

Welist is a list of node-pairs to flip.
If an edge (i,j) is in Welist, then z_{ij} = 1
meaning that we must flip it (add an edge, or label it weak).
Otherwise, z_{ij} = 0.

We need to apply the following change of variables:

x_{ij} = z_{ij} if (i,j) is an edge in A

x_{ij} = 1-z_{ij} if (i,j) is not an edge in A

"""
function flip_set_to_soln(A,Elist,Welist)
    n = size(A,1)
    pairDict = Dict()
    for t = 1:size(Welist,1)
        i = Welist[t,1]
        j = Welist[t,2]
        @assert(i<j)
        pairDict[(i,j)] = t
    end

    soln = zeros(round(Int64,n*(n-1)/2))
    for k = 1:size(Elist,1)
        i = Elist[k,1]
        j = Elist[k,2]

        # Need to determine x_{ij}

        # First, has this edge been flipped?
        # If so, it will be somewhere in pairDict
        flipped = haskey(pairDict,(i,j))

        # Second, is this an edge?
        isedge = A[i,j] == 1

        # If it's an edge and it's flipped, then x_{ij} = 1
        if flipped && isedge
            soln[k] = 1
        end

        # If it's not an edge and it's flipped, then x_{ij} = 0
        if flipped && ~isedge
            soln[k] = 0
        end

        # If it's an edge and it's not flipped, then x_{ij} = 0
        if ~flipped && isedge
            soln[k] = 0
        end

        # If it's not an edge and it's not flipped, then x_{ij} = 1
        if ~flipped && ~isedge
            soln[k] = 1
        end
    end

    return soln

end

"""
Check if a solution vector is feasible for the STC+ ILP .

A is the adjacency matrix for a graph, and
Elist is a linearization of its node pairs, consistent
with the way the CC_STCplus_Gurobi function linearizes node
pairs.
"""
function is_STCplus_feasible(A,Elist,soln)
    n = size(A,1)
    T,W = get_wedges_triangles(A)

    EdgeDict = Dict()
    for k = 1:size(Elist,1)
        i = Elist[k,1]
        j = Elist[k,2]
        @assert(i<j)
        EdgeDict[(i,j)] = k
    end

    # There's a constraint for every open wedge. Check them
    for w in W
        i = w[1]    # this is the center of the wedge
        j = w[2]
        k = w[3]

        # get the indices for each node pair
        # in Elist and soln
        ij = EdgeDict[min(i,j),max(i,j)]
        ik = EdgeDict[min(i,k),max(i,k)]
        jk = EdgeDict[min(j,k),max(j,k)]

        xij = soln[ij]
        xik = soln[ik]
        xjk = soln[jk]

        # We require
        #     xjk <= xik + xjk
        # theoretically.
        # Would be okay with
        #     xjk - xik - xjk <= 1e-8
        # numerically.
        if xjk - xik - xij > 1e-8
            return false
        end
    end
    return true
end

"""
Checks whether a solution to the STC+ LP
satisfies the constraints of the correlation clustering LP.
I.e., are all triangle constraints satisfied?

soln[k] is the variable for node pair Elist[k,:].

"""
function is_CC_feasible(A,Elist,soln)

    # map from edge to edge index
    n = size(A,1)
    D = zeros(n,n)
    for k = 1:size(Elist,1)
        i = Elist[k,1]
        j = Elist[k,2]
        D[i,j] = soln[k]
        D[j,i] = soln[k]
    end
    epsi = 1e-8
    for i = 1:n-2
        for j = i+1:n-1
            a = D[j,i]
            for k = j+1:n
                b = D[k,i]
                c = D[k,j]
                if a-b-c > epsi
                    return false
                end
                if b-a-c > epsi
                    return false
                end

                if c-a-b > epsi
                    return false
                end
            end
        end
    end
    return true
end


"""
Rounds a lower bound for STC+ into a feasible solution
    for a cluster editing problem.

This works for:
    * rounding the CC LP relaxation (2-approx)
    * rounding the STC+ LP relaxation (4-approx)
    * rounding a feasible solution to the STC+ ILP (6-approx)

Input:
    A = adjacency matrix

    Elist = linear ordering of node-pairs

    soln = feasible solution to the CC LP relaxation
            or the STC+ LP relaxation
            or the STC+ ILP

Output:
    Approximation for CC.

    Step 1: Construct derived graph where edges are pairs with x[e] < 1/2
    Step 2: Apply pivot to the new graph
    Step 3: Compute CC objective score
    Step 4: Check how far this is from the lower bound

    pivtimes is the number of times to run the pivot
    method on the derived graph.
"""
function STCplus_to_CC_round(A,Elist,soln,pivtimes = 1)
    n = size(A,1)
    @assert(size(Elist,1) == n*(n-1)/2)

    # Elist is a linearization of all node pairs in A,
    # not just the edges in A.

    # Find all the node pairs Elist[k] that have soln[k] < 1/2
    tol = 1e-8
    keep = findall(x->x<0.5-tol,soln)

    # These all become edges in a new derived graph Anew
    Inew = Elist[keep,1]
    Jnew = Elist[keep,2]
    Anew = sparse(Inew,Jnew,ones(length(Jnew)),n,n)
    Anew = Anew + Anew'

    # Apply pivot on this graph multiple times,
    # returning the best output
    clus = permutation_pivot(Anew)
    cc_obj = check_cc_obj(A,clus)
    for jj = 1:pivtimes
        # Pivot is fast, so we can run it multiple times
        clusnew = permutation_pivot(Anew)
        objnew = check_cc_obj(A,clus) # Check output on original graph
        if objnew < cc_obj
            cc_obj = objnew
            clus = clusnew
        end
    end

    return round(Int64,cc_obj), clus
end

"""
Check complete unweighted correlation clustering objective
"""
function check_cc_obj(A::SparseArrays.SparseMatrixCSC,c)
    n = length(c)
    m = sum(A)/2
    w_volA = n
    lam = 1/2
    obj = n*(n-1)/4 - m/2
    for i = 1:maximum(c)
        S = findall(x->x==i,c)
        AS = A[:,S]
        vol = sum(AS.nzval);
        SAS = AS[S,:]
        edges = sum(SAS.nzval);
        cut = vol-edges
        w_volS = length(S)
        obj += 1/2*(cut - lam*w_volS*(w_volA-w_volS))
    end
    return round(Int64,2*obj)
end

# Run pivot many times
function many_pivot(A,pivtimes = 10)
    clus = permutation_pivot(A)
    cc_obj = check_cc_obj(A,clus)
    for jj = 1:pivtimes
        # Pivot is fast, so we can run it multiple times
        clusnew = permutation_pivot(A)
        objnew = check_cc_obj(A,clus) # Check output on original graph
        if objnew < cc_obj
            cc_obj = objnew
            clus = clusnew
        end
    end
    return cc_obj, clus
end

"""
Rounding for the matching lower bound,
with deterministic rounding procedure.

This is slow and serves just as a proof of concept
"""
function STCplus_to_CC_det_round(A,Flipped)
    n = size(A,1)

    # Apply pivot on this graph multiple times,
    # returning the best output
    clus = deterministic_pivot(A,Flipped)
    cc_obj = check_cc_obj(A,clus)

    return round(Int64,cc_obj), clus
end

"""
Deterministic pivoting for match-flip-pivot,
when given a set of edges to flip

This is slow and serves just as a proof of concept
"""
function deterministic_pivot(A,Flipped)

    n = size(A,1)

    # create a new graph with the flipped edges
    B = copy(A)

    for k = 1:size(Flipped,1)
        i = Flipped[k,1]
        j = Flipped[k,2]
        @assert(i<j)
        # This edge is opposite in graph B
        B[i,j] = 1-A[i,j]
        B[j,i] = 1-A[i,j]
    end
    dropzeros!(B)

    c = zeros(Int64,n)  # cluster indicator vector
    unclus = ones(n)    # indicator for unclustered nodes
    clusnum = 1
    while sum(unclus) > 0
        p, nbs = get_best_pivot(A,B,unclus)
        C = union(p,nbs)
        c[C] .= clusnum
        clusnum += 1
        unclus[C] .= 0
        B[:,C] .= 0
        B[C,:] .= 0
        dropzeros!(B)
    end
    return c
end


"""
Find weights for the best pivot to select in a
    deterministic pivoting strategy for cc

Elist = a linearization of node pairs from A
Flipped = a set of edges that are flipped
"""
function get_best_pivot(A,B,unclus)
    # A is the original graph
    # B is the flipped graph

    n = size(A,1)
    NeighbsB = ConstructAdj(B,n)
    Nums = zeros(n)
    Denoms = zeros(n)

    for k = 1:n
        nbs = NeighbsB[k]
        if unclus[k] == 0
            # this ensures we won't pick it as a pivot
            Nums[k] = 10
        end
        for ii = 1:length(nbs)
            for jj = ii+1:length(nbs)
                i = nbs[ii]
                j = nbs[jj]
                if B[i,j] == 0
                    # Then ijk is an open wedge in B.
                    # It's centered at k

                    # To understand these update, study the determinstic pivoting
                    # requirements in Theorem 3.1 in the work of van Zuylan and Williamson

                    # Updates for node k: ij is in Tk- in B
                    if A[i,j] == 1
                        # Then ij was flipped, which increases denominator
                        Denoms[k] += 1
                    else
                        Nums[k] += 1
                    end

                    # Updates for node j: ik is in Tj+ in B
                    if A[i,k] == 1
                        # then w_{ik}^+ == 1, edge not flipped
                        Nums[j] += 1
                    else
                        Denoms[j] += 1
                    end

                    # Updates for node i: jk is in Tk+ in B
                    if A[j,k] == 1
                        # then w_{jk}^+ == 1, edge not flipped
                        Nums[i] += 1
                    else
                        Denoms[i] += 1
                    end

                    # we visit each open wedge EXACTLY once with this
                    # iteration, so this shouldn't overcount or double
                    # count anything when we are updating scores for
                    # each node
                end
            end
        end
    end
    P = Nums./Denoms
    themin, p = findmin(P)
    # println("$themin $p")
    return p, NeighbsB[p]
end


"""
Chawla et al. 2015 Round.

The 2.06 approximation for rounding the canonical relaxation
of cluster editing.
"""
function chawla_round(A,Elist,soln)

    n = size(A,1)

    a = 0.19
    b = 0.5905

    # Edge list for a new graph to run pivot on
    Is = Vector{Int64}()
    Js = Vector{Int64}()

    for k = 1:size(Elist,1)
        i = Elist[k,1]
        j = Elist[k,2]
        x = soln[k]
        if A[j,i] == 1

            if x < a
                fij = 0
            elseif x >= b
                fij = 1
            else
                fij = ((x-a)/(b-a))^2
            end

            if rand() < 1-fij
                push!(Is,i)
                push!(Js,j)
            end

        else
            fij = x
            if rand() < 1-fij
                push!(Is,i)
                push!(Js,j)
            end

        end
    end
    B = sparse(Is,Js,ones(length(Is)),n,n)
    B = B + B'

    c = permutation_pivot(B)        # permute with respect to B
    obj = check_cc_obj(A,c)      # check output with respect to A

    return c, obj
end


"""
Chawla et al. Round.

The 2.06 approximation for rounding the canonical relaxation
of cluster editing

A slightly tweaked version
"""
function chawla_round_faster(A,ElistA,m,soln)

    n = size(A,1)

    a = 0.19
    b = 0.5905
    bma = b-a
    # Edge list for a new graph to run pivot on
    Is = Vector{Int64}()
    Js = Vector{Int64}()

    k = 0
    # N = round(Int64,n*(n-1)/2)
    # rand_set = rand(N)
    for j=1:n
    for i=1:j-1
        k += 1

        x = soln[k]
        if A[i,j] == 1

            if x < a
                # fij = 0
                push!(Is,i)
                push!(Js,j)

            elseif x >= b
                fij = 1
            else
                fij = ((x-a)/bma)^2
                if rand() < 1-fij
                    push!(Is,i)
                    push!(Js,j)
                end
            end

        else
            fij = x
            if x == 0
                push!(Is,i)
                push!(Js,j)
            elseif x == 1
                continue
            else
                if rand() < 1-fij
                    push!(Is,i)
                    push!(Js,j)
                end
            end

        end
    end
    end
    B = sparse(Is,Js,ones(length(Is)),n,n)
    B = B + B'
    NeighbsB = ConstructAdj(B,n)
    # c = permutation_pivot(B)        # permute with respect to B
    # obj = check_cc_obj(A,c)      # check output with respect to A
    c, Sbin = permutation_pivot_fastest(B,NeighbsB)
    obj = check_cc_obj_faster(ElistA,c,m,Sbin)
    # obj2 = check_cc_obj(A,c)
    # obj3 = check_cc_obj_fastish(A,ElistA,c,m)
    # @show obj, obj2, obj3
    # @assert(obj == obj2)
    return c, obj
end

"""
Run the randomized 2.06 approx many times for standard cc.

numtimes is the number of times to generate the random graph
and run pivot on it
"""
function many_206_faster(A,soln,numtimes = 1)
    Is,Js = findnz(triu(A))
    ElistA = [Is Js]
    m = round(Int64,sum(A)/2)
    n = size(A,1)

    clus, cc_obj = chawla_round_faster(A,ElistA,m,soln)
    for jj = 1:numtimes
        # Pivot is fast, so we can run it multiple times
        clusnew, objnew = chawla_round_faster(A,ElistA,m,soln)
        if objnew < cc_obj
            cc_obj = objnew
            clus = clusnew
        end
    end
    return cc_obj, clus
end


"""
Run the randomized 2.06 approx many times for standard cc.

numtimes is the number of times to generate the random graph
and run pivot on it
"""
function many_206(A,Elist,soln,numtimes = 1)
    clus, cc_obj = chawla_round(A,Elist,soln)
    for jj = 1:numtimes
        # Pivot is fast, so we can run it multiple times
        clusnew, objnew = chawla_round(A,Elist,soln)
        if objnew < cc_obj
            cc_obj = objnew
            clus = clusnew
        end
    end
    return cc_obj, clus
end

"""
Lazy constraints correlation clustering, directly in Gurobi without using JuMP.

This first solves the LP relaxation of STC+, then iteratively (i.e. "lazily") adds in
violated constraints until the full complete unweighted correlation clustering
LP is solved.
"""
function LazyExactCCGurobi(A;outputflag = true,verbose = true,LP = false, timelimit = Inf)
    n = size(A,1)
    W = ones(n,n)
    for i = 1:n
        for j = 1:n
            if A[i,j] == 0
                W[i,j] = -1
            end
        end
    end
    #nvars = div(n*(n-1), 2)
    # build varmap and obj
    starttime = time()
    vmap = -ones(Int,n,n)
    obj = Float64[]
    nvars = 0
    Elist = zeros(Int64,round(Int64,n*(n-1)/2),2)
    for j=1:n
        for i=1:j-1
            nvars += 1
            Elist[nvars,1] = i
            Elist[nvars,2] = j
            vmap[i,j] = nvars-1
            vmap[j,i] = nvars-1
            push!(obj, W[i,j])
        end
    end
    aptr = Ref{Ptr{Cvoid}}()

    if LP
        vtypes = repeat(GRB_CONTINUOUS, nvars)
        ub = ones(nvars)
        err = GRBnewmodel(gurobi_env, aptr, "CCLP", nvars, obj, C_NULL, ub, vtypes, C_NULL)
    else
        vtypes = repeat(GRB_BINARY, nvars)
        err = GRBnewmodel(gurobi_env, aptr, "ExactCC", nvars, obj, C_NULL, C_NULL, vtypes, C_NULL)
    end
    m = aptr[]
    timespent = time()-starttime
    timeremaining = timelimit-timespent
    tlim = max(0,timeremaining)
    GRBsetdblparam(GRBgetenv(m), GRB_DBL_PAR_TIMELIMIT, tlim)
    GRBsetintparam(GRBgetenv(m), "OutputFlag",outputflag)

    try

    cind = Int32[0,0,0]
    cval = Float64[0,0,0]
    for i = 1:n
        NeighbsI = findall(x->x>0,(A[:,i]))
        numNeighbs = size(NeighbsI,1)
        for u = 1:numNeighbs-1
            j = NeighbsI[u]
            for v = u+1:numNeighbs
                k = NeighbsI[v]
                if A[j,k] == 0
                    # Then we have a bad triangle: (i,j), (i,k) \in E
                    # but (j,k) is not an edge, so D(j,k) wants to be 1
                    #assert(i<j<k)
                    #@constraint(m,x[min(j,k),max(j,k)] - x[min(i,k),max(i,k)] - x[min(i,j),max(i,j)] <= 0)
                    cind[1] = vmap[min(j,k),max(j,k)]
                    cind[2] = vmap[min(i,k),max(i,k)]
                    cind[3] = vmap[min(i,j),max(i,j)]
                    cval[1] = 1
                    cval[2] = -1
                    cval[3] = -1
                    error = GRBaddconstr(m, 3, cind, cval, GRB_LESS_EQUAL, 0.0, C_NULL)
                end
            end
        end
    end

    # Find intial first solution
    if verbose
        println("First round of optimization")
    end
    #JuMP.optimize!(m)
    GRBoptimize(m)
    stat = Ref{Int32}(0)
    GRBgetintattr(m, GRB_INT_ATTR_STATUS, stat)
    status = stat[]
    if status == 2
        optimal = true
    else
        optimal = false
    end
    robj = Ref{Float64}(0.0)
    GRBgetdblattr(m, GRB_DBL_ATTR_OBJVAL, robj)

    obj = robj[] + (nvars - sum(A)/2)   # adjust by number of non-edges to get the actual objective

    soln = zeros(nvars)
    GRBgetdblattrarray(m, GRB_DBL_ATTR_X, 0, nvars, soln)
    D = zeros(n,n)
    for j=1:n
        for i=1:j-1
            D[i,j] = soln[vmap[i,j]+1]
            D[j,i] = soln[vmap[i,j]+1]
        end
    end

    while true
        # x will naturally be upper triangular, but 'find_violations'  wants lower triangular
          #D = Matrix(JuMP.value.(x)')

         # Store violating tuples in a vector
         violations = Vector{Tuple{Int,Int,Int}}()

         find_violations!(D,violations)

         # Iterate through and add constraints
         numvi = size(violations,1)
         if verbose
            print("Adding in $numvi violated constraints.")
         end

         # Violations (a,b,c) are ordered so that a < b, and (a,b) is the value
         # that needs to be positive in the inequality:
         # x_ab - x_ac - x_bc <= 0.
         # The other two (b and c) could satisfy either b < c or c <= b.
         # We need to make sure we are only placing constraints in the upper
         # triangular portion of the matrix.
         # We do so by just calling min and max on pairs of nodes
         for v in violations
             #assert(v[1]<v[2])
             cind[1] = vmap[v[1],v[2]]
             cind[2] = vmap[min(v[1],v[3]),max(v[1],v[3])]
             cind[3] = vmap[min(v[2],v[3]),max(v[2],v[3])]
             cval[1] = 1
             cval[2] = -1
             cval[3] = -1
             #@show cind, cval
             error = GRBaddconstr(m, 3, cind, cval, GRB_LESS_EQUAL, 0.0, C_NULL)
             #z = sparsevec(cind.+1, cval, nvars)
             #@show z'*soln
             #@constraint(m,x[v[1],v[2]] - x[min(v[1],v[3]),max(v[1],v[3])] - x[min(v[2],v[3]),max(v[2],v[3])] <= 0)
         end
         if numvi == 0
             if verbose
                println(" Optimal solution found.")
             end
             break
         end
         if verbose
            println(" And re-solving the optimization problem.")
         end

         # update time remaining
         timespent = time()-starttime
         timeremaining = timelimit-timespent
         tlim = max(0,timeremaining)
         GRBsetdblparam(GRBgetenv(m), GRB_DBL_PAR_TIMELIMIT, tlim)
         if verbose
            println("time left until time is up = $timeremaining")
         end
         GRBupdatemodel(m)
         err = GRBoptimize(m)
         stat = Ref{Int32}(0)
         GRBgetintattr(m, GRB_INT_ATTR_STATUS, stat)
         #Status codes: https://www.gurobi.com/documentation/9.5/refman/optimization_status_codes.html#sec:StatusCodes
         status = stat[]
         # println("Status = $status")
         if status == 2
             optimal = true
         else
             optimal = false
         end
         robj = Ref{Float64}(0.0)
         GRBgetdblattr(m, GRB_DBL_ATTR_OBJVAL, robj)

         obj = robj[] + (nvars - sum(A)/2)   # adjust by number of non-edges to get the actual objective
         soln = zeros(nvars)
         GRBgetdblattrarray(m, GRB_DBL_ATTR_X, 0, nvars, soln)
         for j=1:n
             for i=1:j-1
                 D[i,j] = soln[vmap[i,j]+1]
                 D[j,i] = soln[vmap[i,j]+1]
             end
         end

     end
     return optimal, obj, Elist, soln
 finally
     #@show "freemodel here!"
     GRBfreemodel(m)
 end
end

"""
This version of LazyExactCC spits out the STC+ LP relaxation
as a subproblem, since it is solved along the way to solving 
the CC LP relaxation.
"""
function LazyCCandSTCplus(A;outputflag = true,verbose = true,LP = false, timelimit = Inf)
    n = size(A,1)
    W = ones(n,n)
    for i = 1:n
        for j = 1:n
            if A[i,j] == 0
                W[i,j] = -1
            end
        end
    end
    # nvars = div(n*(n-1), 2)
    # build varmap and obj
    starttime = time()
    vmap = -ones(Int,n,n)
    obj = Float64[]
    nvars = 0
    Elist = zeros(Int64,round(Int64,n*(n-1)/2),2)
    for j=1:n
        for i=1:j-1
            nvars += 1
            Elist[nvars,1] = i
            Elist[nvars,2] = j
            vmap[i,j] = nvars-1
            vmap[j,i] = nvars-1
            push!(obj, W[i,j])
        end
    end
    aptr = Ref{Ptr{Cvoid}}()

    if LP
        vtypes = repeat(GRB_CONTINUOUS, nvars)
        ub = ones(nvars)
        err = GRBnewmodel(gurobi_env, aptr, "CCLP", nvars, obj, C_NULL, ub, vtypes, C_NULL)
    else
        vtypes = repeat(GRB_BINARY, nvars)
        err = GRBnewmodel(gurobi_env, aptr, "ExactCC", nvars, obj, C_NULL, C_NULL, vtypes, C_NULL)
    end
    m = aptr[]
    timespent = time()-starttime
    timeremaining = timelimit-timespent
    tlim = max(0,timeremaining)
    GRBsetdblparam(GRBgetenv(m), GRB_DBL_PAR_TIMELIMIT, tlim)
    GRBsetintparam(GRBgetenv(m), "OutputFlag",outputflag)

    try

    cind = Int32[0,0,0]
    cval = Float64[0,0,0]
    for i = 1:n
        NeighbsI = findall(x->x>0,(A[:,i]))
        numNeighbs = size(NeighbsI,1)
        for u = 1:numNeighbs-1
            j = NeighbsI[u]
            for v = u+1:numNeighbs
                k = NeighbsI[v]
                if A[j,k] == 0
                    # Then we have a bad triangle: (i,j), (i,k) \in E
                    # but (j,k) is not an edge, so D(j,k) wants to be 1
                    #assert(i<j<k)
                    #@constraint(m,x[min(j,k),max(j,k)] - x[min(i,k),max(i,k)] - x[min(i,j),max(i,j)] <= 0)
                    cind[1] = vmap[min(j,k),max(j,k)]
                    cind[2] = vmap[min(i,k),max(i,k)]
                    cind[3] = vmap[min(i,j),max(i,j)]
                    cval[1] = 1
                    cval[2] = -1
                    cval[3] = -1
                    error = GRBaddconstr(m, 3, cind, cval, GRB_LESS_EQUAL, 0.0, C_NULL)
                end
            end
        end
    end

    # Find intial first solution
    if verbose
        println("First round of optimization")
    end
    #JuMP.optimize!(m)
    GRBoptimize(m)
    stat = Ref{Int32}(0)
    GRBgetintattr(m, GRB_INT_ATTR_STATUS, stat)
    status = stat[]
    if status == 2
        optimal = true
    else
        optimal = false
    end
    robj = Ref{Float64}(0.0)
    GRBgetdblattr(m, GRB_DBL_ATTR_OBJVAL, robj)

    obj = robj[] + (nvars - sum(A)/2)   # adjust by number of non-edges to get the actual objective

    soln = zeros(nvars)
    GRBgetdblattrarray(m, GRB_DBL_ATTR_X, 0, nvars, soln)
    D = zeros(n,n)
    for j=1:n
        for i=1:j-1
            D[i,j] = soln[vmap[i,j]+1]
            D[j,i] = soln[vmap[i,j]+1]
        end
    end

    stcplus_optimal = optimal
    stcplus_obj = obj
    stcplus_soln = copy(soln)
    stcplus_time = time()-starttime

    while true
        # x will naturally be upper triangular, but 'find_violations'  wants lower triangular
          #D = Matrix(JuMP.value.(x)')

         # Store violating tuples in a vector
         violations = Vector{Tuple{Int,Int,Int}}()

         find_violations!(D,violations)

         # Iterate through and add constraints
         numvi = size(violations,1)
         if verbose
            print("Adding in $numvi violated constraints.")
         end

         # Violations (a,b,c) are ordered so that a < b, and (a,b) is the value
         # that needs to be positive in the inequality:
         # x_ab - x_ac - x_bc <= 0.
         # The other two (b and c) could satisfy either b < c or c <= b.
         # We need to make sure we are only placing constraints in the upper
         # triangular portion of the matrix.
         # We do so by just calling min and max on pairs of nodes
         for v in violations
             #assert(v[1]<v[2])
             cind[1] = vmap[v[1],v[2]]
             cind[2] = vmap[min(v[1],v[3]),max(v[1],v[3])]
             cind[3] = vmap[min(v[2],v[3]),max(v[2],v[3])]
             cval[1] = 1
             cval[2] = -1
             cval[3] = -1
             #@show cind, cval
             error = GRBaddconstr(m, 3, cind, cval, GRB_LESS_EQUAL, 0.0, C_NULL)
             #z = sparsevec(cind.+1, cval, nvars)
             #@show z'*soln
             #@constraint(m,x[v[1],v[2]] - x[min(v[1],v[3]),max(v[1],v[3])] - x[min(v[2],v[3]),max(v[2],v[3])] <= 0)
         end
         if numvi == 0
             if verbose
                println(" Optimal solution found.")
             end
             break
         end
         if verbose
            println(" And re-solving the optimization problem.")
         end

         # update time remaining
         timespent = time()-starttime
         timeremaining = timelimit-timespent
         tlim = max(0,timeremaining)
         GRBsetdblparam(GRBgetenv(m), GRB_DBL_PAR_TIMELIMIT, tlim)
         if verbose
            println("time left until time is up = $timeremaining")
         end
         GRBupdatemodel(m)
         err = GRBoptimize(m)
         stat = Ref{Int32}(0)
         GRBgetintattr(m, GRB_INT_ATTR_STATUS, stat)
         #Status codes: https://www.gurobi.com/documentation/9.5/refman/optimization_status_codes.html#sec:StatusCodes
         status = stat[]
         # println("Status = $status")
         if status == 2
             optimal = true
         else
             optimal = false
         end
         robj = Ref{Float64}(0.0)
         GRBgetdblattr(m, GRB_DBL_ATTR_OBJVAL, robj)

         obj = robj[] + (nvars - sum(A)/2)   # adjust by number of non-edges to get the actual objective
         soln = zeros(nvars)
         GRBgetdblattrarray(m, GRB_DBL_ATTR_X, 0, nvars, soln)
         for j=1:n
             for i=1:j-1
                 D[i,j] = soln[vmap[i,j]+1]
                 D[j,i] = soln[vmap[i,j]+1]
             end
         end

     end
     return stcplus_optimal, stcplus_obj, stcplus_soln, stcplus_time, optimal, obj, Elist, soln
 finally
     #@show "freemodel here!"
     GRBfreemodel(m)
 end
end