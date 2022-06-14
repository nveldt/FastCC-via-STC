## Functions for rounding LP relaxations and bounds for
# cluster deletion
using MatrixNetworks
using Gurobi
using JuMP

"""
Rounds a linear program for cluster deletion into a
    feasible cluster deletion solution.

This works for the canonical LP relaxation, or the 
    weaker STC relaxation, and also works for rounding a 
    feasible solution to the STC ILP.

Input:
    A = adjacency matrix

    Elist = linear ordering of edges

    dist = solution to an LP relaxation (dist[k] corresponds to Elist[k,:])
            This is just 0's and 1's if this is rounding a feasible STC ILP solution.

Output:
    Approximation algorithm for CD. Approximation factor depends on where "dist" comes from.

    Step 1: Delete all edges e with x[e] >= 1/2
    Step 2: Apply pivot to the new graph
    Step 3: Compute CD objective score
    Step 3b: (Sanity check: all clusters are cliques)
    Step 4: Check how far this is from the lower bound
""" 
function CD_LP_round(A,Elist,dist,pivtimes = 1)
    n = size(A,1)
    # Find all the edges that we keep
    tol = 1e-8
    keep = findall(x->x<0.5-tol,dist)
    Inew = Elist[keep,1]
    Jnew = Elist[keep,2]
    Anew = sparse(Inew,Jnew,ones(length(Jnew)),n,n)
    Anew = Anew + Anew'

    # Round attempt 1: 
    sc = scomponents(Anew)
    c1 = sc.map

    if isCDfeasible(A,c1)
        comp = true
        clus = c1
        cd_obj = check_cd_obj(Elist,clus)
    else
        clus = permutation_pivot(Anew)
        comp = false
        cd_obj = check_cd_obj(Elist,clus)
        for jj = 1:pivtimes
            # Pivot is fast, so we can run it multiple times
            clusnew = permutation_pivot(Anew)
            objnew = check_cd_obj(Elist,clus)
            if objnew < cd_obj
                cd_obj = objnew
                clus = clusnew
            end
        end
    end

    return round(Int64,cd_obj), clus, comp
end

# Compute the cluster deletion objective
function check_cd_obj(Elist,clus)
    cd_obj = 0
    for k = 1:size(Elist,1)
        if clus[Elist[k,1]] != clus[Elist[k,2]] 
            cd_obj += 1
        end
    end
    return cd_obj
end
