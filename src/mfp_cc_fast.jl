"""
Find maximal node-pair disjoint set of open wedges in A.

This lower bounds cluster editing.
"""
function maximal_disjoint_openwedge_cc_dict(A)
    n = size(A,1)
    Neighbs = ConstructAdj(A,n)
    Add = Dict()
    Del = Dict()
    for i = 1:n
        N = Neighbs[i]
        for jj = 1:length(N)
            j = N[jj]
            ij1 = max(i,j)
            ij2 = min(i,j)
            if haskey(Del,(ij1,ij2))
                continue
            end

            for kk = jj+1:length(N)
                k = N[kk]
                ik1 = max(i,k)
                ik2 = min(i,k)

                jk1 = max(j,k)
                jk2 = min(j,k)

                if A[k,j] > 0
                    continue
                end
                if haskey(Del,(ik1,ik2)) || haskey(Add,(jk1,jk2))
                    continue

                else
                    # if we reach here: (i,j) and (i,k) are edges, but (j,k) is not
                    # So (i,j,k) is an open wedge centered at i.
                    # Also, neither (i,j) nor (i,k) nor (j,k) already appears in our open wedge set.

                    Del[(ik1,ik2)] = true
                    Del[(ij1,ij2)] = true
                    Add[(jk1,jk2)] = true

                    # Now exit this and start with a new j, as we have already
                    # used edge (i,j) and we don't need to check for more
                    # edges involving (i,j)
                    break
                end
            end
        end
    end
    Edel = collect(keys(Del))
    Eadd = collect(keys(Add))
    bd = (size(Edel,1)+size(Eadd,1))/3   # lower bound on the cluster editing objective
    return Edel, Eadd, round(Int64,bd)
end



"""
Given a set of edges to add and delete, create a
new graph with all these edges flipped
"""
function flip_graph(A,Eadd,Edel)
    n = size(A,1)
    Anew = copy(A)
    for k = 1:size(Edel,1)
        i = Edel[k][1]
        j = Edel[k][2]
        Anew[i,j] = 0
        Anew[j,i] = 0
    end
    dropzeros!(Anew)
    II,JJ = findnz(Anew)
    for k = 1:size(Eadd,1)
        push!(II,Eadd[k][1])
        push!(JJ,Eadd[k][2])
        push!(II,Eadd[k][2])
        push!(JJ,Eadd[k][1])
    end
    Anew = sparse(II,JJ,1,n,n)

    return Anew
end

"""
Rounds a feasible STC+ solution into a solution
    for a cluster editing problem.

Input:
    A = adjacency matrix

    Welist = node pairs to be flipped

Output:
    Approximation for CC.

    Step 1: Construct derived graph which is A with some flips
    Step 2: Apply pivot to the new graph
    Step 3: Compute CC objective score
    Step 4: Check how far this is from the lower bound

    pivtimes is the number of times to run the pivot
    method on the derived graph.

    Also returns the average time for a pivot step and the average
    quality solution
"""
function STCplus_to_CC_round_rand(A,Eadd, Edel,pivtimes = 5)

    tic = time()
    Anew = flip_graph(A,Eadd, Edel)
    setuptime = time()-tic

    # Apply pivot on this graph multiple times,
    # returning the best output
    clus = permutation_pivot(Anew)
    cc_obj = check_cc_obj(A,clus)
    objs = cc_obj
    tic = time()
    for jj = 1:pivtimes
        # Pivot is fast, so we can run it multiple times
        clusnew = permutation_pivot(Anew)
        objnew = check_cc_obj(A,clus) # Check output on original graph
        objs += objnew
        if objnew < cc_obj
            cc_obj = objnew
            clus = clusnew
        end
    end
    totaltime = time()-tic
    avg_time = totaltime/pivtimes + setuptime
    avg_obj = objs/pivtimes
    return round(Int64,cc_obj), clus, avg_obj, avg_time
end


"""
Rounding for the matching lower bound,
with deterministic rounding procedure.
"""
function STCplus_to_CC_round_det(A,Eadd,Edel)
    n = size(A,1)

    clus = deterministic_pivot_fast(A,Eadd,Edel)
    cc_obj = check_cc_obj(A,clus)

    return round(Int64,cc_obj), clus
end

"""
Deterministic pivoting for match-flip-pivot,
when given a set of edges to flip
"""
function deterministic_pivot_fast(A,Eadd,Edel)

    n = size(A,1)

    # create a new graph with the flipped edges
    B = flip_graph(A,Eadd, Edel)

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
