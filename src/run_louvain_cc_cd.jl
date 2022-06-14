# include("../include/lambda_louvain.jl")
include("../include/lam-louv.jl")
include("../src/faster_many_pivot.jl")

"""
Run LambdaLouvain (cluster deletion version) for one graph A.

numtimes = number of times you run Louvain with different random node orderings.

maxits = number of iterations through nodes to allow

If maxits = 1, the code is set up so that you do not enter Phase 2 of the Louvain algorithm.
"""
function run_ll_cd(A,numtimes=10,maxits=1000)
    lam = 0.999999
    n = size(A,1)
    tic = time()
    c, lcc_obj = Many_Louvain(A,ones(n),lam,numtimes,maxits)
    timer = time()-tic

    Is, Js = findnz(triu(A))
    Elist = [Is Js]
    m = sum(A)/2

    mistakes = check_cc_obj_fastish(A,Elist,c,m)

    @assert(isCDfeasible(A,c))
    # @assert(m2 == mistakes)
    # @show m2, mistakes

    return mistakes, timer

end

"""
Run LambdaLouvain (cluster editing version) for one graph A.

numtimes = number of times you run Louvain with different random node orderings.

maxits = number of iterations through nodes to allow
"""
function run_ll_cc(A,numtimes=10,maxits=1000)
    lam = 0.5
    n = size(A,1)
    tic = time()
    c, lcc_obj = Many_Louvain(A,ones(n),lam,numtimes,maxits)
    timer = time()-tic

    Is, Js = findnz(triu(A))
    Elist = [Is Js]
    m = sum(A)/2

    mistakes = check_cc_obj_fastish(A,Elist,c,m)
    # @time mistakes = cc_cd_obj(A,c)
    # @assert(m2 == mistakes)

    return mistakes, timer
    
end
    