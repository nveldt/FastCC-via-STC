using MAT

using SparseArrays

for j = 1:2:199
mr = repeat("0",round(Int64,2-floor(log10(j))))

    A = readgr("exact/exact$mr$j.gr")


    matwrite("exactmat/ex$j.mat", Dict("A"=>A))
end



## Read and check density

for j = 1:2:199
    M = matread("exactmat/ex$j.mat")
    A = M["A"]
    n = size(A,1)
    m = sum(A)/2
    sparsity = round(m/binomial(n,2),digits = 2)
    println("$n $m $sparsity")
end



## Do same for heuristics

for j = 1:2:199
    mr = repeat("0",round(Int64,2-floor(log10(j))))
    
        A = readgr("heur/heur$mr$j.gr")
    
    
        matwrite("heurmat/heur$j.mat", Dict("A"=>A))
end

## read density

for j = 1:2:199
    M = matread("heurmat/heur$j.mat")
    A = M["A"]
    n = size(A,1)
    m = sum(A)/2
    @assert(issymmetric(A))
    sparsity = round(m/binomial(n,2),digits = 2)
    println("$n $m $sparsity")
end
