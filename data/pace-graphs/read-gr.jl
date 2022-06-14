function readgr(filename)
    g = readlines(filename)
    U = Vector{Int64}()
    V = Vector{Int64}()
    for j = 2:length(g)
        edge = parse.(Int64,split(g[j]))
        push!(U,edge[1])
        push!(V,edge[2])
    end
    n = max(maximum(U),maximum(V))

    A = sparse(U,V,1,n,n)
    A = A'+A
    @assert(sum(diag(A)) ==0 )
    return A
end

# mr = repeat("0",round(Int64,2-floor(log10(j))))