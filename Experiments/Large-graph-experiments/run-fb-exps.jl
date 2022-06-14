include("run_mfp_functions.jl")
using MAT

# In order to run these experiments, access to Facebook100 datasets is needed
fbs = readlines("fb-results/Facebook_Sets.txt")
graphs = split(fbs[1])

for i = 1:length(graphs)
    graph = graphs[i]
    F = matread(homedir()*"/data/Facebook100/$graph.mat")
    A = F["A"]
    n = size(A,1)
    m = sum(A)/2
    d = sum(diag(A))
    mx = maximum(A.nzval)
    @assert(issymmetric(A))
    @assert(d == 0)
    if mx != 1
        Is, Js = findnz(A)
        A = sparse(Is,Js,1,n,n)
    end
    num_tri, num_ow = count_wedges_triangles(A)
    println("$graph \t $n \t $m \t $num_ow \t $num_tri")
    graph_stats = [n;m;num_ow;num_tri]
    pivtimes = 100
    runtimes_cc, results_cc = run_matchflippivot_cc(A,pivtimes)
    runtimes_cd, results_cd = run_matchflippivot_cd(A,pivtimes)
    matwrite("fb-results/mfp_results_$(graph).mat", Dict("graph_stats"=>graph_stats,"runtimes_cc"=>runtimes_cc,"results_cc"=>results_cc,"runtimes_cd"=>runtimes_cd,"results_cd"=>results_cd))

    cd_t = sum(runtimes_cd)
    cc_t = sum(runtimes_cc)
    cc_rat = results_cc[4]/results_cc[1]
    cd_rat = results_cd[2]/results_cd[1]
    println("\t \t $cd_t \t $cc_t \t")
    println("\t \t $cd_rat \t $cc_rat \t")

end

