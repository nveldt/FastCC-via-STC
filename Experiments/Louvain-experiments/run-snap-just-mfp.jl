include("run_mfp_functions.jl")
using MAT

snaps = readdir(homedir()*"/data/simple-snap")

for i = 1:length(snaps)
    graph = snaps[i]
    @show i, graph
end

##
snaps = readlines("snap-graph-names.txt")

for i = 1:length(snaps)
    graph = snaps[i]
    pivtimes = 1
    F = matread(homedir()*"/data/simple-snap/$graph")
    A = F["A"]
    n = size(A,1)
    m = round(Int64,sum(A)/2)
    
    num_tri, num_ow = count_wedges_triangles(A)
    println("$graph \t $n \t $m \t $num_ow \t $num_tri")
    graph_stats = [n;m;num_ow;num_tri]
    runtimes_cc, results_cc = run_matchflippivot_cc(A,pivtimes)
    runtimes_cd, results_cd = run_matchflippivot_cd(A,pivtimes)
    matwrite("snap-mfp/mfp_results_$(graph)_1piv.mat", Dict("graph_stats"=>graph_stats,"runtimes_cc"=>runtimes_cc,"results_cc"=>results_cc,"runtimes_cd"=>runtimes_cd,"results_cd"=>results_cd))

    cd_t = sum(runtimes_cd)
    cc_t = sum(runtimes_cc)
    cc_rat = results_cc[4]/results_cc[1]
    cd_rat = results_cd[2]/results_cd[1]
    println("\t \t $cd_t \t $cc_t \t")
    println("\t \t $cd_rat \t $cc_rat \t")

end