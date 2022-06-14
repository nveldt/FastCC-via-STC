include("run_mfp_functions.jl")
using MAT

snaps = readlines("snap-graph-names.txt")

for i = 1:length(snaps)
    graph = snaps[i]

    # Update to point to snap graph folder
    F = matread(homedir()*"/data/simple-snap/$graph")
    A = F["A"]
    n = size(A,1)
    m = round(Int64,sum(A)/2)
   
    numtimes = 1
    maxits = 1
    ll_cc, llcc_time = run_ll_cc(A,numtimes,maxits)
    ll_cd, llcd_time = run_ll_cd(A,numtimes,maxits)
    louv_results = [ll_cc, ll_cd]
    louv_time = [llcc_time, llcd_time]
    println("$graph \t $llcc_time \t $llcd_time")
    matwrite("snap_louvain/louv_$(graph).mat", Dict("louv_results"=>louv_results, "louv_time"=>louv_time))

end