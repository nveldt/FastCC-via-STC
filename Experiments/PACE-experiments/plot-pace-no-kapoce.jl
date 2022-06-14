using Plots

# Let's plot as many heuristic stuff as we can

L = Vector{Float64}()   # louvain upper bounds
Lt = Vector{Float64}()  # louvain runtimes


MFP = Vector{Float64}()   # mfp upper bounds
MFPt = Vector{Float64}()  # mfp runtimes
Bds = Vector{Float64}()
S = Vector{Float64}()
V = Vector{Float64}()
E = Vector{Float64}()

for j = 1:2:199
    M = matread("pace_heuristic/heur$(j)_mfp_louvain.mat")

    graph_stats = M["graph_stats"] # [n;m;num_ow;num_tri]
    runtimes_cc = M["runtimes_cc"] # runtimes = [lb_time; time_rand; time_rand_avg; time_piv]
    results_cc = M["results_cc"]   # [bd; ub_rand_best; ub_rand_avg; ub_piv]
    louv_results = M["louv_results"] # [ll_cc, ll_cd]
    louv_time = M["louv_time"]       # [llcc_time, llcd_time]

    lb = results_cc[1]
    mfp = results_cc[4]
    mfpt = runtimes_cc[4]
    louv = louv_results[1]
    louvt = louv_time[1]

    push!(V, graph_stats[1])
    push!(E, graph_stats[2])
    push!(S, graph_stats[1]+graph_stats[2])

    push!(L,louv)
    push!(Lt,louvt)
    push!(MFP,mfp)
    push!(MFPt,mfpt)
    push!(Bds,lb)

end

## Plots
Lr = L./Bds
Mr = MFP./Bds
p = sortperm(Lr)
x = collect(1:length(L))
s1 = 300
s2 = 300
f = plot(title = "PACE graphs",grid = false,legend = :topleft, xlabel = "Problem instance", ylabel = "Approx Ratio (to MFP lower bd)")
scatter!(f,x,Mr[p],markerstrokewidth = 0, label = "MFP", size = (s1,s2))
scatter!(f,x,Lr[p],markerstrokewidth = 0, label = "Louvain")

## Runtime
x = E
f = plot(title = "PACE graphs", grid = false,legend = false, yaxis = :log10, xlabel = "|E|", ylabel = "Runtime (sec)")
scatter!(f,x,MFPt,markerstrokewidth = 0, label = "MFP", size = (s1,s2))
scatter!(f,x,Lt,markerstrokewidth = 0, label = "Louvain")
