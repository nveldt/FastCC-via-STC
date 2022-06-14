using Plots
using MAT

# Let's plot as many heuristic stuff as we can

L = Vector{Float64}()   # louvain upper bounds
Lt = Vector{Float64}()  # louvain runtimes

LPr = Vector{Float64}()   # LP ratio
LPt = Vector{Float64}()  # LP runtimes

H = Vector{Float64}()   # heuristic upper bounds
Ht = Vector{Float64}()  # heuristic runtimes

MFP = Vector{Float64}()   # mfp upper bounds
MFPt = Vector{Float64}()  # mfp runtimes

LB = Vector{Float64}() # Time to computer lower bound

Bds = Vector{Float64}()
S = Vector{Float64}()
V = Vector{Float64}()
E = Vector{Float64}()

for j = 1:2:57
    M = matread("pace_heuristic/heur$(j)_mfp_louvain.mat")

    graph_stats = M["graph_stats"]   # [n;m;num_ow;num_tri]
    runtimes_cc = M["runtimes_cc"]   # runtimes = [lb_time; time_rand; time_rand_avg; time_piv]
    results_cc = M["results_cc"]     # [bd; ub_rand_best; ub_rand_avg; ub_piv]
    louv_results = M["louv_results"] # [ll_cc, ll_cd]
    louv_time = M["louv_time"]       # [llcc_time, llcd_time]

    # if graph_stats[1] < 100
    #     # Ignore small graphs
    #     continue
    # end
    lb = results_cc[1]
    lb_time = runtimes_cc[1]
    mfp = results_cc[4]
    mfpt = runtimes_cc[4] + lb_time
    louv = louv_results[1]
    louvt = louv_time[1] + lb_time

    push!(V, graph_stats[1])
    push!(E, graph_stats[2])
    push!(S, graph_stats[1]+graph_stats[2])

    push!(L,louv)
    push!(Lt,louvt)
    push!(MFP,mfp)
    push!(MFPt,mfpt)
    push!(Bds,lb)
    push!(LB,lb_time)

    M2 = matread("heur_out/heur$j.mat")
    Heur = M2["mistakes"]
    Heurt = M2["runtime"]
    push!(H,Heur)
    push!(Ht,Heurt + lb_time)

    M3 = matread("pace_heuristic/heur$(j)_stcplus_lp.mat")
    runt = M3["runtimes"]
    res = M3["results"]
    push!(LPr,res[3]/res[1])
    push!(LPt, runt[1]+runt[3])
end

## Plots

Lr = L./Bds
Mr = MFP./Bds
Hr = H./Bds
p = sortperm(Hr)
x = collect(1:length(L))
s1 = 300
s2 = 200
f = plot(title = "PACE graphs",grid = false,legend = :topleft, xlabel = "Problem instance", ylabel = "Approx Ratio") #,foreground_color_legend = nothing) #, background_color_legend = nothing)
scatter!(f,x,Mr[p],markerstrokewidth = 0, label = "MFP", size = (s1,s2),color = :blue)
scatter!(f,x,Lr[p],markerstrokewidth = 0, label = "Louv+MFP",color = :orange)
scatter!(f,x,Hr[p],markerstrokewidth = 0, label = "KaPoCE+MFP",color = :purple,yaxis = [1,2.4])
# scatter!(f,x,LPr[p],markerstrokewidth = 0, label = "LP-STC+",color = :green)
savefig("Figures/pace_ratios.pdf")


## Runtime
x = E
f = plot(title = "PACE graphs", grid = false,legend = false,foreground_color_legend = nothing, background_color_legend = nothing, yaxis = :log10, xaxis = :identity,xlabel = "|E|", ylabel = "Runtime (sec)")
scatter!(f,x,Ht,markerstrokewidth = 0, label = "KaPoCE+MFP",color = :purple)
scatter!(f,x,Lt,markerstrokewidth = 0, label = "Louv+MFP",color = :orange)
scatter!(f,x,MFPt,markerstrokewidth = 0, label = "MFP", size = (s1,s2),color = :blue)
# scatter!(f,x,LB,markerstrokewidth = 0, markershape = :star,color = :green, label = "LB Time")
# scatter!(f,x,LPt,markerstrokewidth = 0, label = "LP-STC+",color = :green)
savefig("Figures/pace_runtimes.pdf")


## Show sizes of graphs where heuristic method is slow
for i = 1:29
    println("$(V[i]) $(E[i]) $(Ht[i])")
end
