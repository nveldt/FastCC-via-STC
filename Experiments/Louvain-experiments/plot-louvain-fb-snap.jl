using Plots

## Read in SNAPS

L = Vector{Float64}()   # louvain upper bounds
Lt = Vector{Float64}()  # louvain runtimes
LB = Vector{Float64}()  # lower bound runtimes


MFP = Vector{Float64}()   # mfp upper bounds
MFPt = Vector{Float64}()  # mfp runtimes
Bds = Vector{Float64}()
S = Vector{Float64}()
V = Vector{Float64}()
E = Vector{Float64}()

tt = "totaltime" # other option is roundtime: plots show only time it takes to round
cdce = "cd"      # choose between cluster editing (ce) or cluster deletion (cd)

snaps = readlines("snap-graph-names.txt")

for i = 1:length(snaps)
    global cdce, tt
    graph = snaps[i]
    if i == 24
        continue
    end
    F = matread("snap-mfp-1piv/mfp_results_$(graph)_1piv.mat")
    rcc = F["results_cc"]   # [cc_bd; ub_rand_best; ub_rand_avg; ub_piv]
    rncc = F["runtimes_cc"] # [lb_time; time_rand; time_rand_avg; time_piv]
    rcd = F["results_cd"]   # [cc_bd, ub, comp], # if comp = true, this means the rounding just took connected components after deleting edges
    rncd = F["runtimes_cd"] # [lb_time; time_round]
    graph_stats = F["graph_stats"] # [n,m,open wedges, triangles]

    M = matread("snap_louvain/louv_$(graph).mat")
    louv_results = M["louv_results"] # [ll_cc, ll_cd]
    louv_time = M["louv_time"]       # [llcc_time, llcd_time]

    if cdce == "ce"

        lb = rcc[1]
        mfp = rcc[4]
        
        lbtime = rncc[1]
        mfpt = rncc[4]

        louv = louv_results[1]
        louvt = louv_time[1]

    else

        cdce = "cd"

        lb = rcd[1]
        mfp = rcd[2]

        lbtime = rncd[1]
        mfpt = rncd[2]

        louv = louv_results[2]
        louvt = louv_time[2]   
    end

    if tt == "totaltime"
        mfpt += lbtime
        louvt += lbtime
    else
        tt == "rountime"
    end

    push!(V, graph_stats[1])
    push!(E, graph_stats[2])
    push!(S, graph_stats[1]+graph_stats[2])

    push!(L,louv)
    push!(Lt,louvt)
    push!(MFP,mfp)
    push!(MFPt,mfpt)
    push!(Bds,lb)
    push!(LB, lbtime)

end

# Plots
Lr = L./Bds
Mr = MFP./Bds
p = sortperm(Lr)
x = collect(1:length(L))
s1 = 300
s2 = 225
f = plot(title = "snap graphs",grid = false,legend = :false, xlabel = "Problem instance", ylabel = "Approx Ratio")
scatter!(f,x,Mr[p],markerstrokewidth = 0, label = "MFP", size = (s1,s2),color = :blue)
scatter!(f,x,Lr[p],markerstrokewidth = 0, label = "Louv+MFP",color = :orange)
savefig("Figures/snap_ratios_$(cdce)_$tt.pdf")


# Runtime
x = E
f = plot(title = "snap graphs", grid = false,legend = :topleft, yaxis = :log10, xaxis = :log10, xlabel = "|E|", ylabel = "Runtime (sec)")
scatter!(f,x,Lt,markerstrokewidth = 0, label = "Louv+MFP",color = :orange)
scatter!(f,x,MFPt,markerstrokewidth = 0, label = "MFP",color = :blue, size = (s1,s2))
if tt == "totaltime"
    scatter!(f,x,LB,markerstrokewidth = 0, label = "LB time",markershape = :star,markersize = 5,color = :green)
end
scatter!(f,x,MFPt,markerstrokewidth = 0, label = "",color = :blue, size = (s1,s2))


savefig("Figures/snap_runtime_$(cdce)_$tt.pdf")


## Now facebook

L = Vector{Float64}()   # louvain upper bounds
Lt = Vector{Float64}()  # louvain runtimes
LB = Vector{Float64}()  # lower bound runtimes


MFP = Vector{Float64}()   # mfp upper bounds
MFPt = Vector{Float64}()  # mfp runtimes
Bds = Vector{Float64}()
S = Vector{Float64}()
V = Vector{Float64}()
E = Vector{Float64}()

fbs = readlines("fb-graph-names.txt")

graphs = split(fbs[1])

cdce = "cd"
tt = "totaltime"

for i = 1:100
    global cdce, tt
    graph = graphs[i]

    F = matread("fb_louvain/$(graph)_fblouv_longround.mat")

    rcc = F["results_cc"]   # [cc_bd; ub_rand_best; ub_rand_avg; ub_piv]
    rncc = F["runtimes_cc"] # [lb_time; time_rand; time_rand_avg; time_piv]
    rcd = F["results_cd"]   # [cc_bd, ub, comp], # if comp = true, this means the rounding just took connected components after deleting edges
    rncd = F["runtimes_cd"] # [lb_time; time_round]
    graph_stats = F["graph_stats"] # [n,m,open wedges, triangles]

    louv_results = F["louv_results"] # [ll_cc, ll_cd]
    louv_time = F["louv_time"]       # [llcc_time, llcd_time]

    if cdce == "ce"

        lb = rcc[1]
        mfp = rcc[4]
        
        lbtime = rncc[1]
        mfpt = rncc[4]

        louv = louv_results[1]
        louvt = louv_time[1]

    else

        cdce = "cd"

        lb = rcd[1]
        mfp = rcd[2]

        lbtime = rncd[1]
        mfpt = rncd[2]

        louv = louv_results[2]
        louvt = louv_time[2]   
    end

    if tt == "totaltime"
        mfpt += lbtime
        louvt += lbtime
    else
        tt == "rountime"
    end

    
    push!(V, graph_stats[1])
    push!(E, graph_stats[2])
    push!(S, graph_stats[1]+graph_stats[2])

    push!(L,louv)
    push!(Lt,louvt)
    push!(MFP,mfp)
    push!(MFPt,mfpt)
    push!(Bds,lb)
    push!(LB,lbtime)

end

# Plots
Lr = L./Bds
Mr = MFP./Bds
p = sortperm(Lr)
x = collect(1:length(L))
s1 = 300
s2 = 225
f = plot(title = "fb100 graphs",grid = false,legend = :false, xlabel = "Problem instance", ylabel = "Approx Ratio")
scatter!(f,x,Mr[p],markerstrokewidth = 0, label = "MFP", size = (s1,s2),color = :blue)
scatter!(f,x,Lr[p],markerstrokewidth = 0, label = "Louv+MFP",color = :orange)
savefig("Figures/fb_ratios_$(cdce)_$(tt)_longround.pdf")


# Runtime
x = E
f = plot(title = "fb100 graphs", grid = false,legend = :topleft, yaxis = :log10, xaxis = :log10, xlabel = "|E|", ylabel = "Runtime (sec)")
scatter!(f,x,Lt,markerstrokewidth = 0, label = "Louv+MFP",color = :orange)
scatter!(f,x,MFPt,markerstrokewidth = 0, label = "MFP",color = :blue, size = (s1,s2))
if tt == "totaltime"
    scatter!(f,x,LB,markerstrokewidth = 0, label = "LB time",markershape = :star,markersize = 5,color = :green)
end
scatter!(f,x,MFPt,markerstrokewidth = 0, label = "",color = :blue, size = (s1,s2))


savefig("Figures/fb_runtime_$(cdce)_$(tt)_longround.pdf")
