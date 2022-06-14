using MAT
snaps = [
"amazon0302.mat"
"amazon0312.mat"
"amazon0505.mat"
"amazon0601.mat"
"ca-AstroPh.mat"
"ca-CondMat.mat"
"ca-GrQc.mat"
"ca-HepPh.mat"
"ca-HepTh.mat"
"cit-HepPh.mat"
"cit-HepTh.mat"
"cit-Patents.mat"
"com-Amazon.mat"
"com-DBLP.mat"
"com-LiveJournal.mat"
"com-Youtube.mat"
"email-Enron.mat"
"email-EuAll.mat"
"loc-Brightkite.mat"
"loc-Gowalla.mat"
"roadNet-CA.mat"
"roadNet-PA.mat"
"roadNet-TX.mat"
"soc-Epinions1.mat"
"soc-LiveJournal1.mat"
"soc-Slashdot0811.mat"
"soc-Slashdot0902.mat"
"web-BerkStan.mat"
"web-Google.mat"
"web-NotreDame.mat"
"web-Stanford.mat"
"wiki-Talk.mat"
"wiki-topcats.mat"]

G = zeros(length(snaps),4)
CC = zeros(length(snaps),3)
CD = zeros(length(snaps),3)
T = zeros(length(snaps),2)

for i = 1:length(snaps)
    graph = snaps[i]
    F = matread("snap-results/mfp_results_simple-$(graph)")
    rcc = F["results_cc"]   # [cd_bd; ub_rand_best; ub_rand_avg; ub_piv]
    rncc = F["runtimes_cc"] # [lb_time; time_rand; time_rand_avg; time_piv]
    rcd = F["results_cd"]   # [cc_bd, ub, comp]
    rncd = F["runtimes_cd"] # [lb_time; time_round]
    stats = F["graph_stats"] # [n,m,open wedges, triangles]
    rat_cd = rcd[2]/rcd[1]
    run_cd = sum(rncd)
    rat_cc = rcc[4]/rcc[1]
    run_cc = rncc[1]+rncc[4]
    G[i,:] = stats
    CC[i,:] = [rcc[1] rcc[4] rat_cc]
    CD[i,:] = [rcd[1] rcd[2] rat_cd]
    T[i,:] = [run_cd run_cc]
end


## Now for FB graphs

f = readlines("fb-results/Facebook_Sets.txt")
fbs = split(f[1])

Gf = zeros(length(fbs),4)
CCf = zeros(length(fbs),3)
CDf = zeros(length(fbs),3)
Tf = zeros(length(fbs),2)

for i = 1:length(fbs)
    graph = fbs[i]
    F = matread("fb-results/mfp_results_$(graph)")
    rcc = F["results_cc"]   # [cd_bd; ub_rand_best; ub_rand_avg; ub_piv]
    rncc = F["runtimes_cc"] # [lb_time; time_rand; time_rand_avg; time_piv]
    rcd = F["results_cd"]   # [cc_bd, ub, comp]
    rncd = F["runtimes_cd"] # [lb_time; time_round]
    stats = F["graph_stats"] # [n,m, open wedges, triangles]
    rat_cd = rcd[2]/rcd[1]
    run_cd = sum(rncd)
    rat_cc = rcc[4]/rcc[1]
    run_cc = rncc[1]+rncc[4]
    Gf[i,:] = stats
    CCf[i,:] = [rcc[1] rcc[4] rat_cc]
    CDf[i,:] = [rcd[1] rcd[2] rat_cd]
    Tf[i,:] = [run_cd run_cc]
end

## Types
#           1       2           3            4        5             6             7           8
typs = ["o-soc", "web", "comm", "road", "prod", "collab", "cit","l-soc"]

snaptype = [5;5;5;5;6;6;6;6;6;7;7;7;5;6;1;1;3;3;8;8;4;4;4;1;1;1;1;2;2;2;2;3;2]

## CD vs CC plot
ms = 6.5
msw = 0
yl = "Clus Deletion approx"
xl = "Clus Editing approx"
P = plot(grid = false,  size = (500,270),legend = :topright, xtickfont=font(10),ytickfont=font(10))
colors = [:magenta, :red, :black, :green, :orange, :purple, :yellow, :cyan]
scatter!(P,CCf[:,3],CDf[:,3], color = :blue, 
        ylabel = yl, xlabel = xl, label = "fb-soc",
        markersize = ms,markerstrokewidth = msw)
for k = [ 8; collect(1:7)]
    J = findall(x->x==k,snaptype)
    # @show snaps[J]
    scatter!(P,CC[J,3],CD[J,3], color = colors[k], 
        ylabel = yl, xlabel = xl, label = typs[k],
        markersize = ms,markerstrokewidth = msw, )
end
J = findall(x->x==8,snaptype)
scatter!(P,CC[J,3],CD[J,3], color = colors[8], #xaxis = [1.8, 2.6],
        ylabel = yl, xlabel = xl, label = "",legendfont=font(9),
        markersize = ms,markerstrokewidth = msw, yaxis = [1.99,2.07])
P
savefig("CDvsCE_long.pdf")

## CD vs CC plot
ms = 6.5
msw = 0
yl = "Clus Deletion approx"
xl = "Clus Editing approx"
P = plot(grid = false,  size = (250,300),legend = :topright)
colors = [:magenta, :red, :black, :green, :orange, :purple, :yellow, :cyan]
scatter!(P,CDf[:,3],CCf[:,3], color = :blue, 
        ylabel = yl, xlabel = xl, label = "fb-soc",
        markersize = ms,markerstrokewidth = msw)
for k = [ 8; collect(1:7)]
    J = findall(x->x==k,snaptype)
    # @show snaps[J]
    scatter!(P,CD[J,3],CC[J,3], color = colors[k], 
        ylabel = yl, xlabel = xl, label = typs[k],
        markersize = ms,markerstrokewidth = msw, )
end
J = findall(x->x==8,snaptype)
scatter!(P,CD[J,3],CC[J,3], color = colors[8], 
        ylabel = yl, xlabel = xl, label = "",legendfont=font(9),
        markersize = ms,markerstrokewidth = msw, xaxis = [1.99,2.15])
P

## Objective score plots
P = plot(grid = false, size = (300,250), legend = false)
ys = CDf[:,2]./Gf[:,2]
ns = Gf[:,1]
Ns = ns .* (ns .- 1) /2
xs = CCf[:,2]./Gf[:,2]
ms = 5
msw = 0
scatter!(P,xs,ys, color = :blue, 
        ylabel = "CD/|E|", xlabel = "UCC/|E|", label = "fb-soc",
        markersize = ms,markerstrokewidth = msw)

ys = CD[:,2]./G[:,2]
xs = CC[:,2]./G[:,2]
for k = 1:8
    J = findall(x->x==k,snaptype)
    scatter!(P,xs[J],ys[J], color = colors[k], 
        label = typs[k], 
        markersize = ms,markerstrokewidth = msw)
end
P

## Lower bound plots
P = plot(grid = false, size = (300,250), legend = false)
ys = CDf[:,1]./Gf[:,2]
xs = CCf[:,1]./Gf[:,2]
ms = 5
msw = 0
scatter!(P,xs,ys, color = :blue, 
        ylabel = "CD/|E|", xlabel = "UCC/|E|", label = "fb-soc",
        markersize = ms,markerstrokewidth = msw)
ys = CD[:,1]./G[:,2]
xs = CC[:,1]./G[:,2]
for k = 1:8
    J = findall(x->x==k,snaptype)
    scatter!(P,ys[J],xs[J], color = colors[k], 
        label = typs[k], 
        markersize = ms,markerstrokewidth = msw)
end
P
savefig("LB_CDvsUCC.pdf")

## Runtime plots for MFP-CD 
ys = Tf[:,1]
xs = Gf[:,2]
shp = :circle
ms = 5
msw = 0
P = plot(grid = false, size = (300,250),legend = false) #, title = "MFP-CD runtime")
scatter!(P,xs,ys, color = :blue, yaxis = :log10, xaxis = :log10,markeralpha = 0.7,
        xlabel = "|E|", ylabel = "runtime (s)", label = "fb-soc",
        markersize = ms,markerstrokewidth = msw, markershape = shp)
ys = T[:,1]
xs = G[:,2]
for k = 1:8
    J = findall(x->x==k,snaptype)
    scatter!(P,xs[J],ys[J], color = colors[k], markeralpha = 0.7,
        label = typs[k], markershape = shp, markersize = ms, 
        markerstrokewidth = msw)
end
P
savefig("Run_MFPCD.pdf")


## Runtime plots for MFP-CC
ys = Tf[:,2]
Ms = Gf[:,2]
ns = Gf[:,1]
Ns = ns .* (ns .- 1) /2
xs = Ms
shp = :circle
ms = 5
msw = 0
P = plot(grid = false,size = (300,250), legend = false) # title = "MFP-CC runtime")
scatter!(P,xs,ys, color = :blue, 
        yaxis = :log10, xaxis = :log10,
        # xaxis = [1,1e6], yaxis = [0,2],
        markeralpha = 0.7,
        xlabel = "|E|", ylabel = "runtime (s)", label = "fb-soc",
        markersize = ms,markerstrokewidth = msw, markershape = shp)
ys = T[:,2]
Ms = G[:,2]
ns = G[:,1]
Ns = ns .* (ns .- 1) /2
xs = Ms
for k = 1:8
    J = findall(x->x==k,snaptype)
    scatter!(P,xs[J],ys[J], color = colors[k], markeralpha = 0.7,
        label = typs[k], markershape = shp, markersize = ms, 
        markerstrokewidth = msw)
end
P
savefig("Run_MFPCC.pdf")



