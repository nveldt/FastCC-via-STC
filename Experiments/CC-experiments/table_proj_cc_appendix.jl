using MAT


graphs = [
    "celegansneuralA";
    "Netscience"; 
    "Erdos991A";
    "Harvard500A";
    "celegansmetabolicA";
    "RogetA";
    "SmaGriA";
    "standard_emailA";
    "standard_polblogsA";
    "standard_ca-GrQc";
    "standard_caHepThA";  
]

psizes = zeros(length(graphs))
for i = 1:length(graphs)
    graph = graphs[i]
    
    # Load graph
    F = matread("../../data/smallgraphs/$graph.mat")
    A = F["A"]
    n = size(A,1)
    m = round(Int64,sum(A)/2)
    psizes[i] = m
end

using Random
p = sortperm(psizes)

## collect
function prun(runtime,dig)
    if runtime > 1
        rt = round(runtime,digits = 1)
    else
        rt = round(runtime, sigdigits = dig)
    end
    return rt
end

## 

for i = 1:length(graphs)

graph = graphs[p[i]]

# Load graph
F = matread("../../data/smallgraphs/$graph.mat")
A = F["A"]
n = size(A,1)
m = round(Int64,sum(A)/2)
dig = 3
dig2 = 4

# STC(plus)LP
M = matread("cc-metricopt-results/$(graph)_stcplus_lp.mat")
times = M["runtimes"] # [lb_time; time_round; time_round_206]
results = M["results"] # [stclp_bd, ub_stclp, ub_stclp_206, status, iscclp]
lb_stclp =  prun(results[1],dig)
run_stclp =  prun(times[1],dig)

# CCLP: lazy LP
M = matread("cc-metricopt-results/$(graph)_lazycc_lp.mat")
time_cclp = M["runtimes"] #[lb_time; time_rd]
res_cclp = M["results"]  # [bd, cclp_round, status]
lb_cclp =  prun(res_cclp[1],dig)
run_cclp = prun(time_cclp[1],dig)

# CCLP: Metric OPT
M = matread("cc-metricopt-results/$(graph)_metricopt_cc_lp.mat")
time_metricopt = M["runtimes"] #[timer; timer_r]
res_metricopt = M["results"]  # [bd, mopt_chaw_round, status, FinalCon]
lb_metricopt =  prun(res_metricopt[1],dig)
run_metricopt = prun(time_metricopt[1],dig)


gname = graph
if graph[end] == 'A'
    gname = graph[1:end-1]
end
if gname[1:4] == "stan"
    gname = gname[10:end]
end
println("\\midrule")
println("\\textsc{$gname}&\$n = $n\$  & LB & $lb_stclp & $lb_cclp & $lb_metricopt \\\\")
println(" & \$m = $m\$ & Run & $run_stclp & $run_cclp  & $run_metricopt \\\\")


end