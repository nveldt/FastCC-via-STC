# Graphs we consider for CC and CD experiments
graphs = [
    "Harvard500A";
    "Erdos991A"; 
    "celegansneuralA";
    "Netscience";
    "celegansmetabolicA";
    "RogetA";
    "SmaGriA";
    "standard_emailA";
    "standard_polblogsA";
    "standard_ca-GrQc";
    "standard_caHepThA";
    "standard_EmailEnronA";
    "standard_condmat2005A";
    "standard_ca-AstroPhA";
    "standard_loc-Brightkite";
    "Caltech36";
    "Reed98";
    "Simmons81";
    "Haverford76";
    "Swarthmore42";
    "Amherst41";
    "Bowdoin47";
    "Rice31";
    "Lehigh96";
]

using MAT

## Printing results for one graph
graphs = ["standard_EmailEnronA", "standard_condmat2005A", "standard_ca-AstroPhA","standard_loc-Brightkite"]
for i = 1:4
graph = graphs[i]

# Load graph
F = matread("../../data/smallgraphs/$graph.mat")
A = F["A"]
n = size(A,1)
m = round(Int64,sum(A)/2)
dig = 3
dig2 = 4

M = matread("cd-experiment-results/$(graph)_mfp.mat")
time_mfp = M["runtimes"]
res_mfp = M["results"]
lb_mfp = res_mfp[1]
ub_mfp = res_mfp[2]
rat_mfp = round(ub_mfp/lb_mfp,sigdigits = dig2)
run_mfp = round(sum(time_mfp),sigdigits = dig)

M = matread("cd-experiment-results/$(graph)_stclp.mat")
time_stclp = M["runtimes"]
res_stclp = M["results"]
lb_stclp = res_stclp[1]
ub_stclp = res_stclp[2]
rat_stclp = round(ub_stclp/lb_stclp,sigdigits = dig2)
runtime_stclp = round(sum(time_stclp),sigdigits = dig)
status_stclp = res_stclp[4] # status--if false, we didn't converge
iscd_stclp = res_stclp[5]   # if true, this is the cd solution too!

status_cdlp = 0.0
lb_cdlp = 0.0
ub_cdlp = 0.0
runtime_cdlp = 0.0
rat_cdlp = 0.0
try
    M = matread("cd-experiment-results/$(graph)_cdlp.mat")
    time_cdlp = M["runtimes"]
    res_cdlp = M["results"]
    lb_cdlp = res_cdlp[1]
    ub_cdlp = res_cdlp[2]
    rat_cdlp = round(ub_cdlp/lb_cdlp,sigdigits = dig2)
    runtime_cdlp = round(sum(time_cdlp),sigdigits = dig)
    status_cdlp = res_cdlp[4] # status--if false, we didn't converge
catch
    status_cdlp = 0.0
end


if status_stclp == 0.0
    lb_stclp = "--"
    ub_stclp = "--"
    rat_stclp = "--"
    runtime_stclp = "--"
    @assert(iscd_stclp == 0.0)
end

if status_cdlp == 0.0
    lb_cdlp = "--"
    ub_cdlp = "--"
    rat_cdlp = "--"
    runtime_cdlp = "--"
end


gname = graph
if graph[end] == 'A'
    gname = graph[1:end-1]
end
if gname[1:4] == "stan"
    gname = gname[10:end]
end

N = n
if n > 10000
    N = round(Int64,n/1000)
    N = "$(N)k"
end
M = m
if m > 10000
    M = round(Int64,m/1000)
    M = "$(M)k"
end

# Print out for LaTeX table
println("\\midrule")
if iscd_stclp == 1.0
println("\\textsc{$gname} & LB & $lb_mfp & $lb_stclp\$^*\$ & $lb_cdlp \\\\")
else
    println("\\textsc{$gname} & LB & $lb_mfp & $lb_stclp & $lb_cdlp  \\\\")
end
println("& UB & $ub_mfp & $ub_stclp & $ub_cdlp  \\\\")
println("\$n = $n\$ & Ratio & $rat_mfp & $rat_stclp & $rat_cdlp \\\\")
println("\$m = $m\$ & Run & $run_mfp & $runtime_stclp & $runtime_cdlp \\\\")



end