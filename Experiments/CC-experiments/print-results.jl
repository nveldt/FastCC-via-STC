using MAT

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
    # "Caltech36";
    # "Reed98";
    # "Simmons81";
    # "Haverford76";
    # "Swarthmore42";
    # "Amherst41";
    # "Bowdoin47";
    # "Rice31";
    # "Lehigh96";
]


## collect
function prun(runtime,dig)
    if runtime > 1
        rt = round(runtime,digits = 1)
    else
        rt = round(runtime, sigdigits = dig)
    end
    return rt
end

for i = 1:length(graphs)
graph = graphs[i]

# Load graph
F = matread("../../data/smallgraphs/$graph.mat")
A = F["A"]
n = size(A,1)
m = round(Int64,sum(A)/2)
dig = 3
dig2 = 4


# MFP-CC
M = matread("cc-experiment-results/$(graph)_mfp_cc_det.mat")

# [lb_time; time_rand; time_rand_avg; time_det; time_piv]
times = M["runtimes"] 

# [bd; ub_rand_best; ub_rand_avg; ub_det; ub_piv]
results = M["results"]   

lb_mfp = round(Int64,results[1])
ub_mfp = round(Int64,results[2])
ub_rand_avg = results[3]
# ub_det = round(Int64,results[4])
ub_dmfp = round(Int64,results[4])

rat_mfp = round(ub_mfp/lb_mfp,sigdigits = dig2)
rat_dmfp = round(ub_dmfp/lb_mfp,sigdigits = dig2)

run_mfp = prun(sum(times[1:2]),dig)
run_dmfp = prun(times[1]+times[4],dig)


# pivot
ub_piv = round(Int64,results[5])
run_piv = prun(times[5], dig)
run_piv_mfp = prun(times[1]+times[5],dig)
rat_piv_mfp = round(ub_piv/lb_mfp,sigdigits = dig2)
rat_piv = "--"

# deterministic rounding
ub_det = round(Int64,results[4])
run_dmfp = prun(sum(times[[1,4]]), dig)
rat_dmfp = round(ub_dmfp/lb_mfp,sigdigits = dig2)

# STC(plus)LP
oom = "--"
lb_stclp = oom
ub_stclp_206 = oom
ub_stclp_rand = oom
rat_stclp_rand = oom
run_stclp_rand = oom
rat_stclp_206 = oom
run_stclp_206 = oom
status_stclp = 1.0
iscc_stclp = 0.0
status_cclp = 1.0
run_cclp = 0.0
lb_cclp = 0.0
ub_cclp = 0.0
rat_cclp = 0.0
try
    M = matread("cc-experiment-results/$(graph)_stcplus_cc_combined.mat")
    times = M["runtimes_stclp"] # [lb_time; time_round; time_round_206]
    results = M["results_stclp"] # [stclp_bd, ub_stclp, ub_stclp_206, status, iscclp]
    lb_stclp = round(results[1],digits = 1)
    ub_stclp_rand = round(Int64,results[2])
    ub_stclp_206 = round(Int64,results[3])

    rat_stclp_rand = round(ub_stclp_rand/lb_stclp,sigdigits = dig2)
    rat_stclp_206 = round(ub_stclp_206/lb_stclp,sigdigits = dig2)

    run_stclp_rand = prun(times[1]+times[2],dig)
    run_stclp_206 = prun(times[1]+times[3],dig)
    status_stclp = results[4] # status--if false, we didn't converge
    iscc_stclp = results[5]   # if true, this is the cc solution too!

    times_cclp = M["runtimes_cclp"] # [lb_time_cclp; time_rd]
    res_cclp = M["results_cclp"] # [lb_cclp, cclp_round, status_cclp]
    lb_cclp = round(res_cclp[1],digits = 1)
    ub_cclp = round(Int64,res_cclp[2])
    rat_cclp = round(ub_cclp/lb_cclp,sigdigits = dig2)
    run_cclp = prun(sum(times_cclp),dig)
    status_cclp = res_cclp[3] # status--if false, we didn't converge
catch
    oom = "--"
    lb_stclp = oom
    ub_stclp_206 = oom
    ub_stclp_rand = oom
    rat_stclp_rand = oom
    run_stclp_rand = oom
    rat_stclp_206 = oom
    run_stclp_206 = oom
    status_stclp = 1.0
    lb_cclp = oom
    ub_cclp= oom
    rat_cclp = oom
    run_cclp = oom
    status_cclp = 1.0
end

oot = "Timeout"
if status_stclp == 0.0
    lb_stclp = oot
    ub_stclp_rand = oot
    rat_stclp_rand = oot
    ub_stclp_206 = oot
    rat_stclp_206 = oot
    @assert(run_stclp_rand > 1800)
    run_stclp_rand = oot
    run_stclp_206 = oot
end

if status_cclp == 0.0
    lb_cclp = oot
    ub_cclp = oot
    rat_cclp = oot
    run_cclp = oot
end

gname = graph
if graph[end] == 'A'
    gname = graph[1:end-1]
end
if gname[1:4] == "stan"
    gname = gname[10:end]
end

println("\\midrule")
if iscc_stclp == 1.0 && run_cclp != oot
println("\\textsc{$gname} & LB & $lb_mfp & & & $lb_stclp\$^*\$ & &$lb_cclp\\\\")
else
    println("\\textsc{$gname} & LB & $lb_mfp & & & $lb_stclp  & &$lb_cclp  \\\\")
end
println("            & UB & $ub_mfp & $ub_dmfp & $ub_piv & $ub_stclp_rand & $ub_stclp_206 & $ub_cclp  \\\\")
println("\$n = $n \$ & Ratio & $rat_mfp & $rat_dmfp & $rat_piv_mfp & $rat_stclp_rand & $rat_stclp_206 & $rat_cclp  \\\\")
println(" \$m = $m\$ & Run & $run_mfp & $run_dmfp & $run_piv_mfp  & $run_stclp_rand & $run_stclp_206 & $run_cclp  \\\\")


end