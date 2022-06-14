include("run_mfp_functions.jl")
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
"com-Orkut.mat"
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



i = 7
graph = snaps[i]
pivtimes = 100
F = matread(homedir()*"/data/simple-snap/simple-$graph")
A = F["A"]
n = size(A,1)
m = round(Int64,sum(A)/2)

num_tri, num_ow = count_wedges_triangles(A)
println("$graph \t $n \t $m \t $num_ow \t $num_tri")
graph_stats = [n;m;num_ow;num_tri]
runtimes_cc, results_cc = run_matchflippivot_cc(A,pivtimes)
runtimes_cd, results_cd = run_matchflippivot_cd(A,pivtimes)

cd_t = sum(runtimes_cd)
cc_t = sum(runtimes_cc)
cc_rat = results_cc[4]/results_cc[1]
cd_rat = results_cd[2]/results_cd[1]
println("\t \t $cd_t \t $cc_t \t")
println("\t \t $cd_rat \t $cc_rat \t")
