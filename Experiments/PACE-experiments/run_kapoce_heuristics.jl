cd(homedir()*"/GitHubRepos/Working/cluster_editing/build/")

using MAT

for j = 43:2:199
    mr = repeat("0",round(Int64,2-floor(log10(j))))     
    mycmd = pipeline(`cat ../instances/heur/heur$mr$j.gr`,`./heuristic --time-limit=600`)

    tic = time()
    a = readchomp(mycmd)
    runtime = time()-tic
    edits = split(a, "\n")
    mistakes = length(edits)
    println("$j $runtime $mistakes")
    matwrite("heur_out/heur$j.mat", Dict("runtime"=>runtime,"mistakes"=>mistakes,"edits"=>edits))
end