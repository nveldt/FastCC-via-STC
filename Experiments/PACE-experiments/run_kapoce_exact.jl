cd(homedir()*"/GitHubRepos/Working/cluster_editing/build/")

using MAT

#Killed after over half hour: 43, 51, 53

for j = 53:2:199
    mr = repeat("0",round(Int64,2-floor(log10(j))))     
    # mycmd = pipeline(`cat ../instances/exact/exact$mr$j.gr`,`./exact '>' exact_out/exact$j.txt`)
    mycmd = pipeline(`cat ../instances/exact/exact$mr$j.gr`,`./exact --time-limit=180`)

    tic = time()
    a = readchomp(mycmd)
    runtime = time()-tic
    edits = split(a, "\n") 
    mistakes = length(edits)
    println("$j $runtime $mistakes")
    matwrite("exact_out/exact$j.mat", Dict("runtime"=>runtime,"mistakes"=>mistakes,"edits"=>edits))
end