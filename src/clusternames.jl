function clusternames(numclusters)
    clusternamesvec = Array{String}(undef,numclusters);

    for i = 1:numclusters
        clusternamesvec[i] = "Cluster " * string(i)
    end

    return clusternamesvec
end
