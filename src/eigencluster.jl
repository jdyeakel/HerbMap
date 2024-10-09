function eigencluster(sp, evecs, n)
    # sp: Array of species labels
    # evecs: Matrix of eigenvectors (columns are eigenvectors)
    # n: Number of eigenvectors to use for clustering

    # Initialize carray with all species indices
    carray = [collect(1:size(evecs, 1))]  # Start with all species in one cluster

    # Optional: Print the ordered list of species according to the second eigenvector
    ordered_indices = sortperm(evecs[:, 2])
    println("Ordered List =", sp[ordered_indices])

    # Loop over the eigenvectors from i = 2 to n
    for i = 2:n
        eigenvec = evecs[:, i]

        # Initialize rarray to hold the new clusters
        rarray = Vector{Vector{Int}}()

        # For each cluster in carray
        for j = 1:length(carray)
            # Get the indices of species in the current cluster
            cluster_indices = carray[j]

            # Find species with positive and negative components in the current eigenvector
            positive_indices = cluster_indices[findall(x -> x > 0, eigenvec[cluster_indices])]
            negative_indices = cluster_indices[findall(x -> x < 0, eigenvec[cluster_indices])]

            # Add non-empty clusters to rarray
            if !isempty(positive_indices)
                push!(rarray, positive_indices)
            end
            if !isempty(negative_indices)
                push!(rarray, negative_indices)
            end
        end
        # Update carray with the new clusters
        carray = copy(rarray)
    end

    # Print the clusters
    for i = 1:length(carray)
        println()
        println("Cluster ", i, "=", sp[carray[i]])
    end

    return carray
end