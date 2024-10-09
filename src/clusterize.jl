function clusterize(category,nsp)

    cluster_labels = zeros(Int, nsp);

    if isa(category, Vector{<:Vector})
        
        for cluster_number in 1:length(category)
            indices = category[cluster_number];
            cluster_labels[indices] .= cluster_number;
        end

    else
        
        unique_cat = unique(category)
        for cluster_number in 1:length(unique_cat)
            indices = findall(x->x==unique_cat[cluster_number],category);
            cluster_labels[indices] .= cluster_number;
        end

    end

    num_clusters = length(unique(cluster_labels));

    palette = ColorScheme(distinguishable_colors(num_clusters, transform=protanopic))

    return cluster_labels, num_clusters, palette

end