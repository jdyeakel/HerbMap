function clusterize(category,nsp)

    cluster_labels = zeros(Int, nsp);
    unique_cat = String[];

    if isa(category, Vector{<:Vector})
        
        unique_cat = ["C$i" for i in 1:length(category)];
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

    cluster_legend = reshape(string.(unique_cat), 1, :);

    return cluster_labels, num_clusters, cluster_legend, palette

end