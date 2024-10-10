# scripts/herb_diffusion_analysis.jl
using Revise

using HerbMap

using LinearAlgebra
using DataFrames 
using CSV 
using UnicodePlots 
using RCall 
using Combinatorics 
using MultivariateStats 
using ProgressMeter 
using SharedArrays 
using Arpack 
using Distances
using Plots
using ColorSchemes
using Distributions


data = load_ungulate_data()

sp = data[!,:species];
nsp = length(sp);

# Grab mass data
mass = data[!,:BM];

# Grab diet data
diet = data[!,:diet];
# clean up spaces
diet = strip.(diet);

# Grab diet2 data
diet2 = data[!,:diet2];
# clean up spaces
diet2 = strip.(diet2);

# Grab family data
family = data[!,:family];
family = strip.(family);

# Grab the list of measures
measures = Matrix(data[:, 7:end]);

# If we want to normalize to certain measures (body mass) we would do it here.
# raise mass to 1/3 because mass scales with the cube of linear dimensions
# SO THIS ASSUMES THESE MEASURES ARE LINEAR DIMENSIONS
# Area measures should be divided by mass^2/3
measures = measures ./ (mass .^(1/3));

# Diffusion function
evalues, evecs = diffusionmap(measures);

# # Build Similarity matrix
# PC = sim_matrix(measures);

# # Construct a laplacian matrix from the PC
# # Note: this matches the original version
# S = laplacian(PC,10);

# #Obtain the eigenvalues and eigenvectors
# ev = eigs(S; nev=10,which=:SR);

# # Subselect the Laplacian eigenvalues
# evalues = ev[1];

# # Subselect the Laplacian eigenvectors
# evecs = ev[2];

# # Assume evalues and evecs are from eigs(S)
# # Exclude the first eigenvalue and eigenvector
# evalues_nonzero = evalues[2:end]
# evecs_nonzero = evecs[:, 2:end]

# Scale eigenvectors by dividing by eigenvalues
# This should serve to emphasize the larger-scale patterns in the data
scaled_evecs = copy(evecs);
for i = 1:length(evalues)
    scaled_evecs[:, i] = evecs[:, i] / evalues[i];
end

# Alternatively, we could emphasize the smaller-scale patterns by multipltying:
scaled_evecs_minor = copy(evecs);
for i = 1:length(evalues)
    scaled_evecs_minor[:, i] = evecs[:, i] * evalues[i];
end

# What is the ranking of the second eigenvector?
ranked = sortperm(evecs[:,2]);
ecluster = eigencluster(collect(1:nsp),evecs,4);


# Visualize full eigenvector landscape with PHATE!
# From Ashkaan: 
# (1) Get eigenvectors from your diffusion map
# (2) Calculate pairwise distances in diffusion space
# (3) Perform principal coordinates analysis on diffusion distances

# NOTE: Violates PCA assumptions to feed in distance matrix
# Instead use Multidimensional Scaling (MDS) - see notes below

# Pairwise distances - Euclidean distance
pdist = Distances.pairwise(Euclidean(),scaled_evecs[:,2:10],dims=1);
pdist_minor = Distances.pairwise(Euclidean(),scaled_evecs_minor[:,2:10],dims=1);

# Jaccard
# Chebyshev
# BrayCurtis
# pcafit = fit(PCA,pdist; maxoutdim=2)
# tpdist_all = MultivariateStats.transform(pcafit, pdist)
# scatterplot(tpdist_all[1,:],tpdist_all[2,:])

#NOTE: we get different results for scaled_evecs vs. evecs

# PCA may not be appropriate to use on distance metrics
# Instead: Multidimensional Scaling (MDS):
	# •	Designed to handle distance or dissimilarity matrices.
	# •	Finds a configuration of points in a low-dimensional space that best preserves the pairwise distances.
	# •	Types of MDS:
	# •	Classical (Metric) MDS: Uses the distances directly.
	# •	Non-metric MDS: Uses the rank order of the distances.

# Perform Classical MDS to reduce the dimensionality to 2D
mds_result = MultivariateStats.classical_mds(pdist, 2);
scatterplot(mds_result[1,:],mds_result[2,:])

mds_result_minor = MultivariateStats.classical_mds(pdist_minor, 2);
scatterplot(mds_result_minor[1,:],mds_result_minor[2,:])

# PLOTTING
# Extract x and y coordinates from mds_result
x_coords = mds_result[1,:];  # First dimension
y_coords = mds_result[2,:];  # Second dimension

# Species names
species_names = sp

# PICK ONE
# Create clusters based on diffusion map and/or independent measures
# Assign cluster numbers based on ecluster
cluster_labels, num_clusters, cluster_legend, palette = clusterize(ecluster,nsp);
cluster_labels, num_clusters, cluster_legend, palette = clusterize(diet,nsp);
cluster_labels, num_clusters, cluster_legend, palette = clusterize(diet2,nsp);
cluster_labels, num_clusters, cluster_legend, palette = clusterize(family,nsp);


# Create the scatter plot - by ecluster
scatter(x_coords, y_coords,
		group = cluster_labels,        # Color points by cluster labels
		palette = palette,             # Use the defined color palette
		xlabel = "Dimension 1",
		ylabel = "Dimension 2",
		title = "Diffusion Map of Herbivore Morphology",
		legend = :outertopright,
		label = cluster_legend,
		markerstrokewidth = 0.5,         # Remove marker stroke
		markerstrokecolor = :black,    # Optional: outline markers
		markersize = 4)                # Optional: adjust marker size



# scatter!(x_coords, y_coords,
# group = dietcluster_labels,        # Color points by cluster labels
# palette = dietpalette,             # Use the defined color palette
# xlabel = "Dimension 1",
# ylabel = "Dimension 2",
# title = "Diffusion Map of Herbivore Morphology",
# legend = :outertopright,
# markerstrokewidth = 0,         # Remove marker stroke
# # markerstrokecolor = :black,    # Optional: outline markers
# markersize = 2)                # Optional: adjust marker size

# Option to include species names as labels - but should label clusters not spp
# This can clutter the plot if there are many species

# Uncomment the following lines to add labels to all points
# Indices of species to annotate
# label_indices = [1, 50, 100]  # Replace with desired indices

# # Create annotations as an array of tuples
# annotations = [ (x_coords[i], y_coords[i], species_names[i]) for i in label_indices ]

# # Create the scatter plot with annotations
# scatter(x_coords, y_coords,
#         group = cluster_labels,
#         palette = palette,
#         xlabel = "Dimension 1",
#         ylabel = "Dimension 2",
#         title = "MDS Plot with Annotations",
#         legend = :outertopright,
#         markersize = 8,
#         markerstrokewidth = 0,
#         annotations = annotations)


# Alternatively, label a subset of points (e.g., cluster centroids or specific species)
