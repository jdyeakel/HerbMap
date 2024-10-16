module HerbMap

using CSV
using DataFrames
using LinearAlgebra 
using ColorSchemes
using Plots
using Arpack 


include("clusternames.jl")
include("eigencluster.jl")
include("eigenclusterold.jl")
include("eigendistance.jl")
include("laplacian.jl")
include("sim_matrix.jl")
include("laplacianold.jl")
include("clusterize.jl")
include("diffusionmap.jl")

function load_ungulate_data()
    path_to_data = joinpath(@__DIR__, "..", "data", "Ungulate.data.dynamic.modeling.csv")
    return CSV.read(path_to_data, DataFrame)
end

function load_primate_data()
    path_to_data = joinpath(@__DIR__, "..", "data", "Primate.dental.data.dynamic.modeling.csv")
    data = CSV.read(path_to_data, DataFrame);
    return data
end

# # NOTE: THIS FUNCTION IS NOT BEING USED
# # Function to convert a dissimilarity (distance) matrix to a Gram matrix using double-centering
# function double_centering(D)
#     n = size(D, 1)
#     J = I - ones(n, n) / n  # Centering matrix
#     B = -0.5 * J * (D .^ 2) * J
#     return B
# end

# # NOTE: THIS FUNCTION IS NOT BEING USED
# # Function to perform Classical MDS
# function classical_mds(D, k)
#     B = double_centering(D)  # Convert distance matrix to Gram matrix
#     gevals, gevecs = eigen(B)  # Eigen decomposition of the Gram matrix
#     # Take the top k eigenvectors corresponding to the largest eigenvalues
#     X = gevecs[:, 1:k] * Diagonal(sqrt.(max.(gevals[1:k], 0)))
#     return X  # The coordinates in k-dimensional space
# end

# function classical_mds(D, k)
#     B = double_centering(D)  # Convert distance matrix to Gram matrix
    
#     # Perform eigen decomposition
#     gevals, gevecs = eigen(B)
    
#     # Take only the real part of the eigenvalues and ensure non-negative before taking the square root
#     real_gevals = real(gevals)
#     X = gevecs[:, 1:k] * Diagonal(sqrt.(max.(real_gevals[1:k], 0.0)))
    
#     return X  # The coordinates in k-dimensional space
# end



# export any necessary functions
export 

load_ungulate_data, 
load_primate_data,

sim_matrix,
laplacian,
eigendistance,
eigencluster,
diffusionmap,

clusterize, 
clusternames,

eigenclusterold,
laplacianold

# double_centering, classical_mds,

end # module HerbMap