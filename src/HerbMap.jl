module HerbMap

using Revise
using CSV
using DataFrames
using LinearAlgebra 


include("clusternames.jl")
include("eigencluster.jl")
include("eigenclusterold.jl")
include("eigendistance.jl")
include("laplacian.jl")
include("sim_matrix.jl")
include("laplacianold.jl")

function load_ungulate_data()
    path_to_data = joinpath(@__DIR__, "..", "data", "Ungulate.data.dynamic.modeling.csv")
    return CSV.read(path_to_data, DataFrame)
end

function load_primate_data()
    path_to_data = joinpath(@__DIR__, "..", "data", "Primate.dental.data.dynamic.modeling.csv")
    return CSV.read(path_to_data, DataFrame)
end

# NOTE: LOOK CAREFULLY AT THIS!
# Function to convert a dissimilarity (distance) matrix to a Gram matrix using double-centering
function double_centering(D)
    n = size(D, 1)
    J = I - ones(n, n) / n  # Centering matrix
    B = -0.5 * J * (D .^ 2) * J
    return B
end

# Function to perform Classical MDS
function classical_mds(D, k)
    B = double_centering(D)  # Convert distance matrix to Gram matrix
    gevals, gevecs = eigen(B)  # Eigen decomposition of the Gram matrix
    # Take the top k eigenvectors corresponding to the largest eigenvalues
    X = gevecs[:, 1:k] * Diagonal(sqrt.(max.(gevals[1:k], 0)))
    return X  # The coordinates in k-dimensional space
end


# export any necessary functions
export clusternames, eigencluster, eigenclusterold, eigendistance, laplacian, laplacianold, sim_matrix, load_ungulate_data, load_primate_data, double_centering, classical_mds

end # module HerbMap