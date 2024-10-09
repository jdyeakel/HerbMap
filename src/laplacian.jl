function laplacian(PC, branchnum)
    nsp = size(PC)[1]  # Number of species (nodes)
    maxdiff = min(branchnum, nsp - 1)  # Max number of connections per node
    S = zeros(Float64, nsp, nsp)  # Initialize Laplacian matrix

    for i = 1:nsp
        val = fill(-Inf, maxdiff)  # Initialize to negative infinity
        lok = zeros(Int64, maxdiff)

        for j = 1:nsp
            if i == j
                continue  # Skip self-comparison
            end

            # Check if current similarity is larger than the smallest in val
            if PC[i, j] > val[maxdiff]
                val[maxdiff] = PC[i, j]
                lok[maxdiff] = j
            end

            # Insertion sort to maintain descending order in val
            for k = 1:maxdiff - 1
                if val[maxdiff + 1 - k] > val[maxdiff - k]
                    # Swap values in val
                    v = val[maxdiff + 1 - k]
                    val[maxdiff + 1 - k] = val[maxdiff - k]
                    val[maxdiff - k] = v
                    # Swap indices in lok
                    l = lok[maxdiff + 1 - k]
                    lok[maxdiff + 1 - k] = lok[maxdiff - k]
                    lok[maxdiff - k] = l
                end
            end
        end

        # Assign negative similarities to off-diagonal entries in S
        S[i, lok] = -PC[i, lok]
        S[lok, i] = -PC[lok, i]  # Ensure symmetry
    end

    # Set diagonal entries to ensure zero row sums
    rowsums = sum(S, dims=2)
    S[diagind(S)] = -rowsums

    return S
end
