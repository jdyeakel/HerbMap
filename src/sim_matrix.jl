function sim_matrix(measures)
    
    nsp = size(measures)[1];
    nmeas = size(measures)[2];


    # Build Pairwise Comparisons
    PC = Array{Float64}(undef,nsp,nsp);

    
    #Build similarity matrix
    # 0 = Similar
    # 1 = dissimilar

    for i = 0:(nsp^2 - 1)
        a = mod(i,nsp) + 1;
        b = Int64(floor(i/nsp)) + 1;

        # Skip calculation if comparing the same species
        if a == b
            PC[a,b] = 1.0; # Set similarity to 1 for identical species
            continue
        end

        # If not the same species, we will go through the measures for each pair and build an index of similarity 
        
        # Initialize accumulators for log ratios and valid measure count
        ct = 0; # Sum of log(min/max) ratios
        ctones = 0; # Count of measures where both species have data
        
        # Loop over each measure to compute similarity
        for j = 1:nmeas

            # Check if both species have data for the current measure
            if !ismissing(measures[a,j]) && !ismissing(measures[b,j])
                
                # Add the log ratio of min to max measure value to ct
                ct += log(minimum([measures[a,j],measures[b,j]])/maximum([measures[a,j],measures[b,j]]));

                ctones += 1; # Increment the count of valid measures
            end
        end
        # Calculate the geometric mean of the min/max ratios
        ctscaled = exp(ct/ctones);
        # Assign the similarity score to the matrix (values between 0 and 1)
        PC[a,b] = Float64(ctscaled); 
    end

    return PC
end
