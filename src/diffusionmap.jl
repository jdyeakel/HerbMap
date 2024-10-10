function diffusionmap(measures)

    # Build Similarity matrix
    PC = sim_matrix(measures);

    # Construct a laplacian matrix from the PC
    # Note: this matches the original version
    S = laplacian(PC,10);

    #Obtain the eigenvalues and eigenvectors
    ev = eigs(S; nev=10,which=:SR);

    # Subselect the Laplacian eigenvalues
    evalues = ev[1];

    # Subselect the Laplacian eigenvectors
    evecs = ev[2];

    return evalues, evecs

end