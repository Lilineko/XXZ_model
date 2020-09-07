using Plots
using LinearAlgebra

function main()
    systemSize = 10
    couplingJ = -1.0
    anisotropy = 0.5 # for anisotropy = 1 we have additional degeneration, thus allowed range is [0, 1)
    magnonInteractions = 1

    # construct the basis (with reorganization into proper subspaces)
    basis = constructBasis(systemSize)

    # calculate the matrix of the Hamiltonian
    blockHamiltonian = constructBlockHamiltonian(basis, systemSize, couplingJ, anisotropy, magnonInteractions)

    # diagonalize all the subspaces
    factorization = diagonalizeHamiltonian(blockHamiltonian)

    # get ground state info
    groundSubspaceIndex, groundStateSubspace = getGroundStateSubspace(factorization)
    groundStateEnergy = groundStateSubspace.values[1] # note: energies are sorted in ascending order
    groundStateVector = groundStateSubspace.vectors[:, 1]

    # apply spin flip to the ground state and calculate overlaps with eigenstates of the model
    overlaps = calculateOverlaps(groundSubspaceIndex, groundStateVector, factorization, basis, systemSize)

    # merge system info
    systemInfo = (groundSubspaceIndex, groundStateEnergy, factorization, overlaps)

    # calculate single spin flip spectral function
    kRange = [n for n in 1:systemSize] .* (2π / systemSize)
    ωRange = [n for n in -2:0.1:6]
    δ = 0.01
    kRange, ωRange, calculateSpectralFunction(kRange, ωRange, δ, systemInfo)
end

kRange, ωRange, A = main()

heatmap(kRange, ωRange, A, clim = (0, 0.01), colorbar = true)
