using Gadfly
using LinearAlgebra

function main()
    systemSize = 6
    couplingJ = -1.0
    anisotropy = 0 # for anisotropy = 1 we have additional degeneration, thus allowed range is [0, 1)
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

    # apply spin flip to the ground state
    flippedSubspaceIndex, flippedState = applySpinFlipUp(4, groundStateVector, groundSubspaceIndex, basis, systemSize)

    # calculate single spin flip spectral function
end

main()

# A = BitArray(undef, 8)
# A[[3,4,7]] .= true;
# A
# circshift(A, -3)

# # Notes:
# 1. If number of lattice sites n = 4k for integer k
#    then both highest and lowest magnetization states belong
#    to the same conserved subspace of the XXZ model.
#    This is not true for n = 4k + 2, i.e. the lowest and
#    highest magnetization states belong to different subspaces.
