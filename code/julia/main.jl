function main()
    # construct the basis (with reorganization into proper subspaces)
    # calculate the matrix of the Hamiltonian
    # diagonalize all the subspaces
    # get the ground state
    # calculate single spin flip spectral function
end

main()

A = BitArray(undef, 8)
A[[3,4,7]] .= true;
A
circshift(A, -3)

# # Notes:
# 1. If number of lattice sites n = 4k for integer k
#    then both highest and lowest magnetization states belong
#    to the same conserved subspace of the XXZ model.
#    This is not true for n = 4k + 2, i.e. the lowest and
#    highest magnetization states belong to different subspaces.
