function getMagnetization(state, systemSize)
    systemSize/2 - sum(digits(state, base = 2))
end

function getStaggeredMagnetization(state, systemSize)
    spinRepresentation = digits(state, base = 2, pad = systemSize) .- 0.5
    filter = [(-1)^(it - 1) for it in 1:systemSize]
    sum(filter .* spinRepresentation)
end

function constructBasis(systemSize)
    conservedSubspaces = Vector{Vector{Int64}}(undef, 4)
    conservedSubspaces[1:4] = [[],[],[],[]]
    for state in 0:(2^systemSize - 1)
        magnetization = getMagnetization(state, systemSize)
        staggeredMagnetization = getStaggeredMagnetization(state, systemSize)
        if (mod(magnetization, 2) == 0) && (mod(staggeredMagnetization, 2) == 0)
            push!(conservedSubspaces[1], state)
        elseif (mod(magnetization, 2) == 0) && (mod(staggeredMagnetization, 2) == 1)
            push!(conservedSubspaces[2], state)
        elseif (mod(magnetization, 2) == 1) && (mod(staggeredMagnetization, 2) == 0)
            push!(conservedSubspaces[3], state)
        else
            push!(conservedSubspaces[4], state)
        end
    end
    conservedSubspaces
end

function applyHamiltonian(state, systemSize, couplingJ, anisotropy, magnonInteractions)
    magnonRepresentation = digits(state, base = 2, pad = systemSize)
    resultingStates = [state]
    resultingTransitions = [0.0]
    V, T = 0, 0
    for siteI in 1:systemSize
        siteJ = mod1(siteI + 1, systemSize)
        V += couplingJ * (0.25 - 0.5 * (magnonRepresentation[siteI] + magnonRepresentation[siteJ]) + magnonInteractions * (magnonRepresentation[siteI] * magnonRepresentation[siteJ]))
        newState = xor(state, 2^(siteI-1) + 2^(siteJ-1))
        if (magnonRepresentation[siteI] + magnonRepresentation[siteJ]) == 1
            T = 0.25 * couplingJ * (1 + anisotropy)
        else
            T = -0.25 * couplingJ * (1 - anisotropy)
        end
        push!(resultingStates, newState)
        push!(resultingTransitions, T)
    end
    resultingTransitions[1] = V
    (resultingStates, resultingTransitions)
end

function constructSubspaceMatrix(subspace, systemSize, couplingJ, anisotropy, magnonInteractions)
    dimensions = length(subspace)
    result = zeros(Float64, dimensions, dimensions)
    if dimensions > 0
        for is in 1:dimensions
            state = subspace[is]
            adjacentStates, coefficients = applyHamiltonian(state, systemSize, couplingJ, anisotropy, magnonInteractions)
            for js in 1:length(adjacentStates)
                result[is, searchsorted(subspace, adjacentStates[js])[1]] += coefficients[js]
            end
        end
    end
    result
end

function constructBlockHamiltonian(basis, systemSize, couplingJ, anisotropy, magnonInteractions)
    result = []
    for subspace in basis
        matrixBlock = constructSubspaceMatrix(subspace, systemSize, couplingJ, anisotropy, magnonInteractions)
        push!(result, matrixBlock)
    end
    result
end

function diagonalizeBlock(matrixBlock)
    eigen(matrixBlock)
end

function diagonalizeHamiltonian(blockHamiltonian)
    result = []
    for matrixBlock in blockHamiltonian
        if length(matrixBlock) > 0
            push!(result, diagonalizeBlock(matrixBlock))
        else
            push!(result, nothing)
        end
    end
    result
end

function getGroundStateSubspace(factorization) # in case of the Heisenberg model there are 2S + 1 degenerated ground states for ferromagnetic coupling constant (this case should be calculated separately, e.g. assumin sponataneous symmetry breaking)
    groundStateEnergy = Inf
    groundStateSubspaceIndex = 0
    for it in 1:length(factorization)
        subspace = factorization[it]
        if subspace != nothing
            subspaceGroundEnergy = subspace.values[1]
            if groundStateEnergy > subspaceGroundEnergy
                groundStateEnergy = subspaceGroundEnergy
                groundStateSubspaceIndex = it
            end
        end
    end
    (groundStateSubspaceIndex, factorization[groundStateSubspaceIndex])
end