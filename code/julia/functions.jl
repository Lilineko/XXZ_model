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

function getRepresentativeState(state, systemSize)
    magnonRepresentation = digits(state, base = 2, pad = systemSize)
    result = state
    magnonRepresentation = circshift(magnonRepresentation, 1)
    newState = sum([magnonRepresentation[k] * 2^(k-1) for k in 1:systemSize])
    while newState != state
        result = min(result, newState)
        magnonRepresentation = circshift(magnonRepresentation, 1)
        newState = sum([magnonRepresentation[k] * 2^(k-1) for k in 1:systemSize])
    end
    result
end

function getPeriod(state, systemSize)
    magnonRepresentation = digits(state, base = 2, pad = systemSize)
    magnonRepresentation = circshift(magnonRepresentation, 1)
    newState = sum([magnonRepresentation[k] * 2^(k-1) for k in 1:systemSize])
    result = 1
    while newState != state
        magnonRepresentation = circshift(magnonRepresentation, 1)
        newState = sum([magnonRepresentation[k] * 2^(k-1) for k in 1:systemSize])
        result += 1
    end
    result
end

function findContainingMomentumSubspaces(period, systemSize)
    if mod(systemSize, period) != 0
        println("Wrong Periodicity")
        return nothing
    end
    step = div(systemSize, period)
    result = [1]
    index = 1 + step
    while index <= systemSize
        push!(result, index)
        index += step
    end
    result
end

function constructMomentumBasis(systemSize)
    basis = constructBasis(systemSize)
    result = [[[] for _ in basis] for _ in 1:systemSize]
    for is in 1:length(basis)
        subspace = basis[is]
        for state in subspace
            representative = getRepresentativeState(state, systemSize)
            if length(searchsorted(result[1][is], representative)) == 0
                p = getPeriod(representative, systemSize)
                containingMomentumSubspaces = findContainingMomentumSubspaces(p, systemSize)
                for it in containingMomentumSubspaces
                    push!(result[it][is], representative)
                end
            end
        end
    end
    result
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

function getTranslationalFamily(state, systemSize)
    period = getPeriod(state, systemSize)
    magnonRepresentation = digits(state, base = 2, pad = systemSize)
    result = Vector{Int64}(undef, period)
    for it in 1:period
        newMember = circshift(magnonRepresentation, it-1)
        result[it] = sum([newMember[k] * 2^(k-1) for k in 1:systemSize])
    end
    result
end

function getTranslationsToRepresentative(state, systemSize)
    representative = getRepresentativeState(state, systemSize)
    result = 0
    stateRepresentation = digits(state, base = 2, pad = systemSize)
    representativeRepresentation = digits(representative, base = 2, pad = systemSize)
    while stateRepresentation != representativeRepresentation
        representativeRepresentation = circshift(representativeRepresentation, 1)
        result += 1
    end
    result
end

# TODO: write proper algoryth mor phase calculation or better momentum hamiltonian creation
function getPhase(momentum, state, adjacentState, systemSize)
    stateShift = getTranslationsToRepresentative(state, systemSize)
    adjacentStateShift = getTranslationsToRepresentative(adjacentState, systemSize)
    exp(im * (momentum - 1) * (stateShift - adjacentStateShift) / systemSize)
end

function constructMomentumSubspaceMatrix(momentum, subspace, systemSize, couplingJ, anisotropy, magnonInteractions)
    dimensions = length(subspace)
    result = zeros(ComplexF64, dimensions, dimensions)
    if dimensions > 0
        for is in 1:dimensions
            representativeState = subspace[is]
            statesFamily = getTranslationalFamily(representativeState, systemSize)
            for it in 1:length(statesFamily)
                state = statesFamily[it]
                adjacentStates, coefficients = applyHamiltonian(state, systemSize, couplingJ, anisotropy, magnonInteractions)
                for js in 1:length(adjacentStates)
                    adjacentState = adjacentStates[js]
                    indices = searchsorted(subspace, adjacentState)
                    if length(indices) != 0
                        phase = getPhase(momentum, state, adjacentState, systemSize)
                        result[is, indices[1]] += phase * coefficients[js]
                    end
                end
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

function constructMomentumBlockHamiltonian(momentumBasis, systemSize, couplingJ, anisotropy, magnonInteractions)
    result = [[] for _ in 1:systemSize]
    for momentum in 1:systemSize
        for subspace in momentumBasis[momentum]
            matrixBlock = constructMomentumSubspaceMatrix(momentum, subspace, systemSize, couplingJ, anisotropy, magnonInteractions)
            push!(result[momentum], matrixBlock)
        end
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
