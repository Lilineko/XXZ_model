function getMagnetization(state, systemSize)
    sum(digits(state, base = 2)) - systemSize/2
end

function getMagnonRepresentation(state, systemSize)
    digits(state, base = 2, pad = systemSize)
end

function getStaggeredMagnetization(state, systemSize)
    spinRepresentation = getMagnonRepresentation(state, systemSize) .- 0.5
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

function getStateIndex(representation, systemSize)
    sum([representation[k] * 2^(k-1) for k in 1:systemSize])
end

function getRepresentativeState(state, systemSize)
    magnonRepresentation = getMagnonRepresentation(state, systemSize)
    result = state
    magnonRepresentation = circshift(magnonRepresentation, 1)
    newState = sum([magnonRepresentation[k] * 2^(k-1) for k in 1:systemSize])
    while newState != state
        result = min(result, newState)
        magnonRepresentation = circshift(magnonRepresentation, 1)
        newState = getStateIndex(magnonRepresentation, systemSize)
    end
    result
end

# function getPeriod(state, systemSize)
#     magnonRepresentation = getMagnonRepresentation(state, systemSize)
#     magnonRepresentation = circshift(magnonRepresentation, 1)
#     newState = sum([magnonRepresentation[k] * 2^(k-1) for k in 1:systemSize])
#     result = 1
#     while newState != state
#         magnonRepresentation = circshift(magnonRepresentation, 1)
#         newState = getStateIndex(magnonRepresentation, systemSize)
#         result += 1
#     end
#     result
# end

# function findContainingMomentumSubspaces(period, systemSize)
#     if mod(systemSize, period) != 0
#         println("Wrong Periodicity")
#         return nothing
#     end
#     step = div(systemSize, period)
#     result = [1]
#     index = 1 + step
#     while index <= systemSize
#         push!(result, index)
#         index += step
#     end
#     result
# end

# function constructMomentumBasis(systemSize)
#     basis = constructBasis(systemSize)
#     result = [[[] for _ in basis] for _ in 1:systemSize]
#     for is in 1:length(basis)
#         subspace = basis[is]
#         for state in subspace
#             representative = getRepresentativeState(state, systemSize)
#             if length(searchsorted(result[1][is], representative)) == 0
#                 p = getPeriod(representative, systemSize)
#                 containingMomentumSubspaces = findContainingMomentumSubspaces(p, systemSize)
#                 for it in containingMomentumSubspaces
#                     push!(result[it][is], representative)
#                 end
#             end
#         end
#     end
#     result
# end

function applyHamiltonian(state, systemSize, couplingJ, anisotropy, magnonInteractions)
    magnonRepresentation = getMagnonRepresentation(state, systemSize)
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

# function getTranslationalFamily(state, systemSize)
#     magnonRepresentation = getMagnonRepresentation(state, systemSize)
#     result = [state]
#     newRepresentation = circshift(magnonRepresentation, 1)
#     while newRepresentation != magnonRepresentation
#         push!(result, getStateIndex(newRepresentation, systemSize))
#         newRepresentation = circshift(newRepresentation, 1)
#     end
#     result
# end

# function getTranslationsToRepresentative(state, representative, systemSize)
#     result = 0
#     stateRepresentation = getMagnonRepresentation(state, systemSize)
#     representativeRepresentation = digits(representative, base = 2, pad = systemSize)
#     while stateRepresentation != representativeRepresentation
#         representativeRepresentation = circshift(representativeRepresentation, 1)
#         result += 1
#     end
#     result
# end

# function getPhase(momentum, adjacentState, adjacentRepresentative, systemSize)
#     adjacentShift = getTranslationsToRepresentative(adjacentState, systemSize)
#     exp(2π * (momentum - 1) * im * adjacentShift / systemSize)
# end

# function constructMomentumSubspaceMatrix(momentum, subspace, systemSize, couplingJ, anisotropy, magnonInteractions)
#     dimensions = length(subspace)
#     result = zeros(ComplexF64, dimensions, dimensions)
#     if dimensions > 0
#         for is in 1:dimensions
#             representativeState = subspace[is]
#             statePeriod = getPeriod(representativeState, systemSize)
#             magnonRepresentation = getMagnonRepresentation(representativeState, systemSize)
#             for r in 0:(systemSize-1)
#                 state = getStateIndex(magnonRepresentation, systemSize)
#                 adjacentStates, coefficients = applyHamiltonian(state, systemSize, couplingJ, anisotropy, magnonInteractions)
#                 for js in 1:length(adjacentStates)
#                     adjacentState = adjacentStates[js]
#                     adjacentRepresentative = getRepresentativeState(adjacentState, systemSize)
#                     indices = searchsorted(subspace, adjacentRepresentative)
#                     if length(indices) != 0
#                         adjacentPeriod = getPeriod(adjacentState, systemSize)
#                         phase = getPhase(momentum, adjacentState, adjacentRepresentative, systemSize)
#                         result[is, indices[1]] += sqrt(adjacentPeriod / statePeriod) * phase * coefficients[js]
#                     end
#                 end
#                 magnonRepresentation = circshift(magnonRepresentation, 1)
#             end
#         end
#     end
#     result
# end

function constructBlockHamiltonian(basis, systemSize, couplingJ, anisotropy, magnonInteractions)
    result = []
    for subspace in basis
        matrixBlock = constructSubspaceMatrix(subspace, systemSize, couplingJ, anisotropy, magnonInteractions)
        push!(result, matrixBlock)
    end
    result
end

# function constructMomentumBlockHamiltonian(momentumBasis, systemSize, couplingJ, anisotropy, magnonInteractions)
#     result = [[] for _ in 1:systemSize]
#     for momentum in 1:systemSize
#         for subspace in momentumBasis[momentum]
#             matrixBlock = constructMomentumSubspaceMatrix(momentum, subspace, systemSize, couplingJ, anisotropy, magnonInteractions)
#             push!(result[momentum], matrixBlock)
#         end
#     end
#     result
# end

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
        if subspace !== nothing
            subspaceGroundEnergy = subspace.values[1]
            if groundStateEnergy > subspaceGroundEnergy
                groundStateEnergy = subspaceGroundEnergy
                groundStateSubspaceIndex = it
            end
        end
    end
    (groundStateSubspaceIndex, factorization[groundStateSubspaceIndex])
end

function applySpinFlipUp(position, state, subspace, basis, systemSize)
    dimensions = length(state)
    newSubspace = 5 - subspace
    result = zeros(Float64, dimensions) # we can do this because both subspaces are of the same size
    for it in 1:dimensions
        coefficient = state[it]
        magnonRepresentation = getMagnonRepresentation(basis[subspace][it], systemSize)
        if magnonRepresentation[position] == 0
            magnonRepresentation[position] = 1
            newState = getStateIndex(magnonRepresentation, systemSize)
            index = searchsorted(basis[newSubspace], newState)[1]
            result[index] += coefficient
        end
    end
    (newSubspace, result)
end

function calculateOverlaps(groundSubspaceIndex, groundStateVector, factorization, basis, systemSize)
    flipSubspaceIndex = 5 - groundSubspaceIndex
    numberOfEigenstates = length(factorization[flipSubspaceIndex].values)
    result = Array{Float64, 2}(undef, numberOfEigenstates, systemSize) # we know vectors have real coefficient thus overlaps must be real too
    for position in 1:systemSize
        flippedSubspaceIndex, flippedState = applySpinFlipUp(position, groundStateVector, groundSubspaceIndex, basis, systemSize)
        for index in 1:numberOfEigenstates
            result[index, position] = sum(factorization[flippedSubspaceIndex].vectors[:, index] .* flippedState)
        end
    end
    result
end

function calculateGreensFunctionValue(k, ω, groundSubspaceIndex, groundStateEnergy, factorization, overlaps)
    numberOfEigenstates, systemSize = size(overlaps)
    flipSubspaceIndex = 5 - groundSubspaceIndex
    eigenValues = factorization[flipSubspaceIndex].values
    denominators = Vector{ComplexF64}(undef, numberOfEigenstates)
    for index in 1:numberOfEigenstates
        denominators[index] = ω - eigenValues[index] + groundStateEnergy
    end
    result = 0
    for it in 1:systemSize
        ri = it - 1
        rightOverlaps = overlaps[:, it]
        for jt in 1:systemSize
            rj = jt - 1
            leftOverlaps = overlaps[:, jt]
            result += exp(-im * k * (ri - rj)) * sum(rightOverlaps .* leftOverlaps ./ denominators)
        end
    end
    result
end

function calculateSpectralFunction(kRange, ωRange, δ, systemInfo)
    groundSubspaceIndex, groundStateEnergy, factorization, overlaps = systemInfo
    dimensions = (length(kRange), length(ωRange))
    result = Array{Float64, 2}(undef, dimensions[1], dimensions[2])
    for it in 1:dimensions[1]
        for jt in 1:dimensions[2]
            result[it, jt] = -imag(calculateGreensFunctionValue(kRange[it], ωRange[jt] + δ*im, groundSubspaceIndex, groundStateEnergy, factorization, overlaps)) / π
        end
    end
    result
end
