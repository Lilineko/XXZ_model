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

function mapToEigenBasis(state, eigenBasis)
    dimensions = length(state)
    result = zeros(Float64, dimensions)
    for it in 1:dimensions
        result[it] = sum(state .* eigenBasis[:, it])
    end
    result
end

function calculateOverlaps(groundSubspaceIndex, groundStateVector, factorization, basis, systemSize)
    flippedSubspaceIndex = 5 - groundSubspaceIndex
    numberOfEigenstates = length(factorization[flippedSubspaceIndex].values)
    result = Array{Float64, 2}(undef, numberOfEigenstates, systemSize) # we know vectors have real coefficient thus overlaps must be real too
    for position in 1:systemSize
        _, flippedState = applySpinFlipUp(position, groundStateVector, groundSubspaceIndex, basis, systemSize)
        result[:, position] .= mapToEigenBasis(flippedState, factorization[flippedSubspaceIndex].vectors)
    end
    result
end

function getSystemInfo(systemSize = 2, magnonInteractions = 1,  anisotropy = 0.0, couplingJ = -1.0, info = "")
    # construct the basis (with reorganization into proper subspaces)
    basis = constructBasis(systemSize)
    if info == "basis"
        return basis
    end

    # calculate the matrix of the Hamiltonian
    blockHamiltonian = constructBlockHamiltonian(basis, systemSize, couplingJ, anisotropy, magnonInteractions)
    if info == "matrix"
        return blockHamiltonian
    end

    # diagonalize all the subspaces
    factorization = diagonalizeHamiltonian(blockHamiltonian)
    if info == "eigen"
        return factorization
    end

    # get ground state info
    groundSubspaceIndex, groundStateSubspace = getGroundStateSubspace(factorization)
    groundStateEnergy = groundStateSubspace.values[1] # note: energies are sorted in ascending order
    groundStateVector = groundStateSubspace.vectors[:, 1]
    if info == "ground"
        return (groundSubspaceIndex, groundStateSubspace)
    end

    # apply spin flip to the ground state and calculate overlaps with eigenstates of the model
    overlaps = calculateOverlaps(groundSubspaceIndex, groundStateVector, factorization, basis, systemSize)
    if info == "overlap"
        return overlaps
    end

    # merge system info
    (groundSubspaceIndex, groundStateEnergy, factorization, overlaps)
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
            numerators = leftOverlaps .* rightOverlaps
            result += exp(-im * k * (ri - rj)) * sum(numerators ./ denominators)
        end
    end
    result / systemSize
end

# function calculateSpectralFunction(kRange, ωRange, δ, systemInfo)
#     groundSubspaceIndex, groundStateEnergy, factorization, overlaps = systemInfo
#     dimensions = (length(kRange), length(ωRange))
#     result = Array{Float64, 2}(undef, dimensions[1], dimensions[2])
#     println()
#     println(" > Calculating Spectral Function : ")
#     progress = 0
#     nSteps = ceil(Int64, dimensions[2] / Threads.nthreads())
#     Threads.@threads for jt in 1:dimensions[2]
#         if Threads.threadid() == 1
#             progress += 1
#             print(" --> Evaluating Step : ", progress, "/", nSteps)
#         end
#         @simd for it in 1:dimensions[1]
#             result[it, jt] = (-1.0 / π) * imag(calculateGreensFunctionValue(kRange[it], ωRange[jt] + δ*im, groundSubspaceIndex, groundStateEnergy, factorization, overlaps))
#         end
#         if Threads.threadid() == 1
#             print("\r")
#         end
#     end
#     println(" --> Evaluating Step : ", nSteps, "/", nSteps)
#     println(" --> Evaluation Complete!")
#     result
# end

function calculateSpectralFunction(kRange, ωRange, δ, systemInfo)
    groundSubspaceIndex, groundStateEnergy, factorization, overlaps = systemInfo
    dimensions = (length(kRange), length(ωRange))
    result = Array{Float64, 2}(undef, dimensions[1], dimensions[2])
    println()
    println(" > Calculating Spectral Function : ")
    progress = 0
    nSteps = dimensions[2]
    for jt in 1:nSteps
        progress += 1
        print(" --> Evaluating Step : ", progress, "/", nSteps)
        for it in 1:dimensions[1]
            result[it, jt] = (-1.0 / π) * imag(calculateGreensFunctionValue(kRange[it], ωRange[jt] + δ*im, groundSubspaceIndex, groundStateEnergy, factorization, overlaps))
        end
        print("\r")
    end
    println(" --> Evaluating Step : ", nSteps, "/", nSteps)
    println(" --> Evaluation Complete!")
    result
end

function getNextFigureName(path)
    fileNames = readdir(path)
    index = 0;
    for name in fileNames
        if name[1:3] == "fig"
            index = max(index, parse(Int, name[4:6]))
        end
    end
    string(path, "fig", lpad(index + 1, 3, "0"), ".png")
end
