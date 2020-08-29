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
    V = 0
    for siteI in 1:systemSize
        siteJ = mod1(siteI + 1, systemSize)
        V += couplingJ * (0.25 - 0.5 * (magnonRepresentation[siteI] + magnonRepresentation[siteJ]) + magnonInteractions * (magnonRepresentation[siteI] * magnonRepresentation[siteJ]))
        tempStates = []
        tempTransitions = []
        if (magnonRepresentation[siteI] + magnonRepresentation[siteJ]) == 1
            newState = xor(state, 2^(siteI-1) + 2^(siteJ-1))
            T = 0.25 * couplingJ * (1 + anisotropy)
            push!(tempStates, newState)
            push!(tempTransitions, T)
        else
            newState = xor(state, 2^(siteI-1) + 2^(siteJ-1))
            T = -0.25 * couplingJ * (1 - anisotropy)
            push!(tempStates, newState)
            push!(tempTransitions, T)
        end
        append!(resultingStates, tempStates)
        append!(resultingTransitions, tempTransitions)
    end
    resultingTransitions[1] = V
    result = (resultingStates, resultingTransitions)
end

function constructSubspaceMatrix(subspace, systemSize, couplingJ, anisotropy, magnonInteractions)
    dimensions = length(subspace)
    result = zeros(Float64, dimensions, dimensions)
    if dimensions > 0
        for is in 1:dimensions
            state = subspace[is]
            adjacentStates, coefficients = applyHamiltonian(state, systemSize, couplingJ, anisotropy, magnonInteractions)
            for js in 1:length(adjacentStates)
                result[is, searchsorted(subspace, adjacentStates[js])[1]] = coefficients[js]
            end
        end
    end
    result
end
