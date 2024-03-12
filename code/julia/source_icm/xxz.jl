module XXZ

using OrderedCollections
using SparseArrays
using KrylovKit
using JSON

"""
`struct System` immutable structure for system input parameters:
# Fields
*   `size::Int64`: number of lattice sites
*   `momentum::Int64`: index indicating which momentum sector should be calculated. Index `momentum` corresponds to momentum `k = 2π * momentum / size`.
*   `magnetization::Int64`: index in range `0:floor(size/2)` indicating which magnetization sector should be calculated. Magnetization is taken without sign.
*   `coupling::Float64`: value of the coupling constant J. For ferromagnet `coupling < 0`, while for antiferromagnet `coupling > 0`.
*   `interaction::FLoat64`: parameter for scaling magnon-magnon interactions. For pure Heisenberg model `interaction = 1.0`.
"""
struct System
    size::Int64
    parity::Int64
    momentum::Int64
    coupling::Float64
    anisotropy::Float64
    interaction::Float64
end

"`Basis === OrderedDict{Int64, Int64}`"
Basis = OrderedDict{Int64, Int64}

"""
`mutable struct LinearCombination`: structure for storing result of operators action on states belonging to `basis::Basis`
# Fields
*   `state::Int64`: magnon configuration in binary representation written as decimal number
*   `coefficient::Vector{Complex{Float64}}`: coeffcient multiplying state in the linear combination
"""
mutable struct LinearCombination
    state::Vector{Int64}
    coefficient::Vector{Complex{Float64}}
end

"`Model === SparseMatrixCSC{Complex{Float64},Int64}`"
Model = SparseMatrixCSC{Complex{Float64},Int64}

"""
    run()

Run model diagonalization procedure.
"""
function run(input::Union{Missing, OrderedDict} = missing; howmany = 1, factor = true)
    system::System = if input === missing
        readInput()
    else
        System(
            input["system size"],
            input["parity sector"],
            input["momentum sector"],
            input["coupling constant"],
            input["magnon interaction"],
            input["anisotropy"]
        )
    end
    basis::Basis = makeBasis(system)
    model::Model = makeModel(basis, system)
    if factor
        @time factorization = factorize(model, howmany = howmany, kryldim = max(30, howmany))
        return system, basis, factorization
    else
        return system, basis, model
    end
end

"Read `input.json` file and return `System` structure with input data. It requires `input.json` file to be located in the current working directory."
function readInput()::System
    path = "./code/julia/02_impl/"
    input = JSON.parsefile(path * "input.json", use_mmap = false) # use_mmap = false is a workaroud to ensure one can change JSON file without restarting Julia
    return System(
        input["system size"],
        input["parity sector"],
        input["momentum sector"],
        input["coupling constant"],
        input["magnon interaction"],
        input["anisotropy"]
    )
end

"""
    makeBasis(system::System) -> Basis

Return `Basis === Dict{Int64, Int64}` dictionary representing basis of given magnetization and momentum sector specified by `system.magnetization` and `system.momentum` respectively. Each index in `Basis` corresponds to state value, and each value in `Basis` corresponds to position in the basis.
"""
function makeBasis(system::System)::Basis
    if isodd(system.size)
        error("Requested 'system.size' is odd! Only even sizes are supported.")
    end
    if system.parity != 0 && system.parity != 1
        error("Wrong magnetization parity sector in the input file! Possible subspaces are 0 or 1.")
    end
    ### comment: there are two possible parity sectors
    ### assigned based on magnetization mod 2 of the state
    ### from this staggered magnetization mod 2 can be infered
    ### for given lattice size (see table in pdf file)
    ### comment 2: different total magnetization values are mixed

    ### potrzebujemy zdiagonalizować obie podprzestrzenie
    ### by znaleźć GS i zrobić spin flip
    ### iteracje po parzystości magnetyzacji mozna zrobić zewnętrznie
    ### zatem bazę tworzymy wybierając podzbiór stanów,
    ### które mają określoną parzystość

    ### dodatkowo chcemy wprowadzić przestrzeń pędową
    ### czyli potrzebujemy tylko reprezentantów z każdej klasy
    ### abstrakcji utworzonej przez relację translacji

    ### ponieważ rozważamy przypadek ferromagnetyczny,
    ### to nie dokonujemy rotacji podsieci by wprowadzić
    ### magnony, zatem translacja w przestrzeni magnonów
    ### jest tożsama z translacją w przestrzeni spinów

    ### pytanie: czy istenieje możlwiość wyboru reprezentantów
    ### szybsza niż iteracja po całej dostępnej przestrzeni stanów?

    ### initialize the basis
    basis::Basis = Basis()

    index = 0
    ### iterate over all states
    for state in 0:(2^system.size - 1)

        ### check if state belongs to proper parity sector
        if hasParity(state, system)

            ### check if state belongs to requested momentum subspace
            if hasMomentum(state, system)
                ### pick representative state
                repState = getRepresentative(state, system)

                ### check if repState is already included
                if get(basis, repState, nothing) === nothing
                    ### if not then add repState to basis
                    push!(basis, repState => (index += 1))
                end
            end
        end
    end

    return basis
end

# """
#     sublatticeRotation(state::Int64, mask::Int64) -> Int64
#
# Reverse bits according to mask.
# """
# function sublatticeRotation(state::Int64, mask::Int64)::Int64
#     return xor(state, mask)
# end

"""
    hasParity(state::Int64, system::System) -> Bool

Return `True` if `state` belongs to `system.parity` subspace or returns `False` otherwise.
"""
function hasParity(state::Int64, system::System)::Bool
    return mod(sum(digits(state, base = 2, pad = system.size)), 2) == system.parity
end

"""
    hasMomentum(state::Int64, system::System) -> Bool

Return `True` if `state` belongs to `system.momentum` subspace or returns `False` otherwise.
"""
function hasMomentum(state::Int64, system::System)::Bool
    ### system.size must divide system.momentum times periodicity
    return rem(system.momentum * getPeriodicity(state, system), system.size) == 0
end

"""
    getPeriodicity(state::Int64, system::System) -> Int64

Return periodicity of the `state` within given `system` parameters.
"""
function getPeriodicity(state::Int64, system::System)::Int64
    ### initialize some constants for faster evaluation
    l::Int = system.size
    highestBit::Int = 1 << (l - 1)
    highestValue::Int = (1 << l) - 1

    ### initialize state translation
    stateTranslation = state

    ### periodicity is smallest positive number of translations
    ### that transform state onto itself
    periodicity::Int64 = 1
    while state != (stateTranslation = bitmov(stateTranslation, l, false, hb = highestBit, hv = highestValue))
        periodicity += 1
    end

    return periodicity
end

"""
    getRepresentative(state::Int64, system::System) -> Int64

Return representative state `repState::Int64` of `state::Int64` within given `system::System` parameters.
`repState::Int64` is a cyclic translation of `state::Int64` that has smallest possible value.
"""
function getRepresentative(state::Int64, system::System)::Int64
    ### initialize some constants for faster evaluation
    l::Int = system.size
    highestBit::Int = 1 << (l - 1)
    highestValue::Int = (1 << l) - 1

    ### initialize result as state
    result = state

    ### loop over translations of state
    newState = state
    while state != (newState = bitmov(newState, l, false, hb = highestBit, hv = highestValue))
        ### replace returned state if translation has lower value
        if result > newState
            result = newState
        end
    end

    ### returned value is so called representative state
    ### i.e. it is the same state as initial one (with respect
    ### to cyclic translation) but it has smallest possible
    ### value within the familly of translations of initial state
    return result
end

"""
    getStateInfo(state::Int64, system::System) -> (Int64, Int64)

Combine `hasMomentum` and `getRepresentative` and `getPeriodicity` in one more efficient function. In addition returns distance between `state` and its represantative state calculated in number of translations needed to transform one onto another.
"""
function getStateInfo(state::Int64, system::System)::Tuple{Bool, Int64, Int64, Int64}
    ### initialize some constants for faster evaluation
    l::Int = system.size
    highestBit::Int = 1 << (l - 1)
    highestValue::Int = (1 << l) - 1

    ### initialize representative as state
    representative = state

    ### loop over translations of state
    newState = state
    distance::Int64 = 0
    periodicity::Int64 = 1
    while state != (newState = bitmov(newState, l, false, hb = highestBit, hv = highestValue))
        if representative > newState
            representative = newState
            distance = periodicity
        end
        periodicity += 1
    end

    hasMomentum = rem(system.momentum * periodicity, system.size) == 0
    return (hasMomentum, representative, periodicity, distance)
end

"""
    bitmov(s::Int, l::Int, f::Bool = false; hb::Int = 1 << (l - 1), hv::Int = (1 << l) - 1) -> Int

Cyclic bit shift for calculationg bit translations with periodic boundary conditions.

# Arguments
*   `s::Int` - value which binary representation will be shifted
*   `l::Int` - size of the cycle (total number of bits in the cycle)
*   `f::Bool` - `true ->` move forward (~mult by `2`); `false ->` backward (~div by `2`); ~mult and ~div mean multiplication and division with preservance of periodic boundary conditions within cycle size `l`
*   `hb::Int` [optional] - highest bit (for speed up put the value of 2^(l-1))
*   `hv::Int` [optional] - highest value (for speed put (2^l)-1)
"""
@inline bitmov(s::Int, l::Int, f::Bool = false; hb::Int = 1 << (l - 1), hv::Int = (1 << l) - 1) = f ? 2s - div(s, hb) * hv : div(s, 2) + rem(s, 2) * hb

"""
    act(operator::Function, state::Int64, basis::Basis, system::System) -> LinearCombination

Apply `operator` to `state` belonging to `basis` and returns `LinearCombination  === Dict{Int64, Complex{Float64}}` representing states with their coefficients.
"""
function act(operator::Function, state::Int64, basis::Basis, system::System)::LinearCombination
    return operator(state, basis, system)
end

"""
    hamiltonian(state::Int64, basis::Basis, system::System) -> LinearCombination

Apply Hamiltonian to `state` written in hole-magnon momentum `basis` obtained for input `system` parameters. Returns `LinearCombination` structure representing resulting states with their coefficients.
"""
function hamiltonian(state::Int64, basis::Basis, system::System)::LinearCombination
    ### initialize result as empty linear combination
    result = LinearCombination(fill(state, system.size + 1), zeros(Complex{Float64}, system.size + 1))

    ### check if initial state belongs to basis
    if haskey(basis, state)
        ### initialize ik for faster exponent calculations
        ik::Complex{Float64} = 2.0 * pi * im * system.momentum / system.size

        ### calculate state periodicity
        periodicity = getPeriodicity(state, system)

        ### loop over lattice sites
        for i in 1:system.size
            j = mod1(i + 1, system.size)

            ### get bit value at i and j bit positions
            iValue, jValue = (1 << (i - 1)), (1 << (j - 1))

            ## work out off-diagonal coeffcients
            iBit, jBit = div(state & iValue, iValue), div(state & jValue, jValue)

            ### magnon single and pair hopping
            newState = xor(state, iValue + jValue)
            ### remember we applied the sublattice rotation creating the basis
            ### this is why we flip spins when they are the same
            ### it corresponds to creating pair of magnons in two nearest
            ### empty sites or annihilating pair of neighboring magnons

            ### get info about state after spin flip
            hasMomentum, repState, repPeriodicity, distance = getStateInfo(newState, system)

            ### check if it belongs to correct momentum subspace
            ### if it does not, then it will cancel out with other terms
            ### after summing over all the sites
            if hasMomentum
                ### take proper tau (for pair or single hopping)
                tau = ifelse(iBit == jBit, system.anisotropy - 1.0, system.anisotropy + 1.0)

                ### calculate matrix coefficient
                coefficient = 0.25 * tau * exp(ik * distance) * sqrt(periodicity / repPeriodicity)

                ### create a new entry in linear combination
                ### and set its corresponding coeffcient
                result.state[i + 1] = repState
                result.coefficient[i + 1] = coefficient
            end

            ## work out diagonal coefficient
            ## comment: bits represent magnons (0 -> no magnon, 1 -> magnon present)
            # iBit, jBit = div(state & iValue, iValue), div(state & jValue, jValue)
            result.coefficient[1] += 0.25 - 0.5 * (iBit + jBit) + system.interaction * iBit * jBit
        end

        ### multiply the result by coupling constant
        result.coefficient .*= system.coupling
    end

    ### return resulting linear combination
    return result
end

"""
    makeModel(basis::Basis, system::System) -> Model

Calculate sparse matrix of the `Model` Hamiltonian. Returns `Model === Array{Complex{Float64},2}`.
"""
function makeModel(basis::Basis, system::System)::Model
    subspaceSize = length(basis)
    linearCombinationLength = system.size + 1
    I = Vector{Int64}(undef, linearCombinationLength * subspaceSize)
    J = Vector{Int64}(undef, linearCombinationLength * subspaceSize)
    V = Vector{Complex{Float64}}(undef, linearCombinationLength * subspaceSize)
    for (state, index) in basis
        linearCombination::LinearCombination = act(hamiltonian, state, basis, system)
        for it in 1:linearCombinationLength
            I[(index - 1) * linearCombinationLength + it] = index
            J[(index - 1) * linearCombinationLength + it] = basis[linearCombination.state[it]]
            V[(index - 1) * linearCombinationLength + it] = linearCombination.coefficient[it]
        end
    end
    return dropzeros!(sparse(I, J, V, subspaceSize, subspaceSize, +))
end

"""
    factorize(model::Model, [howmany = 1, which = :SR])

Compute eigenvalues (by default with smallest real part) and their corresponding eigenvectors.
"""
function factorize(model::Model; howmany = 1, which = :SR, kryldim)
    if length(model) != 0
        return eigsolve(model, howmany, which, ishermitian = true, krylovdim = max(30, howmany))
    else
        return (missing, missing, missing)
    end
end

# """
#     saveResult(factorization)
#
# Take output of the `factorize(model::Model)` and write file with convergance info, norm of ritz resudual for smallest eigenvalue, smallest eigenvalue and its corresponding eigenvector.
# """
# function saveResult(factorization)
#     vals, vecs, info = factorization
#     file = open("result.txt", "w")
#     write(file, string("Converged:", "\n", info.converged, "\n\n"))
#     write(file, string("Norm of Residual:", "\n", info.normres[1], "\n\n"))
#     write(file, string("Eigenvalue:", "\n", vals[1], "\n\n"))
#     write(file, string("Eigenvector:", "\n"))
#     for coeff in vecs[1]
#         write(file, string(coeff, "\n"))
#     end
#     close(file)
# end
#
end
