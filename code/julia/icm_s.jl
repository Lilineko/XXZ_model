using OrderedCollections
using LinearAlgebra
using DelimitedFiles
using JSON

### Needed due to a bug in libopenblas64_.dll in Julia-1.5.0
if Sys.iswindows()
    LinearAlgebra.BLAS.set_num_threads(1)
end

dir = Dict(
    "xxz" => "./02_impl/"
)

include(dir["xxz"] * "xxz.jl")
include(dir["xxz"] * "spectral.jl")

struct SPC
    size::Int64
    anisotropy::Float64
    interaction::Float64
    parity::Int64
    coupling::Float64
    momentum::Vector{Float64}
    energy::Vector{Float64}
    spectrum::Array{Float64, 2}
end

struct Parameters
    n::Int64
    α::Float64
    β::Float64
    J::Float64
end

function run_s(parameters::Parameters)
    ### System Parameters
    n = parameters.n
    α = parameters.α
    β = parameters.β
    J = parameters.J


    iDelta = 0.05im
    ωRange = collect(-3:0.005:7)

    kRange = (2 / n) .* collect(0:n) ### k range in π
    ### ground state momentum depends on system size: 4k vs 4k + 2
    q = ifelse(mod(n, 4) == 0, 0, div(n, 2))
    ### gound state parity depends on system size: 4k vs 4k + 2
    parity = ifelse(mod(n, 4) == 0, 0, 1)

    ### XXZ Ground State
    input = OrderedDict(
        "system size" => n,
        "parity sector" => parity,
        "momentum sector" => q,
        "coupling constant" => J,
        "magnon interaction" => β,
        "anisotropy" => α
    )

    println()
    @time gsSystem, gsBasis, gsFactorization = Main.XXZ.run(input)

    vals, vecs, info = gsFactorization
    GSE, GSV = vals[1], vecs[1]

    println("\n", info)

    ### Spectral Function of a single spin flip
    spectrum = Array{Float64, 2}(undef, length(ωRange), gsSystem.size + 1)
    for k in 0:gsSystem.size
        println("Evaluating k = ", k, " in 0:", gsSystem.size)
        p = mod(k + q, gsSystem.size)

        input = OrderedDict(
            "system size" => n,
            "parity sector" => mod(parity + 1, 2),
            "momentum sector" => p,
            "coupling constant" => J,
            "magnon interaction" => β,
            "anisotropy" => α
        )

        @time xxzSystem, xxzBasis, xxzModel = Main.XXZ.run(input, factor = false)

        initialState = getInitialState(GSV, gsBasis, gsSystem, xxzBasis, xxzSystem)

        ### lanczos method for spectral function
        @time spectrum[:, k + 1] .= Main.SpectralFunction.run(ωRange .+ GSE, iDelta, initialState, xxzModel)

        println()
    end

    return SPC(n, α, β, parity, J, kRange, ωRange, spectrum)
end

function getInitialState(GSV, gsBasis, gsSystem, xxzBasis, xxzSystem)::Vector{Complex{Float64}}
    l::Int = gsSystem.size

    q = gsSystem.momentum
    p = xxzSystem.momentum

    ### not proved
    k = mod(p - q, gsSystem.size)

    ### spin to flip, 0: down -> up, 1: up -> down
    flipSpin = 0

    function getNewState(state, flipSpin, R)
        newState = state
        if flipSpin == 0
            newState = newState | (1 << R)
        else
            newState = newState & ~(1 << R)
        end
        return newState
    end


    ### calculate initial state
    initialState = zeros(Complex{Float64}, length(xxzBasis))
    for (state, index) in gsBasis
        statePeriodicity = Main.XXZ.getPeriodicity(state, gsSystem)
        for R in 0:(statePeriodicity-1)
            phase = exp(2π * im * k * R / l)
            newState = getNewState(state, flipSpin, R)
            if (newState != state)
                hasMomentum, repState, periodicity, distance = Main.XXZ.getStateInfo(newState, xxzSystem)
                if hasMomentum
                    periodicity = Main.XXZ.getPeriodicity(repState, gsSystem)
                    coefficient = GSV[index] * phase * sqrt(periodicity * statePeriodicity) / l
                    initialState[xxzBasis[repState]] += coefficient
                end
            end
        end
    end
    return initialState
end

function saveData(data::Vector{SPC})
    spcData = Vector{OrderedDict{String, Union{Int64, Float64, Vector{Float64}, Array{Float64, 2}}}}(undef, length(data))

    for it in 1:length(data)
        spcData[it] = OrderedDict(
            "size" => data[it].size,
            "anisotropy" => data[it].anisotropy,
            "interaction" => data[it].interaction,
            "parity" => data[it].parity,
            "coupling" => data[it].coupling,
            "momentum" => data[it].momentum,
            "energy" => data[it].energy,
            "spectrum" => data[it].spectrum
        )
    end

    tail = prod(["_" * replace(arg, ":" => "-") for arg in ARGS])

    file = open(string("../../data/spc", tail, ".json"), "w")
    JSON.print(file, spcData, 1)
    close(file)

    return nothing
end

### System Parameters
J = -0.4
αRange = [1.0]
nRange = [n for n in 16:2:16]
βRange = [Β for Β in 1.0:1.0:1.0]

if length(ARGS) > 0
    J = eval(Meta.parse(ARGS[1]))
end
if length(ARGS) > 1
    nRange = [n for n in eval(Meta.parse(ARGS[2]))]
end
if length(ARGS) > 2
    αRange = [x for x in eval(Meta.parse(ARGS[3]))]
end
if length(ARGS) > 3
    βRange = [x for x in eval(Meta.parse(ARGS[4]))]
end

parameters = Vector{Parameters}()
for n in nRange
    for α in αRange
        for β in βRange
            push!(parameters, Parameters(n, α, β, J))
        end
    end
end

data = Vector{SPC}(undef, length(parameters))
for it in 1:length(parameters)
    data[it] = run_s(parameters[it])
end

saveData(data)
