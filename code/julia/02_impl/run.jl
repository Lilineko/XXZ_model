using OrderedCollections
using JSON
using BenchmarkTools

# if !any(names(Main, imported = true) .== :Heisenberg)
    include("xxz.jl")
# end

# n = 4
# J = 1.0
# β = 1.0
#
# μ = 1 + div(n, 2)
#
# E = Inf
# V = []
# k = 0
# for momentum in 0:(n-1)
#     input = OrderedDict(
#         "system size" => n,
#         "momentum sector" => momentum,
#         "magnetization sector" => μ,
#         "coupling constant" => J,
#         "magnon interaction" => β
#     )
#     print(">>> ")
#     @time system, bs = Main.Heisenberg.run(input)
#     vals, vecs, info = fc
#     if info !== missing
#         if info.converged == 0
#             println("Warning! Calculation did not converge.")
#             println()
#         end
#         println("Energy: ", vals[1])
#         if vals[1] < E
#             global E = vals[1]
#             global v = vecs[1]
#             global k = 2 * momentum / n
#         end
#     end
# end
# println("\n", "Ground state found for momentum: ", k, " π")
# println("Ground state energy: ", E)


@time sys, bs, fc = Main.XXZ.run()
