using Plots
using Plots.PlotMeasures
using LinearAlgebra
using ExcelFiles

include("functions.jl")

function main(info) #, systemSize)
    # parameters of the model
    systemSize = 14
    magnonInteractions = 1.0
    anisotropy = 0.0 # allowed range is [0, 1); for anisotropy = 1 we have additional degeneration
    couplingJ = -1.0

    # prepare and diagonalize the model
    println("\n", " > Diagonalization - Running...")
    @time systemInfo = getSystemInfo(systemSize, magnonInteractions, anisotropy, couplingJ, info)
    println(" > Diagonalization - Complete!")

    # calculate single spin flip spectral function
    kRange = [n for n in 0:systemSize] .* (2π / systemSize)
    ωRange = [n for n in -0.25:0.005:3.25]
    δ = 0.02
    A = calculateSpectralFunction(kRange, ωRange, δ, systemInfo)

    # return figure
    # return systemInfo
    return (systemSize, kRange, ωRange, A)
end

# @time figure = main("")
# save(getNextFigureName(string(pwd(), "/figures/")), figure)


# # calculation of mean number of magnons
# mean_n = []
# lRange = 2:2:10
# for L in lRange
#     @time gsInfo = main("", L)
#     s = gsInfo[1]
#     gs = gsInfo[3][s].vectors[:,1]
#
#     @time bsInfo = main("basis", L)
#     bs = bsInfo[s]
#
#     push!(mean_n, sum([sum(getMagnonRepresentation(bs[it], L)) * gs[it]^2 for it in 1:length(gs)]) / L)
# end
#
# display(scatter(1.0 ./ [x for x in lRange] , mean_n, xlim = (0, 0.55), ylim = (0, 0.1)))
#
# mean_n


# calculations for linear plots
# @time systemSize, kRange, ωRange, A = main("S+-")
# @time systemSize, kRange, ωRange, B = main("Szz)


# plot the result (linear)
figure = plot()
scale = systemSize / 10.0
kticks = (kRange .* systemSize) / (2pi)
kticklabels = round.(kRange ./ pi, digits = 2)
for it in 1:length(kRange)
    figure = plot!(ωRange, scale * (A[it, :] .- B[it, :]) .+ (it-1),
                    color = :black, legend = :none, ylim = (-0.05 * systemSize, 1.4 * systemSize),
                    w = 2, frame = :box, fill = ((it-1), 0.2, :black),
                    yticks = (kticks, kticklabels), ylabel = "k / π", xlabel = "ω",
                    ann = [
                    (2.3, 2pi * 2.8, text("XY : β = 0", 20)),
                    # (0.0, 2pi * 2.0, text("LSW", 12)),
                    # (0.25, 2pi * 1.8, text("LSW+", 12)),
                    # (0.55, 2pi * 2.1, text("ED", 12))
                    ],
                    xtickfontsize = 14, ytickfontsize = 14,
                    xguidefontsize = 16, yguidefontsize = 16,
                    left_margin=8mm, bottom_margin = 6mm,
                    size = (800, 600)
                    )
end
# qRange = [q for q in 0:(2pi / 200):2pi]
# dispersion = bogoliubov(0.14, qRange, 0.0, -1.0)
# figure = plot!(dispersion, (qRange .* systemSize) / (2pi))
display(figure)

# save(getNextFigureName(string(pwd(), "/figures/")), figure)


# plot the result (heatmap)
figure = heatmap(kRange, ωRange, transpose(A .- B), clim = (0, 1), colorbar = true);
display(figure)
