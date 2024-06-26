using Plots
using Plots.PlotMeasures
using LinearAlgebra
using ExcelFiles
using DelimitedFiles

include("functions.jl")

function main(info) #, systemSize)
    # parameters of the model
    systemSize = 10
    magnonInteractions = 1.0
    anisotropy = 0 # allowed range is [0, 1); for anisotropy = 1 we have additional degeneration
    couplingJ = -1.0

    # prepare and diagonalize the model
    println("\n", " > Diagonalization - Running...")
    @time systemInfo = getSystemInfo(systemSize, magnonInteractions, anisotropy, couplingJ, info)
    println(" > Diagonalization - Complete!")

    # calculate single spin flip spectral function
    kRange = [n for n in 0:systemSize] .* (2π / systemSize)
    ωRange = [n for n in -0.25:0.005:2.75]
    δ = 0.02
    A = calculateSpectralFunction(kRange, ωRange, δ, systemInfo)

    return (systemSize, kRange, ωRange, A)
end

# calculations for linear plots
@time systemSize, kRange, ωRange, A = main("S+-")
@time systemSize, kRange, ωRange, B = main("Szz")

# systemSize, kRange, ωRange, A = readSpectralFunction("./01_impl/results/A_L=16_m=1.0.txt")
# systemSize, kRange, ωRange, B = readSpectralFunction("./01_impl/results/B_L=16_m=1.0.txt")

# plot the result (linear)
figure = plot()
scale = systemSize / 12.0
kticks = (kRange .* systemSize) / (2pi)
kticklabels = round.(kRange ./ pi, digits = 2)
for it in 1:length(kRange)
    figure = plot!(ωRange, scale * (A[it, :] .- B[it, :]) .+ (it-1),
    # figure = plot!(ωRange, scale * (B[it, :]) .+ (it-1),
        color = :black, legend = :none, ylim = (-0.05 * systemSize, 1.4 * systemSize),
        w = 2, frame = :box, fill = ((it-1), 0.2, :black),
        yticks = (kticks, kticklabels), ylabel = "k / π", xlabel = "ω",
        ann = [
        # (2.0, 2pi * 3.2, text("XY, \"S+- - Szz\" : β = 1", 20)),
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

# save(getNextFigureName(string(pwd(), "/")), figure)


# plot the result (heatmap)
# mapfig = heatmap(kRange, ωRange, transpose(A .- B), clim = (0, 1), colorbar = true);
# display(mapfig)

save(getNextFigureName(string(pwd(), "/")), figure)
