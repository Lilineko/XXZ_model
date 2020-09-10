using Plots
using LinearAlgebra
using ExcelFiles


function main()
    # parameters of the model
    systemSize = 12
    magnonInteractions = 0.0
    anisotropy = 0.0 # for anisotropy = 1 we have additional degeneration, thus allowed range is [0, 1)
    couplingJ = -1.0

    # prepare and diagonalize the model
    systemInfo = getSystemInfo(systemSize, magnonInteractions, anisotropy, couplingJ, "")

    # calculate single spin flip spectral function
    kRange = [n for n in 0:systemSize] .* (2π / systemSize)
    ωRange = [n for n in -0.25:0.01:2.75]
    δ = 0.02
    A = calculateSpectralFunction(kRange, ωRange, δ, systemInfo)

    # plot the result
    figure = heatmap(kRange, ωRange, transpose(A), clim = (0, 1), colorbar = true)
    display(figure)
    return figure
end

@time figure = main()

# save(string(pwd(), "/figures/", "fig.png"), figure)
