using Plots
using LinearAlgebra
using ExcelFiles

function main()
    # parameters of the model
    systemSize = 12
    magnonInteractions = 0.0
    anisotropy = 0.25 # allowed range is [0, 1); for anisotropy = 1 we have additional degeneration
    couplingJ = -1.0

    # prepare and diagonalize the model

    println("\n", " > Diagonalization - Running...")
    @time systemInfo = getSystemInfo(systemSize, magnonInteractions, anisotropy, couplingJ, "")
    println(" > Diagonalization - Complete!")
    # calculate single spin flip spectral function
    kRange = [n for n in 0:systemSize] .* (2π / systemSize)
    ωRange = [n for n in -0.25:0.01:2.75]
    δ = 0.01
    A = calculateSpectralFunction(kRange, ωRange, δ, systemInfo)

    # plot the result
    figure = heatmap(kRange, ωRange, transpose(A), clim = (0, 1), colorbar = true)
    display(figure)
    # return figure
end

@time main()
# @time save(string(pwd(), "/figures/", "fig000.png"), main())
