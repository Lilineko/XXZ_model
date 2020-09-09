using Plots
using LinearAlgebra

# parameters of the model
systemSize = 12
magnonInteractions = 0.0
anisotropy = 0.0 # for anisotropy = 1 we have additional degeneration, thus allowed range is [0, 1)
couplingJ = -1.0

# prepare and diagonalize the model
systemInfo = getSystemInfo(systemSize, magnonInteractions, anisotropy, couplingJ)

# calculate single spin flip spectral function
kRange = [n for n in 0:systemSize] .* (2π / systemSize)
ωRange = [n for n in -0.5:0.01:2.5]
δ = 0.05
A = calculateSpectralFunction(kRange, ωRange, δ, systemInfo)

# plot(ωRange, transpose(A))
heatmap(kRange, ωRange, transpose(A), clim = (0, 2), colorbar = true)
