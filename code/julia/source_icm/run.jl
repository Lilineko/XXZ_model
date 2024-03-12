using OrderedCollections
using JSON
using BenchmarkTools

include("xxz.jl")

@time sys, bs, fc = Main.XXZ.run()
