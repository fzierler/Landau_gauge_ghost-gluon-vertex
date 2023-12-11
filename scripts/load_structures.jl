using Pkg; Pkg.activate("."); Pkg.instantiate()
using UnquenchingQCD
using JLD

pathtodata = "/home/fabian/Documents/daten/UnquenchingQCD"
file1 = joinpath(pathtodata,"SO9/Nf=1.1/prop/parameters.jld")
file2 = joinpath(pathtodata,"SO9/Nf=1.1/prop/theory.jld")
file3 = joinpath(pathtodata,"SO9/Nf=6.35/vertex/parametersGGV.jld")

_load_propagator_parameters(file1)
_load_theory(file2)
_load_vertex_parameters(file3)
