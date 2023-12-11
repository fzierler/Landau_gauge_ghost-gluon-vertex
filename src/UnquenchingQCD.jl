module UnquenchingQCD

using JLD
using Plots
using ProgressMeter
using Parameters
using KahanSummation
using LoopVectorization
using Quadmath
import FastGaussQuadrature: gausslegendre, gausschebyshev
using Base.Threads

# DSE specifics
include("Structures.jl")
include("IterateVertex.jl")
include("IteratePropagators.jl")
include("IntegrateDiagrams.jl")
include("SaveLoad.jl")
include("Scale.jl")
# integration
include("Integrate.jl")
# interpolation
include("FindIndex.jl")
include("Splines.jl")
include("Multilinear.jl")
include("Chebyshev.jl")

export theory
export loadGGV, interpolateGGV, parameterGGV, iterate_vertex, saveGGV
export parameter, gausslegendre, load_hybrid, quad_kernels, iterate_propagators
export savedata, IterateSystem

export _load_vertex_parameters, _load_theory, _load_propagator_parameters

end # module
