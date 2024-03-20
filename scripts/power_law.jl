using Pkg; Pkg.activate(".")
using UnquenchingQCD
using Parameters
using JLD
using Plots; gr()
include("IterationHelper.jl")

pathtodata = "./input"

function estimate_ρ(loadpath)
    GIR, ρ1, ρ2, plt = extrapolate_GIR(loadpath,IR=-4.75)
    GIR, ρ3, ρ4, plt = extrapolate_GIR(loadpath,IR=-4.50)
    GIR, ρ5, ρ6, plt = extrapolate_GIR(loadpath,IR=-4.25)
    GIR, ρ7, ρ8, plt = extrapolate_GIR(loadpath,IR=-4.0)
    ρs = (ρ1,ρ2,ρ3,ρ4,ρ5,ρ6,ρ7,ρ8)

    ρ  = (maximum(ρs) + minimum(ρs))/2
    Δρ = (maximum(ρs) - minimum(ρs))/2
    return ρ, Δρ
end
function estimate_from_two_paths(loadpath1,loadpath2)
    GIR, ρ1, ρ2, plt = extrapolate_GIR(loadpath1,IR=-4.75)
    GIR, ρ3, ρ4, plt = extrapolate_GIR(loadpath1,IR=-4.50)
    GIR, ρ5, ρ6, plt = extrapolate_GIR(loadpath1,IR=-4.25)
    GIR, ρ7, ρ8, plt = extrapolate_GIR(loadpath1,IR=-4.0)
    GIR, ρ9, ρA, plt = extrapolate_GIR(loadpath2,IR=-4.75)
    GIR, ρB, ρC, plt = extrapolate_GIR(loadpath2,IR=-4.50)
    GIR, ρD, ρE, plt = extrapolate_GIR(loadpath2,IR=-4.25)
    GIR, ρF, ρG, plt = extrapolate_GIR(loadpath2,IR=-4.0)
    ρs = (ρ1,ρ2,ρ3,ρ4,ρ5,ρ6,ρ7,ρ8,ρ9,ρA,ρB,ρC,ρD,ρE,ρF,ρG)

    ρ  = (maximum(ρs) + minimum(ρs))/2
    Δρ = (maximum(ρs) - minimum(ρs))/2
    display(plt)
    return ρ, Δρ
end
    
group = "SU2"
Nf = "2.75"

loadpath1 = joinpath(pathtodata,"$group/Nf=$Nf/prop_scaling")
loadpath2 = joinpath(pathtodata,"$group/Nf=$Nf/prop_overshoot")
ρ, Δρ = estimate_from_two_paths(loadpath1,loadpath2)
@show ρ, Δρ

loadpath1 = joinpath(pathtodata,"$group/Nf=$Nf/vertex_scaling")
loadpath2 = joinpath(pathtodata,"$group/Nf=$Nf/vertex_overshoot")
ρ, Δρ = estimate_from_two_paths(loadpath1,loadpath2)
@show ρ, Δρ
