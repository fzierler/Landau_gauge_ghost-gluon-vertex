using Pkg; Pkg.activate(".")
using Plots
using JLD
using UnquenchingQCD
using Parameters
using LaTeXStrings
using Measures 
pgfplotsx(tickfontsize=14,labelfontsize=16,legendfontsize=15,titlefontsize=20,framestyle=:box,lw=3,size=(600,400),color_palette=palette(:Dark2_3))

function correct_vertex_path(pathtodata,gauge_group,Nfdir)
    ggvpath = joinpath(pathtodata,gauge_group,Nfdir,"vertex","GGV.jld")
    if !isfile(ggvpath)
        ggvpath = joinpath(pathtodata,gauge_group,Nfdir,"vertex_scaling","GGV.jld")
    end
    return ggvpath
end
function correct_prop_path(pathtodata,gauge_group,Nfdir)
    ggvpath = joinpath(pathtodata,gauge_group,Nfdir,"prop","hybrid.jld")
    if !isfile(ggvpath) 
        #ggvpath = joinpath(pathtodata,gauge_group,Nfdir,"prop_scaling","hybrid.jld")
    end
    return ggvpath
end
function _load_correct_vertex(pathtodata,gauge_group,Nfdir)
    ggvpath = correct_vertex_path(pathtodata,gauge_group,Nfdir)
    x,y,cosθ,GGV = loadGGV(ggvpath)
    return x,y,cosθ,GGV
end
function _vertex_special_1d(GGV,cosθ)
    θindex = findmin(abs.(cosθ.-1/2))[2]
    N = first(size(GGV))
    GGV_van_x = zeros(N)
    GGV_van_y = zeros(N)
    GGV_equalmomenta = zeros(N)
    for i in 1:N
        GGV_van_x[i] = GGV[1,i,end]
        GGV_van_y[i] = GGV[i,1,end]
        GGV_equalmomenta[i] = GGV[i,i,θindex]
    end
    return GGV_van_x, GGV_van_y, GGV_equalmomenta
end
function vanishing_vertex(path, gauge_group)
    Nf = Float64[]
    A0 = Float64[]
    for Nfdir in readdir(joinpath(path,gauge_group))
        append!(Nf, parse(Float64,last(split(Nfdir,"="))))
        A000 = NaN
        try 
            x,y,cosθ,GGV = _load_correct_vertex(path,gauge_group,Nfdir)
            GGV_van_x, GGV_van_y, GGV_equalmomenta = _vertex_special_1d(GGV,cosθ)
            A000 = GGV_equalmomenta[1]
        catch
        end
        append!(A0, A000)
    end
    return (Nf,A0)
end
function vanishing_momentum_mass_function(path, gauge_group; with_vertex = true)
    Nf = Float64[]
    M0 = Float64[]
    for Nfdir in readdir(joinpath(path,gauge_group))
        append!(Nf, parse(Float64,last(split(Nfdir,"="))))
        if with_vertex
            ggvpath = correct_vertex_path(path,gauge_group,Nfdir)
        else
            ggvpath = correct_prop_path(path,gauge_group,Nfdir)
        end
        if isfile(ggvpath)
            data = load(joinpath(dirname(ggvpath),"hybrid.jld")) 
            append!(M0, first(data["B"] ./ data["A"]))
        else
            @warn "$Nfdir: no vertex results"
            append!(M0, NaN)
        end
    end
    return (Nf,M0)
end

path = "input"
gauge_group = "SU3"


plt = plot(title = "Correlation functions at vanishing momenta: $gauge_group",xlabel=L"N_f")

Nf, M0 = vanishing_momentum_mass_function(path, gauge_group; with_vertex = false)
plot!(plt,Nf[1:8]  ,M0[1:8]  ,markershape=:rect,linestyle=:dash,color=:black,ms=5,label=L"$M(p^2=0)$ (propagators only)")
plot!(plt,Nf[9:end],M0[9:end],markershape=:rect,linestyle=:dash,color=:black,ms=5,label="")

Nf, M0 = vanishing_momentum_mass_function(path, gauge_group)
plot!(plt,Nf[1:8]  ,M0[1:8]  ,markershape=:circle,color=:green,ms=5,label=L"$M(p^2=0)$ (dynamic vertex)")
plot!(plt,Nf[9:end],M0[9:end],markershape=:circle,color=:green,ms=5,label="")

plot!(plt,ylims=(-0.05,0.5),legend=:topright)
savefig("plots/M0_vs_Nf_decoupling.pdf")