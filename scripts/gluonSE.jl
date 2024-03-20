using Pkg; Pkg.activate(".")
using UnquenchingQCD
using Parameters
using JLD
using Plots
using LaTeXStrings

pathtodata = "./input"

size = (400,300)
pgfplotsx(tickfontsize=11,labelfontsize=15,legendfontsize=11,titlefontsize=15,framestyle=:box,lw=2,size=size)
pgfplotsx(color_palette=palette(:Dark2_3))

function gluon_contribution(dirProp)
    pathYM  = joinpath(dirProp,"hybrid.jld")
    pathdat  = joinpath(dirProp,"data.jld")
    dQCD = load(joinpath(dirProp,"theory.jld"))
    dPRM = load(joinpath(dirProp,"parameters.jld"))
    #########################################################
    QCD   = _load_theory(joinpath(dirProp,"theory.jld"))
    param = _load_propagator_parameters(joinpath(dirProp,"parameters.jld"))
    
    @unpack ϵ, Λ, nCheby,μ = param
    data = load(pathdat)
    Z2 = data["Z2"]
    Zm = data["Zm"]
    Z3 = data["Z3"]
    Z3_tilde = data["Z3_tilde"]
    RCs = UnquenchingQCD.RenormalizationConstants(Z2,Zm,Z3,Z3_tilde)
    x  = gausslegendre(1000,ϵ/5,5Λ)[1]
    A,B,p2,coeffG,coeffZ = load_hybrid(pathYM,x,nCheby,ϵ,Λ,µ)
    Quad,Kernels = quad_kernels(p2,param,QCD)
    interYM = UnquenchingQCD.interpolate_ym(p2,coeffG,coeffZ,QCD,Quad,param)
    interQ = UnquenchingQCD.interpolate_quark(x,A,B,Quad,QCD)
    @unpack k2, q2, k2_inv, q2_inv = Quad
    UnquenchingQCD.interpolate_ym!(k2,k2_inv,interYM.Gy,interYM.Zy,p2,coeffG,coeffZ,QCD,param)
    UnquenchingQCD.interpolate_ym!(q2,q2_inv,interYM.Gz,interYM.Zz,p2,coeffG,coeffZ,QCD,param)
    ghost, gluon, quark =  UnquenchingQCD.integrate(p2,coeffG,coeffZ,QCD,Kernels,Quad,param,interYM,interQ,RCs)[2:end]
    return p2, ghost, gluon, quark, Z2, dQCD["Nf"], Z3
end
defaultcolor(i) = get_color_palette(:auto, plot_color(:white))[i]
plotSE(dir;kws...) = plotSE!(plot(),dir;kws...)
function plotSE!(pltR,dir;ylims,label=true)
    p2, ghost, gluon, quark, Z2, Nf, Z3 = gluon_contribution(dir)

    Nfquark = Nf.*quark
    plot!(pltR,p2,-gluon./ghost,label=L"$+\tiny{\tilde{\Pi}_{\rm gluon}}(p^2)$",color=defaultcolor(1),ls=:dash)
    plot!(pltR,p2,ghost./ghost,label=L"$-\tiny{\tilde{\Pi}_{\rm ghost}}(p^2) \equiv 1$",color=:black,ls=:black,lw=1)
    plot!(pltR,p2,-Nfquark./ghost,label=L" $+\tiny{\tilde{\Pi}_{\rm quark}}(p^2)$",color=defaultcolor(2),ls=:dot)
    plot!(pltR,p2,-(gluon .+ Nfquark )./ghost,label=L"$\tiny{\tilde{\Pi}_{\rm gluon}}(p^2) + \tiny{\tilde{\Pi}_{\rm quark}}(p^2)$",color=defaultcolor(3))

    plot!(pltR,xscale=:log10,legend_columns=2,legendfonthalign=:left)
    plot!(pltR,yticks=min(ylims...):0.2:max(ylims...))

    ylims!(pltR,ylims)
    xlims!(pltR,(10^-5,10^4))

    # annotate number of flavours 
    Nfstring = isinteger(round(Nf,digits=1)) ? string(Int(Nf)) : string(round(Nf,digits=1))
    annotate!(pltR,10^-3.5,ylims[2]*0.84,L"N_f = %$Nfstring")
    if label 
        log_tickval = collect(-5:1:5)
        plot!(pltR,xticks=(exp10.(log_tickval), (x->LaTeXString("\$ 10 ^ {$x} \$")).(log_tickval)))
    else
        log_tickval = collect(-5:1:5)
        plot!(pltR,xticks=(exp10.(log_tickval), (x->"").(log_tickval)))
    end
    return pltR
end

dir2 = joinpath(pathtodata,"SU3","Nf=0.1","prop")
dir3 = joinpath(pathtodata,"SU3","Nf=0.8","prop")
dir4 = joinpath(pathtodata,"SU3","Nf=1.0","prop")

plt2 = plotSE(dir2;ylims=(-0.4,1.4))
plt3 = plotSE(dir3;ylims=(-0.4,1.4))
plt4 = plotSE(dir4;ylims=(-0.4,1.4))

plot!(plt2,legend=:left,xlabel=L"$p^2$ [GeV$^2$]",minorgrid=true)
plot!(plt3,legend=:left,xlabel=L"$p^2$ [GeV$^2$]",minorgrid=true)
plot!(plt4,legend=:left,xlabel=L"$p^2$ [GeV$^2$]",minorgrid=true)

l = @layout [a; b; c]
plt = plot(plt2, plt3, plt4, layout = l, size=(size[1],3*size[2]))
savefig(plt,"plots/figure4_gluon_self_energy.pdf")