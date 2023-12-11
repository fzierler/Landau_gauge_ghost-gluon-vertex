using Pkg; Pkg.activate(".")
using Plots
using JLD
using UnquenchingQCD
using Parameters
using LaTeXStrings
using Measures 

size=(600,400)
pgfplotsx(tickfontsize=14,labelfontsize=16,legendfontsize=15,titlefontsize=20,framestyle=:box,lw=3,size=size)

pathtodata = "./input"
isdir("plots") || mkdir("plots")

pgfplotsx(color_palette=reverse(palette(:Paired_6)))
xlab = L"$p^2$ [GeV$^2$]"
plt1 = plot(title=L"Quark dressing function $A_f(p^2)$",xlabel=xlab,xscale=:log10,legend=:topright,bottom_margin=8mm)
plt2 = plot(title=L"Quark mass function $M_f(p^2)$ [GeV]",xlabel=xlab,xscale=:log10,legend=:topright,bottom_margin=8mm)
plt3 = plot(title=L"Ghost dressing $G(p^2)$",xlabel=xlab,xscale=:log10,legend=:topright,bottom_margin=8mm)
plt4 = plot(title=L"Gluon propagator $D(p^2)$",xlabel=xlab,xscale=:log10,legend=:topright,bottom_margin=8mm)
plt5 = plot(title=L"Coupling $\alpha(p^2)$",xlabel=xlab,xscale=:log10,legend=:topright,bottom_margin=8mm)
plt6 = plot(title=L"Gluon dressing $Z(p^2)$",xlabel=xlab,xscale=:log10,legend=:topleft,bottom_margin=8mm)
pgfplotsx(color_palette=palette(:Dark2_3))
pltGGV1 = plot(legend=:topright,title=L"Orthogonal momenta: $A(p^2,p^2,$cos$\theta=0)$",  xlabel=xlab,bottom_margin=8mm)
pltGGV2 = plot(legend=:topright,title=L"Vanishing ghost momentum: $A(p^2,0,$cos$\theta)$",xlabel=xlab,bottom_margin=8mm)
pltGGV3 = plot(legend=:topright,title=L"Vanishing gluon momentum: $A(0,p^2,$cos$\theta)$",xlabel=xlab,bottom_margin=8mm)

pgfplotsx(color_palette=reverse(palette(:Paired_6)))
for Nf in [0,2,5], type in ["vertex","prop"]
    if Nf == 5
        type = type*"_scaling"
    end
    qcd = load(joinpath(pathtodata,"SU3","Nf=$Nf",type,"theory.jld"))
    data = load(joinpath(pathtodata,"SU3","Nf=$Nf",type,"hybrid.jld"))
    param = load(joinpath(pathtodata,"SU3","Nf=$Nf",type,"parameters.jld"))
    x, p2, A, B, coeffG, coeffZ =data["x"], data["p2"], data["A"], data["B"], data["coeffG"], data["coeffZ"]
    p = exp10.(log10(minimum(p2)):0.1:log10(maximum(p2)))
    G = exp.(UnquenchingQCD.clenshaw(p,coeffG,param["ϵ"],param["Λ"]))
    Z = exp.(UnquenchingQCD.clenshaw(p,coeffZ,param["ϵ"],param["Λ"]))
    ls = ifelse(contains(type,"vertex"),:solid,:dash)
    db = ifelse(contains(type,"prop"),"bare","dynamic")
    ticks = exp10.(-5:1:4)
    l =L"$N_f=%$Nf$ (%$db)"
    xlims = (10^-5, 5*10^4)
    plot!(plt1,x,A,label=l,ls=ls,xlims=xlims,xticks=ticks)
    plot!(plt2,x,B./A,label=l,ls=ls,xlims=xlims)
    plot!(plt3,p,G,label=l,ls=ls,xlims=xlims)
    plot!(plt4,p,Z./p,label=l,ls=ls,yaxis=:log10,xlims=xlims)
    plot!(plt5,p,G.*G.*Z.*(qcd["g2"])/(4π),label=l,ls=ls,xlims=xlims)
    plot!(plt6,p,Z,label=l,ls=ls,xlims=xlims)
end
for plt in [plt1,plt2,plt3,plt4,plt5,plt6]
    log_tickval = collect(-5:2:5)
    plot!(plt,xticks=(exp10.(log_tickval), (x->LaTeXString("\$ 10 ^ {$x} \$")).(log_tickval)),minorgrid=true)
end

l = @layout [a; b; c; d]
pltPROP = plot(plt1, plt2, plt6, plt3, layout = l,  size=(size[1],4*size[2]))
savefig(pltPROP,"plots/figure2_props.pdf")

pgfplotsx(color_palette=palette(:Dark2_3))
for (j,Nf) in enumerate([0,2,5])
    type = Nf == 5 ? "vertex_scaling" : "vertex"
    ggvpath = joinpath(pathtodata,"SU3","Nf=$Nf",type,"GGV.jld")
    x,y,cosθ,GGV = loadGGV(ggvpath)
    θindex = findmin(abs.(cosθ.-1/2))[2]
    GGV_equalmomenta = zeros(length(x))
    GGV_van_x = zeros(length(x))
    GGV_van_y = zeros(length(x))
    for i in 1:length(x)
        GGV_equalmomenta[i] = GGV[i,i,θindex]
        GGV_van_x[i] = GGV[1,i,end]
        GGV_van_y[i] = GGV[i,1,end]
    end
    linestyles = (:solid,:dash,:dashdot)
    plot!(pltGGV1,x,GGV_equalmomenta,xaxis=:log10,label=L"N_f=%$Nf",ls=linestyles[j])
    plot!(pltGGV2,x,GGV_van_x,xaxis=:log10,label=L"N_f=%$Nf",ls=linestyles[j])
    plot!(pltGGV3,x,GGV_van_y,xaxis=:log10,label=L"N_f=%$Nf",ls=linestyles[j])
end
for plt in [pltGGV1,pltGGV2,pltGGV3]
    log_tickval = collect(-7:2:7)
    plot!(plt,xticks=(exp10.(log_tickval), (x->LaTeXString("\$ 10 ^ {$x} \$")).(log_tickval)),minorgrid=true)
end

l = @layout [a; b; c]
pltGGV = plot(pltGGV1, pltGGV2, pltGGV3, layout = l, size=(size[1],3*size[2]))
savefig(pltGGV,"plots/figure3_vertex.pdf")
pltPROP
pltGGV