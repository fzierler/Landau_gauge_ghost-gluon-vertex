using Pkg; Pkg.activate(".")
using Plots
using JLD
using UnquenchingQCD
using Parameters
using LaTeXStrings
using Measures 
pgfplotsx(tickfontsize=14,labelfontsize=16,legendfontsize=15,titlefontsize=20,framestyle=:box,lw=3,size=(600,400),color_palette=palette(:Dark2_3))

_update_logticks!(plt,log_tickval) = plot!(plt,xticks=(exp10.(log_tickval), (x->LaTeXString("\$ 10 ^ {$x} \$")).(log_tickval)),minorgrid=true)
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
function _load_correct_vertex(pathtodata,gauge_group,Nf)
    ggvpath = joinpath(pathtodata,gauge_group,"Nf=$Nf","vertex_scaling","GGV.jld")
    if !isfile(ggvpath) 
        ggvpath = joinpath(pathtodata,gauge_group,"Nf=$Nf","vertex","GGV.jld")
    end
    x,y,cosθ,GGV = loadGGV(ggvpath)
    return x,y,cosθ,GGV
end
function plot_vertex_1d(pathtodata,groups,Ncs,Nfs)   
    linestyles = (:solid,:dash,:dashdot,:dot,:dashdotdot)
    xlab = L"$p^2$ [GeV$^2$]"
    
    pltGGV1 = plot(legend=:outerright,title=L"Orthogonal momenta: $A(p^2,p^2,$cos$\theta=0)$",  xlabel=xlab,bottom_margin=8mm)
    pltGGV2 = plot(legend=:outerright,title=L"Vanishing ghost momentum: $A(p^2,0,$cos$\theta)$",xlabel=xlab,bottom_margin=8mm)
    pltGGV3 = plot(legend=:outerright,title=L"Vanishing gluon momentum: $A(0,p^2,$cos$\theta)$",xlabel=xlab,bottom_margin=8mm)
    
    for (j,Nf) in enumerate(Nfs)
        gauge_group = groups[j]*"$(Ncs[j])"
        gauge_group_fmt = replace(gauge_group,last(gauge_group)=>"($(last(gauge_group)))")
        
        x,y,cosθ,GGV = _load_correct_vertex(pathtodata,gauge_group,Nf)
        GGV_van_x, GGV_van_y, GGV_equalmomenta = _vertex_special_1d(GGV,cosθ)
        
        plot!(pltGGV1,x,GGV_equalmomenta,label=L"%$gauge_group_fmt : N_f=%$Nf",ls=linestyles[j])
        plot!(pltGGV2,x,GGV_van_x,label=L"%$gauge_group_fmt : N_f=%$Nf",ls=linestyles[j])
        plot!(pltGGV3,x,GGV_van_y,label=L"%$gauge_group_fmt : N_f=%$Nf",ls=linestyles[j])
    end
    xscale = :log10
    plot!(pltGGV1;xaxis=xscale)
    plot!(pltGGV2;xaxis=xscale)
    plot!(pltGGV3;xaxis=xscale)
    _update_logticks!.([pltGGV1,pltGGV2,pltGGV3],Ref(collect(-7:2:7)))
    return pltGGV1,pltGGV2,pltGGV3
end

pathtodata = "./input"

groups = ["Sp","Sp","Sp"]
Ncs = [4,4,4]
Nfs = ["0.0","2.0","4.50"]
plt1, plt2, plt3 = plot_vertex_1d(pathtodata,groups,Ncs,Nfs)
savefig(plt1,"plots/vertex_Sp4_Nf_dependence.pdf")

groups = ["SO","SO","SO"]
Ncs = [8,8,8]
Nfs = ["0.0","2.0","5.68"]
plt1, plt2, plt3 = plot_vertex_1d(pathtodata,groups,Ncs,Nfs)
savefig(plt1,"plots/vertex_SO8_Nf_dependence.pdf")

groups = ["SU","Sp","SO"]
Ncs = [3,4,8]
Nfs = ["2.0","2.0","2.0"]
plt1, plt2, plt3 = plot_vertex_1d(pathtodata,groups,Ncs,Nfs)
savefig(plt1,"plots/vertex_Nf2_gauge_groups.pdf")

groups = ["SU","SU","SU","SU"]
Ncs = [2,3,4,5]
Nfs = ["2.0","2.0","2.0","2.0"]
plt1, plt2, plt3 = plot_vertex_1d(pathtodata,groups,Ncs,Nfs)
#plot!(plt1,legend=:topleft)
savefig(plt1,"plots/vertex_SUN_largeN.pdf")

groups = ["SU","SU","SU","SU"]
Ncs = [2,3,4,5]
Nfs = ["2.0","3.0","4.0","5.0"]
plt1, plt2, plt3 = plot_vertex_1d(pathtodata,groups,Ncs,Nfs)
#plot!(plt1,legend=:topleft)
savefig(plt1,"plots/vertex_SUN_largeN_Nf_over_Nc_fixed.pdf")
