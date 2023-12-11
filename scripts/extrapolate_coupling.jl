using LsqFit, JLD
using Plots
include("IterationHelper.jl")
pathtodata = "/home/zierler_fabian/Documents/daten/UnquenchingQCD"
pathtodata = "/home/fabian/Dokumente/daten/UnquenchingQCD"

Nf = Float64[]
g2 = Float64[]
strings = ["0.0","0.5","1.0","1.1","1.15","1.35","1.4","1.5","1.8"] #SU(5)
strings = ["0","0.5","0.8","0.9","1.1","1.2","1.5"] #SU(4)
strings = ["0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1"] #SU(3)
strings = ["0","0.1","0.2","0.4","0.6","0.8"] #SU(2)
strings = ["0","0.5","0.6","0.7","0.8","0.9","1.0","1.5"] #Sp(4)
strings = ["0","0.8","0.9","1.1","1.2","1.5"] #Sp(6)
strings = ["0","0.5","0.6","0.8","1.0","1.4","1.5","2.0"] #Sp(8)
strings = ["0.0","0.3","0.4","0.6","1.0"] #SO(6)
strings = ["0.0","0.1","0.5","0.8","0.9","1.5"] #SO(7)
strings = ["0.0","0.5","0.6","0.65","0.85","0.9","1.0","1.5"] #SO(8)
strings = ["0.0","0.5","0.7","0.8","1.0","1.1","1.2","1.5"] #SO(8)

for i in strings
    qcd = load(joinpath(pathtodata,"SO9","Nf=$i/prop/theory.jld"))
    push!(g2,qcd["g2"])
    push!(Nf,parse(Float64,i))
end
#@. model(x, p) = p[1] + p[2]*x + x*x*p[3]
#fit = curve_fit(model, Nf[3:end], g2[3:end], [1.0,1.0,1.0])
#@info "param = $(fit.param)"

#@show evalpoly(3,fit.param)
#Nfs = 0.01:0.01:3
#scatter(Nf,g2)
#plot!(Nfs,g2_a1_GIR15_SU4.(Nfs))
#plot!(Nfs,g2_a1_GIR15_Sp6.(Nfs))
