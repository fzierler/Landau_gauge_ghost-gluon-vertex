function savedata(path,x,A,B,p2,coeffG,coeffZ,Z2,Zm,Z3,Z3_tilde,QCD::theory,param::parameter;plot=true)
    ispath(path) || mkpath(path) # create that folder
    # interpolate G, Z from Chebyshev nodes for plotting
    p = exp10.(range(log10(param.ϵ),length=param.nCheby*100,stop=log10(param.Λ)))
    pQ= exp10.(range(log10(minimum(x)),length=param.nCheby*100,stop=log10(maximum(x))))
    G = exp.(clenshaw(p2,coeffG,param.ϵ,param.Λ))
    Z = exp.(clenshaw(p2,coeffZ,param.ϵ,param.Λ))
    μ  = renormalization_point(p2,G,Z,QCD)
    fπ = pion_decay(x,A,B,Z2,QCD.Nc)
    Gint = exp.(clenshaw(p,coeffG,param.ϵ,param.Λ))
    Zint = exp.(clenshaw(p,coeffZ,param.ϵ,param.Λ))
    save(joinpath(path,"hybrid.jld"),"x",x,"A",A,"B",B,"p2",p2,"coeffG",coeffG,"coeffZ",coeffZ,"Z2",Z2,"fπ",fπ)
    save(joinpath(path,"data.jld"),"x",x,"A",A,"B",B,"p2",p2,"G",G,"Z",Z,"Z2",Z2,"fπ",fπ,"µ",µ,"Z3",Z3,"Z3_tilde",Z3_tilde,"Zm",Zm)
    save_structure(joinpath(path,"theory.jld"),QCD)
    save_structure(joinpath(path,"parameters.jld"),param)
    if plot
        plotpath = joinpath(path,"Plots")
        ispath(plotpath) || mkpath(plotpath)
        plot_alpha(p,Gint,Zint,QCD.g2,plotpath)
        plot_ym(p,Gint,Zint,plotpath)
        plot_ym_UV(p,Gint,Zint,QCD,plotpath)
        plot_ym_IR(p,Gint,Zint,QCD,plotpath)
        plot_quark(x,A,B,plotpath,"Quark.pdf")
    end
end
function plotdata(path,x,A,B,p2,coeffG,coeffZ,ϵ,Λ,nCheby,QCD)
    p = exp10.(range(log10(ϵ),length=nCheby*100,stop=log10(Λ)))
    pQ= exp10.(range(log10(minimum(x)),length=nCheby*100,stop=log10(maximum(x))))
    G = exp.(clenshaw(p2,coeffG,ϵ,Λ))
    Z = exp.(clenshaw(p2,coeffZ,ϵ,Λ))
    Gint = exp.(clenshaw(p,coeffG,ϵ,Λ))
    Zint = exp.(clenshaw(p,coeffZ,ϵ,Λ))
    # create plots
    plotpath = joinpath(path,"Plots")
    ispath(plotpath) || mkpath(plotpath)
    plot_alpha(p,Gint,Zint,QCD.g2,plotpath)
    plot_ym(p,Gint,Zint,plotpath)
    plot_ym_UV(p,Gint,Zint,QCD,plotpath)
    plot_ym_IR(p,Gint,Zint,QCD,plotpath)
    plot_quark(x,A,B,plotpath,"Quark.pdf")
end
function redoplots(plotpath)
    d = load(joinpath(plotpath,"hybrid.jld"))
    p = load(joinpath(plotpath,"parameters.jld"))
    q = load(joinpath(plotpath,"theory.jld"))
    A, B, Nf = d["A"],d["B"],q["Nf"]
    x,p2,coeffG,coeffZ = d["x"],d["p2"],d["coeffG"],d["coeffZ"]
    ϵ,Λ,nCheby = p["ϵ"],p["Λ"],Int(p["nCheby"])
    QCD = _load_theory(joinpath(plotpath,"theory.jld"))
    plotdata(plotpath,x,A,B,p2,coeffG,coeffZ,ϵ,Λ,nCheby,QCD)
    try
        v = load(joinpath(plotpath,"GGV.jld"))
        plotGGV(plotpath,v["x"],v["y"],v["costh"],v["GGV"])
    catch
        nothing
    end
end
function save_structure(path,structure)
    type = typeof(structure)
    jldopen(path, "w") do file
        for i in 1:nfields(structure)
            val = getfield(structure,fieldname(type,i))
            str = string(fieldname(type,i))
            typeof(val) ∈ (Symbol,String) ? write(file, str, string(val)) : write(file, str, Float64(val))
        end
    end
end
function load_hybrid(path,x,n,ϵ,Λ,mu)
    p2 = chebyshevnodes(n,ϵ,Λ)
    p2 = append!(p2,mu)
    d  = load(path)
    prm = load(joinpath(dirname(path),"parameters.jld"))
    x0, p, coeffG, coeffZ = d["x"], d["p2"], d["coeffG"], d["coeffZ"]
    A, B = d["A"], d["B"]
    G,Z = zero(p2), zero(p2)
    QCD = _load_theory(joinpath(dirname(path),"theory.jld"))
    param = parameter(ϵ=prm["ϵ"],Λ=prm["Λ"])
    interpolate_ym!(p2,G,Z,p,coeffG,coeffZ,QCD,param)
    coeffG0 = chebycoeff(log.(G[1:end-1]))
    coeffZ0 = chebycoeff(log.(Z[1:end-1]))
    A0 = splines(x,x0,A)
    B0 = splines(x,x0,B)
    return A0, B0, p2, coeffG0, coeffZ0
end
function plot_alpha(p2,G,Z,g2,filename)
    n = length(p2)
    plot(p2[1:n-1],G[1:n-1].*G[1:n-1].*Z[1:n-1]*(g2/(4*π)))
    plot!(title="alpha(p^2)",xaxis=:log10,xticks=10.0 .^(-6:6))
    savefig(joinpath(filename,"Coupling.pdf"))
end
function plot_ym(p2,G,Z,path)
    n = length(p2)
    plot(p2[1:n-1],G[1:n-1],label="G(p^2)")
    plot!(p2[1:n-1],Z[1:n-1],label="Z(p^2)",xaxis=:log10)
    savefig(joinpath(path,"YM.pdf"))
    plot(p2[1:n-1],Z[1:n-1]./p2[1:n-1],label="Z(p^2)/p^2",yaxis=:log10,xaxis=:log10)
    savefig(joinpath(path,"GluonProp.pdf"))
end
function plot_quark(p2,A,B,filename,name="Quark.pdf")
    n = length(p2)
    plot(p2[1:n-1],(B ./ A)[1:n-1],label="M(p^2)",xaxis=:log10)
    savefig(joinpath(filename,"M"*name))
    plot!(p2[1:n-1],B[1:n-1],label="B(p^2)",xaxis=:log10)
    plot(p2[1:n-1],(1 ./ A)[1:n-1],label="1/A(p^2)",xaxis=:log10)
    savefig(joinpath(filename,"A"*name))
    savefig(joinpath(filename,name))
end
function plot_ym_UV(p2,G,Z,QCD,filename)
    UV,n  = findmax(p2)
    logx  = log10(UV)
    x     = exp10.(collect(logx:0.01:(logx+3)))
    β0 = (11QCD.CA-4*QCD.T*QCD.Nf)/3
    ω  = β0*QCD.g2*G[n]^2*Z[n]/(16*pi^2)
    Zextr = @. Z[n]*(ω*log(x/UV) + 1)^QCD.γ
    Gextr = @. G[n]*(ω*log(x/UV) + 1)^QCD.δ
    plot(p2[end-10^3:end-1],Z[end-10^3:end-1],label="Z")
    plot!(p2[end-10^3:end-1],G[end-10^3:end-1],label="G")
    plot!(x,Zextr,label="")
    plot!(x,Gextr,label="",xaxis=:log10,yaxis=:log10)
    savefig(joinpath(filename,"UV.pdf"))
end
function plot_ym_IR(p2,G,Z,QCD,filename)
    p2min,n= findmin(p2)
    logx   = log10(p2min)
    IR  = exp10.(collect((logx-3):0.01:logx))
    ZIR = @. Z[n]*(IR/p2min)^QCD.κGl
    GIR = @. G[n]*(IR/p2min)^QCD.κGh
    plot(p2,Z./p2,label="Z(p^2)/p^2")
    plot!(p2,G,label="G(p^2)")
    plot!(IR,ZIR./IR,label="")
    plot!(IR,GIR,label="")
    plot!(xaxis=:log10)
    savefig(joinpath(filename,"IR.pdf"))
    plot!(xaxis=:log10,yaxis=:log10)
    savefig(joinpath(filename,"IRloglog.pdf"))
end
#################################
# GGV
#################################
function saveGGV(x,y,cosθ,GGV,paramGGV,QCD,path)
    file = joinpath(path,"GGV.jld")
    ispath(path) || mkpath(path)
    save(file,"x",x,"y",y,"costh",cosθ,"GGV",GGV)
    save_structure(joinpath(path,"parametersGGV.jld"),paramGGV)
    save_structure(joinpath(path,"theory.jld"),QCD)
    plotGGV(path,x,y,cosθ,GGV,"GGV")
end
function saveGGV(x,y,cosθ,GGV,p2,G,Z,paramGGV,QCD,path)
    file = joinpath(path,"GGV.jld")
    ispath(path) || mkpath(path)
    save(file,"x",x,"y",y,"costh",cosθ,"GGV",GGV,"p2",p2,"G",G,"Z",Z)
    save_structure(joinpath(path,"parametersGGV.jld"),paramGGV)
    save_structure(joinpath(path,"theory.jld"),QCD)
    plotGGV(path,x,y,cosθ,GGV,"GGV")
end
function loadGGV(pathGGV)
    d = load(pathGGV); T = Float64
    x, y, cosθ, GGV = d["x"], d["y"], d["costh"], d["GGV"]
    return x::Array{T,1}, y::Array{T,1}, cosθ::Array{T,1}, GGV::Array{T,3}
end
function plotGGV(savepath,x,y,cosθ,GGV,name="GGV")
    θind = div(length(cosθ),2)
    wireframe(log10.(x),log10.(y),GGV[:,:,θind])
    plot!( xlabel="log10(x)", ylabel="log10(y)")
    plot!( title = "cosθ = $(cosθ[θind]), ghost: x [GeV^2], gluon: y [GeV^2]")
    savefig(joinpath(savepath,"Plots","$(name).pdf"))
    # now plot symmetric momentum configuration
    @assert isequal(x,y)
    θindex = findmin(abs.(cosθ.-1/2))[2]
    GGV_equalmomenta = zeros(length(x))
    for i in 1:length(x)
        GGV_equalmomenta[i] = GGV[i,i,θindex]
    end
    plot(x,GGV_equalmomenta,xaxis=:log10)
    savefig(joinpath(savepath,"Plots",name*"symm.pdf"))
end
function _load_propagator_parameters(file;kws...)
    dict  = load(file)
    # overwrite parameters if the are provided by keyword arguments
    for i in kws
        dict[String(i.first)] = i.second
    end
    param = parameter(
        #renormalization
        µ=dict["μ"],
        ZUV=dict["ZUV"], # value of dressing at gluon subtraction point
        GIR=dict["GIR"], # value of the ghost dressing at the lowest external momentum
        # iteration
        maxiterNR=dict["maxiterNR"], # maximal number of iterations using Newton's method
        linesearch=dict["linesearch"], # maximal number of steps in line search
        epsilonNR=dict["epsilonNR"], #convergence criterion for Newton's method
        quarkrelax=dict["quarkrelax"], # relaxation parameter for quark iteration
        h=dict["h"], #step size for forward differentiation in Jacobian
        propiter=dict["propiter"],
        # integration and interpolation
        nr=dict["nr"],  # radial integration points
        na=dict["na"],  # angular integration points
        nrQ=dict["nrQ"], # radial integration points
        naQ=dict["naQ"], # angular integration points
        nCheby=dict["nCheby"], # number of Chebyshev polynomials for expansion
        nQ=dict["nQ"], # external momenta for quark integration
        # external momentum range and integeration range
        IR=dict["IR"], # IR cutoff
        UV=dict["UV"], # UV cutoff
        ϵ=dict["ϵ"], # smallest external momentum
        Λ=dict["Λ"], # highest external momentum
        # quark iteration
        maxiterQ=dict["maxiterQ"],
        ϵQ=dict["ϵQ"],
        ϵtotal=dict["ϵtotal"],
        splinesQ=dict["splinesQ"],
        # 3-Gluon Vertex Model
        GluonModel=Symbol(dict["GluonModel"]),
        hIR=dict["hIR"],
        Λ3G=dict["Λ3G"],
    )
    return param
end
function _load_theory(file;kws...)
    dict = load(file)
    # overwrite parameters if the are provided by keyword arguments
    for i in kws
        dict[String(i.first)] = i.second
    end
    theo = theory(
        Nc = dict["Nc"],
        Nf = dict["Nf"],
        g2 = dict["g2"],
        mQ = dict["mQ"],
        # theory and representation
        # IR power law coefficient
        κGl = dict["κGl"],
        κGh = dict["κGh"],
        # Casimirs 
        CA = dict["CA"],
        CF = dict["CF"],
        T  = dict["T"],
        # general γ aus Muta 3.4.18
        γ = dict["γ"],
        δ = dict["δ"],
        # exponents for three-gluon vertex model
        a = dict["a"], # scaling - for decoupling choose a = -1
        α = dict["α"],
        β = dict["β"]
    )
    return theo 
end
function _load_vertex_parameters(file;kws...)
    dict = load(file)
    # overwrite parameters if the are provided by keyword arguments
    for i in kws
        dict[String(i.first)] = i.second
    end
    # construct structure
    param = parameterGGV(
        nr = dict["nr"],
        na1 = dict["na1"],
        na2 = dict["na2"],
        IR = dict["IR"],
        UV = dict["UV"],
        maxiter = dict["maxiter"],
        converged = dict["converged"],
        sysitermax = dict["sysitermax"],
        eps_sys = dict["eps_sys"],
        rel = dict["rel"],
        # spacing of external momenta
        Δxy = dict["Δxy"],
        Δcosθ = dict["Δcosθ"],
        logϵ = dict["logϵ"],
        logΛ = dict["logΛ"],
    )
    return param    
end