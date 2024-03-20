using Parameters
using LsqFit, JLD, UnquenchingQCD, Plots
function mainPropagators(dirProp,savepath;Δg2=0,ΔNf=0,mQ=0,GIR=15,relax=1,κGl=1,κGh=0)
    pathYM  = joinpath(dirProp,"hybrid.jld")
    dataQCD = load(joinpath(dirProp,"theory.jld"))
    dataPRM = load(joinpath(dirProp,"parameters.jld"))
    #############################################################
    # further changes in the parameters can be provided through
    # keyword arguements to the loading function
    #############################################################
    g2 = dataQCD["g2"] + Δg2
    Nf = dataQCD["Nf"] + ΔNf
    QCD   = _load_theory(joinpath(dirProp,"theory.jld");g2,Nf,mQ)
    param = _load_propagator_parameters(joinpath(dirProp,"parameters.jld");GIR,relax,κGl,κGh)
    nr  = max(param.nr,1000)
    nrQ = max(param.nrQ,1000)
    nQ  = max(param.nQ,1000)
    nCheby = max(param.nCheby,40)
    param  = _load_propagator_parameters(joinpath(dirProp,"parameters.jld");GIR,relax,κGl,κGh,nr,nrQ,nCheby,nQ)
    log_structure(QCD)
    log_structure(param)
    #########################################################
    x  = gausslegendre(param.nQ,param.ϵ/10,10*param.Λ)[1]
    A0,B0,p2,coeffG0,coeffZ0 = load_hybrid(pathYM,x,param.nCheby,param.ϵ,param.Λ,param.µ)
    Quad,Kernels = quad_kernels(p2,param,QCD)
    A,B,coeffG,coeffZ,RenConstants,success = iterate_propagators(x,A0,B0,p2,coeffG0,coeffZ0,QCD,param,Quad,Kernels)
    if success
        @unpack Z2,Zm,Z3,Z3_tilde  = RenConstants
        savedata(savepath,x,A,B,p2,coeffG,coeffZ,Z2,Zm,Z3,Z3_tilde,QCD,param)
        G = exp.(UnquenchingQCD.clenshaw(p2[1:end-1],coeffG,param.ϵ,param.Λ))
        Z = exp.(UnquenchingQCD.clenshaw(p2[1:end-1],coeffZ,param.ϵ,param.Λ))
        α = @. QCD.g2*G^2*Z/(4*π)
        data  = load(joinpath(savepath,"data.jld"))
        fπ, µ = data["fπ"], data["µ"]
        M0 = B[1]/A[1]
        @info "fπ, µ, M0 = $fπ, $µ, $M0"
        αmax, p2max_ind = findmax(α) 
        @info "maximal coupling: α(p²=$(p2[p2max_ind]))=$αmax"
        return fπ, M0
    end
    return nothing, nothing
end
function mainSystem(dirsys,savepath;vertexbare=false,Δg2=0,ΔNf=0,mQ=0,GIR=15,relax=1,quarkrelax=1,κGl=1,κGh=0)
    pathYM  = joinpath(dirsys,"hybrid.jld")
    pathGGV = joinpath(dirsys,"GGV.jld")
    dataQCD = load(joinpath(dirsys,"theory.jld"))
    dataPRM = load(joinpath(dirsys,"parameters.jld"))
    #############################################################
    # further changes in the parameters can be provided through
    # keyword arguements to the loading function
    #############################################################
    g2 = dataQCD["g2"] + Δg2
    Nf = dataQCD["Nf"] + ΔNf
    QCD = _load_theory(joinpath(dirsys,"theory.jld");g2,Nf,mQ)
    paramYM = _load_propagator_parameters(joinpath(dirsys,"parameters.jld");GIR,relax=quarkrelax,κGl,κGh,propiter=1)
    nr  = max(paramYM.nr,1000)
    nrQ = max(paramYM.nrQ,1000)
    nQ  = max(paramYM.nQ,1000)
    nCheby = max(paramYM.nCheby,40)
    paramYM = _load_propagator_parameters(joinpath(dirsys,"parameters.jld");GIR,relax=quarkrelax,κGl,κGh,nr,nrQ,nCheby,nQ,propiter=1)
    log_structure(QCD)
    log_structure(paramYM)
    #########################################################
    xQ = gausslegendre(paramYM.nQ,paramYM.ϵ/10,10*paramYM.Λ)[1]
    A0,B0,p2,coeffG0,coeffZ0 = load_hybrid(pathYM,xQ,paramYM.nCheby,paramYM.ϵ,paramYM.Λ,paramYM.µ)
    if vertexbare
        paramGGV = parameterGGV(rel=relax,nr=240,Δxy=0.25,Δcosθ=0.222)
        log_structure(paramGGV)        
        x, y, cosθ, GGV0 = defaultGGV(logϵ=paramGGV.logϵ,logΛ=paramGGV.logΛ,Δxy=paramGGV.Δxy,Δcosθ=paramGGV.Δcosθ)
    else
        paramGGV = _load_vertex_parameters(joinpath(dirsys,"parametersGGV.jld");rel=relax)
        nr  = max(paramGGV.nr,240)
        Δxy = min(paramGGV.Δxy,0.25)
        paramGGV = _load_vertex_parameters(joinpath(dirsys,"parametersGGV.jld");rel=relax,Δxy,nr)
        log_structure(paramGGV)        
        x, y, cosθ, GGV0 = defaultGGV(logϵ=paramGGV.logϵ,logΛ=paramGGV.logΛ,Δxy=paramGGV.Δxy,Δcosθ=paramGGV.Δcosθ)
        x0, y0, cosθ0, GGVload = loadGGV(pathGGV)
        GGV0 = interpolateGGV(x,y,cosθ,GGVload,x0,y0,cosθ0)
    end
    A,B,coeffG,coeffZ,RenConstants,x,y,cosθ,GGV,success = IterateSystem(xQ,A0,B0,p2,coeffG0,coeffZ0,x,y,cosθ,GGV0,QCD,paramYM,paramGGV)
    if true
        @unpack Z2,Zm,Z3,Z3_tilde  = RenConstants
        savedata(savepath,xQ,A,B,p2,coeffG,coeffZ,Z2,Zm,Z3,Z3_tilde,QCD,paramYM)
        saveGGV(x,y,cosθ,GGV,paramGGV,QCD,savepath)
        G = exp.(UnquenchingQCD.clenshaw(p2[1:end-1],coeffG,paramYM.ϵ,paramYM.Λ))
        Z = exp.(UnquenchingQCD.clenshaw(p2[1:end-1],coeffZ,paramYM.ϵ,paramYM.Λ))
        α = @. QCD.g2*G^2*Z/(4*π)
        data  = load(joinpath(savepath,"data.jld"))
        fπ, µ = data["fπ"], data["µ"]
        M0 = B[1]/A[1]
        @info "fπ, µ, M0 = $fπ, $µ, $M0"
        αmax, p2max_ind = findmax(α) 
        @info "maximal coupling: α(p²=$(p2[p2max_ind]))=$αmax"
        return fπ, M0
    end
    return nothing, nothing
end

function defaultGGV(;logϵ=-7.01,logΛ=8.0,Δxy=0.10,Δcosθ=0.3333)
    x    = exp10.(collect(logϵ:Δxy:logΛ))
    y    = exp10.(collect(logϵ:Δxy:logΛ))
    cosθ = collect(-0.9999:Δcosθ:0.9999)
    GGV  = ones(length(x),length(y),length(cosθ))
    return x, y, cosθ, GGV
end
# rewrite in terms of fπ_target
function SetScaleProp(dirProp,savepath;mQ=0,GIR=15,target=0.075,Δtarget=0.001)
    fπ, M0 = mainPropagators(dirProp,savepath;Δg2=0,mQ,GIR)
    fπ > target + Δtarget && (s = -1)
    fπ < target - Δtarget && (s = +1)
    (target - Δtarget < fπ < target + Δtarget) && return nothing
    while !(target - Δtarget < fπ < target + Δtarget)
        step = find_step_size(fπ;target,Δtarget)
        @info "fπ, step = $fπ, $step"
        fπ, M0 = mainPropagators(savepath,savepath;Δg2=step,mQ,GIR)
        s == -1 && fπ < target - Δtarget && break
        s == +1 && fπ > target + Δtarget && break
    end
end
function SetScalePropM0(dirProp,savepath;mQ=0,GIR=15,target=0.5,Δtarget=0.01)
    fπ, M0 = mainPropagators(dirProp,savepath;Δg2=0,mQ,GIR)
    M0 > target + Δtarget && (s = -1)
    M0 < target - Δtarget && (s = +1)
    (target - Δtarget < M0 < target + Δtarget) && return nothing
    while !(target - Δtarget < M0 < target + Δtarget)
        step = find_step_size(M0;target,Δtarget)
        @info "M0, step = $M0, $step"
        fπ, M0 = mainPropagators(savepath,savepath;Δg2=step,mQ,GIR)
        s == -1 && M0 < target + Δtarget && break
        s == +1 && M0 > target - Δtarget && break
    end
end
function SetScaleSystem(dirsys,savepath;mQ=0,GIR=15,target=0.075,Δtarget=0.01)
    fπ = mainSystem(dirsys,savepath;vertexbare=false,Δg2=0,mQ,GIR)
    fπ > target + Δtarget && (s = -1)
    fπ < target - Δtarget && (s = +1)
    target - Δtarget < fπ < target + Δtarget && return nothing
    while !(target - Δtarget < fπ < target + Δtarget)
        step = find_step_size(fπ;target,Δtarget)
        @info "fπ, step = $fπ, $step"
        fπ   = mainSystem(savepath,savepath;vertexbare=false,Δg2=step,mQ,GIR)
        s == -1 && fπ < target + Δtarget && break
        s == +1 && fπ > target - Δtarget && break
    end
end
function find_step_size(x;target=0.35,Δtarget=0.01)
    if x < target*0.3
        return 0.020
    elseif x < target*0.5
        return 0.010
    elseif x < target*0.9
        return 0.005
    elseif x < target - Δtarget
        return 0.005
    elseif x > target*2.0
        return -0.01
    elseif x > target*1.1
        return -0.005
    elseif x > target + Δtarget
        return -0.002
    end
end
function extrapolate_GIR(loadpath;IR=-4.5)
    data = load(joinpath(loadpath,"hybrid.jld"))
    param = load(joinpath(loadpath,"parameters.jld"))
    x, p2, coeffG, coeffZ =data["x"], data["p2"], data["coeffG"], data["coeffZ"]
    p = exp10.(log10(minimum(p2)):0.01:log10(maximum(p2)))
    pfit = exp10.(IR:0.01:0)
    G0 = exp.(UnquenchingQCD.clenshaw(p,coeffG,param["ϵ"],param["Λ"]))
    Z0 = exp.(UnquenchingQCD.clenshaw(p,coeffZ,param["ϵ"],param["Λ"]))
    G = exp.(UnquenchingQCD.clenshaw(pfit,coeffG,param["ϵ"],param["Λ"]))
    Z = exp.(UnquenchingQCD.clenshaw(pfit,coeffZ,param["ϵ"],param["Λ"]))

    @. model1(x, p) = p[1]*x^-(p[2])
    @. model2(x, p) = p[1]*x^(2p[2])
    fit1 = curve_fit(model1, pfit, G, ones(2))
    fit2 = curve_fit(model2, pfit, Z, ones(2))
    GIR = model1(param["ϵ"],fit1.param)
    ρ1 = fit1.param[2]
    ρ2 = fit2.param[2]

    # create a plot for visualisation
    plt = plot(p,G0,xscale=:log10,yscale=:log10,lw=2,label="G (ρ=$(round(ρ1,digits=3)))")
    plot!(plt, pfit,G,xscale=:log10,yscale=:log10,lw=3,label="")
    scatter!(plt, [param["ϵ"]],[GIR],label="GIR=$(round(GIR,digits=2))")

    return GIR, ρ1, ρ2, plt
end


g2_a1_GIR15_SU2(Nf) = 4.0925510186590730 - Nf*2.12804544164290330 + Nf*Nf*1.390676478553773700
g2_a1_GIR15_SU3(Nf) = 2.0209668357753214 - Nf*0.31430595000991635 + Nf*Nf*0.162924128645046900
g2_a1_GIR15_SU4(Nf) = 1.4306566023886830 - Nf*0.15610000021861350 + Nf*Nf*0.059000000048131274
g2_a1_GIR15_SU5(Nf) = 1.1184907300259854 - Nf*0.10615754639695112 + Nf*Nf*0.031133763653489136

g2_a1_GIR15_Sp4(Nf) = 2.1557529477944610 - Nf*0.49710095331420906 + Nf*Nf*0.2282571498415587
g2_a1_GIR15_Sp6(Nf) = 1.4169066023494337 - Nf*0.16660000018182833 + Nf*Nf*0.06600000004003292
g2_a1_GIR15_Sp8(Nf) = 1.1147226605521388 - Nf*0.12134571997565849 + Nf*Nf*0.035209024552546306

g2_a1_GIR15_SO9(Nf) = 0.7660499999994401 - Nf*0.045299999999475295 + Nf*Nf*0.02299999999988449
g2_a1_GIR15_SO8(Nf) = 0.9692000000007700 - Nf*0.136657142858056750 + Nf*Nf*0.051714285714518166
g2_a1_GIR15_SO7(Nf) = 1.1165560853520240 - Nf*0.147662574817760460 + Nf*Nf*0.0727063677245942
g2_a1_GIR15_SO6(Nf) = 1.3998999980021245 - Nf*0.275799997775127000 + Nf*Nf*0.14999999945445422
g2_a1_GIR15_SO5(Nf) = 1.7004285680412565 - Nf*0.223428566696477540 + Nf*Nf*0.2171428555799307
g2_a1_GIR15_SO4(Nf) = 2.9763762082616880 - Nf*1.682573057107197000 + Nf*Nf*1.1730931749470026
g2_a1_GIR15_SO3(Nf) = 3.7907999975381452 + Nf*2.076000003381764700
