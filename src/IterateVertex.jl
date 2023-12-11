function IterateSystem(xQ,A0,B0,p2,coeffG0,coeffZ0,x,y,cosθ,GGV0,QCD,paramYM,paramGGV)
    @unpack ϵ, Λ = paramYM
    savetmp=true
    for iter in 1:paramGGV.sysitermax
        Quad,Kernels = quad_kernels(p2,paramYM,QCD)
        kernelGGVmod!(x,y,cosθ,GGV0,Kernels,Quad.p2,Quad.k2,Quad.q2)
        A,B,coeffG,coeffZ,RCs,success = iterate_propagators(xQ,A0,B0,p2,coeffG0,coeffZ0,QCD,paramYM,Quad,Kernels;savetmp)
        savetmp = success
        G = exp.(clenshaw(p2,coeffG,ϵ,Λ))
        Z = exp.(clenshaw(p2,coeffZ,ϵ,Λ))
        G0 = exp.(clenshaw(p2,coeffG0,ϵ,Λ))
        Z0 = exp.(clenshaw(p2,coeffZ0,ϵ,Λ))
        GGV = iterate_vertex(x,y,cosθ,GGV0,QCD,paramGGV,p2,G,Z)
        saveGGV(x,y,cosθ,GGV,paramGGV,QCD,"tmp")
        diff  = diff_w(A,A0) + diff_w(B,B0) + diff_w(G,G0)  + diff_w(Z,Z0)
        diff += diff_w(GGV,GGV0)
        @info "macrocycle $iter : diff = $diff"
        if diff < paramGGV.eps_sys
            return A,B,coeffG,coeffZ,RCs,x,y,cosθ,GGV,success
        end
        A0 .= A
        B0 .= B
        coeffG0 .= coeffG
        coeffZ0 .= coeffZ
        GGV0 .= GGV
        if iter == paramGGV.sysitermax
            return A,B,coeffG,coeffZ,RCs,x,y,cosθ,GGV,success
        end
    end
end
function interpolateGGV(x,y,cosθ,GGV0,x0,y0,cosθ0)
    lx, ly, lcos = length(x), length(y), length(cosθ)
    xint = zeros(lx,ly,lcos)
    yint = zeros(lx,ly,lcos)
    cosθint = zeros(lx,ly,lcos)
    for i in 1:lx, j in 1:ly, k in 1:lcos
        xint[i,j,k] = x[i]
        yint[i,j,k] = y[j]
        cosθint[i,j,k] = cosθ[k]
    end
    GGV = trilinear(xint,yint,cosθint,x0,y0,cosθ0,GGV0)
    return GGV
end
###############################################################################
# GGV corrections
###############################################################################
function kernelGGVmod!(x,y,cosθ,GGV,Kernel,p2,k2,q2)
    na, nr, ext = size(q2)
    x1 = q2
    x2 = zero(q2)
    y1 = zero(q2)
    y2 = q2
    cosθ1 = zero(q2)
    cosθ2 = zero(q2)
    for k in 1:na, j in 1:size(k2)[1], i in 1:ext
        y1[k,j,i] = k2[j,i]
        x2[k,j,i] = p2[i]
        cosθ1[k,j,i] = (p2[i] - q2[k,j,i] - k2[j,i])/(2*sqrt(p2[i]*k2[j,i]))
        cosθ2[k,j,i] = (k2[j,i] - q2[k,j,i] - p2[i])/(2*sqrt(p2[i]*q2[k,j,i]))
    end
    A1 = trilinear(x1,y1,cosθ1,x,y,cosθ,GGV)
    A2 = trilinear(x2,y2,cosθ2,x,y,cosθ,GGV)
    @. Kernel.Ghost = Kernel.Ghost * A1
    @. Kernel.GhostLoop = Kernel.GhostLoop * A2
end
# x : antighost-momentum squared
# y : gluon momentum squared
function iterate_vertex(x,y,cosθ,GGV0,QCD,param,p,G,Z)
    @info "start vertex iteration"
    #second array to compare changes between iteration steps
    GGV,newGGV = copy(GGV0),copy(GGV0)
    @unpack nr,na1,na2,IR,UV = param
    p2, Gp2, Zp2 = permute_ym(p,G,Z,QCD)
    # preallocate arrays for calculation
    Inter  = [ InterpolationPointsGGV(zeros(nr),zeros(na2,nr),zeros(na1,na2,nr),zeros(na2,nr),zeros(na1,na2,nr),zeros(na2),zeros(na2,nr),zeros(na1,na2,nr),zeros(na2,nr)) for i in 1:nthreads()]
    Kernel = [ KernelsGGV(zeros(na1,na2,nr),zeros(na1,na2,nr)) for i in 1:nthreads()]
    IntGGV = [ InterpolatedGGV(ones(na2,nr),ones(na1,na2,nr),ones(na2,nr),ones(na1,na2,nr),ones(nr),ones(na2,nr),ones(na1,na2,nr),ones(nr),ones(na2,nr),ones(na1,na2,nr)) for i in 1:nthreads()]
    Index  = [ IndexEqualGGV(zeros(na1,na2,nr),zeros(na1,na2,nr),zeros(na1,na2,nr),zeros(na2,nr),zeros(na2,nr),zeros(na2,nr),zeros(na2,nr),zeros(nr),zeros(na2),zeros(nr),zeros(na2,nr),zeros(na1,na2,nr)) for i in 1:nthreads()]
    for iter in 1:param.maxiter
        threaded_sweep!(newGGV,param,QCD,p2,Gp2,Zp2,x,y,cosθ,GGV,Inter,Kernel,IntGGV,Index)
        diff = maximum(2abs.(GGV-newGGV))
        @info "iteration $iter : max(|reldiff|)=$diff"
        rel  = diff > 0.05 ? param.rel : one(param.rel)
        if diff < param.converged
            @info "convergence"
            return @. newGGV*rel + GGV*(1-rel)
        end
        @. GGV = newGGV*rel + GGV*(1-rel)
    end
    @info "not converging"
    return GGV
end
function threaded_sweep!(newGGV,param,QCD,p2,Gp2,Zp2,x,y,cosθ,GGV,Inter,Kernel,IntGGV,Index)
    p = Progress(length(x)*length(y)*length(cosθ),1,"Vertex sweep: ",20)
    indices = [(i,j,k) for i in 1:length(x), j in 1:length(y), k in 1:length(cosθ) ]
    @threads :static for ind in 1:length(indices)
        i,j,k = indices[ind]
        id    = threadid()
        newGGV[i,j,k] = iterate_point(i,j,k,param,QCD,p2,Gp2,Zp2,x,y,cosθ,GGV,Inter[id],Kernel[id],IntGGV[id],Index[id])
        next!(p)
    end
end
function iterate_point(i,j,k,param,QCD,p2,G,Z,x,y,cosθ,A,Inter,Kernel,IntGGV,Index)
    nr,na1,na2,IR,UV = param.nr,param.na1,param.na2,param.IR,param.UV
    Quad = nodes_weights(x[i],y[j],nr,na1,na2,IR,UV)
    xᵢ,yⱼ,cosθₖ = x[i],y[j],cosθ[k]
    # load Quad values
    radW, θ₁W, θ₂W = Quad.RadWeights, Quad.θ₁Weights, Quad.θ₂Weights
    nr,na1,na2,IR,UV = param.nr,param.na1,param.na2,param.IR,param.UV
    # precalculate kernels for integration
    setkernels!(xᵢ,yⱼ,cosθₖ,Quad,Inter,Kernel)
    # interpolate propagators
    @unpack GI1, GI2, GI3, ZI1, ZI2, ZI3 = IntGGV
    @unpack indp2I1, indp2I2, indp2I3 = Index
    interpolate_ym_ggv!(Inter.I1, p2, G, Z, QCD, indp2I1, GI1, ZI1)
    interpolate_ym_ggv!(Inter.I2, p2, G, Z, QCD, indp2I2, GI2, ZI2)
    interpolate_ym_ggv!(Inter.I3, p2, G, Z, QCD, indp2I3, GI3, ZI3)
    # interpolate GhostGluonVertex
    multilinearGGV!(i,Inter,y,x,cosθ,A,IntGGV,Index)
    #perform iteration
    pdotk     = cosθₖ*sqrt(xᵢ*yⱼ)
    prefactor = QCD.CA*QCD.g2/((xᵢ*yⱼ - pdotk^2)*(2*(2*pi)^3))
    abelianang1,nabelianang1 = zeros(na1),zeros(na1)
    abelianang2,nabelianang2 = zeros(na2),zeros(na2)
    abelianrad,nabelianrad = zeros(nr),zeros(nr)
    @unpack AI1, AI2, AI3, AI4 = IntGGV
    @unpack Abelian, NonAbelian = Kernel
    @inbounds for m in 1:nr
        for n in 1:na2
            for l in 1:na1
                abelianang1[l] = θ₁W[l]*Abelian[l,n,m]*GI3[l,n,m]*AI2[l,n,m]
                nabelianang1[l] = θ₁W[l]*NonAbelian[l,n,m]*ZI3[l,n,m]*AI4[l,n,m]
            end
            abelianang2[n] = sum(abelianang1)*θ₂W[n]*GI1[m]*ZI2[n,m]*AI1[n,m]
            nabelianang2[n] = sum(nabelianang1)*θ₂W[n]*ZI1[m]*GI2[n,m]*AI3[n,m]
        end
        abelianrad[m] = radW[m]*sum(abelianang2)
        nabelianrad[m] = radW[m]*sum(nabelianang2)
    end
    A = 1 + prefactor*(sum(abelianrad) + sum(nabelianrad))
    return A
end
function nodes_weights(xᵢ,yⱼ,nr,na1,na2,IR,UV)
    cosθ₁, θ₁Weights = gausslegendre(na1)
    cosθ₂, θ₂Weights = gausschebyshev(na2,2)
    if xᵢ == yⱼ
        nr = div(nr,2)
        a, Jacobian = zeros(2nr), zeros(2nr)
        a[1:nr], Jacobian[1:nr] = gausslegendre(nr,xᵢ,UV)
        a[nr+1:end], Jacobian[nr+1:end] = gausslegendre(nr,IR,xᵢ)
    else
        nr = div(nr,3)
        a, Jacobian = zeros(3nr), zeros(3nr)
        splita, splitb = minmax(xᵢ,yⱼ)
        a[1:nr], Jacobian[1:nr] = gausslegendre(nr,splitb,UV)
        a[nr+1:2nr], Jacobian[nr+1:2nr] = gausslegendre(nr,splita,splitb)
        a[2nr+1:end], Jacobian[2nr+1:end] = gausslegendre(nr,IR,splita)
    end
    RadWeights = Jacobian .* a  #Jacobian from d⁴r
    Quad = QuadratureGGV(a,cosθ₂,cosθ₁,RadWeights, θ₁Weights, θ₂Weights)
    return Quad
end
function setkernels!(x,y,cosθ,Quad,Inter,Kernel)
    nr,na1,na2 = length(Quad.a),length(Quad.cosθ₁),length(Quad.cosθ₂)
    a, cosθ₂, cosθ₁ = Quad.a, Quad.cosθ₂, Quad.cosθ₁
    pk = cosθ*sqrt(x*y)
    xpk = x + pk
    sinθ = sqrt(1-cosθ^2)
    factor3 = x*y - pk^2
    Inter.I1 .= a
    Inter.I6 .= -cosθ₂
    @inbounds for i in 1:nr
        aᵢ = a[i]
        rootay = sqrt(y*aᵢ)
        for k in 1:na2
            pr = sqrt(x*aᵢ)*cosθ₂[k]
            aprx = muladd(-2,pr,aᵢ+x)
            div1 = 2*aᵢ*aprx*aprx
            div2 = aᵢ*aᵢ*aprx
            sinθ₂ = sqrt(1 - cosθ₂[k]^2)
            prpkpr = pr*(pk + pr)
            sum1 = aᵢ*muladd(-aᵢ,xpk,prpkpr)
            h1 = muladd(-2,xpk,pr)
            h2 = muladd(aᵢ,h1,prpkpr)
            h3 = -aᵢ*(xpk-pr)
            h4 = aᵢ*(aᵢ-pr-pk)
            h5 = pk + pr - aᵢ
            Inter.I2[k,i] = aprx
            Inter.I4[k,i] = (pr - x)/sqrt(x*aprx)
            for j in 1:na1
                kr   = rootay*muladd(sinθ*sinθ₂,cosθ₁[j],cosθ₂[k]*cosθ)
                akry = muladd(2,kr,aᵢ+y)
                sum2 = y*muladd(kr,pr,h3) + kr*muladd(kr,pr,h2)
                factor1 = prpkpr + muladd(kr,(pr - x),-aᵢ*xpk)
                factor2 = muladd(kr,pk,-pr*y)
                factor4 = muladd(kr,muladd(kr,(2pr - aᵢ), -pk*aᵢ),y*muladd(kr,pr,h4))
                Kernel.NonAbelian[j,k,i] = muladd(muladd(kr,pk,-pr*y),(sum1+sum2),factor3*factor4)/(div2*akry*akry)
                Kernel.Abelian[j,k,i] = factor1*factor2/(div1*akry)
                Inter.I3[j,k,i] = akry
                Inter.I5[j,k,i] = (h5 - kr)/sqrt(akry*aprx)
            end
        end
    end
end
function interpolate_ym_ggv!(x,p2,Gp2,Zp2,QCD,ind,Gx,Zx)
    # p2 needs to be already sorted (Gp2 and Zp2 accordingly)
    @assert issorted(p2) "error: array p2 not sorted"
    IR, UV = p2[1], p2[end]
    GIR, ZIR = Gp2[1], Zp2[1]
    GUV, ZUV = Gp2[end], Zp2[end]
    # find inidces
    findindex!(x,p2,ind)
    splines2functions2!(Gx,Zx,x,p2,Gp2,Zp2,ind)
    # Perform extrapolation
    β0 = (11QCD.CA-4*QCD.T*QCD.Nf)/3
    ω  = β0*QCD.g2*GUV^2*ZUV/(16*pi^2)
    for i in 1:length(x)
        if ind[i] == length(p2)
            Gx[i] = GUV*(ω*log(x[i]/UV)+1)^QCD.δ
            Zx[i] = ZUV*(ω*log(x[i]/UV)+1)^QCD.γ
        end
        if ind[i] == 0
            Gx[i] = GIR*(x[i]/IR)^QCD.κGh
            Zx[i] = ZIR*(x[i]/IR)^QCD.κGl
        end
    end
end
function multilinearGGV!(xpos,Inter,ks,ps,cosθs,A,IntGGV,Index)
    ks!=ps && @warn "error: ghost and gluon momenta must be the same"
    na1,na2,nr=size(Inter.I3)
    # Load precalculated interpolation points
    @unpack newI1,newI2,newI6,I1,I2,I3,I4,I5,I6 = Inter
    @unpack xind1,yind1,zind,newIndI1,newIndI6,indI1,indI2,indI4,indI6 = Index
    # Find indices in reference
    findindex!(I1,ks,indI1)
    findindex!(I2,ks,indI2)
    findindex!(I3,ps,yind1)
    findindex!(I4,cosθs,indI4)
    findindex!(I5,cosθs,zind)
    findindex!(I6,cosθs,indI6)
    # create arrays so that arguments of multilinear methods are of same size
    for k in 1:na2, i in 1:nr
        newI6[k,i] = I6[k]
        newIndI6[k,i] = indI6[k]
    end
    for j in 1:na2, i in 1:nr
        newI1[j,i] = I1[i]
        newIndI1[j,i] = indI1[i]
    end
    @inbounds for k in 1:na2, i in 1:nr
        xind = indI2[k,i]
        nI2 = I2[k,i]
        for j in 1:na1
            xind1[j,k,i] = xind
            newI2[j,k,i] = nI2
        end
    end
    trilinear!(IntGGV.AI2,newI2,I3,I5,ks,ps,cosθs,A,xind1,yind1,zind)
    trilinear!(IntGGV.AI4,I3,newI2,I5,ks,ps,cosθs,A,yind1,xind1,zind)
    bilinear!(IntGGV.AI1,I2,I4,ks,cosθs,A[:,xpos,:],indI2,indI4)
    bilinear!(IntGGV.AI3,newI1,newI6,ks,cosθs,A[:,xpos,:],newIndI1,newIndI6)
end
function permute_ym(p,G,Z,QCD)
    p2, Gp2, Zp2 = deepcopy(p), deepcopy(G), deepcopy(Z)
    # optionally extrapolate YM propagators into UV
    UVind = findmax(p2)[2]
    UV, GUV, ZUV = p2[UVind], Gp2[UVind], Zp2[UVind]
    β0 = (11QCD.CA-4*QCD.T*QCD.Nf)/3
    ω  = β0*QCD.g2*GUV^2*ZUV/(16*pi^2)
    # perform UV extrapolation once
    p2UV = exp10.(collect(range(log10(maximum(p2)),stop=15.0,length=25)))
    p2UV = p2UV[2:end]
    G2UV = @. GUV*(ω*log(p2UV/UV) + 1)^QCD.δ
    Z2UV = @. ZUV*(ω*log(p2UV/UV) + 1)^QCD.γ
    append!(p2,p2UV)
    append!(Gp2,G2UV)
    append!(Zp2,Z2UV)
    #sort YM data
    perm = sortperm(p2)
    permute!(p2,perm)
    permute!(Gp2,perm)
    permute!(Zp2,perm)
    return p2, Gp2, Zp2
end
