function integrate(p2,coeffG,coeffZ,QCD,Kernels,Quad,param,IYM,IQ,RCs)
    Z3,Z3t = RCs.Z3, RCs.Z3_tilde
    @unpack α, β = QCD
    @unpack Gz, Zz, Gy, Zy = IYM
    @unpack Az, Bz, Ay, By = IQ
    @unpack k2, q2 = Quad
    @unpack Ghost,GhostLoop,GluonLoop,QuarkLoop = Kernels
    angW,radW = Quad.AngWeights,Quad.RadWeights
    gh = ghost_selfenergy(Ghost,Zz,Gz,Zy,Gy,angW,radW)
    ghloop = ghost_loop(GhostLoop,Gz,Gy,angW,radW)
    quloop = quark_loop(QuarkLoop,q2,Gz,Az,Bz,k2,Gy,Ay,By,angW,radW)
    if param.GluonModel == :HS
        glloop = gluon_loop_HS(GluonLoop,IYM,Quad)
    elseif param.GluonModel== :CSF
        glloop = gluon_loop(GluonLoop,Zz,Gz,Zy,Gy,angW,radW,α,β)
    end
    return gh, ghloop, glloop, quloop
end
function ghost_selfenergy(kernel,Zz,Gz,Zy,Gy,wA,wR)
    na, nr, ext = size(Zz)
    ghost, radial, angular = zeros(ext), zeros(nr), zeros(na)
    @inbounds for i in 1:ext
        for j in 1:nr
            for k in 1:na
                angular[k] = wA[k] * Gz[k,j,i] * kernel[k,j,i]
            end
            radial[j] = wR[j,i]* Zy[j,i] * sum_kbn(angular)
        end
        ghost[i] = sum_kbn(radial)
    end
    return ghost
end
function ghost_loop(kernel,Gz,Gy,wA,wR)
    na, nr, ext = size(Gz)
    ghostloop, radial, angular = zeros(ext), zeros(nr), zeros(na)
    @inbounds for i in 1:ext
        for j in 1:nr
            for k in 1:na
                angular[k] = wA[k]*Gz[k,j,i]*kernel[k,j,i]
            end
            radial[j] = wR[j,i]*sum_kbn(angular)*Gy[j,i]
        end
        ghostloop[i] = sum_kbn(radial)
    end
    return ghostloop
end
function gluon_loop(kernel,Zz,Gz,Zy,Gy,wA,wR,α,β)
    na, nr, ext = size(Zz)
    gluonloop, radial, angular = zeros(ext), zeros(nr), zeros(na)
    @inbounds for i in 1:ext
        for j in 1:nr
            @turbo for k in 1:na
                Γ = (Zy[j,i]*Zz[k,j,i])^β*(Gy[j,i]*Gz[k,j,i])^α
                angular[k] = Γ*wA[k]*Zz[k,j,i]*kernel[k,j,i]
            end
            radial[j] = sum_kbn(angular)*wR[j,i]*Zy[j,i]
        end
        gluonloop[i] = sum_kbn(radial)
    end
    return gluonloop
end
function ThreeGluonModelDressing!(Quad,IYM)
    na, nr, ext = size(Quad.q2)
    @unpack p2, damping = Quad
    @unpack Gs3, Gs6, Zs3 = IYM
    @inbounds for i in 1:ext, j in 1:nr, k in 1:na
        DIR = damping[k,j,i]*Gs6[k,j,i]^3
        DUV = Gs3[k,j,i]*Zs3[k,j,i]
        IYM.VertexHS[k,j,i] = DUV*(DIR+DUV)
    end
end
f3g(s,Λ3g) = Λ3g/(s+Λ3g)
function gluon_loop_HS(kernel,IYM,Quad)
    na, nr, ext = size(Quad.q2)
    wA, wR = Quad.AngWeights, Quad.RadWeights
    @unpack Zz,Zy,VertexHS = IYM
    gluonloop, radial, angular = zeros(ext), zeros(nr), zeros(na)
    ThreeGluonModelDressing!(Quad,IYM)
    @inbounds for i in 1:ext
        for j in 1:nr
            for k in 1:na
                angular[k] = VertexHS[k,j,i]*wA[k]*Zz[k,j,i]*kernel[k,j,i]
            end
            radial[j] = sum_kbn(angular)*Zy[j,i]*wR[j,i]
        end
        gluonloop[i] = sum_kbn(radial)
    end
    return gluonloop
end
function quark_loop(KQL,z,Gz,Az,Bz,y,Gy,Ay,By,wA,wR)
    na, nr, ext = size(Gz)
    quarkloop, radial, angular = zeros(ext), zeros(nr), zeros(na)
    @inbounds for i in 1:ext
        for j in 1:nr
            @turbo for k in 1:na
                σz = Az[k,j,i]*(z[k,j,i]+(Bz[k,j,i]/Az[k,j,i])^2)
                Γ = Gz[k,j,i]*Gy[j,i]*(Ay[j,i]+Az[k,j,i])/2
                angular[k] = wA[k]*Γ*KQL[k,j,i]/σz
            end
            σy = Ay[j,i]*(y[j,i]+(By[j,i]/Ay[j,i])^2)
            radial[j] = sum_kbn(angular)*wR[j,i]/σy
        end
        quarkloop[i] = sum_kbn(radial)
    end
    return quarkloop
end
#-------------------------------------------
# Quark DSE
#-------------------------------------------
function precalc_quark(x,p2,k2,u,coeffG,coeffZ,param,QCD)
    @unpack ϵ, Λ = param
    na, nr, n = size(u)[1],size(k2)[1], length(x)
    q2, q2ind = zeros(na,nr,n),zeros(Int,(na,nr,n))
    kernel = zeros(na,nr,n)
    Gq2,Zq2 = zeros(na,nr,n),zeros(na,nr,n)
    p2ind = zeros(Int,(na,nr,n))
    @threads :static for i in 1:n
        xᵢ = Float128(x[i])
        for j in 1:nr
            y = Float128(k2[j,i])
            a = (xᵢ+y)/2
            b = (xᵢ-y)*(xᵢ-y)/2
            rtxy = 2*sqrt(xᵢ)*sqrt(y)
            for k in 1:na
                z = xᵢ+y-rtxy*u[k]
                q2[k,j,i] = Float64(z)
                kernel[k,j,i] = Float64((a+b/z-z)/xᵢ)
            end
        end
    end
    findindex!(q2,x,q2ind)
    G0 = exp.(clenshaw(p2,coeffG,ϵ,Λ))
    Z0 = exp.(clenshaw(p2,coeffZ,ϵ,Λ))
    p,G,Z = permute_ym(p2,G0,Z0,QCD)
    findindex!(q2,p,p2ind)
    sdG = second_derivative(p,G)
    sdZ = second_derivative(p,Z)
    splines2functions2sd!(Gq2,Zq2,q2,p,G,Z,p2ind,sdG,sdZ)
    return q2, q2ind, Gq2, Zq2, kernel
end
function quark_integral(x,k2,wR,wA,A,Ak2,Bk2,q2,q2ind,Gq2,Zq2,kernel,Aq2,QCD)
    na,nr,n = size(q2)
    @unpack CF, g2 = QCD
    tmpA = [zeros(na) for _ in 1:nthreads()]
    tmpB = [zeros(na) for _ in 1:nthreads()]
    radA = [zeros(nr) for _ in 1:nthreads()]
    radB = [zeros(nr) for _ in 1:nthreads()]
    ΣA,ΣB = zeros(n),zeros(n)
    sdA = second_derivative(x,A)
    splines_threaded!(Aq2,q2,x,A,q2ind,sdA)
    pref = g2 * (2*pi)^(-3) * CF / 2
    @threads :static for i in 1:n
        id = threadid()
        for j in 1:nr
            @turbo for k in 1:na
                Γ = Gq2[k,j,i]*Gq2[k,j,i]*Zq2[k,j,i]*(Ak2[j,i] + Aq2[k,j,i])
                tmpB[id][k] = wA[k]*Γ*(k2[j,i]/q2[k,j,i])
                tmpA[id][k] = wA[k]*Γ*(k2[j,i]/q2[k,j,i])*kernel[k,j,i]
            end
            σw = wR[j,i] / (k2[j,i]*Ak2[j,i]^2 + Bk2[j,i]^2)
            radA[id][j] = sum_kbn(tmpA[id])*  Ak2[j,i] * σw
            radB[id][j] = sum_kbn(tmpB[id])* 3Bk2[j,i] * σw
        end
        ΣA[i] = sum_kbn(radA[id])*pref
        ΣB[i] = sum_kbn(radB[id])*pref
    end
    return ΣA,ΣB
end
#-------------------------------------------
# Quark DSE (modified momentum routing)
#-------------------------------------------
function precalc_quark_v2(x,p2,k2,u,coeffG,coeffZ,param,QCD)
    @unpack ϵ, Λ = param
    na, nr, n = size(u)[1],size(k2)[1], length(x)
    q2 = zeros(eltype(x),(na,nr,n))
    Gk2,Zk2 = similar(k2),similar(k2)
    KA = similar(q2)
    @threads :static for i in 1:n
        for j in 1:nr
            rtxy = 2*sqrt(x[i])*sqrt(k2[j,i])
            for k in 1:na
                q2[k,j,i] = x[i]+k2[j,i]-rtxy*u[k]
                KA[k,j,i] = ((x[i]+q2[k,j,i])/2+(x[i]-q2[k,j,i])*(x[i]-q2[k,j,i])/2/k2[j,i]-k2[j,i])/x[i]
            end
        end
    end
    G = exp.(clenshaw(p2,coeffG,ϵ,Λ))
    Z = exp.(clenshaw(p2,coeffZ,ϵ,Λ))
    Gk2 = exp.(clenshaw(k2,coeffG,ϵ,Λ))
    Zk2 = exp.(clenshaw(k2,coeffZ,ϵ,Λ))
    extrapolate_ym!(Gk2,Zk2,k2,p2,G,Z,QCD)
    return q2, KA, Gk2, Zk2
end
function quark_integral_v2(k2,wR,wA,Ak2,q2,Gk2,Zk2,Aq2,Bq2,KA,QCD)
    na,nr,n = size(q2)
    ΣA,ΣB = zeros(n),zeros(n)
    tmpA = [zeros(na) for _ in 1:nthreads()]
    tmpB = [zeros(na) for _ in 1:nthreads()]
    radA = [zeros(nr) for _ in 1:nthreads()]
    radB = [zeros(nr) for _ in 1:nthreads()]
    pref = QCD.g2 * (2*pi)^(-3) * QCD.CF
    @threads :static for i in 1:n
        id = threadid()
        for j in 1:nr
            # radial weight and jacobian factor
            R = wR[j,i]*k2[j,i]
            # structure below only depends on contractions
            Γ = Gk2[j,i]*Gk2[j,i]*Zk2[j,i]*Ak2[j,i]/k2[j,i]
            @turbo for k in 1:na
                σ = inv(q2[k,j,i]*Aq2[k,j,i]^2 + Bq2[k,j,i]^2)
                S = R*wA[k]*Γ*σ
                tmpA[id][k] = S*KA[k,j,i]*Aq2[k,j,i]
                tmpB[id][k] = S*3Bq2[k,j,i]
            end
            radA[id][j] = sum_kbn(tmpA[id])
            radB[id][j] = sum_kbn(tmpB[id])
        end
        ΣA[i] = sum_kbn(radA[id])*pref
        ΣB[i] = sum_kbn(radB[id])*pref
    end
    return ΣA,ΣB
end
