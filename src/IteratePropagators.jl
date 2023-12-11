function diff_w(f,f0)
    sum1 = zero(eltype(f))
    sum2 = zero(eltype(f))
    for i in eachindex(f,f0)
        sum1 += abs(f[i]+f0[i])
        sum2 += abs(f[i]-f0[i])
    end
    return (sum1 < 10^-10) ? sum2 : min(2sum2/sum1,sum2)
end
function iterate_propagators(x,A0,B0,p2,coeffG0,coeffZ0,QCD,param,Quad,Kernels;savetmp=true)
    @unpack ϵ, Λ, propiter, quarkrelax = param
    G0 = exp.(clenshaw(p2,coeffG0,ϵ,Λ))
    Z0 = exp.(clenshaw(p2,coeffZ0,ϵ,Λ))
    Z2 = one(eltype(x))
    Zm = one(eltype(x))
    Z3,Z3_tilde = one(eltype(x)),one(eltype(x))
    RCs = RenormalizationConstants(Z2,Zm,Z3,Z3_tilde)
    diff = +Inf
    success = false
    for iter in 1:propiter
        quark_first = true
        GC.gc()
        if quark_first
            A,B = iterate_quark(x,p2,A0,B0,coeffG0,coeffZ0,QCD,param,RCs)
            t = diff > 10param.ϵtotal ? quarkrelax : 1.0
            @. A = A*t+A0*(1-t)
            @. B = B*t+B0*(1-t)
            GC.gc()
            cG,cZ,success = newton_ym(x,A,B,p2,coeffG0,coeffZ0,RCs,QCD,param,Quad,Kernels)
        else
            cG,cZ,success = newton_ym(x,A0,B0,p2,coeffG0,coeffZ0,RCs,QCD,param,Quad,Kernels)
            GC.gc()
            A,B = iterate_quark(x,p2,A0,B0,cG,cZ,QCD,param,RCs)
            t = diff > 10param.ϵtotal ? quarkrelax : 1.0
            @. A = A*t+A0*(1-t)
            @. B = B*t+B0*(1-t)
        end
        #t = diff > 10param.ϵtotal ? quarkrelax : 1.0
        #@. cG = cG*t+coeffG0*(1-t)
        #@. cZ = cZ*t+coeffZ0*(1-t)
        G = exp.(clenshaw(p2,cG,ϵ,Λ))
        Z = exp.(clenshaw(p2,cZ,ϵ,Λ))
        diff   = diff_w(A,A0) + diff_w(B,B0) + diff_w(G,G0) + diff_w(Z,Z0)
        diffYM = diff_w(G,G0) + diff_w(Z,Z0)
        diffQ  = diff_w(A,A0) + diff_w(B,B0)
        coeffG0 .= cG
        coeffZ0 .= cZ
        G0 .= G
        Z0 .= Z
        A0 .= A
        B0 .= B
        @unpack Z2,Zm,Z3,Z3_tilde  = RCs
        if savetmp
            savedata("tmp",x,A,B,p2,cG,cZ,Z2,Zm,Z3,Z3_tilde,QCD,param)
        end
        @info "cycle $iter :  diff = $diff (YM: $diffYM) (Q: $diffQ)"
        @info "M(0) = $(B[1]/A[1])"
        if diff<param.ϵtotal && iter > 2
            return A,B,cG,cZ,RCs,success
        end
    end
    @info "no convergence"
    success = false
    return A0,B0,coeffG0,coeffZ0,RCs,success
end
#####################################################
# quark propagator iteration
#####################################################
function kernels_quark(x,param)
    n = length(x)
    @unpack IR, UV = param
    na,nr = param.naQ, div(param.nrQ,2)
    u,wA = gausschebyshev(na,2)
    k2,wR = zeros(2nr,n),zeros(2nr,n)
    for i in 1:n
        k2[1:nr,i],wR[1:nr,i] = gausslegendre(nr,IR,x[i])
        k2[nr+1:end,i],wR[nr+1:end,i] = gausslegendre(nr,x[i],10UV)
    end
    return k2,u,wA,wR
end
function iterate_quark(x,p2,A0,B0,coeffG,coeffZ,QCD,param,RCs;renormalize=true)
    A,B = copy(A0), copy(B0)
    RenA, RenB = similar(x), similar(x)
    k2,u,wA,wR = kernels_quark(x,param)
    q2,KA,Gk2,Zk2 = precalc_quark_v2(x,p2,k2,u,coeffG,coeffZ,param,QCD)
    prog = ProgressThresh(param.ϵQ,1,"quark iteration:")
    # findindex
    k2ind = findindex(k2,x)
    q2ind = findindex(q2,x)
    Ak2 = zero(k2)
    Bk2 = zero(k2)
    Aq2 = zero(q2)
    Bq2 = zero(q2)
    for iter in 1:param.maxiterQ
        # interpolate quark dressing
        # use splines
        if param.splinesQ
            splines2functions!(Ak2,Bk2,k2,x,A,B,k2ind)
            splines2functions!(Aq2,Bq2,q2,x,A,B,q2ind)
        else
            # Alternateively use linear interpolation
            monolinear_thread!(Ak2,k2,x,A,k2ind)
            monolinear_thread!(Aq2,q2,x,A,q2ind)
            monolinear_thread!(Bq2,q2,x,B,q2ind)
        end
        # perform quark integral
        ΣA,ΣB = quark_integral_v2(k2,wR,wA,Ak2,q2,Gk2,Zk2,Aq2,Bq2,KA,QCD)
        # get renormalization constant 
        Z2 = RCs.Z2
        if renormalize
            s = findmin(abs.(x.-param.µ))[2]
            Z2 = 1/(1 + ΣA[s])
            ZmZ2m = QCD.mQ + Z2*ΣB[s]
            iszero(QCD.mQ) ? Zm=1.0 : Zm = ZmZ2m/(Z2*QCD.mQ)
            RCs.Z2 = Z2
            RCs.Zm = Zm
        end
        # update quark dressing functions
        @. RenA = Z2*(1 + ΣA)
        @. RenB = QCD.mQ + Z2*(ΣB-ΣB[s])
        diff = diff_w(A,RenA) + diff_w(B,RenB)
        update!(prog, diff)
        if diff < param.ϵQ
            finish!(prog)
            return RenA,RenB
        end
        @. A = RenA
        @. B = RenB
    end
    return A,B
end
function iterate_quark_complex(x,p2,A0,B0,coeffG,coeffZ,Ak2,QCD,param,RCs)
    A,B = copy(A0), copy(B0)
    RenA, RenB = similar(x), similar(x)
    k2,u,wA,wR = kernels_quark(x,param)
    q2,KA,Gk2,Zk2 = precalc_quark_v2(x,p2,k2,u,coeffG,coeffZ,param,QCD)
    prog = ProgressThresh(param.ϵQ,1,"quark iteration:")
    for iter in 1:param.maxiterQ
        # interpolate quark dressing
        # Ak2 = monolinear(k2,x,A)
        Aq2 = monolinear(q2,x,A)
        Bq2 = monolinear(q2,x,B)
        # perform quark integral
        ΣA,ΣB = quark_integral_v2(k2,wR,wA,Ak2,q2,Gk2,Zk2,Aq2,Bq2,KA,QCD)
        # get renormalization constant 
        Z2 = RCs.Z2
        Zm = RCs.Zm
        # update quark dressing functions
        @. RenA = Z2*(1 + ΣA)
        @. RenB = QCD.mQ + Z2*ΣB-(Zm*Z2-1)*QCD.mQ
        diff = diff_w(A,RenA) + diff_w(B,RenB)
        update!(prog, diff)
        if diff < param.ϵQ
            finish!(prog)
            return RenA,RenB
        end
        @. A = RenA
        @. B = RenB
    end
    return A,B
end
#####################################################
# Yang-Mills propagator iteration
#####################################################
function quad_kernels(p2,param,QCD)
    @unpack na, nr, IR, UV, ϵ, Λ = param
    k2,q2,rad,ang = quadrature(p2,nr,na,IR,UV)
    k2_inv = log_mapping_inv(k2,ϵ,Λ)
    q2_inv = log_mapping_inv(q2,ϵ,Λ)
    if param.GluonModel == :HS
        @unpack hIR, Λ3G = param
        ext = length(p2)
        s6, s3, f = zeros(na,nr,ext), zeros(na,nr,ext), zeros(na,nr,ext)
        x,y,z = p2, k2, q2
        @inbounds for i in 1:ext, j in 1:nr, k in 1:na
            s6[k,j,i] = x[i]+y[j,i]+z[k,j,i]
            s3[k,j,i] = s6[k,j,i]/2
            f[k,j,i] = hIR*(f3g(x[i],Λ3G)*f3g(y[j,i],Λ3G)*f3g(z[k,j,i],Λ3G))^4
        end
        s3_inv = log_mapping_inv(s3,ϵ,Λ)
        s6_inv = log_mapping_inv(s6,ϵ,Λ)
    else
        # otherwise the objects are still around but allocate minimally
        s3,s3_inv = zeros(1,1,1),zeros(1,1,1)
        s6,s6_inv = zeros(1,1,1),zeros(1,1,1)
        f = zeros(1,1,1)
    end
    Quad = QuadratureProp(p2,k2,q2,rad,ang,k2_inv,q2_inv,s3,s6,s3_inv,s6_inv,f)
    Ghost,GhostLoop,GluonLoop,QuarkLoop = kernel_momenta(p2,k2,q2,QCD)
    Kernels = KernelsProp(Ghost,GhostLoop,GluonLoop,QuarkLoop)
    return Quad, Kernels
end
function quadrature(x,nr,na,IR,UV)
    ext = length(x)
    k2,wR = zeros(nr,ext), zeros(nr,ext)
    u, wA = gausschebyshev(na,2)
    q2 = zeros(na,nr,ext)
    for i in 1:ext
        k2[:,i], wR[:,i] = gausslegendre(nr,IR,UV)
    end
    for l in 1:ext, j in 1:nr
        x0,y0 = Float128(x[l]), Float128(k2[j,l])
        rtxy0 = 2sqrt(x0*y0)
        for i in 1:na
            q2[i,j,l] = Float64(x0 + y0 + u[i]*rtxy0)
        end
    end
    return k2, q2, wR, wA
end
function kernel_momenta(p2,k2,q2,QCD)
    @unpack CA,T,Nf,g2 = QCD
    na, nr, ext = size(q2)
    GhostSE = zeros(na,nr,ext)
    GhostLoop= zeros(na,nr,ext)
    GluonLoop = zeros(na,nr,ext)
    QuarkLoop= zeros(na,nr,ext)
    f  = g2 * (2*pi)^(-3) * CA
    f2 = g2 * (2*pi)^(-3) * T
    @inbounds for i in 1:ext
        @threads :static for j in 1:nr
            for k in 1:na
                x,y,z = Float128(p2[i]), Float128(k2[j,i]), Float128(q2[k,j,i])
                GhostSE[k,j,i] = Float64((x^2+(y-z)^2-2x*(y+z))/(4x*y*z))
                GhostLoop[k,j,i] = Float64(-(x^2+(y-z)^2-2x*(y+z))/(12x^2*z))
                GluonLoop[k,j,i] = Float64((x^4+8x^3*(y+z)+x^2*(-18y^2-32y*z-
                    18z^2)+(y-z)^2*(y^2+10y*z+z^2)+2x*(y+z)*(4y^2-20y*z+4z^2))/
                    (24x^2*y*z^2)+(5/4+2*T*Nf/CA)/x)
                QuarkLoop[k,j,i] = Float64(2y*(2-(y+z)/x-((y-z)/x)^2)/(3))
            end
        end
    end
    return f*GhostSE, f*GhostLoop, f*GluonLoop, f2*QuarkLoop
end
nr_error(nodesG,nodesZ) = sqrt(sum(nodesG.^2)) + sqrt(sum(nodesZ.^2))
function newton_ym(x,A,B,p2,coeffG,coeffZ,RCs,QCD,param,Quad,Kernels)
    maxiter = param.maxiterNR
    println()
    @info "start propagator iteration:"
    interYM = [interpolate_ym(p2,coeffG,coeffZ,QCD,Quad,param) for i in 1:nthreads()+1]
    interQ = interpolate_quark(x,A,B,Quad,QCD)
    nodesG,nodesZ = residual(p2,RCs,coeffG,coeffZ,QCD,param,Kernels,Quad,interYM[1],interQ)
    error = nr_error(nodesG,nodesZ)
    @info string("initial error:",error)
    success = false
    for j in 1:maxiter
        t0 = time()
        #J = jacobian_ym(p2,coeffG,coeffZ,nodesG,nodesZ,RCs,QCD,param,Kernels,Quad,interYM,interQ)
        #J = jacobian_ym_central(p2,coeffG,coeffZ,RCs,QCD,param,Kernels,Quad,interYM,interQ)
        J = jacobian_ym_five_point(p2,coeffG,coeffZ,RCs,QCD,param,Kernels,Quad,interYM,interQ,param.h)
        found, error, cG, cZ, nG, nZ = linesearch_ym(p2,coeffG,coeffZ,nodesG,nodesZ,error,RCs,J,param,Kernels,QCD,Quad,interYM[1],interQ,t0,j)
        coeffG .= cG
        coeffZ .= cZ
        nodesG .= nG
        nodesZ .= nZ
        if error < param.epsilonNR
            @info "Propagators : convergence!"
            success = true
            return coeffG, coeffZ, success
        end
        if !found
            @info "Error: Line search not succesfull"
            success = false
            return coeffG, coeffZ, success
        end
    end
    success = false
    @info "Propagators : No convergence!"
    return coeffG, coeffZ, success
end
scalecoeff(coeff;s=1) = max(abs(coeff)/s,1)*sign(coeff)
function jacobian_ym(p2,coeffG,coeffZ,nodesG,nodesZ,RCs,QCD,param,Kernels,Quad,interYM,interQ,h)
    @unpack nCheby = param
    @unpack k2, q2, k2_inv, q2_inv = Quad
    @unpack Gy, Zy, Gz, Zz = interYM[1]
    interpolate_ym!(k2,k2_inv,Gy,Zy,p2,coeffG,coeffZ,QCD,param)
    interpolate_ym!(q2,q2_inv,Gz,Zz,p2,coeffG,coeffZ,QCD,param)
    J = zeros(2nCheby,2nCheby)
    @threads :static for i in 1:nCheby
        id = threadid()+1
        coeffGmod, coeffZmod = deepcopy(coeffG), deepcopy(coeffZ)
        scaleG = scalecoeff(coeffG[i])
        scaleZ = scalecoeff(coeffZ[i])
        hG = scaleG*h
        hZ = scaleZ*h
        coeffGmod[i] += hG
        coeffZmod[i] += hZ
        copy_interYM!(interYM[id], interYM[1])
        nodesG_modG, nodesZ_modG = residual(p2,deepcopy(RCs),coeffGmod,coeffZ,QCD,param,Kernels,Quad,interYM[id],interQ;skipZ=true)
        copy_interYM!(interYM[id], interYM[1])
        nodesG_modZ, nodesZ_modZ = residual(p2,deepcopy(RCs),coeffG,coeffZmod,QCD,param,Kernels,Quad,interYM[id],interQ;skipG=true)
        for j in 1:nCheby
            J[j,i] = ( nodesG_modG[j] - nodesG[j] )/hG
            J[j,i+nCheby] = ( nodesG_modZ[j] - nodesG[j] )/hZ
            J[j+nCheby,i] = ( nodesZ_modG[j] - nodesZ[j] )/hG
            J[j+nCheby,i+nCheby] = ( nodesZ_modZ[j] - nodesZ[j] )/hZ
        end
    end
    return J
end
function copy_interYM!(A, B)
    A.Gy .= B.Gy
    A.Zy .= B.Zy
    A.Gz .= B.Gz
    A.Zz .= B.Zz
end
function jacobian_ym_central(p2,coeffG,coeffZ,RCs,QCD,param,Kernels,Quad,interYM,interQ,h)
    @unpack nCheby = param
    J = zeros(2nCheby,2nCheby)
    @unpack k2, q2, k2_inv, q2_inv = Quad
    @unpack Gy, Zy, Gz, Zz = interYM[1]
    interpolate_ym!(k2,k2_inv,Gy,Zy,p2,coeffG,coeffZ,QCD,param)
    interpolate_ym!(q2,q2_inv,Gz,Zz,p2,coeffG,coeffZ,QCD,param)
    @threads :static for i in 1:nCheby
        id = threadid()+1
        # vary coefficients and calculate modified g and z functions
        scaleG = scalecoeff(coeffG[i])
        scaleZ = scalecoeff(coeffZ[i])
        hG = scaleG*h
        hZ = scaleZ*h
        coeffGmod, coeffZmod = deepcopy(coeffG), deepcopy(coeffZ)
        coeffGmod[i] += hG
        coeffZmod[i] += hZ
        coeffGmod2, coeffZmod2 = deepcopy(coeffG), deepcopy(coeffZ)
        coeffGmod2[i] -= hG
        coeffZmod2[i] -= hZ
        copy_interYM!(interYM[id], interYM[1])
        nodesG_modG, nodesZ_modG = residual(p2,deepcopy(RCs),coeffGmod,coeffZ,QCD,param,Kernels,Quad,interYM[id],interQ,skipZ=true)
        copy_interYM!(interYM[id], interYM[1])
        nodesG_modZ, nodesZ_modZ = residual(p2,deepcopy(RCs),coeffG,coeffZmod,QCD,param,Kernels,Quad,interYM[id],interQ,skipG=true)
        copy_interYM!(interYM[id], interYM[1])
        nodesG_modG2, nodesZ_modG2 = residual(p2,deepcopy(RCs),coeffGmod2,coeffZ,QCD,param,Kernels,Quad,interYM[id],interQ,skipZ=true)
        copy_interYM!(interYM[id], interYM[1])
        nodesG_modZ2, nodesZ_modZ2 = residual(p2,deepcopy(RCs),coeffG,coeffZmod2,QCD,param,Kernels,Quad,interYM[id],interQ,skipG=true)
        for j in 1:nCheby
            J[j,i] = ( nodesG_modG[j] - nodesG_modG2[j] )/(2hG)
            J[j,i+nCheby] = ( nodesG_modZ[j] - nodesG_modZ2[j] )/(2hZ)
            J[j+nCheby,i] = ( nodesZ_modG[j] - nodesZ_modG2[j] )/(2hG)
            J[j+nCheby,i+nCheby]=( nodesZ_modZ[j] - nodesZ_modZ2[j] )/(2hZ)
        end
    end
    return J
end
function jacobian_ym_five_point(p2,coeffG,coeffZ,RCs,QCD,param,Kernels,Quad,interYM,interQ,h)
    @unpack nCheby = param
    J = zeros(2nCheby,2nCheby)
    @unpack k2, q2, k2_inv, q2_inv = Quad
    @unpack Gy, Zy, Gz, Zz = interYM[1]
    interpolate_ym!(k2,k2_inv,Gy,Zy,p2,coeffG,coeffZ,QCD,param)
    interpolate_ym!(q2,q2_inv,Gz,Zz,p2,coeffG,coeffZ,QCD,param)
    #p = Progress(nCheby,5,"Calculating Jacobian",25)
    @threads :static for i in 1:nCheby
        id = threadid()+1
        # vary coefficients and calculate modified g and z functions
        scaleG = scalecoeff(coeffG[i])
        scaleZ = scalecoeff(coeffZ[i])
        hG = scaleG*h
        hZ = scaleZ*h
        coeffGp, coeffZp = deepcopy(coeffG), deepcopy(coeffZ)
        coeffGm, coeffZm = deepcopy(coeffG), deepcopy(coeffZ)
        coeffG2p, coeffZ2p = deepcopy(coeffG), deepcopy(coeffZ)
        coeffG2m, coeffZ2m = deepcopy(coeffG), deepcopy(coeffZ)
        coeffGp[i] += hG
        coeffGm[i] -= hG
        coeffZp[i] += hZ
        coeffZm[i] -= hZ
        coeffG2p[i] += 2hG
        coeffG2m[i] -= 2hG
        coeffZ2p[i] += 2hZ
        coeffZ2m[i] -= 2hZ
        copy_interYM!(interYM[id], interYM[1])
        nodesG_Gp, nodesZ_Gp = residual(p2,deepcopy(RCs),coeffGp,coeffZ,QCD,param,Kernels,Quad,interYM[id],interQ,skipZ=true)
        copy_interYM!(interYM[id], interYM[1])
        nodesG_Gm, nodesZ_Gm = residual(p2,deepcopy(RCs),coeffGm,coeffZ,QCD,param,Kernels,Quad,interYM[id],interQ,skipZ=true)
        copy_interYM!(interYM[id], interYM[1])
        nodesG_Zp, nodesZ_Zp = residual(p2,deepcopy(RCs),coeffG,coeffZp,QCD,param,Kernels,Quad,interYM[id],interQ,skipG=true)
        copy_interYM!(interYM[id], interYM[1])
        nodesG_Zm, nodesZ_Zm = residual(p2,deepcopy(RCs),coeffG,coeffZm,QCD,param,Kernels,Quad,interYM[id],interQ,skipG=true)

        copy_interYM!(interYM[id], interYM[1])
        nodesG_G2p, nodesZ_G2p = residual(p2,deepcopy(RCs),coeffG2p,coeffZ,QCD,param,Kernels,Quad,interYM[id],interQ,skipZ=true)
        copy_interYM!(interYM[id], interYM[1])
        nodesG_G2m, nodesZ_G2m = residual(p2,deepcopy(RCs),coeffG2m,coeffZ,QCD,param,Kernels,Quad,interYM[id],interQ,skipZ=true)
        copy_interYM!(interYM[id], interYM[1])
        nodesG_Z2p, nodesZ_Z2p = residual(p2,deepcopy(RCs),coeffG,coeffZ2p,QCD,param,Kernels,Quad,interYM[id],interQ,skipG=true)
        copy_interYM!(interYM[id], interYM[1])
        nodesG_Z2m, nodesZ_Z2m = residual(p2,deepcopy(RCs),coeffG,coeffZ2m,QCD,param,Kernels,Quad,interYM[id],interQ,skipG=true)
        for j in 1:nCheby
            J[j,i] =                sum_kbn((-nodesG_G2p[j],+8nodesG_Gp[j],-8nodesG_Gm[j],+nodesG_G2m[j]))/(12hG)
            J[j,i+nCheby] =         sum_kbn((-nodesG_Z2p[j],+8nodesG_Zp[j],-8nodesG_Zm[j],+nodesG_Z2m[j]))/(12hZ)
            J[j+nCheby,i] =         sum_kbn((-nodesZ_G2p[j],+8nodesZ_Gp[j],-8nodesZ_Gm[j],+nodesZ_G2m[j]))/(12hG)
            J[j+nCheby,i+nCheby] =  sum_kbn((-nodesZ_Z2p[j],+8nodesZ_Zp[j],-8nodesZ_Zm[j],+nodesZ_Z2m[j]))/(12hZ)
        end
    end
    return J
end
function residual(p2,RCs,coeffG,coeffZ,QCD,param,Kernels,Quad,interYM,interQ;kws...)
    @unpack ϵ, Λ = param
    @unpack Gz, Gy, Zz, Zy = interYM
    @unpack k2, q2, k2_inv, q2_inv = Quad
    G = exp.(clenshaw(p2,coeffG,ϵ,Λ))
    Z = exp.(clenshaw(p2,coeffZ,ϵ,Λ))
    interpolate_ym!(k2,k2_inv,Gy,Zy,p2,coeffG,coeffZ,QCD,param;kws...)
    interpolate_ym!(q2,q2_inv,Gz,Zz,p2,coeffG,coeffZ,QCD,param;kws...)
    if param.GluonModel == :HS
        @unpack Gs3, Gs6, Zs3 = interYM
        @unpack s3, s6, s3_inv, s6_inv = Quad
        dummy = ones(1,1,1)
        interpolate_ym!(s3,s3_inv,Gs3,Zs3,p2,coeffG,coeffZ,QCD,param,powG=QCD.α,powZ=QCD.β)
        interpolate_ym!(s6,s6_inv,Gs6,dummy,p2,coeffG,coeffZ,QCD,param,skipZ=true)
    end
    renG,renZ = renormalize(p2,coeffG,coeffZ,RCs,QCD,Kernels,Quad,param,interYM,interQ)
    nodesG, nodesZ = G - renG, Z - renZ
    return nodesG[1:end-1], nodesZ[1:end-1]
end
function renormalize(p2,coeffG,coeffZ,RCs,QCD,Kernels,Quad,param,interYM,interQ)
    n = length(p2)
    s = findmin(p2)[2]
    ghost, ghostloop, gluonloop, quarkloop = integrate(p2,coeffG,coeffZ,QCD,Kernels,Quad,param,interYM,interQ,RCs)
    gluon = @. gluonloop + ghostloop +  QCD.Nf*RCs.Z2*quarkloop
    renG  = @. ((1.0/param.GIR - ghost[s]) + ghost)^(-1)
    renZ  = @. ((1.0/param.ZUV - gluon[n]) + gluon)^(-1)
    Z̃3 = 1.0/param.GIR - ghost[s]
    Z3 = 1.0/param.ZUV - gluon[n]
    RCs.Z3,RCs.Z3_tilde = Z3, Z̃3
    return renG, renZ
end
#####################################################
# interpolation
#####################################################
function interpolate_ym(p2,coeffG,coeffZ,QCD,Quad,param)
    @unpack k2, q2, k2_inv, q2_inv = Quad
    Gk2,Zk2,Gq2,Zq2 = zero(k2),zero(k2),zero(q2),zero(q2)
    if param.GluonModel == :HS
        Gs3,Zs3,Gs6,vertex = zero(q2),zero(q2),zero(q2),zero(q2)
    else
        Gs3,Zs3,Gs6,vertex = zeros(1,1,1),zeros(1,1,1),zeros(1,1,1),zeros(1,1,1)
    end
    interpolate_ym!(k2,k2_inv,Gk2,Zk2,p2,coeffG,coeffZ,QCD,param)
    interpolate_ym!(q2,q2_inv,Gq2,Zq2,p2,coeffG,coeffZ,QCD,param)
    if param.GluonModel == :HS
        dummy = ones(1,1,1)
        @unpack s3,s3_inv,s6,s6_inv = Quad
        interpolate_ym!(s3,s3_inv,Gs3,Zs3,p2,coeffG,coeffZ,QCD,param,powG=QCD.α,powZ=QCD.β)
        interpolate_ym!(s6,s6_inv,Gs6,dummy,p2,coeffG,coeffZ,QCD,param,skipZ=true)
    end
    interYM = InterpolatedYM(Gk2, Zk2, Gq2, Zq2, Gs3, Zs3, Gs6, vertex)
    return interYM
end
function interpolate_quark(x,A,B,Quad,QCD)
    k2,q2 = Quad.k2,Quad.q2
    Ak2 = splines(k2,x,A)
    Bk2 = splines(k2,x,B)
    Aq2 = splines(q2,x,A)
    Bq2 = splines(q2,x,B)
    I = InterpolatedQuark(Ak2, Bk2, Aq2, Bq2)
end
function interpolate_ym!(y,Gy,Zy,x,coeffG,coeffZ,QCD,param;kws...)
    y_inv = log_mapping_inv(y,param.ϵ,param.Λ)
    interpolate_ym!(y,y_inv,Gy,Zy,x,coeffG,coeffZ,QCD,param;kws...)
end
function interpolate_ym!(y,y_inv,Gy,Zy,x,coeffG,coeffZ,QCD,param;powG=1,powZ=1,skipG=false,skipZ=false)
    x̃ = collect(extrema(x))
    G = exp.(clenshaw(x̃,powG*coeffG,param.ϵ,param.Λ))
    Z = exp.(clenshaw(x̃,powZ*coeffZ,param.ϵ,param.Λ))
    interpolate_ym!(y,y_inv,Gy,Zy,x,coeffG,coeffZ,x̃,G,Z,QCD,param;powG=powG,powZ=powZ,skipG=skipG,skipZ=skipZ)
end
function interpolate_ym!(y,y_inv,Gy,Zy,x,coeffG,coeffZ,x̃,G,Z,QCD,param;powG=1,powZ=1,skipG=false,skipZ=false)
    skipG || clenshaw!(Gy,y_inv,powG*coeffG)
    skipZ || clenshaw!(Zy,y_inv,powZ*coeffZ)
    skipG || (@turbo @. Gy = exp(Gy))
    skipZ || (@turbo @. Zy = exp(Zy))
    extrapolate_ym!(Gy,Zy,y,x̃,G,Z,QCD;powG=powG,powZ=powZ,skipG=skipG,skipZ=skipZ)
end
function extrapolate_ym!(Gy,Zy,y,x,G,Z,QCD;powG=1,powZ=1,skipG=false,skipZ=false)
    @unpack δ, γ, κGh, κGl, Nc, Nf, g2, CA, T = QCD
    IRind, UVind  = findmin(x)[2],findmax(x)[2]
    GUV, ZUV, xUV = G[UVind], Z[UVind], x[UVind]
    GIR, ZIR, xIR = G[IRind], Z[IRind], x[IRind]
    β0 = (11CA-4*T*Nf)/3
    ω  = β0*g2*GUV^(2/powG)*ZUV^(1/powZ)/(16*pi^2)
    @inbounds for i in eachindex(y)
        if y[i] > xUV
            base  = ω*log(y[i]/xUV) + one(xUV)
            skipG || (Gy[i] = GUV*base^(δ*powG))
            skipZ || (Zy[i] = ZUV*base^(γ*powZ))
        end
    end
    @inbounds for i in eachindex(y)
        if y[i] < xIR
            skipG || (Gy[i] = GIR*(y[i]/xIR)^(κGh*powG))
            skipZ || (Zy[i] = ZIR*(y[i]/xIR)^(κGl*powZ))
        end
    end
end
###############################################################################
# line search
###############################################################################
function linesearch_ym(p2,coeffG,coeffZ,nodesG,nodesZ,error,RCs,jacobian,param,Kernels,QCD,Quad,InterYM,InterQ,t0,iter)
    jinvv = jacobian\vcat(nodesG,nodesZ)
    coeff = vcat(coeffG,coeffZ)
    ϕ0, ϕ1, ϕ2 = 0.0, 0.0, 0.0 # store these values for more involved
    λ1, λ2, λ3 = 1.0 , 0.0, 0.0 # line search algorithms
    for i in 1:param.linesearch
        if i > 2 && !isnan(ϕ0) && !isinf(ϕ0)
            λ = parabola_search(ϕ0, ϕ1, ϕ2, λ1, λ2, λ3, 2, 10.0)
        else
            λ = (1/2)^(i-1)
        end
        ncoeff = coeff - λ*jinvv
        cG, cZ = ncoeff[1:param.nCheby],ncoeff[param.nCheby+1:2*param.nCheby]
        nG, nZ = residual(p2,RCs,cG,cZ,QCD,param,Kernels,Quad,InterYM,InterQ)
        nerr = nr_error(nG,nZ)
        # update parabola
        (i == 1) && (ϕ0 = nerr)
        ϕ1,ϕ2 = ϕ2,nerr
        λ2,λ3 = λ3,λ
        # end of parabola search
        α = 10^-4
        if nerr < error*(1.0 - α) # sufficient decrease
            # TODO: check if even smaller λ leads to smaller error
            @info string("line search #",lpad(i,2)," : error: ",rpad(nerr,26),"\t (elapsed time: ",round(time()-t0,digits=2)," s) \t ($(lpad(iter,4))/$(param.maxiterNR))")
            return true, nerr, cG, cZ, nG, nZ
        end
    end
    return false, error, coeffG, coeffZ, nodesG, nodesZ
end
# for now: assume that λ0 = 1
function parabola_search(ϕ0, ϕ2, ϕ3, λ0, λ2, λ3, a, b)
    curvature = (λ2*(ϕ3-ϕ0) - λ3*(ϕ2-ϕ0))/(λ3*λ2*(λ3-λ2))
    if curvature > 0
        mp0  = ((λ3/λ2)*(ϕ2-ϕ0) + (λ2/λ3)*(ϕ0-ϕ3))/(λ2 - λ3)
        λmin = mp0/curvature
        if λmin < λ3/b
            return  λ3/b
        elseif λmin > λ3/a
            return λ3/a
        else
            return λmin
        end
    else
        return 2*λ3/(a + b)
    end
end
