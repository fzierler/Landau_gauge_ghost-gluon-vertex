function Deriv(x::Array{Float64,1},y::Array{Float64,1})
    n = length(x)
    perm = sortperm(x)
    x,y = permute!(x,perm), permute!(y,perm)
    yprime = zeros(n)
    for i in 1:n
        if i == 1
            yprime[i] = (y[i+1] - y[i])/(x[i+1] - x[i])
        elseif i == n
            yprime[i] = (y[i] - y[i-1])/(x[i] - x[i-1])
        else
            yprime[i] = (y[i+1] - y[i-1])/(x[i+1] - x[i-1])
        end
    end
    return yprime
end
function pion_decay(p20::Array{Float64,1},A0::Array{Float64,1},B0::Array{Float64,1},Z2::Real,Nc::Real)
    p2, A, B = copy(p20), copy(A0), copy(B0)
    perm = sortperm(p2)
    p2, A ,B =permute!(p2,perm), permute!(A,perm), permute!(B,perm)
    n = length(p2)
    # TODO Casimir
    prefactor = Nc*Z2/(4*pi^2)
    M = B ./ A
    derivM = Deriv(p2,M)
    integrand = prefactor * (M./A) .* p2 .* ( M .+ 0.5*p2.*derivM) ./ (p2 .+ M.^2).^2
    decayconstant2 = 0.0
    for i in 1:n-1
        decayconstant2 += 0.5*(integrand[i]+integrand[i+1])*(p2[i+1]-p2[i])
    end
    return sqrt(decayconstant2)
end
function renormalization_point(p2::Array{Float64,1},G0::Array{Float64,1},Z0::Array{Float64,1},QCD)
    # find highest value of p2 for which G^2*Z = 1 i.e. find root of G²Z-1=0
    x,G,Z = permute_ym(p2,G0,Z0,QCD)
    y = copy(G.*G.*Z .- 1.0)
    # sort x and y
    # find highest [x0,x1] interval in x where y changes sign
    x0,y0 = NaN, NaN
    x1,y1 = NaN, NaN
    for i in length(x):-1:2
        if sign(y[i]) != sign(y[i-1])
            x0, y0 = x[i-1], y[i-1]
            x1, y1 = x[i], y[i]
            break
        end
    end
    highestroot = y0*(x0-x1)/(y1-y0) + x0
    return highestroot
end
