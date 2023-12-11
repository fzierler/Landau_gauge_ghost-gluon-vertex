chebyshevnodes(n)     = gausschebyshev(n,1)[1]
chebyshevnodes(n,a,b) = gausschebyshev(n,a,b,order=1)[1]
function chebyshev(x,n)
    p2 = one(x)
    p1 = copy(x)
    for i in 1:n
        p3 = p2
        p2 = p1
        p1 = muladd(2x,p2,-p3)
    end
    return p2
end
# use cubic spline inerpolation to determine y at the Cheby nodes
chebycoeff(x,y,n,a,b) =  chebycoeff(splines(chebyshevnodes(n,a,b),x,y))
# assumes x-values are transformed Chebyshev nodes
function chebycoeff(y)
    n     = length(y)
    nodes = chebyshevnodes(n)
    coeff = zero(y)
    for i in 1:n
        coeff[i] = (2/n)*sum(y.*chebyshev.(nodes,i-1))
    end
    coeff[1] /= 2
    return coeff
end
function clenshaw(x,coeff,a,b)
    ret = zero(x)
    clenshaw!(ret,x,coeff,a,b)
    return ret
end
function clenshaw!(ret,x,coeff,a,b)
    y = log_mapping_inv(x,a,b)
    return clenshaw!(ret,y,coeff)
end
function clenshaw!(ret,x,coeff)
    @turbo for j in eachindex(ret)
        ret[j] = clenshaw(x[j],coeff)
    end
end
function clenshaw(x,coeff)
    m   = length(coeff)
    d   = zero(x)
    dd  = zero(x)
    ret = zero(x)
    @inbounds for i in m:-1:2
        ret    = muladd(x,2d,coeff[i]-dd)
        dd,d   = d,ret
    end
    @inbounds ret = muladd(x,d,coeff[1]-dd)
    return ret
end
