#####################################################
# broaden scope of functions from FastGaussQuadrature
#####################################################
function gausslegendre(n,a,b)
    x,w = gausslegendre(n)
    A = -log(a*b)/log(b/a)
    B = 2/log(b/a)
    @. w = w*exp((x-A)/B)/B
    @. x = exp((x-A)/B)
    return x,w
end
function gausschebyshev(n,a,b;order=2)
    x,w = gausschebyshev(n,order)
    A = -log(a*b)/log(b/a)
    B = 2/log(b/a)
    @. w = w*exp((x-A)/B)/B
    @. x = exp((x-A)/B)
    return x,w
end
function log_mapping_inv(x,a,b)
    log_mapping_inv!(zero(x),x,a,b)
end
function log_mapping_inv!(ret,x,a,b)
    A = -log(a*b)/log(b/a)
    B = 2/log(b/a)
    @turbo @. ret = A + B*log(x)
end
