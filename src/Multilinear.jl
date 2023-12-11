function monolinear(x,x0,y0)
    xind = findindex(x,x0)
    monolinear(x,x0,y0,xind)
end
function monolinear(x,x0,y0,xind)
    int  = zero(x)
    monolinear!(int,x,x0,y0,xind)
end
function monolinear!(int,x,x0,y0,xind)
    lx0  = length(x0)
    @inbounds for i in eachindex(x)
        xi = xind[i]
        xi == 0 ? xi = 1 : xi == lx0 ? xi -= 1 : nothing
        t = (x[i]-x0[xi])/(x0[xi+1]-x0[xi])
        t < 0.0 ? t = 0.0 : t > 1.0 ? t = 1.0 : nothing
        int[i] = (1-t)*y0[xi] + t*y0[xi+1]
    end
    return int
end
function monolinear_thread!(int,x,x0,y0,xind)
    lx0  = length(x0)
    @threads :static for i in eachindex(x)
        xi = xind[i]
        xi == 0 ? xi = 1 : xi == lx0 ? xi -= 1 : nothing
        t = (x[i]-x0[xi])/(x0[xi+1]-x0[xi])
        t < 0.0 ? t = 0.0 : t > 1.0 ? t = 1.0 : nothing
        int[i] = (1-t)*y0[xi] + t*y0[xi+1]
    end
    return int
end
function bilinear(x,y,x0,y0,data,xind,yind)
    inter = zeros(eltype(data),size(x))
    bilinear!(inter,x,y,x0,y0,data, xind, yind)
    return inter
end
function bilinear(x,y,x0,y0,data)
    xind = findindex(x,x0)
    yind = findindex(y,y0)
    return bilinear(x,y,x0,y0,data, xind, yind)
end
function bilinear!(inter,x,y,x0,y0,data,xind,yind)
    lx0, ly0 = length(x0), length(y0)
    @inbounds for i in eachindex(x)
        xi, yi = xind[i], yind[i]
        xi == 0 ? xi = 1 : xi == lx0 ? xi -= 1 : nothing
        yi == 0 ? yi = 1 : yi == ly0 ? yi -= 1 : nothing
        u = (y[i]-y0[yi])/(y0[yi+1]-y0[yi])
        t = (x[i]-x0[xi])/(x0[xi+1]-x0[xi])
        u < 0.0 ? u = 0.0 : u > 1.0 ? u = 1.0 : nothing
        t < 0.0 ? t = 0.0 : t > 1.0 ? t = 1.0 : nothing
        d1 = data[xi,  yi]
        d2 = data[xi+1,yi]
        d3 = data[xi+1,yi+1]
        d4 = data[xi  ,yi+1]
        inter[i] = (1-u)*((1-t)*d1 + t*d2) + u*(t*d3 + (1-t)*d4)
    end
end
function trilinear(x,y,z,x0,y0,z0,data,xind,yind,zind)
    ret = zero(x)
    trilinear!(ret,x,y,z,x0,y0,z0,data,xind,yind,zind)
    return ret
end
function trilinear(x,y,z,x0,y0,z0,data)
    xind = findindex(x,x0)
    yind = findindex(y,y0)
    zind = findindex(z,z0)
    return trilinear(x,y,z,x0,y0,z0,data,xind,yind,zind)
end
function trilinear!(ret,x,y,z,x0,y0,z0,data,xind,yind,zind)
    lx0, ly0, lz0 = length(x0), length(y0), length(z0)
    T = eltype(x)
    @inbounds for i in eachindex(x)
        xi, yi, zi = xind[i], yind[i], zind[i]
        xi == 0 ? xi = 1 : xi == lx0 ? xi = lx0-1 : nothing
        yi == 0 ? yi = 1 : yi == ly0 ? yi = ly0-1 : nothing
        zi == 0 ? zi = 1 : zi == lz0 ? zi = lz0-1 : nothing
        v = (z[i]-z0[zi])/(z0[zi+1]-z0[zi])
        u = (y[i]-y0[yi])/(y0[yi+1]-y0[yi])
        t = (x[i]-x0[xi])/(x0[xi+1]-x0[xi])
        # constant extrapolation (for linear extrapolation remove next 3 lines)
        v < zero(T) ? v = zero(T) : v > one(T) ? v = one(T) : nothing
        u < zero(T) ? u = zero(T) : u > one(T) ? u = one(T) : nothing
        t < zero(T) ? t = zero(T) : t > one(T) ? t = one(T) : nothing
        d1 = data[xi,  yi,  zi]
        d2 = data[xi+1,yi,  zi]
        d3 = data[xi+1,yi+1,zi]
        d4 = data[xi,  yi+1,zi]
        d5 = data[xi,  yi,  zi+1]
        d6 = data[xi+1,yi,  zi+1]
        d8 = data[xi,  yi+1,zi+1]
        d7 = data[xi+1,yi+1,zi+1]
        um1, tm1, vm1 = (one(T) - u), (one(T) - t), (one(T) - v)
        sum1 = muladd(um1,muladd(tm1,d1,t*d2),u*muladd(t,d3,tm1*d4))
        sum2 = muladd(um1,muladd(tm1,d5,t*d6),u*muladd(t,d7,tm1*d8))
        ret[i] = muladd(sum1,vm1,sum2*v)
    end
end
