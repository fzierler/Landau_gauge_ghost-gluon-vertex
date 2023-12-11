function solve_tridiagonal(a,b,c,r)
    N = size(a)[1]
    T = eltype(a)
    @inbounds begin
        beta,rho,z = zeros(T,N),zeros(T,N),zeros(T,N)
        beta[1],rho[1] = b[1],r[1]
        for i in 2:N
            beta[i] = b[i]-(a[i]*c[i-1])/beta[i-1]
            rho[i] = r[i]-(a[i]*rho[i-1])/beta[i-1]
        end
        z[N] = rho[N]/beta[N]
        for i in 1:N-1
            z[N-i] = (rho[N-i]-c[N-i]*z[N-i+1])/beta[N-i]
        end
    end
    return z
end
function second_derivative(x,y)
    N = size(x)[1]
    T = eltype(x)
    a,b,c,r = zeros(T,N),zeros(T,N),zeros(T,N),zeros(T,N)
    @inbounds begin
        for i in 2:(N-1)
            r[i] = (y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1])
            a[i] = (x[i]-x[i-1])/6
            b[i] = (x[i+1]-x[i-1])/3
            c[i] = (x[i+1]-x[i])/6
        end
        #Zusaetzliche Koffizienten durch Wahl: y2_N = y2_(N-1);y2_0 = y2_1
        a[1],a[N] = 0,-1
        b[1],b[N] = 1,1
        c[1],c[N] = -1,0
        r[1],r[N] = 0,0
    end
    ret = solve_tridiagonal(a,b,c,r)
    return ret
end
function splines(x,xdata0,ydata0)
    xdata, ydata = deepcopy(xdata0),deepcopy(ydata0)
    #Sortiere Datenpunkte
    permutation = sortperm(xdata)
    permute!(xdata,permutation)
    permute!(ydata,permutation)
    #Suche passendes Intervall mit xdata[i] < x < xdata[i+1]
    index = findindex(x,xdata)
    return splines(x,xdata,ydata,index)
end
# assumes that datset is sorted
function splines(x,xdata,ydata,index)
    interpolatedPoints = zeros(size(x))
    return splines!(interpolatedPoints,x,xdata,ydata,index)
end
# precompute second derivat
function splines!(points,x,xdata,ydata,index)
    secondderiv = second_derivative(xdata,ydata)
    return splines!(points,x,xdata,ydata,index,secondderiv)
end
function splines!(points,x,xdata,ydata,index,y2)
    l = length(x)
    ldata = length(xdata)
    #Bilde 2.Ableitungen
    @inbounds for t in 1:l
        ind = index[t]
        if ind == 0
            points[t] = ydata[1]
        elseif ind == ldata
            points[t] = ydata[ldata]
        else
            Δx = xdata[ind+1]-xdata[ind]
            A = (xdata[ind+1]-x[t])/Δx
            B = (x[t]-xdata[ind])/Δx
            C = muladd(A,A,-1)*Δx*Δx/6
            D = muladd(B,B,-1)*Δx*Δx/6
            points[t] = A*muladd(C,y2[ind],ydata[ind])
            points[t] += B*muladd(D,y2[ind+1],ydata[ind+1])
        end
    end
    return points
end
function splines_threaded!(points,x,xdata,ydata,index,y2)
    l = length(x)
    ldata = length(xdata)
    #Bilde 2.Ableitungen
    @threads :static for t in 1:l
        @inbounds begin
            ind = index[t]
            if ind == 0
                points[t] = ydata[1]
            elseif ind == ldata
                points[t] = ydata[ldata]
            else
                Δx = xdata[ind+1]-xdata[ind]
                A = (xdata[ind+1]-x[t])/Δx
                B = (x[t]-xdata[ind])/Δx
                C = muladd(A,A,-1)*Δx*Δx/6
                D = muladd(B,B,-1)*Δx*Δx/6
                points[t] = A*muladd(C,y2[ind],ydata[ind])
                points[t] += B*muladd(D,y2[ind+1],ydata[ind+1])
            end
        end
    end
    return points
end
function splines2functions!(PointsA,PointsB,x,xdata,ydataA,ydataB,index)
    l = length(x)
    ldata = length(xdata)
    y2A = second_derivative(xdata,ydataA)
    y2B = second_derivative(xdata,ydataB)
    @inbounds @threads :static for t in 1:l
        ind = index[t]
        if ind == 0
            PointsA[t] = ydataA[1]
            PointsB[t] = ydataB[1]
        elseif ind == ldata
            PointsA[t] = ydataA[ldata]
            PointsB[t] = ydataB[ldata]
        else
            A = (xdata[ind+1]-x[t])/(xdata[ind+1]-xdata[ind])
            B = (x[t]-xdata[ind])/(xdata[ind+1]-xdata[ind])
            C = (1/6)*(A^3 - A)*(xdata[ind+1]-xdata[ind])^2
            D = (1/6)*(B^3 - B)*(xdata[ind+1]-xdata[ind])^2
            PointsA[t] = A*ydataA[ind] + B*ydataA[ind+1] + C*y2A[ind]+D*y2A[ind+1]
            PointsB[t] = A*ydataB[ind] + B*ydataB[ind+1] + C*y2B[ind]+D*y2B[ind+1]
        end
    end
end
function splines2functions2!(PointsA,PointsB,x,xdata,ydataA,ydataB,index)
    y2A = second_derivative(xdata,ydataA)
    y2B = second_derivative(xdata,ydataB)
    splines2functions2sd!(PointsA,PointsB,x,xdata,ydataA,ydataB,index,y2A,y2B)
end
function splines2functions2sd!(PointsA,PointsB,x,xdata,ydataA,ydataB,index,y2A,y2B)
    l = length(x)
    ldata = length(xdata)
    @inbounds for t in 1:l
        ind = index[t]
        if ind == 0
            PointsA[t] = ydataA[1]
            PointsB[t] = ydataB[1]
        elseif ind == ldata
            PointsA[t] = ydataA[ldata]
            PointsB[t] = ydataB[ldata]
        else
            Δx = xdata[ind+1]-xdata[ind]
            A = (xdata[ind+1]-x[t])/Δx
            B = (x[t]-xdata[ind])/Δx
            C = muladd(A,A,-1)*Δx*Δx/6
            D = muladd(B,B,-1)*Δx*Δx/6
            PointsA[t] = A*muladd(C,y2A[ind],ydataA[ind]) + B*muladd(D,y2A[ind+1],ydataA[ind+1])
            PointsB[t] = A*muladd(C,y2B[ind],ydataB[ind]) + B*muladd(D,y2B[ind+1],ydataB[ind+1])
        end
    end
end
