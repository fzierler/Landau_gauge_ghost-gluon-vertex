function findindex!(x,xdata,index)
    findindex_hunt!(x,xdata,index)
end
function findindex(x,xdata)
    index = zeros(Int,size(x))
    findindex!(x,xdata,index)
    return index
end
function findindex_hunt!(x,xdata,index)
    n = length(xdata)
    m = length(x)
    lower = 1
    xend = xdata[end]
    xone = xdata[1]
    for t in 1:m
        xt = x[t]
        # check if extrapolation is needed
        if xt >= xend
            index[t] = n
            lower = n-1
            continue
        elseif xt <= xone
            index[t] = 0
            lower = 1
            continue
        end
        # check if current intervall is correct
        if xdata[lower] <= xt < xdata[lower + 1]
            index[t] = lower
            continue
        end
        # check wether if we need to search upper or lower intervall
        if xt < xdata[lower]
            upper = lower
            lower = 1
            up    = false
        else
            lower = lower + 1
            upper = n
            up    = true
        end
        # perform hunt
        inc = 1
        if up == false
            while true
                lower  = upper - inc
                if lower < 1
                    lower = 1
                    break
                elseif xt >= xdata[lower]
                    break
                else
                    upper = lower
                    inc *= 2
                end
            end
        else
            while true
                upper = lower + inc
                if upper > n - 1
                    upper = n
                    break
                elseif xt < xdata[upper]
                    break
                else
                    lower = upper
                    inc *= 2
                end
            end
        end
        #now perform bisection
        (xt == xdata[lower]) && (upper=lower)
        while (upper - lower) > 1
            # calculate midpoint, div() is integer divison
            mid = div(upper + lower,2)
            if xt < xdata[mid]
                upper = mid
            else
                lower = mid
                (xt == xdata[lower]) && break
            end
        end
        index[t] = lower
    end
end

function findindex_naive!(x,xdata,index)
    for t in eachindex(x)
        if x[t] > xdata[end]
            index[t] = length(xdata)
        end
        for i in eachindex(xdata)
            if x[t] < xdata[i]
                index[t] = i - 1
                break
            end
        end
    end
end
function findindex_bisect!(x,xdata,index)
    n = length(xdata)
    m = length(x)
    for t in 1:m
        if x[t] > xdata[end]
            index[t] = n
            continue
        elseif x[t] < xdata[1]
            index[t] = 0
            continue
        end
        lower = 1   # initial lower bound for bisection
        upper = n   # initial upper bound for bisection
        while (upper - lower) > 1
            # calculate midpoint, div() is integer divison
            mid = div(upper + lower,2)
            if x[t] < xdata[mid]
                upper = mid
            else
                lower = mid
            end
        end
        index[t] = lower
    end
end
