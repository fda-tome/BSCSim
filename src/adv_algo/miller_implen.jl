using SpecialFunctions
function besselN(N, z)
    cosz = cos(z)
    sinz = sin(z)
    frac = (4 * N^2 - 1) / (8 * z)
    nfrac = (4 * (N+1)^2 - 1) / (8 * z)
    if N % 2 == 0 
        return (-sinz + cosz +  frac * (cosz + sinz)) / (-cosz - sinz + nfrac * (-sinz + cosz))
    else
        return (sinz + cosz +  frac * (-cosz + sinz)) / (cosz - sinz + nfrac * (sinz + cosz))
    end
end


function besselj(orderLim, arg, tol)
    N = tol * orderLim
    println(N)
    values = Vector{typeof(arg)}(undef, N+2)    
    values[N+1] = 1
    values[N+2] = 0
    for i in N+1:-1:2
        values[i-1] = (2 * (i-1) / arg) * values[i] - values[i+1]
    end 
    norm = values[1]
    for i in 3:2:N+2
        norm += 2 * values[i]
    end
    norm = 1/norm
    values = values[1:orderLim+1]
    return norm.*values
end

function shortbess!(arr::Matrix{BigFloat}, dimy, tol)
    j = [1, 2, 6]
    for k in 1:dimy    
        values1 = besselj(5, (2 * big(pi) * big(sind(big(1))) * big(k/10000)) / big(0.00154), tol)
        values2 = besselj(5, (2 * big(pi) * big(sind(big(10))) * big(k/10000)) / big(0.00154), tol)
        values3 = besselj(5, (2 * big(pi) * big(sind(big(40))) * big(k/10000)) / big(0.00154), tol)
        for i in 0:2
            arr[1 + 3*i, k] = values1[j[i+1]]
            arr[2 + 3*i, k] = values2[j[i+1]]
            arr[3 + 3*i, k] = values3[j[i+1]]
        end
    end
end

function longbess!(arr::Matrix{BigFloat}, dimy, tol)
    j=[1, 2, 6]
    for k in 1:dimy
        values1 = besselj(5, (2 * big(pi) * big(sind(big(1))) * k ) / big(0.00154), tol)
        values2 = besselj(5, (2 * big(pi) * big(sind(big(10))) * k) / big(0.00154), tol)
        values3 = besselj(5, (2 * big(pi) * big(sind(big(40))) * k) / big(0.00154), tol)
        for i in 0:2
            arr[1 + 3*i, k] = values1[j[i+1]]
            arr[2 + 3*i, k] = values2[j[i+1]]
            arr[3 + 3*i, k] = values3[j[i+1]]
        end
    end
end


function compare!(func::Function,  arrmat::Matrix{BigFloat}, arrjl::Matrix{BigFloat}, dimx::Int, dimy::Int, csvname::String, tol)
    func(arrjl, dimy, tol)
    arrmat = readdlm(csvname, ',', BigFloat)
    reshape(arrmat, dimx, dimy)
    err = zeros(BigFloat, dimx, dimy)
    for i in 1:dimx
        for j in 1:dimy
            err[i,j] = big(arrjl[i,j] - arrmat[i,j])
        end
    end
    println(100 * findmax(err)[1] / arrmat[findmax(err)[2]])
    println(100 * findmin(err)[1] / arrmat[findmin(err)[2]])
end










