using SpecialFunctions
using LegendrePolynomials
using DelimitedFiles

function shortsphebess!(arr::Matrix{BigFloat}, dimy)
    for j in 1:25    
        for i in 1:10000
            arr[j, i] = big(sphericalbesselj(j-1, (2 * (pi) * (i / 10000)) / (0.00154)))
            arr[j, i] = big(sphericalbesselj(j-1, (2 * (pi) * (i / 10000)) / (0.00154)))
            arr[j, i] = big(sphericalbesselj(j-1, (2 * (pi) * (i / 10000)) / (0.00154)))
        end
    end
end


function shortbess!(arr::Matrix{BigFloat}, dimy)
    j = [0, 1, 5]
    for k in 1:dimy    
        for i in 0:2
            arr[1 + 3*i, k] = big(besselj(j[i+1], (2 * big(pi) * big(sind(big(1))) * big(k / 10000)) / big(0.00154)))
            arr[2 + 3*i, k] = big(besselj(j[i+1], (2 * big(pi) * big(sind(big(10))) * big(k / 10000)) / big(0.00154)))
            arr[3 + 3*i, k] = big(besselj(j[i+1], (2 * big(pi) * big(sind(big(40))) * big(k / 10000)) / big(0.00154)))
        end
    end
end

function longsphebess!(arr::Matrix{BigFloat}, dimy)
    for j in 1:25    
        for i in 1:100000
            arr[j, i] = big(sphericalbesselj(j-1, (2 * (pi) * (i)) / (0.00154)))
            arr[j, i] = big(sphericalbesselj(j-1, (2 * (pi) * (i)) / (0.00154)))
            arr[j, i] = big(sphericalbesselj(j-1, (2 * (pi) * (i)) / (0.00154)))
        end
    end
end


function longbess!(arr::Matrix{BigFloat}, dimy)
    j=[0, 1, 5]
    for k in 1:dimy
        for i in 0:2
            arr[1 + 3*i, k] = big(besselj(j[i+1], (2 * big(pi) * big(sind(big(1))) * k ) / big(0.00154)))
            arr[2 + 3*i, k] = big(besselj(j[i+1], (2 * big(pi) * big(sind(big(10))) * k) / big(0.00154)))
            arr[3 + 3*i, k] = big(besselj(j[i+1], (2 * big(pi) * big(sind(big(40))) * k) / big(0.00154)))
        end
    end
end

function legenpol!(arr::Matrix{BigFloat}, dimy)
    for j in 1:25
        for i in 1:j
            for k in 1:1000
                arr[j, 1000 * (i - 1) + k] = big(Plm(big((k - 500) / 500), j - 1, i - 1 ))
            end
        end
    end
end

function compare!(func::Function,  arrmat::Matrix{BigFloat}, arrjl::Matrix{BigFloat}, dimx::Int, dimy::Int, csvname::String)
    func(arrjl, dimy)
    arrmat = readdlm(csvname, ',')
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


#=arrmat = Matrix{BigFloat}(undef, 25, 10000)
arrjl = Matrix{BigFloat}(undef, 25, 10000)
compare!(shortsphebess!, arrmat, arrjl, 25, 10000, "sphe_bess_curto_mat.csv")
arrmat = Matrix{BigFloat}(undef, 25, 100000)
arrjl = Matrix{BigFloat}(undef, 25, 100000)
compare!(longsphebess!, arrmat, arrjl, 25, 100000, "sphe_bess_longo_mat.csv")=#
#=arrmat = Matrix{BigFloat}(undef, 9, 10000)
arrjl = Matrix{BigFloat}(undef, 9, 10000)
compare!(shortbess!, arrmat, arrjl, 9 , 10000, "bess_curto_mat.csv")
arrmat = Matrix{BigFloat}(undef, 9, 100000)
arrjl = Matrix{BigFloat}(undef, 9, 100000)
compare!(longbess!, arrmat, arrjl, 9, 100000, "bess_longo_mat.csv")
=#
arrmat = zeros(BigFloat, 25, 25000)
arrjl = zeros(BigFloat, 25, 25000)
compare!(legenpol!, arrmat, arrjl, 25, 25000, "pol_leg_mat.csv")



        

