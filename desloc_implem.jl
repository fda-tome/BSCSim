using Dagger
using SpecialFunctions
using LegendrePolynomials
using Plots
using Integrals
using Plots.PlotMeasures
using LaTeXStrings

global fact = zeros(typeof(big(1)*big(1)), 10000)
fact[1] = 1

function factorial(x::T) where T<:BigInt
    if x == 0
        return 1
    end
    if fact[x] == zero(BigInt)
        i::typeof(x*x) = x
        while fact[i] == zero(typeof(x*x))
            i-=1
        end
        i+=1
        while i != x + 1
            fact[i] = fact[i-1] * i
            i+=1
        end
    end
    return fact[x]
end

function SpecialFunctions.besselj(orderLim, arg, tol)
    @assert tol>orderLim "Tolerance must be greater than the order to be calculated"
    N = tol
    values = Vector{typeof(arg)}(undef, N+2)    
    values[N+1] = 1
    values[N+2] = 0
    norm = 0.0
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

function bsccalc(n, m, axiconang)
    fract = (im ^ (n - m)) * (2 * n + 1) * factorial(big(n - m)) / factorial(big(n+m))
    special = Plm(cos(axiconang), n, m)
    return fract * special 
end

function completebeam(x, y, z, axang, ord, amp, k, x0, y0, z0)
    beam = []
    for zi in z, j in y, i in x        
        push!(beam, partialwavexp(amp, axang, ord, k, 0, 0, 0, i - x0, j - y0, zi - z0))
    end
    return reshape(beam, (size(x)[1], size(y)[1], size(z)[1]))
end

function besbeam(psiamp, axiconang, k, order, rho, phi, z) 
    kz = k * cos(axiconang)
    krho = k * sin(axiconang)
    return psiamp * besselj(order, krho * rho) * cis(order * phi) * cis(kz * z)
end

function partialwavexp(psiamp, axiconang, order, k, r, theta, phi, x0, y0, z0)
    kr = k * r
    krho = k * sin(axiconang)
    kz = k * cos(axiconang)
    phi0 = Base.atan(y0,x0)
    rho0 = (x0^2+y0^2)^(1/2)
    nmax = Int64(ceil(kr + (405 / 100) * (kr^(1/3)) + 2))
    psi = 0
    besselm = []
    cisvalm = []
    for m in -nmax:nmax
        push!(cisvalm,  cis(-(m - order) * phi0) * cis(-kz * z0))
        push!(besselm, besselj(m - order, Float64(krho * rho0)))
    end
    for n in 0:nmax
        spher = sphericalbesselj(Int64(n), Float64(kr))
        for m in -n:n
            BSC = bsccalc(n, m, axiconang)
            psi += cisvalm[m + nmax + 1] * besselm[m + nmax + 1] * BSC * spher * Plm(cos(theta), n, m) * cis(m * phi)
        end
    end
    return psi * psiamp
end

function aq(f, q, l)
    g(z,p) = f(z) * cis(- 2 * pi * q * z / l) 
    return (1/l)* solve(IntegralProblem(g, 0.0, l), QuadGKJL()).u
end

function fw(f, n, q, l, k, x, y, z)
    test = []
    axang = []
    w = 2 * pi * 2.5e6
    c = 340
    for j in z
        fwpart = 0
        for i in -n:n
            kzq = q + 2 * pi * i / l
            if ((w^2 / c^2) - kzq^2) < 0
                println("erro krho")
            end
            if w / kzq ≤ 0
                println("erro")
            end
            amp = aq(f, i, l) 
            fwpart += partialwavexp(amp, Base.acos(q/k + (2 * pi * i) / (k * l)), 0, k, 0, 0, 0, 0, 0, j)
            push!(axang, rad2deg(Base.acos(q/k + (2 * pi * i) / (k * l))))
        end
        push!(test, fwpart)
    end
    println(findmax(axang))
    println(findmin(axang))
    return test
end

function pfw(tf, n, q, l, c, f, x, y, z, x0, y0, z0)
    test = []
    axang = []
    lambda = c / f
    k = 2 * pi / lambda
    w = 2 * pi * f
    q = q * k
    space = []
    for i in -n:n
        kzq = q + 2 * pi * i / l
        if ((w^2 / c^2) - kzq^2) < 0
            println("erro krho")
        end
        if w / kzq ≤ 0
            println("erro")
        end
        amp = aq(tf, i, l) 
        push!(space, Dagger.@spawn completebeam(x, y, z, Base.acos(q/k + (2 * pi * i) / (k * l)), 0, amp, k, x0, y0, z0))
        push!(axang, rad2deg(Base.acos(q/k + (2 * pi * i) / (k * l))))
    end
    println("Max and Min angles")
    println(findmax(axang))
    println(findmin(axang))
    println("Max and Min radius")
    println(l*sqrt((k/(k*cos(deg2rad(findmax(axang)[1]))))^2 - 1))
    return Dagger.@spawn reduce(+, fetch.(space))
end

function completebeamorig(x, y, z, axang, ord, amp, x0, y0, z0)
    beam = Vector{BigFloat}()
    for k in z, j in y, i in x
        rho =  sqrt(big(i)^2 + big(j)^2)
        rho0 = sqrt(x0^2 + y0^2)
        phi0 = big(atan(y0, x0))
        phi = big(atan(j, i))
        push!(beam, abs(besbeam(amp, big(axang), ord, sqrt(rho^2 + rho0^2 - 2 * rho * rho0 * cos(phi - phi0)), phi, k))^2)
    end
    return reshape(beam, (size(x)[1], size(y)[1], size(z)[1]))
end

function f1(z)
    l = 0.1
    l1 = l / 10
    l2 = 3 * l / 10
    l3 = 4 * l / 10
    l4 = 6 * l / 10
    l5 = 7 * l / 10
    l6 = 9 * l / 10
    if l1 ≤ z ≤ l2
        return -4 * (z - l1) * (z - l2) / (l2-l1) ^ 2
    end
    if l3 ≤ z ≤ l4
        return 1
    end
    if l5 ≤ z ≤ l6
        return -4 * (z - l5) * (z - l6) / (l6-l5) ^ 2
    end
    return 0
end

function f2(z)
    l = 0.1
    if 0.2 * l ≤ z ≤ 0.8 * l
        return 1
    else
        return 0
    end
end

function f3(z)
    l = 0.1
    if 0.2 * l ≤ z ≤ 0.8 * l
        return exp(-50*((z - 0.5 * l) / l)^2)
    else
        return 0
    end
end

function f4(z)
    l = 0.1
    if 0.2 * l ≤ z ≤ 0.8 * l
        return exp(-50*((z - 0.5 * l) / l)^4) * cos(6 * pi * z / l)
    else
        return 0
    end
end






