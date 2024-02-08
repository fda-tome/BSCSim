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


function bsccalc(n, m, axiconang)
    fract = (im ^ (n - m)) * (2 * n + 1) * factorial(big(n - m)) / factorial(big(n+m))
    special = Plm(cos(axiconang), n, m)
    return fract * special 
end

function completebeam(x, y, z, axang, ord, amp, k, x0, y0, z0)
    beam = Vector{BigFloat}()
    for zi in z, j in y, i in x        
        push!(beam, abs(partialwavexp(amp, big(axang), ord, k, 0, 0, 0, i - x0, j - y0, zi - z0))^2)
    end
    return reshape(beam, (size(x)[1], size(y)[1], size(z)[1]))
end

function besbeam(psiamp, axiconang, k, order, rho, phi, z) 
    #k = 1000 * 2 * big(pi) / big(1.54) #fixed wavelength 1.54mm
    kz = k * cos(axiconang)
    krho = k * sin(axiconang)
    return big(psiamp * besselj(order, krho * rho) * cis(order * phi) * cis(kz * z))
end

function partialwavexp(psiamp, axiconang, order, k, r, theta, phi, x0, y0, z0)
    kr = k * r
    krho = k * sin(axiconang)
    kz = k * cos(axiconang)
    phi0 = Base.atan(y0,x0)
    rho0 = (x0^2+y0^2)^(1/2)
    nmax = Int64(ceil(kr + (big(405) / 100) * (kr^(1/3)) + 2))
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
    g(z,p) = f(z) * ℯ ^ (-2 * pi * q * z * im / l) 
    return (1/l)* solve(IntegralProblem(g, 0, l), HCubatureJL()).u
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
            #push!(test, completebeam(x, y, z, Base.acos(q/k + (2 * pi * i) / (k * l)), 0, amp, k, 0, 0, 0))
            push!(axang, rad2deg(Base.acos(q/k + (2 * pi * i) / (k * l))))
        end
        push!(test, fwpart)
    end
    println(findmax(axang))
    println(findmin(axang))
    return test
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

function makeplt(figname, xlabel, ylabel, x, y, z, amp, x0, y0, z0)
    axang = [deg2rad(1), deg2rad(10), deg2rad(40)]
    ord = [0, 1, 5]
    mag = 10
    alf = ['a','b', 'c', 'd', 'e', 'f', 'g', 'h', 'i' ]
    alfcount = 0
    htmaps = Vector{Any}()
    for axiconang in axang
        for order in ord
            println(order)
            alfcount+=1
            beam = completebeam(x, y, z, axiconang, order, amp, x0, y0, z0)
            beam = reshape(beam, (size(x)[1],size(y)[1]))
            push!(htmaps, heatmap(1000*x, 1000*y, beam, xlabel = latexstring(xlabel), ylabel= latexstring(ylabel), title = latexstring('(' * alf[alfcount] * ')') , titlefontsize = 16, xlabelfontsize = 12, ylabelfontsize = 12, xtickfontsize=9, ytickfontsize=9, left_margin=5mm, bottom_margin = 5mm))
        end
        x = x / mag
        x0 = x0 / mag
        y0 = y0 /mag
        y = y / mag
        mag = 4
    end
    savefig(plot(htmaps..., layout = (3, 3), size =(1200,900)), figname)
end

function makeplterr(figname, xlabel, ylabel, x, y, z, amp, x0, y0, z0)
    axang = big.([deg2rad(1)])
    ord = [0, 1, 5]
    mag = 10
    alf = ['a','b', 'c', 'd', 'e', 'f', 'g', 'h', 'i' ]
    alfcount = 0
    htmaps = Vector{Any}()
    for axiconang in axang
        for order in ord
            alfcount+=1
            beam1 = completebeam(x, y, z, axiconang, order, amp, x0, y0, z0)
            beam1 = reshape(beam1, (size(x)[1],size(y)[1]))
            push!(htmaps, heatmap(1000*x, 1000*y, beam1, xlabel = latexstring(xlabel), ylabel= latexstring(ylabel), title = latexstring('(' * alf[alfcount] * ')') , titlefontsize = 16, xlabelfontsize = 12, ylabelfontsize = 12, xtickfontsize=9, ytickfontsize=9, left_margin=5mm, bottom_margin = 5mm))
            alfcount+=1
            beam2 = completebeamorig(x, y, z, axiconang, order, amp, -x0, -y0, -z0)
            beam2 = reshape(beam2, (size(x)[1],size(y)[1]))
            push!(htmaps, heatmap(1000*x, 1000*y, beam2, xlabel = latexstring(xlabel), ylabel= latexstring(ylabel), title = latexstring('(' * alf[alfcount] * ')') , titlefontsize = 16, xlabelfontsize = 12, ylabelfontsize = 12, xtickfontsize=9, ytickfontsize=9, left_margin=5mm, bottom_margin = 5mm))
            err = abs.(beam1-beam2)
            err = err ./ beam2
            err = 100 * err
            err = log10.(err)
            alfcount+=1
            push!(htmaps, heatmap(1000*x, 1000*y, err, xlabel = latexstring(xlabel), ylabel= latexstring(ylabel), title = latexstring('(' * alf[alfcount] * ')') , titlefontsize = 16, xlabelfontsize = 12, ylabelfontsize = 12, xtickfontsize=9, ytickfontsize=9, left_margin=5mm, bottom_margin = 5mm))
        end
        x = x / mag
        y = y / mag
        mag = 4
    end
    savefig(plot(htmaps..., layout = (3, 3), size =(1200,900)), figname)
end

function test(z)
    l = 0.05
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




