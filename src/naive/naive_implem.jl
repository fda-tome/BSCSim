using Bessels
using SpecialFunctions
using AssociatedLegendrePolynomials
using Plots



sphericalbesselj(nu, x::T) where T = sqrt(T(pi)/(2 * x)) * Bessels.besselj(nu + one(T)/2, x)

function generalPlm(l, m, x)
    if m < 0
        return ((-1)^m) * (factorial(big(l-m)) / factorial(big(l+m))) * Plm(l, -m, x)
    else
        return Plm(l, m, x)
    end
end

function besbeam(psiamp, axiconang, order, rho, phi, z) 
    k = 1000 * 2 * big(pi) / big(1.54) #fixed wavelength 1.54mm
    kz = k * big(cosd(axiconang))
    krho = k * big(sind(axiconang))
    return big(psiamp * besselj(order, krho * rho) * cis(order * phi) * cis(kz * z))
end

function bsccalc(n, m, axiconang, order, krho, kz, phi0, z0, rho0)
    fract = (im ^ (n - m)) * (2 * n + 1) * factorial(big(n - m)) / factorial(big(n+m))
    special = SpecialFunctions.besselj(m - order, krho * rho0) * generalPlm(n, m, cos(axiconang))
    cisval = cis(-(m + order) * phi0) * cis(-kz * z0)
    return fract * special * cisval
end

function partialwavexp(psiamp, axiconang, order, r, theta, phi)
    k = 1000 * 2 * big(pi) / 1.54
    kr = k * r
    krho = k * sin(axiconang)
    kz = k * cos(axiconang)
    phi0 = z0 = rho0 = 0
    nmax = Int64(ceil(kr + (big(405) / 100) * (kr^(1/3)) + 2))
    psi = 0
    for n in 0:nmax
        spher = sphebesselj(Int64(n), Float64(kr))
        for m in -n:n
            BSC = bsccalc(n, m, axiconang, order, krho, kz, phi0, z0, rho0)
            psi += BSC * spher * generalPlm(n, m, cos(theta)) * cis(m * phi)
        end
    end

    return psi * psiamp
end


function makeplotbes(nx, ny, ry, rx, psiamp)
    z = 0
    x = range(-rx, rx, nx)
    y = range(-ry, ry, ny)
    axang = [deg2rad(1), deg2rad(10), deg2rad(40)]
    ord = [0, 1, 5]
    htmaps = Vector{Any}()
    for order in ord, axiconang in axang
        beamaux = Vector{BigFloat}()
        for j in y, i in x
            rho =  sqrt(big(i)^2 + big(j)^2)
            if j > 0
                phi = big(acos(big(i) / rho))
            else
                if i > 0
                    phi = 2 * pi - big(acos(big(i) / rho))
                else
                    phi = pi +  big(acos(-big(i) / rho))
                end
            end
            push!(beamaux,abs(besbeam(psiamp, big(axiconang), order, rho, phi, z)))
        end
        push!(htmaps, heatmap(x, y,  reshape(beamaux, (nx, ny))))
    end
    savefig(plot(htmaps..., layout = (3, 3),  xtickfontsize = 7, yfontsize = 7), "plotorg.png")
end

function makepltpartial(nx, ny, ry, rx, psiamp)
    z = 0
    x = range(-rx, rx, nx)
    y = range(-ry, ry, ny)
    axang = [deg2rad(1), deg2rad(10), deg2rad(40)]
    ord = [0, 1, 5]
    htmaps = Vector{Any}()
    for order in ord, axiconang in axang
        beamaux = Vector{BigFloat}()
        for j in y, i in x
            rho =  sqrt(big(i)^2 + big(j)^2)
            r = sqrt(big(i)^2 + big(j)^2 + z^2)
            if j > 0
                phi = big(acos(big(i) / rho))
            else
                if i > 0
                    phi = 2 * pi - big(acos(big(i) / rho))
                else
                    phi = pi +  big(acos(-big(i) / rho))
                end
            end
            theta = acos(z)
            push!(beamaux,abs(partialwavexp(psiamp, big(axiconang), order, r, theta, phi)))
        end
        push!(htmaps, heatmap(x, y,  reshape(beamaux, (nx, ny))))
    end
    savefig(plot(htmaps..., layout = (3, 3),  xtickfontsize = 7, yfontsize = 7), "plotpart.png")
end

makepltpartial(200, 200, 20, 20, 1)
makeplotbes(200, 200, 20, 20, 1)







