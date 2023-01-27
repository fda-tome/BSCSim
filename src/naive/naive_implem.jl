using SpecialFunctions
using AssociatedLegendrePolynomials
using Plots

function besbeam(psiamp, axiconang, order, rho, phi, z) 
    k = 1000 * 2 * big(pi) / big(1.54) #fixed wavelength 1.54mm
    kz = k * big(cosd(axiconang))
    krho = k * big(sind(axiconang))
    return big(psiamp * besselj(order, krho * rho) * exp(order * im * phi) * exp(kz * im * z))
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
    savefig(plot(htmaps..., layout = (3, 3)), "plot.png")
end

makeplotbes(200, 200, 20, 20, 1)


