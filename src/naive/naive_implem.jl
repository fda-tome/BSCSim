using SpecialFunctions
using AssociatedLegendrePolynomials
using Plots
setprecision(167)


function besbeam(psiamp, axiconang, order, rho, phi, z) 
    k = big(2 * pi) / big(1.54/1000) #fixed wavelength 1.54mm
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
    htmaps = [(heatmap([(abs(besbeam(psiamp, big(axiconang), order, big(sqrt(i^2 + j^2)), big(acos(i / big(sqrt(i^2 + j^2)))), z))) for i in x, j in y])) for order in ord, axiconang in axang]
    gr()
    display(plot(vec(htmaps)..., layout = (3, 3)))
end

makeplotbes(200, 200, 20, 20, 1)


