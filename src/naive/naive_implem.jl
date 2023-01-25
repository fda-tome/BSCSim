using SpecialFunctions
using AssociatedLegendrePolynomials
setprecision(167)


function besbeam(psiamp, axiconang, order, rho, phi, z) 
    k = (2 * pi) / big(1.54/1000) #fixed wavelength 1.54mm
    kz = k * cosd(axiconang)
    krho = k * sind(axiconang)
    return amp * besselj(ord, krho * rho) * exp(ord * im * phi) * exp(kz * im * z)
end

function makeplot
 

        
