export legendre, spheregrids, sphericalmode, sphericalmode_complex, fouriermode, plot_sphere, plot_frequency
export zonal_modes, meridional_modes, latitude_interior_grid

# associated legendre polynomial 
import GSL.sf_legendre_sphPlm
function legendre(m,n,X)
    [sf_legendre_sphPlm(n,m,x) for x in X]
end

# returns meshgrid for our canonical coordinate system on the sphere;
# takes resolution in latitude, longitude as input.
function spheregrids(Nλ::Int, Nφ::Int)
    dλ, dφ = 2*π / Nλ, π/Nφ
    return [λ for λ in -π:dλ:(π-dλ/2), _ in 1:Nφ], [φ for _ in 1:Nλ, φ in π*((0:Nφ-1) .+ 0.5)/Nφ]
end

# spherical to cartesian; should match MATLAB spec
# https://www.mathworks.com/help/matlab/ref/sph2cart.html
function sph2cart(Gλ, Gφ, Gr)
    return Gr .* cos(Gλ) .* cos(Gφ),
    Gr .* sin(Gλ) .* cos(Gφ),
    Gr .* sin(Gφ)
end

# cartesian to spherical; should match MATLAB spec
# https://www.mathworks.com/help/matlab/ref/cart2sph.html
function cart2sph(Gx, Gy, Gz)
    return atan2(Gy, Gx),
    atan2(Gz, sqrt(Gx.^2 + Gy.^2)),
    sqrt(Gx .^2 + Gy .^2 + Gz.^2)
end

# returns the n-m spherical mode on the sphere,
# m = λ wavenumber, longitude
# n = θ wavenumber, latitude
function sphericalmode(m,n)
    return (λ, θ) -> legendre(m,n,cos(θ)) .* cos(m*λ)
end

function sphericalmode_complex(m,n)
    return (λ, θ) -> legendre(m,n,cos(θ)) .* exp(1.0im*m*λ)
end

# returns the n-m fourier mode on the sphere,
# m = λ wavenumber, longitude
# n = θ wavenumber, latitude
function fouriermode(m, n)
    if m==0 
        return (λ, φ) -> cos(n*φ)
    elseif m%2 == 0
        return (λ, φ) -> sin(φ) .* sin(n*φ) .* cos(m*λ )
    else
        return (λ, φ) -> sin(n*φ) .* cos(m*λ)
    end
end

# FIXME
import Plots

# plots the 2D matrix G on our canonical sphere
function plot_sphere(G, title::String = "")
    # Plots.heatmap(longitude_grid(G), latitude_interior_grid(G), G, title=title)
    Plots.heatmap(G, title=title)
end

# ditto, but in frequency space, so no contour
# TODO make something prettier
function plot_frequency(Gf, title::String = "")
    # M, Ms = zonal_modes(Gf)
    # N, Ns0, Ns = meridional_modes(Gf)

    # Plots.heatmap(Ms, Ns, Gf, title=title)
    Plots.heatmap(Gf, title=title)
end

# returns maximal wavenumber and iterable of
# wavenumbers
function zonal_modes(A)
    M = Int(round(size(A,1)/2))
    Ms = [0:M-1 ; -M:-1]
    return M, Ms
end

function meridional_modes(A)
    N = size(A,2)
    Ns0 = 0:N-1
    Ns = 1:N
    return N, Ns0, Ns
end

function longitude_grid(A)
    X = size(A,1)
    0:2*pi/X:(2*pi-1/X)
end

function latitude_interior_grid(A)
    N = size(A, 2)
    # Φs = [ π*(j+0.5)/N for j in 0:N-1]
    Φs = 0.5*π/N:π/N:π*(1-0.5*1/N)
end

function average_sphere_spectral(Gf)
    sum( Gf[1, 1:2:end] .* (1 - ((0:2:size(Gf,2)-1).^2))) / size(Gf, 1) # cause we didn't normalize the longitude fft
end
