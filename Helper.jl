export legendre, spheregrids, sphericalmode, fouriermode, sphereplot
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
  # sin(θ) here (instead of expected cos(θ)) agrees with Muraki's,
  # but as far as I can tell, Muraki's spherical harmonics are not right.
  # test case (m,n) = (2,2)
  # Seemingly, Muraki takes n = 2*n, which is strange
  return (λ, θ) -> legendre(m,n,cos(θ)) .* cos(m*λ)
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

# plots the 2D matrix G on our canonical sphere
using PyPlot: contourf, ColorMap, clf
function sphereplot(G)
  clf()
  M, N = size(G,1), size(G,2)
  contourf(spheregrids(M, N)..., G, cmap=ColorMap("inferno"))
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

function latitude_interior_grid(A)
    N = size(A, 2)
    Φs = [ π*(j+0.5)/N for j in 0:N-1]
end
