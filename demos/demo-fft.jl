using DoubleFourierSphere
using Plots

M = 4
N = 3

Gλ, Gφ = spheregrids(64,32)
U = fouriermode(M,N)(Gλ, Gφ)

M, Ms = zonal_modes(U)
N, Ns0, Ns = meridional_modes(U)

Uf = fft_sphere(U)

# display

pyplot()

p1 = contour(1:64, 1:32, U')
p2 = contour(Ms, Ns, abs(Uf'))

plot(p1, p2)

savefig("demo-fft.svg")
