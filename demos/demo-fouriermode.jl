using Helper

print("Longitudinal wavenumber m: ")
M = parse(Int, readline())

print("Latitudinal wavenumber n: ")
N = parse(Int, readline())

# small coordinate shift: Muraki uses
# φ=0 at equator, vs. Cheong who takes
# it at the south pole.
Gλ, Gφ = spheregrids(128,64)
Gφc = Gφ + π/2
sphereplot(fouriermode(M,N)(Gλ, Gφc))
