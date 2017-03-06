using Helper
using FFTsphere

M = 4
N = 4

Gλ, Gφ = spheregrids(128,64)
Gφc = Gφ + π/2

U = fouriermode(M,N)(Gλ, Gφc)
Uf = fftsphere(U)
Utest = ifftsphere(Uf) 

# check idempotency
println(maximum(abs(U-Utest)))
