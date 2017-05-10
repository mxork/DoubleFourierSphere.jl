using DoubleFourierSphere

M = 4
N = 3

Gλ, Gφ = spheregrids(64,32)
Gφc = Gφ + π/2

U = fouriermode(M,N)(Gλ, Gφc)
Uf = fftsphere(U)
Utest = ifftsphere(Uf) 

# check idempotency
println(maximum(abs(U-Utest)))
