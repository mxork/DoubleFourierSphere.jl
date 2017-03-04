using Helper
using FFTsphere
using DiffMatrices

# poissson equation

M = 4
N = 4

Λ, Φ = spheregrids(128,64)
Θ = Φ + π/2

G = fouriermode(M,N)(Λ, Θ) #forcing function
Gf = fftsphere(G)

Uf = zeros(Gf) # same size, type

# n here is meridional wave resolution
n = size(Gf, 2)

# m=0
D, A = DAzero(n)
Uf[1, :] = D \ A * Gf[1, :] # TODO ON THIS LINE ******

# m odd
for mi in 2:2:size(Gf,1)
  m = mi-1
  D, A = DAzero(n, m)
  Uf[mi, :] = Dz \ Az * Gf[mi, :]
end

# m even
for mi in 3:2:size(Gf,1)

end




U =



Utest = ifftsphere(Uf) 

# check idempotency
println(maximum(abs(U-Utest)))
