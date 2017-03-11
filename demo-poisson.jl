using DoubleFourierSphere

# poisson equation
M = 1
N = 1

Λ, Φ = spheregrids(64,32)
Θ = Φ + π/2

G = sphericalmode(M,N)(Λ, Φ) #forcing function
Gf = fftsphere(G)

Uf = zeros(Gf) # same size, type

# n here is meridional wave resolution
n = size(Gf, 2)

# m=0
#D, A = DAzero(n)
#Uf[1, :] = D \ A * Gf[1, :] # TODO ON THIS LINE, singularity exception ******

# m odd
for mi in 2:2:size(Gf,1)
  m = mi-1
  D, A = DAodd(n, m)
  Uf[mi, :] = D \ A * Gf[mi, :]
end

# m even
for mi in 3:2:size(Gf,1)
  m = mi-1
  D, A = DAeven(n, m)
  Uf[mi, :] = D \ A * Gf[mi, :]
end

U = G # not actually true, but some scalar multiple since spheremode is eigenfunction
Utest = ifftsphere(Uf) 

# check idempotency
println(maximum(abs(U-Utest)))
