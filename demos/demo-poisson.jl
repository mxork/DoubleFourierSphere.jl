using DoubleFourierSphere

# poisson equation
M = 64 
N = 32

# wave numbers
km = 0
kn = 0

Λ, Φ = spheregrids(M,N)
Θ = Φ + π/2

G = sphericalmode(km,kn)(Λ, Φ) #forcing function
Gf = fftsphere(G)

Uf = zeros(Gf) # same size, type

# m=0
D, A = DAzero(N)
Uf[1, :] = D \ A * Gf[1, :] # TODO ON THIS LINE, singularity exception ******
Uf[1,1] = 0; # manual surgery to fix constant term.

# m odd
for mi in 2:2:size(Gf,1)
  m = mi-1
  D, A = DAodd(N, m)
  Uf[mi, :] = D \ A * Gf[mi, :]
end

# m even
for mi in 3:2:size(Gf,1)
  m = mi-1
  D, A = DAeven(N, m)
  Uf[mi, :] = D \ A * Gf[mi, :]
end

U = G # not actually true, but some scalar multiple since spheremode is eigenfunction
Utest = ifftsphere(Uf) 

# check idempotency
println(maximum(abs(U-Utest)))
