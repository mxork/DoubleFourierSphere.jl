using DoubleFourierSphere

M = 128
N = 64

Λ, Φ = spheregrids(M,N)

# for a variety of modes, check inversion property of 
# laplace - inverse laplace
for kn in 1:5, km in 0:kn
    U = sphericalmode(km,kn)(Λ, Φ) #forcing function
    G = laplace_sphere(U)
    Utest = laplace_sphere_inv(G)

    if norm(U - Utest) > 1e-10
        println("Inversion test fails on km: $km, kn: $kn, and $M x $N grid.")
        @assert norm(U - Utest) < 1e-10
    end
end
