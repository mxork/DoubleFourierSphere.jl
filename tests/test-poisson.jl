using DoubleFourierSphere

X = 64
Y = 32

Λ, Φ = spheregrids(X,Y)

# for a variety of modes, check inversion property of 
# laplace - inverse laplace
# note the arbitrary constant we compensate for
for kn in 1:5, km in 0:kn
    U = sphericalmode(km,kn)(Λ, Φ) #forcing function
    G = laplace_sphere(U)
    Utest = laplace_sphere_inv(G)

    u00 = mean(U - Utest)

    if norm((U - Utest) .- u00) > 1e-10
        println("Inversion test fails on km: $km, kn: $kn, and $M x $N grid.")
        @assert norm(U - Utest) < 1e-10
    end
end
