using DoubleFourierSphere

M = 64
N = 32

Λ, Φ = spheregrids(M,N)

# check eigenfunction-ness for non-constant modes
for kn in 1:5, km in 0:kn
    U = sphericalmode(km,kn)(Λ, Φ) #forcing function
    G = laplace_sphere(U)

    scale = norm(U) / norm(G)

    Gmod = scale * G

    # resolve the sign ambiguity in a crude, crude way
    if norm(U - Gmod) > 1e-10 && norm(U + Gmod) > 1e-10
        println("Eigenfunction test fails on km: $km, kn: $kn, and $M x $N grid.")
        @assert norm(U - Gmod) < 1e-10 || norm(U + Gmod) < 1e-10
    end
end
