using DoubleFourierSphere

X = 64
Y = 32

Λ, Φ = spheregrids(X,Y)

# check eigenfunction-ness for non-constant modes
for kn in 1:5, km in 0:kn
    U = sphericalmode(km,kn)(Λ, Φ) #forcing function
    G = laplace(U)

    scale = norm(U) / norm(G)

    Gmod = scale * G

    # resolve the sign ambiguity in a crude, crude way
    if norm(U - Gmod) > 1e-10 && norm(U + Gmod) > 1e-10
        println("Eigenfunction test fails on km: $km, kn: $kn, and $M x $N grid.")
        @assert norm(U - Gmod) < 1e-10 || norm(U + Gmod) < 1e-10
    end
end
