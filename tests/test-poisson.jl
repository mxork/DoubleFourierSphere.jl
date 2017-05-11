using DoubleFourierSphere
using PyPlot

# somewhat convoluted check to make sure G looks like U
# (using eigenfunction)
function isScalarMultipleIsh(G, U)
    scalar, set = 0.0, false

    for j in 1:size(G,2), i in 1:size(G,1)
        if abs(G[i,j]) > 1e-10 
                scalar = G[i,j] / U[i,j]
                break
        end
    end

    norm( G - scalar*U ) < 1e-10
end

# poisson equation
M = 512
N = 256

# wave numbers
Λ, Φ = spheregrids(M,N)

G, U = 0, 0

# currently breaking on km=0, AND (kn, km) = (4,4)
for kn in 1:5, km in 0:kn
    G = sphericalmode(km,kn)(Λ, Φ) #forcing function
    U = invertPoisson(G)

    if !isScalarMultipleIsh(G, U)
        println("kn: $kn, km: $km")
        plot( U[10, :] )
        break
        # @assert isScalarMultipleIsh(G, U)
    end
end
