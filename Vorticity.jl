using DoubleFourierSphere

# goal: solve vorticity on sphere using fts as above

# Z = vorticity grid
# F = coriolis forcing
function dζ(Z, F)
    # some helper cruft
    M = Int(round(size(Z,1)/2))
    Ms = [0:M-1 ; -M:-1]

    N = size(Z,2)
    Ns0 = 0:N-1
    Ns = 1:N

    # the meat begins
    H = Z + F
    Zf = fftsphere(Z)
    Ψf = laplace_sphere_inv_spectral(Zf)

    # Vf = im*Ψf
    Vf = 1.0im*diagm(Ms)*Ψf
    Uf = zeros(Vf)
    for j in 1:size(Uf,2), i in 1:size(Uf, 1)
        # index fixing
        m = Ms[i]
        n = m==0 ? Ns0[j] : Ns[j]

        Uf[i,j] = iseven(m) && m !=0 ? n*()/2
    end

end
