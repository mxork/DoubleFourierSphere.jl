using DoubleFourierSphere

# time step. Assume constant coriolis (F)
# forward euler, baby
function timestep(Z0, F, Δt, n_steps)
    for i in 1:n_steps
        Z0 += real(dζ(Z0, F)*Δt)
    end

    Z0
end

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
    Zf = fftsphere(Z)
    Ff = fftsphere(F) # Terrible notation, but consistent
    Hf = Zf + Ff
    Ψf = laplace_sphere_inv_spectral(Zf)

    # Vf = im*Ψf
    Vf = 1.0im*diagm(Ms)*Ψf

    Uf = zeros(Vf)
    A = sinφdφ(N)
    B = sinφdφ_even_not_zero(N)
    for mi in 1:size(Ψf, 1)
        m = Ms[mi]
        Uf[mi, :] = iseven(m) && m != 0 ?
                    B * Ψf[mi,:] :
                    A * Ψf[mi,:]
    end

    # similar for η derivatives
    Hλf = 1.0im*diagm(Ms)*Hf
    Hφf = zeros(Hλf)
    for mi in 1:size(Hf, 1)
        m = Ms[mi]
        # extra -ve sign intentional; see defs
        # otherwise, same matrices
        Hφf[mi, :] = iseven(m) && m != 0 ?
                    B * Hf[mi,:] :
                    A * Hf[mi,:]
    end

    # now we have the derived quantites, pop back
    # into normal space and do the PW multiplication
    # these can be FASTER
    U, V = ifftsphere(Uf), ifftsphere(Vf)
    Hλ, Hφ = ifftsphere(Hλf), ifftsphere(Hφf)

    # §3-8, Cheong
    Φs = Complex128[ π*(j+0.5)/N for j in 0:N-1]
    S = diagm( 1 ./ sin(Φs) )
    X = (U .* Hλ) *S*S
    Y = (V .* Hφ) *S

    Xf, Yf = fftsphere(X), fftsphere(Y)

    # and here is dζ/dt (spectral, actually)
    ifftsphere(-Xf - Yf)
end

# returns the matrix corresponding to the sinφdφ
# operator in frequency space
# actually, the -ve sinφdφ operator
# TODO truncate Nth mode entirely since out of range?
# ANSWER Yes, and the N=0 mode is 0, so don't worry about it
function sinφdφ(N)
    return [ j==i+1 ?  (i+1)/2:
             j==i-1 && i!=N ? -(i-1)/2:
                             0
        for i in 1:N, j in 1:N]
end

function sinφdφ_even_not_zero(N)
    return [ j==i+1 ?  i/2 :
             j==i-1 && i!=N ? -i/2 :
             0
        for i in 1:N, j in 1:N]
end
