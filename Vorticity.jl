using DoubleFourierSphere

# time step. Assume constant coriolis (F)
# forward euler, baby
function timestep(Z0, F, Δt, n_steps)
    for i in 1:n_steps
        Z0 += real(dζ(Z0, F)*Δt)
    end

    Z0
end

function plan_dζf!(Zf)
    # ALLOCATION
    M, Ms = zonal_modes(Zf)
    N, Ns0, Ns = meridional_modes(Zf)
    Mscale = diagm(sparse(Ms))

    F_into! =  plan_fft_sphere!(Zf)
    Fi_into! = plan_ifft_sphere!(Zf)
    laplace_inv_into! = plan_laplace_inv!(Zf)
    sinφdφ_into! = plan_sinφdφ!(Zf)
    XY_into! = plan_calculate_XY!(Zf, F_into!, Fi_into!)

    Uf = similar(Zf)
    Vf = similar(Zf)

    Hf = similar(Zf)
    Ψf = similar(Zf)

    Hλf = similar(Zf)
    Hφf = similar(Zf)

    Xf = similar(Zf)
    Yf = similar(Zf)

    lat_trunc_mask = latitude_truncation_mask(Zf) 

    # throw in a computation for the CFL number
    # somewhere in here TODO
    function (dζf, Zf, Ff)
        Hf[:] = Zf + Ff
        laplace_inv_into!(Ψf, Zf)

        A_mul_B!(Hλf, Mscale, Hf)
        sinφdφ_into!(Hφf, Hf)
        Hφf *= -1.0

        A_mul_B!(Vf, Mscale, Ψf)
        sinφdφ_into!(Uf, Ψf)

        XY_into!(Xf, Yf, Uf, Vf, Hλf, Hφf)

        dζf[:] = -Xf -Yf
        dζf .*= lat_trunc_mask
    end
end

# goal: solve vorticity on sphere using fts as above

# Z = vorticity grid
# F = coriolis forcing
function dζ(Z, F)

    # MEAT
    Zf = fftsphere(Z)
    Ff = fftsphere(F) # Terrible notation, but consistent
    Hf = Zf + Ff
    Ψf = laplace_sphere_inv_spectral(Zf)

    # calculate spectral coefficients for
    # U, V and η derivatives

    # U,V
    A_mul_B!(Vf, Mscale, Ψf)
    Vf[:,:] *= 1.0im

    sinφdφ_into!(Uf, Ψf)

    # η
    A_mul_B!(Hλf, Mscale, Hf)
    Hλf[:,:] *= 1.0im

    sinφdφ_into!(Hφf, Hf)
    Hφf *= -1.0 # conventional -ve


    # now we have the derived quantites, pop back
    # into normal space and do the PW multiplication
    XY_into!(Xf, Yf, Uf, Vf, Hλf, Hφf)

    dζf[:,:] = -Xf -Yf

    # filter out the higher frequencies near the poles
    # §3.2
    dζf .*= lat_trunc_mask

    # and here is dζ/dt (spectral, actually)
    ifftsphere(-Xf - Yf)
end

function plan_sinφdφ!(Ψf)
    M, Ms = zonal_modes(Ψf)

    S = sinφdφ(size(Ψf, 2))
    T = sinφdφ_even_not_zero(size(Ψf, 2))

    function (Uf, Ψf)
        for mi in 1:size(Ψf, 1)
            m = Ms[mi]
            if isodd(m) || m == 0
                Uf[mi, :] = S * Ψf[mi,:]
            else
                Uf[mi, :] = T * Ψf[mi,:]
            end
        end

        Uf
    end
end

# this is a little specific, but oh well
# pass in Uf, and a plan for FFT/IFFT
function plan_calculate_XY!(Uf, F!, Fi!)
    U = similar(Uf)
    V = similar(Uf)
    Hλ = similar(Uf)
    Hφ = similar(Uf)

    # alias
    X = U
    Y = V

    N = size(U, 2)
    Φs = Complex128[ π*(j+0.5)/N for j in 0:N-1]
    S = diagm( 1 ./ sin(Φs) )

    function (Xf, Yf, Uf, Vf, Hλf, Hφf)
        Fi!(U, Uf)
        Fi!(V, Vf)
        Fi!(Hλ, Hλf)
        Fi!(Hφ, Hφf)

        X .*= Hλ; X *= S; X *= S
        Y .*= Hφ; Y *= S

        F!(Xf, X)
        F!(Yf, Y)

        Xf, Yf
    end
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


function latitude_truncation_mask(A)
    M, Ms = zonal_modes(A)
    Φs = latitude_interior_grid(A)
    upper_limit = (φ) -> min(M, 6+(M-6)*sin(φ))

    return [ abs(m) > upper_limit(φ)? 0 : 1
             for m in Ms, φ in Φs]
end
