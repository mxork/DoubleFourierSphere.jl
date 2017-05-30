function gradient(G)
    Gf = fft_sphere(G)

    # apply sinφ grad in space and undo sinφ afties
    Gλf, Gφf = similar(Gf), similar(Gf)

    M, Ms = zonal_modes(Gf)
    Gλf[:] = 1.0im * diagm(Ms) * Gf
    (plan_sinφdφ!(Gf))(Gφf, Gf)

    Gλm = similar(Gf)
    ift_latitude!(Gλm,Gλf)
    Gλm .*= latitude_truncation_mask(Gλm)
    Gλ = ifft(Gλm, 1)

    Gφm = similar(Gf)
    ift_latitude!(Gφm,Gφf)
    Gφm .*= latitude_truncation_mask(Gφm)
    Gφ = ifft(Gφm, 1)

    # now account for sinφ
    Φ = latitude_interior_grid(G)
    pole_scale = 1 ./ sin(Φ)

    for λi in 1:size(Gλ, 1)
        Gλ[λi, :] .*= pole_scale
    end

    Gx = Gλ

    for λi in 1:size(Gφ, 1)
        Gφ[λi, :] .*= pole_scale
    end

    Gy = Gφ

    Gx, Gy
end

function plan_gradient!(G::Array{Complex128, 2})
    Gf = similar(G)

    Gλf = similar(Gf)
    Gφf = similar(Gf)

    Gλm = similar(Gf)
    Gφm = similar(Gf)

    F! = plan_fft_sphere!(G)

    # Fiφ! = plan_ifft_latitude!(Gf)
    Fiφ! = plan_ift_latitude!(Gf)
    Fiλ! = plan_ifft_longitude!(Gf)

    M, Ms = zonal_modes(G)
    Mscale = 1.0im * sparse(diagm(Ms))

    N, Ns0, Ns = meridional_modes(G)

    Sφ! = plan_sinφdφ!(Gf)

    lat_trunc_mask = latitude_truncation_mask(Gλm)

    pole_scale = 1 ./ sin( latitude_interior_grid(G))

    function (Gx, Gy, G)
        F!(Gf, G)

        # try taking out a big chunk
        for mi in 1:size(Gf,1)
            m = Ms[mi]
            if abs(m) > div(M,2)
                Gf[mi, :] = 0
            end
        end

        for ni in 1:size(Gf,2)
            n = Ns[ni]
            if abs(n) > div(N,2)
                Gf[:, ni] = 0
            end
        end

        Gλf[:] = Mscale * Gf
        Sφ!(Gφf, Gf)

        for ni in 1:size(Gf,2)
            n = Ns[ni]
            if abs(n) > div(N,2)
                Gφf[:, ni] = 0
            end
        end

        # ift_latitude!(Gλm, Gλf) # ah well FIXME
        # ift_latitude!(Gφm, Gφf)
        Fiφ!(Gλm, Gλf)
        Fiφ!(Gφm, Gφf)

        Gλm .*= lat_trunc_mask
        Gφm .*= lat_trunc_mask

        Fiλ!(Gx, Gλm)
        Fiλ!(Gy, Gφm)

        # FIXME
        for λi in 1:size(Gx, 1)
            Gx[λi, :] .*= pole_scale
            Gy[λi, :] .*= pole_scale
        end

        Gx, Gy
    end
end

# plans an application of the sinφdφ to a matrix shaped like Ψf
function plan_sinφdφ!(Ψf)
    M, Ms = zonal_modes(Ψf)

    S0 = sparse(sinφdφ_zero(size(Ψf, 2)))
    S = sparse(sinφdφ_odd(size(Ψf, 2)))
    T = sparse(sinφdφ_even_not_zero(size(Ψf, 2)))

    function (Uf, Ψf)
        for mi in 1:size(Ψf, 1)
            m = Ms[mi]
            if m==0
                Uf[mi, :] = S0 * Ψf[mi,:]
            elseif isodd(m) 
                Uf[mi, :] = S * Ψf[mi,:]
            elseif iseven(m)
                Uf[mi, :] = T * Ψf[mi,:]
            end
        end

        Uf
    end
end

# returns the matrix corresponding to the sinφdφ
# operator in frequency space
function sinφdφ_zero(N)
    return [ j==i+1 ?  -(i+1)/2:
             j==i-1 && i!=N ? (i-1)/2:
             0
             for i in 0:N-1, j in 0:N-1]
end

function sinφdφ_odd(N)
    return [ j==i+1 ?  -(i+1)/2:
             j==i-1 && i!=N ? (i-1)/2:
             0
             for i in 1:N, j in 1:N]
end

function sinφdφ_even_not_zero(N)
    return [ j==i+1 ?  -i/2 :
             j==i-1 && i!=N ? i/2 :
             0
             for i in 1:N, j in 1:N]
end

# Binary matrix (1s and 0s) corresponding to a truncation of high zonal
# frequencies near poles. The exact point of truncation is up for massage.
function latitude_truncation_mask(A)
    M, Ms = zonal_modes(A)
    Φs = latitude_interior_grid(A)
    upper_limit = (φ) -> min(M, 6+(M-6)*sin(φ))

    return [ abs(m) > upper_limit(φ) ? 0 : 1
             for m in Ms, φ in Φs]
end

function iter_step(Gnxt, G, Vx, Vy, Gx, Gy, GRAD)
    Gnxt[:] = G[:]
    GRAD(Gx, Gy, G)

    while norm(Gnxt - G - (Vx.*Gx + Vy.*Gy)) > 1e-10
        Gnxt[:] = G + (Vx.*Gx + Vy.*Gy)
        GRAD(Gx, Gy, Gnxt)
    end

    G[:] = Gnxt[:]
end
