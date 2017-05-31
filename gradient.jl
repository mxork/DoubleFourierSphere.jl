function gradient(G)
    Gf = fft_sphere(G)

    Gλf = similar(Gf)
    Gφf = similar(Gf)

    (plan_sinφ_gradient_spectral!(Gf))(Gλf, Gφf, Gf)

    Gλ = ift_sphere(Gλf)
    Gφ = ift_sphere(Gφf)

    pole_scale = spdiagm(sin(latitude_interior_grid(Gf)))

    Gλ *= pole_scale
    Gφ *= pole_scale

    Gλ, Gφ
end

# still need to divide after this guy
function plan_sinφ_gradient_spectral!(Uf::Array{Complex128, 2})
    M, Ms = zonal_modes(Uf)
    N, _ ,_ = meridional_modes(Uf)

    Dz = sinφdφ_zero!(spzeros(N, N))
    Do = sinφdφ_odd!(spzeros(N, N))
    De = sinφdφ_even!(spzeros(N, N))

    Dλ = dλ(Uf)

    function (Uλf, Uφf, Uf)
        Uλf[:] = Dλ * Uf

        for mi in 1:size(Uf, 1)
            m = Ms[mi]
            D = m==0     ? Dz :
                isodd(m) ? Do :
                           De

            Uφf[mi, :] = D * Uf[mi, :]
        end

        Uλf, Uφf
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


# Binary matrix (1s and 0s) corresponding to a truncation of high zonal
# frequencies near poles. The exact point of truncation is up for massage.
function latitude_truncation_mask(A)
    M, Ms = zonal_modes(A)
    Φs = latitude_interior_grid(A)
    upper_limit = (φ) -> min(M, 6+(M-6)*sin(φ))

    return [ abs(m) > upper_limit(φ) ? 0 : 1
             for m in Ms, φ in Φs]
end
