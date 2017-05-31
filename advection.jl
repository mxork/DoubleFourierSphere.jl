# advect scalar field Uf under velocity field Vxf, Vyf
function plan_advection_spectral!(Uf)
    M, Ms = zonal_modes(Uf)
    N, _ ,_ = meridional_modes(Uf)

    # Az = sin2_zero!(spzeros(N, N))
    # Ao = sin2_odd!(spzeros(N, N))
    # Ae = sin2_even!(spzeros(N, N))

    Dz = sinφdφ_zero!(spzeros(N, N))
    Do = sinφdφ_odd!(spzeros(N, N))
    De = sinφdφ_even!(spzeros(N, N))

    Dλ = dλ(Uf)

    pole_scale = spdiagm(sin(latitude_interior_grid(Uf)))

    Uλf = similar(Uf)
    Uφf = similar(Uf)

    F! = plan_fft_sphere!(Uf)
    Fiφ! = plan_ift_latitude!(Uf)
    Fiλ! = plan_ifft_longitude!(Uf)
    # truncation would go here TODO

    Uλm = similar(Uf) 
    Uφm = similar(Uf)

    lat_trunc_mask = latitude_truncation_mask(Uf)

    Uλ = similar(Uf) 
    Uφ = similar(Uf)

    dU = similar(Uf)
    dUf = similar(Uf)

    function (Uf, Vx, Vy, dt)
        Uλf[:] = Dλ * Uf

        for mi in 1:size(Uf, 1)
            m = Ms[mi]
            D = m==0     ? Dz :
                isodd(m) ? Do :
                           De

            Uφf[mi, :] = D * Uf[mi, :]
        end

        # and the top meridional modes are garbage
        Uλf[:, end] = 0
        Uφf[:, end] = 0

        # let's do this explicit
        Fiφ!(Uλm, Uλf)
        Fiφ!(Uφm, Uφf)

        Uλm .*= lat_trunc_mask
        Uφm .*= lat_trunc_mask

        Fiλ!(Uλ, Uλm)
        Fiλ!(Uφ, Uφm)

        Uλ .*= Vx
        Uφ .*= Vy

        dU = (-Uλ - Uφ) * pole_scale
        F!(dUf, dU)

        Uf .+= dt.*dUf
        Uf
    end
end
