# advect scalar field Uf under velocity field Vxf, Vyf
function plan_advection_spectral!(Uf)
    M, Ms = zonal_modes(Uf)
    N, _ ,_ = meridional_modes(Uf)

    Az = sin2_zero!(spzeros(N, N))
    Ao = sin2_odd!(spzeros(N, N))
    Ae = sin2_even!(spzeros(N, N))

    Dz = sinφdφ_zero!(spzeros(N, N))
    Do = sinφdφ_odd!(spzeros(N, N))
    De = sinφdφ_even!(spzeros(N, N))

    Dλ = dλ(Uf)

    pole_scale = spdiagm(sin(latitude_interior_grid(Uf)))

    Uλf = similar(Uf)
    Uφf = similar(Uf)
    U_working = similar(Uf)

    F! = plan_fft_sphere!(Uf)
    Fi! = plan_ift_sphere!(Uf)
    # truncation would go here TODO

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

        # let's do this explicit
        Fi!(Uλ, Uλf)
        Fi!(Uφ, Uφf)

        Uλ .*= Vx
        Uφ .*= Vy

        dU = (-Uλ - Uφ) * pole_scale
        F!(dUf, dU)

        Uf .+= dt.*dUf
        Uf
    end
end
