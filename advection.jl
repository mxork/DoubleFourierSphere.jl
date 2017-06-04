# advect scalar field Uf under velocity field Vxf, Vyf
function plan_advection_spectral!(Uf)
    advect_delta! = plan_advection_spectral_delta!(Uf)
    advect_delta_implicit! = plan_advection_spectral_implicit_delta!(Uf)

    dUf = similar(Uf)
    dUf_next = similar(Uf)

    function (Uf, Vx, Vy, dt)
        advect_delta!(dUf, Uf, Vx, Vy, dt)
        advect_delta_implicit!(dUf_next, Uf, Vx, Vy, dt)

        Uf .+= 0.4*dUf + 0.6*dUf_next
        Uf
    end
end

function plan_advection_spectral_delta!(Uf)
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

    function (dUf, Uf, Vx, Vy, dt)
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

        dUf .*= dt
        dUf
    end
end

#fixed point
function plan_advection_spectral_implicit_delta!(Uf)
    Uf_next = similar(Uf)
    advect_delta! = plan_advection_spectral_delta!(Uf)

    function (dUf_next, Uf, Vx, Vy, dt)
        Uf_next[:] = Uf
        advect_delta!(dUf_next, Uf_next, Vx, Vy, dt)

        num_iter = 0
        while norm(Uf_next - Uf - dUf_next) > 1e-12
            θ = max((0.98)^num_iter, 0.1)
            Uf_next[:] = θ*(Uf + dUf_next) + (1-θ)*Uf_next #relaxation
            advect_delta!(dUf_next, Uf_next, Vx, Vy, dt)
            num_iter += 1
            if num_iter > 200
                println("not converging")
                @assert num_iter < 200
            end
        end
        println("num_iter: $num_iter")

        dUf_next
    end
end

# secant method
# FIXME
function implicit_solve!(f, F1, F0, x1, x0)
    @assert x1 != x0

    F0[:] = f(x0)
    F1[:] = f(x1)

    while( norm(F1) > 1e-12 ) 
        x0[:] = x1 - F1 .* (x1 - x0) ./ (F1 - F0) 
        for i in 1:length(x0)
            if isnan(x0[i])
                x0[i] = x1[i] #cheat
            end
        end
        F0[:] = f(x0)

        x0, x1 = x1, x0 # swap
        F0, F1 = F1, F0 # swap
    end

    x1
end
