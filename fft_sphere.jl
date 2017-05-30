export fft_sphere, ifft_sphere, plan_fft_sphere!, plan_ifft_sphere!

function fft_sphere(U)
    U = convert(Array{Complex128,2}, U)
    Uf = similar(U)
    Fs! = plan_fft_sphere!(U)
    Fs!(Uf, U)

    # hack
    M, Ms = zonal_modes(Uf)
    N, Ns0, Ns = meridional_modes(Uf)

    for ni in 1:size(Uf,2), mi in 1:size(Uf,1)
        m, n = Ms[mi], mi==1 ? Ns0[ni] : Ns[ni]
        if abs(m) >= div(M,2) || abs(n) >= div(N,2)
            Uf[mi, ni] = 0
        end
    end
    Uf
end

function ifft_sphere(Uf::Array{Complex128, 2})
    Fs! = plan_ifft_sphere!(Uf)
    U = similar(Uf)
    Fs!(U, Uf)
    U
end

function plan_fft_sphere!(U::Array{Complex128, 2})
    Um = similar(U)
    Fλ! = plan_fft_longitude!(U)
    Fφ! = plan_fft_latitude!(Um)

    function (Uf, U)
        Fλ!(Um, U)
        Fφ!(Uf, Um)
        Uf
    end
end

function plan_ifft_sphere!(Uf::Array{Complex128, 2})
    Um = similar(Uf)
    Fiφ! = plan_ifft_latitude!(Uf)
    Fiλ! = plan_ifft_longitude!(Um)

    function (U, Uf)
        Fiφ!(Um, Uf)
        Fiλ!(U, Um)
        U
    end
end

function plan_fft_longitude!(U::Array{Complex128,2})
    F = plan_fft(U,1)

    function (Um, U)
        A_mul_B!(Um, F, U)
    end
end

function plan_ifft_longitude!(Um::Array{Complex128,2})
    Fi = plan_ifft(Um,1)

    function (U, Um)
        A_mul_B!(U, Fi, Um)
        U
    end
end

function plan_fft_latitude!(Um::Array{Complex128, 2})
    M, Ms = zonal_modes(Um)

    # we're offset by half a Δφ
    # F(f(x-z))(k) = e^ikz F(f(x)), missing a 2 somewhere
    Nφ = size(Um,2)

    # two different sets of frequencies → different shifts
    Ns0 = 0:Nφ-1
    interior_grid_shift0 = exp( -1.0im * Ns0 * (π/Nφ/2))

    Ns = 1:Nφ
    interior_grid_shift = exp( -1.0im * Ns * (π/Nφ/2) + 1.0im*π/2)

    pole_scale = 1./ sin(latitude_interior_grid(Um))

    # FFts
    v = [ Um[1, :] ; -Um[1, end:-1:1] ]
    w = similar(v)
    Fs = plan_fft(v)

    return function (Uf, Um)
        # m zero
        # cosine transform
        v[:] = [ Um[1, :] ; Um[1, end:-1:1] ] #even

        A_mul_B!(w, Fs, v)

        Uf[1, :] = w[1:Nφ]
        Uf[1, :] .*= interior_grid_shift0

        # sine transform
        for mi in 2:size(Um, 1)
            m = Ms[mi]
            #pre multiply the m even modes

            v[1:size(Um,2)] = Um[mi, :]

            if iseven(m)
              v[1:size(Um,2)] .*= pole_scale
            end

            v[size(Um,2)+1:end] = -v[size(Um,2):-1:1]

            A_mul_B!(w, Fs, v)

            Uf[mi, :] = w[2:Nφ+1]      # drop the zeroth freq
            Uf[mi, :] .*= interior_grid_shift
        end

    end
end

function plan_ifft_latitude!(Uf::Array{Complex128, 2})
    M, Ms = zonal_modes(Uf)
    Nφ = size(Uf, 2)

    # unshifts: argument is -ve of shift
    Ns0 = 0:Nφ-1
    interior_grid_unshift0 = exp( 1.0im * Ns0 * (π/Nφ/2))

    Ns = 1:Nφ
    interior_grid_unshift = exp( +1.0im * Ns * (π/Nφ/2) - 1.0im*π/2)

    #
    pole_scale = sin(latitude_interior_grid(Uf))

    # plan fts
    v = [ 0.0 ; Uf[1, :]... ; conj(Uf[1, end:-1:2])... ]
    w = similar(v)
    u = view(w, 1:Nφ) # this is an undoubled view into w
    u2 = view(w, Nφ+1:2*Nφ)

    # INPLACE
    Fsi = plan_ifft(v)

    function (Um, Uf)
        # start going backwards
        # zero
        u[:] = Uf[1,:]
        u[:] .*= interior_grid_unshift0
        v[:] = [ u[:] ; 0.0 ; conj(u[end:-1:2]) ] #even

        A_mul_B!(w, Fsi, v)
        Um[1, :] = u[:]

        for mi in 2:size(Um, 1)
            m = Ms[mi]

            u[:] = Uf[mi,:]
            u2[:] = Uf[mi,:]
            u[:] .*= interior_grid_unshift
            v[:] = [ 0.0 ; u[1:end] ; conj(u[end-1:-1:1]) ]

            A_mul_B!(w, Fsi, v)

            # undo change of variable
            if iseven(m)
              u[:] .*= pole_scale
            end

            Um[mi, :] = u[:]
        end

        Um
    end
end

# just do the latitude transform
function ft_latitude(Um)
    # full transform output
    Nλ, Nφ = size(Um)
    Uf = zeros(Complex128, Nλ, Nφ)

    # interior grid for latitude
    Φs = Complex128[ π*(j+0.5)/Nφ for j in 0:Nφ-1]
    basis = zeros(Complex128, Nφ)

    # m zero
    m = 0
    mi = m +1
    for n in 0:Nφ-1
        ni = n+1
        basis = cos(n*Φs)
        b = (n == 0 ? 1 : 2)

        Uf[mi,ni] = b * sum(Um[mi, :] .* basis) / Nφ
    end

    # m odd
    for m in 1:2:Nλ-1
        mi = m+1
        for n in 1:Nφ
            ni = n # no offset, since n=0 is trivial
            basis = sin(n*Φs)
            c = (n == Nφ ? 1 : 2)

            Uf[mi,ni] = c* sum(Um[mi, :] .* basis) / Nφ
        end
    end

    # m even
    for m in 2:2:Nλ-1
        mi = m+1

        for n in 1:Nφ
            ni = n
            basis = sin(n*Φs)
            c = (n == Nφ ? 1 : 2)

            Uf[mi,ni] = c* sum( (Um[mi, :] ./ sin(Φs)) .* basis) / Nφ
        end
    end

    Uf
end


# Naive implementation of the transform
function ft_sphere(U)
    # quantities?
    Nλ, Nφ = size(U)

    # fft out the longitude; Um(φ)
    Um = fft(U, 1) 

    # full transform output
    Uf = zeros(Complex128, Nλ, Nφ)

    # interior grid for latitude
    Φs = Complex128[ π*(j+0.5)/Nφ for j in 0:Nφ-1]
    basis = zeros(Complex128, Nφ)

    # m zero
    m = 0
    mi = m +1
    for n in 0:Nφ-1
        ni = n+1
        basis = cos(n*Φs)
        b = (n == 0 ? 1 : 2)

        Uf[mi,ni] = b * sum(Um[mi, :] .* basis) / Nφ
    end

    # m odd
    for m in 1:2:Nλ-1
        mi = m+1
        for n in 1:Nφ
            ni = n # no offset, since n=0 is trivial
            basis = sin(n*Φs)
            c = (n == Nφ ? 1 : 2)

            Uf[mi,ni] = c* sum(Um[mi, :] .* basis) / Nφ
        end
    end

    # m even
    for m in 2:2:Nλ-1
        mi = m+1

        for n in 1:Nφ
            ni = n
            basis = sin(n*Φs)
            c = (n == Nφ ? 1 : 2)

            Uf[mi,ni] = c* sum( (Um[mi, :] ./ sin(Φs)) .* basis) / Nφ
        end
    end

    Uf
end

function ift_sphere(Uf)
    Nλ, Nφ = size(Uf)
    Um = zeros(Uf)

    Ns = 1:Nφ
    Ns0 = 0:(Nφ-1)
    Φs = Complex128[ π*(j+0.5)/Nφ for j in 0:Nφ-1]

    # cheong basis

    M, Ms = zonal_modes(Uf)

    # literally the worst way of grouping this
    for mi in 1:size(Uf,1)
        m = Ms[mi]
        for φi in 1:Nφ
            φ = Φs[φi]

            if m == 0
                Um[mi, φi] = sum( Uf[mi, :] .* cos( Ns0 * φ ) )
            elseif  isodd(m)
                Um[mi, φi] = sum( Uf[mi, :] .* sin( Ns * φ ) )
            else
                Um[mi, φi] = sum( Uf[mi, :] .* sin( Ns * φ ) )
            end
        end

        if iseven(m) && m != 0
            Um[mi, :] .*= sin(Φs) # undo the CoV
        end
    end

    # ifft meriodonal
    U = ifft(Um, 1)
end

function plan_ift_latitude!(Uf)
    Nλ, Nφ = size(Uf)
    Ns = 1:Nφ
    Ns0 = 0:(Nφ-1)
    Φs = latitude_interior_grid(Uf)
    pole_scale = sin(Φs)
    M, Ms = zonal_modes(Uf)

    C = [ cos(Ns0 * φ) for φ in Φs]
    S = [ sin(Ns * φ) for φ in Φs]

    function (Um, Uf)
        for mi in 1:size(Uf,1)
            m = Ms[mi]

            for φi in 1:Nφ
                φ = Φs[φi]

                if m == 0
                    Um[mi, φi] = sum( Uf[mi, :] .* C[φi] )
                elseif  isodd(m)
                    Um[mi, φi] = sum( Uf[mi, :] .* S[φi] )
                else
                    Um[mi, φi] = sum( Uf[mi, :] .* S[φi] )
                end
            end

            if iseven(m) && m != 0
                Um[mi, :] .*= sin(Φs) # undo the CoV
            end
        end

        Um ./= Nφ
    end
end

function ift_latitude!(Um, Uf)
    Nλ, Nφ = size(Uf)

    Ns = 1:Nφ
    Ns0 = 0:(Nφ-1)
    Φs = Complex128[ π*(j+0.5)/Nφ for j in 0:Nφ-1]

    # cheong basis

    M, Ms = zonal_modes(Uf)

    # literally the worst way of grouping this
    for mi in 1:size(Uf,1)
        m = Ms[mi]

        for φi in 1:Nφ
            φ = Φs[φi]

            if m == 0
                Um[mi, φi] = sum( Uf[mi, :] .* cos( Ns0 * φ ) )
            elseif  isodd(m)
                Um[mi, φi] = sum( Uf[mi, :] .* sin( Ns * φ ) )
            else
                Um[mi, φi] = sum( Uf[mi, :] .* sin( Ns * φ ) )
            end
        end

        if iseven(m) && m != 0
            Um[mi, :] .*= sin(Φs) # undo the CoV
        end
    end

    Um ./= Nφ
end
