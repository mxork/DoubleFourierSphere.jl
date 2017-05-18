export fftsphere, ifftsphere, plan_fft_sphere, plan_ifft_sphere
export modes_from_indices, indices_from_modes

#### Naive implementation

# can we can still FFT through λ direction,
# but we implement the naive algorithm for
# computing the fourier coefficients

# even the for loops here are probably silly;
# could rewrite if I knew the conventions better.

# note the lack of an extra "f"
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

    # m zero
    m = 0
    mi = m+1
    for φi in 1:Nφ
        φ = Φs[φi]
        Um[mi, φi] = sum( Uf[mi, :] .* cos( Ns0 * φ ) )
    end

    # m odd
    for m in 1:2:Nλ-1
        mi = m+1
        for φi in 1:Nφ
            φ = Φs[φi]
            Um[mi, φi] = sum( Uf[mi, :] .* sin( Ns * φ ) )
        end
    end

    # m even
    for m in 2:2:Nλ-1
        mi = m+1
        for φi in 1:Nφ
            φ = Φs[φi]
            Um[mi, φi] = sum( Uf[mi, :] .* sin( Ns * φ ) )
        end

        Um[mi, :] = Um[mi, :] .* sin(Φs) # undo the transform
    end

    # ifft meriodonal
    U = ifft(Um, 1)
end

# export the slow ones for now;
# no "as" for export renaming in Julia
fftsphere = ft_sphere
ifftsphere = ift_sphere

# TODO typeswitch on U real/complex below here

# only uses dims of U here
# FIXME using an oddreflection doubles at least one of the modes, so need to compensate for that
function plan_fft_sphere(U)
    Um = Array{Complex128}(size(U))

    # we're offset by half a Δφ
    # F(f(x-z))(k) = e^ikz F(f(x)), missing a 2 somewhere
    Nφ = size(U,2)

    # two different sets of frequencies → different shifts
    Ns0 = 0:Nφ-1
    interior_grid_shift0 = exp( -1.0im * Ns0 * (π/Nφ/2))

    Ns = 1:Nφ
    interior_grid_shift = exp( -1.0im * Ns * (π/Nφ/2) + 1.0im*π/2)

    # interior grid pole condition scale for even modes
    Φs = Complex128[ π*(j+0.5)/Nφ for j in 0:Nφ-1]
    pole_scale = 1 ./ sin(Φs)

    # pre plan ffts
    F = plan_fft(U,1)

    # v is just a slot for shoving in an odd reflection # OPTIMIZE POINT
    v = [ Um[1, :] ; -Um[1, end:-1:1] ]
    w = similar(v)
    Fs = plan_fft(v)

    # returns a function which fourier sphere transforms
    # U (shadows the closed one) into Uf
    return function (Uf, U)
        # longitudinal
        # Um[:,:] = F*U
        A_mul_B!(Um, F, U) 

        # m zero
        # Uf[1, :] = C * Um[1, :]
        v[:] = [ Um[1, :] ; Um[1, end:-1:1] ] #even

        A_mul_B!(w, Fs, v)

        Uf[1, :] = w[1:Nφ]
        Uf[1, :] .*= interior_grid_shift0

        # m odd
        for mi in 2:2:size(U,1)
            v[:] = [ Um[mi, :] ; -Um[mi, end:-1:1] ]

            A_mul_B!(w, Fs, v)

            Uf[mi, :] = w[2:Nφ+1]      # drop the zeroth freq
            Uf[mi, :] .*= interior_grid_shift
        end

        # m even
        for mi in 3:2:size(U,1)
            Um[mi, :] .*= pole_scale
            v[:] = [ Um[mi, :] ; -Um[mi, end:-1:1] ]

            A_mul_B!(w, Fs, v)

            Uf[mi, :] = w[2:Nφ+1]
            Uf[mi, :] .*= interior_grid_shift
        end
    end
end

# in an ideal world, this would use the coupled inverse
# procedures of the forward FT, but 'cause we're doing other
# things, just split it off
# I should remove u from this code. it is confusing
function plan_ifft_sphere(Uf)
    Um = Array{Complex128}(size(Uf))

    Nφ = size(Uf, 2)

    # unshifts: argument is -ve of shift
    Ns0 = 0:Nφ-1
    interior_grid_unshift0 = exp( 1.0im * Ns0 * (π/Nφ/2))

    Ns = 1:Nφ
    interior_grid_unshift = exp( +1.0im * Ns * (π/Nφ/2) - 1.0im*π/2)

    #
    Φs = Complex128[ π*(j+0.5)/Nφ for j in 0:Nφ-1]
    pole_scale = sin(Φs)

    # plan fts
    Fi = plan_ifft(Um, 1)
    v = [ 0.0 ; Uf[1, :]... ; Uf[1, end:-1:2]... ]
    w = similar(v)
    u = view(w, 1:Nφ) # this is an undoubled view into w

    # INPLACE
    Fsi = plan_ifft(v)

    # swapped, using v as workspace
    function (U, Uf)
        # start going backwards

        # zero
        u[:] = Uf[1,:]
        u[:] .*= interior_grid_unshift0
        v[:] = [ u[:] ; 0.0 ; u[end:-1:2] ] #even

        A_mul_B!(w, Fsi, v)

        Um[1, :] = u[:]

        # odd
        for mi in 2:2:size(Um,1)
            u[:] = Uf[mi,:]
            u[:] .*= interior_grid_unshift
            v[:] = [ 0.0 ; u[:] ; -u[end-1:-1:1] ]

            A_mul_B!(w, Fsi, v)

            Um[mi, :] = u[:]
        end

        # even
        for mi in 3:2:size(Um,1)
            u[:] = Uf[mi,:]
            u[:] .*= interior_grid_unshift
            v[:] = [ 0.0 ; u[:] ; -u[end-1:-1:1] ]

            A_mul_B!(w, Fsi, v)

            u[:] .*= pole_scale
            Um[mi, :] = u[:]
        end

        # longitudinal
        A_mul_B!(U, Fi, Um)
    end
end

# following returns even/odd periodic completion of input vector
function evenrefl(vec)
    [vec[:] ; flipdim(vec,1)]
end

function oddrefl(vec)
    [vec[:] ; -flipdim(vec,1)]
end

# modes go from -M to M
# and from 0/1 : N-1/N
function modes_from_indices(M, N)
    return function (i,j)
      m = i <= M+1 ? i-1 : i-1-2*M
      n = m==0 ? j-1 : j
      return m, n
    end
end

function indices_from_modes(M,N)
    return function (m,n)
        i = m>=0 ? m+1 : 2*M + m
        j = m==0 ? n+1 : n

        return i,j
    end
end
