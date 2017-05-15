export fftsphere, ifftsphere

# returns m-n fourier coefficients of U,
# looks good
function fft_sphere(U)
    Nλ, Nφ = convert(Int, round(size(U, 1)/2)), size(U, 2)

    Um = conj( fft(U, 1) ) # fft out the longitude; Um(φ)
    Uf = zeros(Complex128, Nλ, Nφ)

    # TODO understand wtf, from Muraki
    # pretty sure this is all just for dealing with interior grid, plus some other
    # ops which got rolled in
    # F(f(x-z))(k) = e^ikx F(f(x)), missing a 2 somewhere
    Ts = (π/2) * (0:Nφ-1) / Nφ # first quarter circle?
    Vshiftcos = exp( -1.0im * Ts )
    Vshiftsin = -Vshiftcos * exp(1.0im*π/2 * (1 - 1/Nφ)) # take this on faith
    Vsin = -1./sin((π/Nφ/2)*(1:2:2*Nφ-1));

    # m=0; read this as taking an even reflection,
    # then taking only the first half of the frequencies
    Uf[1,:] = real( fft( evenrefl(Um[1,:]) )[1:Nφ] .* Vshiftcos )

    # m even
      # TODO make sense of the +1 shift of the frequencies,
      # probably something with aliasing... but why is 0th entry bogus?
      # => likely just to get it to line up nicely with the definition of Vshiftsin
    for j in 3:2:Nλ
        Um[j, :] = Um[j, :] .* Vsin
        Uf[j, :] = fft( oddrefl(Um[j, :]) )[2:Nφ+1] .* Vshiftsin
    end

    # m odd
    for j in 2:2:Nλ
        Uf[j, :] = fft( oddrefl(Um[j,:]) )[2:Nφ+1] .* Vshiftsin
    end

    Uf
end

# returns λ-φ values of U
function ifft_sphere(Uf)
    Nλ, Nφ = size(Uf, 1), size(Uf, 2)

    Vunshiftcos    =  exp((1.0im*pi/Nφ/2)*(0:2*Nφ-1))
    Vunshiftsin    =  1.0im*Vunshiftcos
    Vunsin         =  -sin((pi/Nφ/2)*(1:2:2*Nφ-1))

    Um = zeros(Complex128, size(Uf))
    U = zeros(Float64, 2*Nλ, Nφ)

    # TODO why do we 0 out the (N/2)th frequency but not the 0th? Answer:
    Um[1, :] = ifft(  [ Uf[1, :] ; 0.0 ; -Uf[1, Nφ:-1:2] ] .* Vunshiftcos )[1:Nφ]

    # m even
    for j in 3:2:Nλ
        Um[j, :] = ifft( [0.0 ; Uf[j, :] ; Uf[j, Nφ-1:-1:1]] .* Vunshiftsin )[1:Nφ] .* Vunsin
    end

    # m odd
    for j in 2:2:Nλ
        Um[j, :] = ifft( [0.0 ; Uf[j, :] ; Uf[j, Nφ-1:-1:1]] .* Vunshiftsin )[1:Nφ]
    end

    # conjugate flips here because of impl difference
    U = real( ifft([conj(Um) ; zeros(Nφ)' ; Um[Nλ:-1:2,:]], 1) )
end

# following returns even/odd periodic completion of input vector
function evenrefl(vec)
    [vec[:] ; flipdim(vec,1)]
end

function oddrefl(vec)
    [vec[:] ; -flipdim(vec,1)]
end

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
