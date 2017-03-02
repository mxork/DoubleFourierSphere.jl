#TODO brainsweat
module FFTsphere

export fftsphere, ifftsphere

# returns m-n fourier coefficients of U,
# looks good
function fftsphere(U)
    Nλ, Nφ = convert(Int, round(size(U, 1)/2)), size(U, 2)

    Um = conj( fft(U, 1) ) # fft out the longitude; Um(φ)
    Uf = zeros(Complex128, Nλ, Nφ)

    # TODO understand wtf, from Muraki
    Ts = (π/2) * (0:Nφ-1) / Nφ # first quarter circle?
    Vshiftcos = exp( -1.0im * Ts ) / (2*Nλ *Nφ) # denom just normalization, undone later
    Vshiftsin = -Vshiftcos * exp(1.0im*π/2 * (1 - 1/Nφ))
    Vshiftcos[1] /= 2
    Vsin = -1./sin((π/Nφ/2)*(1:2:2*Nφ-1));

    # m=0; read this as taking an even reflection,
    # then taking only the first half of the frequencies
    Uf[1,:] = real( fft( evenrefl(Um[1,:]) )[1:Nφ] .* Vshiftcos )

    # m even; read this as comping for latitude, taking even reflection,
    # then ...
    # Muraki has a tricky with the -ve sign rolled into the Vsin
    # making the reflection work out nicely
    # TODO make sense of the +1 shift of the frequencies,
    # probably something with aliasing... but why is 0th entry bogus?
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
function ifftsphere(Uf)
    Nλ, Nφ = size(Uf, 1), size(Uf, 2)

    Vunshiftcos    =  exp((1.0im*pi/Nφ/2)*(0:2*Nφ-1))*(2*Nφ*Nλ)
    Vunshiftsin    =  1.0im*Vunshiftcos
    Vunshiftcos[1] *= 2
    Vunsin         =  -sin((pi/Nφ/2)*(1:2:2*Nφ-1))

    Um = zeros(Complex128, size(Uf))
    U = zeros(Float64, 2*Nλ, Nφ)

    # TODO why do we 0 out the (N/2)th frequency but not the 0th?
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

end
