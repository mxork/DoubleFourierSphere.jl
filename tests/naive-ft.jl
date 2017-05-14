#### Naive implementation

# can we can still FFT through λ direction,
# but we implement the naive algorithm for
# computing the fourier coefficients

# note the lack of an extra "f"

function ft_sphere(U)
    Nλ, Nφ = convert(Int, round(size(U, 1)/2)), size(U, 2)

    Um = conj( fft(U, 1) ) # fft out the longitude; Um(φ)
    Uf = zeros(Complex128, Nλ, Nφ)

    # interior grid
    Φs = Complex128[ π*(j+0.5)/Nφ for j in 0:Nφ-1]
    basis = zeros(Complex128, Nφ)

    # zero
    m = 0
    mi = m +1
    for n in 0:Nφ-1
        ni = n+1
        basis = cos(n*Φs)
        b = (n == 0 ? 1 : 2)

        Uf[mi,ni] = b * sum(Um[mi, :] .* basis) / Nφ
    end

    # odd
    for m in 1:2:Nλ-1
        mi = m+1
        for n in 1:Nφ
            ni = n
            basis = sin(n*Φs)
            c = (n == Nφ ? 1 : 2)

            Uf[mi,ni] = c* sum(Um[mi, :] .* basis) / Nφ
        end
    end

    # even
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

end
