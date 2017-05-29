using DoubleFourierSphere

X=128
Y=64
Gl, Gp = spheregrids(X, Y)

# idempotence
for  mode in [fouriermode, sphericalmode_complex]
    for M in 0:div(X,2)-1
        for N in 0:M
            U = mode(N,M)(Gl, Gp)
            Utest = ifft_sphere(fft_sphere(U))

            if maximum(abs(U-Utest)) >= 1e-10
                modes = mode ==fouriermode ? "fourier" : "spherical"
                println("Idempotence failing on m=$M, n=$N, $modes")
                @assert maximum(abs(U-Utest)) < 1e-10
            end
        end
	  end
end
