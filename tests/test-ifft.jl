using DoubleFourierSphere

X=32
Y=16
Gl, Gp = spheregrids(X, Y)

# idempotence
for  mode in [fouriermode, sphericalmode]
	for M in 0:3
		for N in 0:M
            U = mode(N,M)(Gl, Gp)
            Utest = ifft_sphere(fft_sphere(U))
            @assert maximum(abs(U-Utest)) < 1e-14
		end
	end
end
