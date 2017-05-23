using DoubleFourierSphere

X=128
Y=64
Gl, Gp = spheregrids(X, Y)

#check mode shows up in the right place
for M in 0:3
    for N in 0:3
        # these modes are degenerate
        if M!=0 && N==0
            continue
        end
        # slight not lining up
        Mi = M+1
        Ni = M == 0 ? N+1 : N

        U = fouriermode(M,N)(Gl, Gp)
        Uf = fft_sphere(U)

        # most of the energy should be here
        # note the conjugate could also show up
        if ! (abs(Uf[Mi, Ni]) == maximum(abs(Uf)))
            println("Energy concentration test failing on n=$N, m=$M")
            @assert !abs(Uf[Mi, Ni]) == maximum(abs(Uf))
        end
    end
end

