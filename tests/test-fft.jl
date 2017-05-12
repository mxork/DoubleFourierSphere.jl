    using DoubleFourierSphere

X=128
Y=64
Gl, Gp = spheregrids(X, Y)

# check mode shows up in the right plac
for M in 0:3
    for N in 0:3
        U = fouriermode(M,N)(Gl, Gp)
        Uf = fftsphere(U)

        println("M: $M, N: $N")
        # most of the energy should be here
        @show abs(Uf[M+1, N+1]) == maximum(abs(Uf))

        Uf[M+1, N+1] = 0
        @show norm(Uf) < 1.0
        println()
    end
end
