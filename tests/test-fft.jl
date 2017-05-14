using DoubleFourierSphere

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
        Uf = fftsphere(U)

        # most of the energy should be here
        @assert abs(Uf[Mi, Ni]) == maximum(abs(Uf))

        Uf[Mi, Ni] = 0
        @assert norm(Uf) < 1.0
    end
end

