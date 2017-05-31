# wave the scalar field U, in frequency space
# using fixed time step dt, with Uf0 = Uf at t-1
function plan_wave_spectral!(Uf, dt)
    M, Ms = zonal_modes(Uf)
    N, _ ,_ = meridional_modes(Uf)

    Ds = Array{ SparseMatrixCSC{Float64, Int64}, 1 }( size(Ms) )
    As = Array{ SparseMatrixCSC{Float64, Int64}, 1 }( size(Ms) )
    Cs = Array{ SparseArrays.UMFPACK.UmfpackLU, 1 }( size(Ms) )

    for mi in 1:size(Ms,1)
        m = Ms[mi]

        Ds[mi] = spzeros(N, N)
        As[mi] = spzeros(N, N)
        D, A = Ds[mi], As[mi]

        if m == 0
            DAzero!(D, A)
        elseif isodd(m)
            DAodd!(D, A, m)
        else
            DAeven!(D, A, m)
        end

        Cs[mi] = lufact( A - dt*dt*D )
    end

    U_working = similar(Uf)

    function (Uf, Uf0)
        U_working[:] = 2*Uf - Uf0

        for mi in 1:size(Uf,1)
            Uf0[mi, :] = Cs[mi] \ (As[mi]*U_working[mi, :])
        end

        # this is a weird convenction CONFUSING
        Uf0, Uf
    end
end
