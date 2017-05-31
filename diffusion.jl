# diffuse the field U in fourier space, with
# timestep dt
function plan_diffuse_spectral!(Uf, dt)
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

        Cs[mi] = lufact( A - dt*D )
    end

    function (Uf)
        for mi in 1:size(Uf,1)
            Uf[mi, :] = Cs[mi]\(As[mi]*Uf[mi, :])
        end
    end
end
