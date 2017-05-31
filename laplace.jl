export laplace, laplace_inv, plan_laplace_inv!

# computes the laplacian of U over a sphere
function laplace(U)
    Ufbig = fft_sphere(U)
    Uf = trunc_modes!( Array{Complex128,2}(div(size(U,1),2), div(size(U,2), 2) ), Ufbig)
    # Uf = fft_sphere(U)
    Gf = zeros(Uf)

    D = zeros( size(Uf,2), size(Uf,2))
    A = zeros( size(Uf,2), size(Uf,2))

    # m zero
    DAzero!(D, A)
    Gf[1, :] = A \ (D * Uf[1, :] )

    M, Ms = zonal_modes(Gf)

    for mi in 1:size(Gf,1)
        m = Ms[mi]
        if isodd(m)
          DAodd!(D, A, m)
        else
          DAeven!(D, A, m)
        end
        Gf[mi, :] = A \ (D * Uf[mi, :])
    end

    Gfbig = expand_modes!(similar(Ufbig), Gf)
    ift_sphere(Gfbig)
    # ift_sphere(Gfbig)
end

function laplace_inv(G)
    Gf = fft_sphere(G)

    Uf = laplace_inv_spectral(Gf)

    ift_sphere(Uf)
end

# like laplace_inv but from spectral to spectral
function laplace_inv_spectral(Gf)
    Uf = similar(Gf)
    (plan_laplace_inv!(Gf))(Uf, Gf)
end

# these need reworking - all the DA functions
# throw a lot of garbage
# ugly as hell

function plan_laplace!(G)
    F! = plan_fft_sphere!(G)
    Fiφ! = plan_ift_latitude!(G) # FIXME
    Fiλ! = plan_ifft_longitude!(G)


    Gf = similar(G)
    Gfs =  Array{Complex128,2}(div(size(G,1),2), div(size(G,2), 2)+2) # +2 for shift

    Lf! = plan_laplace_spectral!(Gfs)

    Ufs = similar(Gfs)
    Uf = similar(Gf)
    Um = similar(Gf)

    lat_trunc_mask = latitude_truncation_mask(Uf)

    function (U, G)
        F!(Gf, G)
        trunc_modes!(Gfs, Gf)
        Lf!(Ufs, Gfs)
        Ufs[:, end-1:end] = 0 # zero the top +2 mode
        expand_modes!(Uf, Ufs)
        Fiφ!(Um, Uf)
        Um .*= lat_trunc_mask
        Fiλ!(U, Um)
    end
end

function plan_laplace_spectral!(Gf)
    M, Ms = zonal_modes(Gf)
    N, Ns0, Ns = meridional_modes(Gf)

    # gonna be backsolving D, so LU it
    Ds = Array{ SparseMatrixCSC{Float64, Int64}, 1 }( size(Ms) )
    As = Array{ SparseArrays.UMFPACK.UmfpackLU, 1 }( size(Ms) )

    D = Array{Float64, 2}(N, N)
    A = similar(D)

    for mi in 1:size(Ms,1)
        m = Ms[mi]

        if m == 0
            DAzero!(D, A)
            D[1,1] = 1.0
        elseif isodd(m)
            DAodd!(D, A, m)
        else
            DAeven!(D, A, m)
        end

        # pretreat the matrices a little to make 'em faster
        Ds[mi] = sparse(D)
        As[mi] = lufact(sparse(A))
    end

    function (Uf, Gf)
        # everything else
        for mi in 1:size(Ms,1)
            D, A = Ds[mi], As[mi]
            Uf[mi, :] = A \ (D * Gf[mi, :])
        end
        Uf[:, end] = 0
        Uf
    end
end

function plan_laplace_spectral_inv!(Gf)
    M, Ms = zonal_modes(Gf)
    N, Ns0, Ns = meridional_modes(Gf)

    # precompute all the matrices; uses
    # a lot of space, but we're going to be using them
    # often

    # gonna be backsolving D, so LU it
    Ds = Array{ SparseArrays.UMFPACK.UmfpackLU, 1 }( size(Ms) )
    As = Array{ SparseMatrixCSC{Float64, Int64}, 1 }( size(Ms) )

    D = Array{Float64, 2}(N, N)
    A = similar(D)

    for mi in 1:size(Ms,1)
        m = Ms[mi]

        if m == 0
            DAzero!(D, A)
            D[1,1] = 1.0
        elseif isodd(m)
            DAodd!(D, A, m)
        else
            DAeven!(D, A, m)
        end

        # pretreat the matrices a little to make 'em faster
        Ds[mi] = lufact(sparse(D))
        As[mi] = sparse(A)
    end

    function (Uf, Gf)
        # zero is a special case
        # FIXME this back solve agrees, right? can't slice on an LU-sparse
        Uf[1, :] = Ds[1] \ (As[1] * Gf[1, :] )
        Uf[1,1] = 0

        # everything else
        for mi in 2:size(Ms,1)
            D, A = Ds[mi], As[mi]
            Uf[mi, :] = D \ (A * Gf[mi, :])
        end

        Uf
    end
end
