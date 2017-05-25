export laplace, laplace_inv, plan_laplace_inv!

# computes the laplacian of U over a sphere
function laplace(U)
    Uf = fftsphere(U)
    Gf = zeros(Uf)

    # m zero
    D, A = DAzero( size(Uf) )
    Gf[1, :] = A \ D * Uf[1, :] 

    M = convert(Int64, round(size(Gf,1)/2))
    Ms = [0:M-1 ; -M:-1]

    # m odd
    for mi in 2:2:size(Gf,1)
        m = Ms[mi]
        DAodd!(D, A, m)
        Gf[mi, :] = A \ D * Uf[mi, :]
    end

    # m even
    for mi in 3:2:size(Gf,1)
        m = Ms[mi]
        DAeven!(D, A, m)
        Gf[mi, :] = A \ D * Uf[mi, :]
    end

    ifftsphere(Gf)
end

function laplace_inv(G)
    Gf = fft_sphere(G)

    Uf = laplace_inv_spectral(Gf)

    ifft_sphere(Uf)
end

# like laplace_inv but from spectral to spectral
function laplace_inv_spectral(Gf)
    Uf = similar(Gf)
    (plan_laplace_inv!(Gf))(Uf, Gf)
end

# these need reworking - all the DA functions
# throw a lot of garbage
# ugly as hell

function plan_laplace_inv!(Gf)
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


# Linear systems corresponding to laplacian in frequency space. 
# Infers the size of the target vector from input arguments, writes result to
# arguments.

# TODO consider using a sparse matrix set up. just call sparse(A), sparse(D)
# TODO consider factoring these into true tridiagonals, of halfsize (would have)
# to hack the stride of the matrices to get the right entries w/o a copy.

function DAzero!(D, A)
    M = 0 # by definition
    @assert size(D) == size(A)

    # some bounds trickery cuz offset here
    for ji in 1:size(D,2), ii in 1:size(D,1)
        i, j = ii-1, ji-1

        D[ii,ji] =
            i==j   ? -(i^2 + 2*M^2)/2 :
            j==i-2 ? (i-1)*(i-2)/4 :
            j==i+2 ? (i+1)*(i+2)/4 :
            0

        D[1,1] = 0 # from comment in Cheong
        D[3,1] = 0 # ditto

        A[ii,ji] =
            i==j   ? 1/2 :
            j==i-2 ? -1/4 :
            j==i+2 ? -1/4 :
            0
    end

    # from comment in cheong, §2.4
    A[2,2] = 1/4
    A[3,1] = -1/2

    # two more equalities given in Cheong
    A[1,3] = -1/4
    D[1,3] = 1/2
end

function DAodd!(D, A, M)
    @assert isodd(M)
    @assert size(D) == size(A)

    # apparently, inner loop should be first in Julia
    # @inbounds
    for j in 1:size(D,2), i in 1:size(D,1)
        D[i,j] =
            i==j   ? -(i^2 + 2*M^2)/2 :
            j==i-2 ? (i-1)*(i-2)/4 :
            j==i+2 ? (i+1)*(i+2)/4 :
            0

        A[i,j] =
            i==j   ? 1/2 :
            j==i-2 ? -1/4 :
            j==i+2 ? -1/4 :
            0
    end

    A[1,1] = 3/4;
end

function DAeven!(D, A, M)
    @assert iseven(M)
    @assert size(D) == size(A)

    # @inbounds
    for j in 1:size(D,2), i in 1:size(D,1)
        D[i,j] =
            i==j   ? -(i^2 + 2*M^2)/2 :
            j==i-2 ? i*(i-1)/4 :
            j==i+2 ? i*(i+1)/4 :
            0

        A[i,j] =
            i==j   ? 1/2 :
            j==i-2 ? -1/4 :
            j==i+2 ? -1/4 :
            0
    end

    A[1,1] = 3/4;
end

