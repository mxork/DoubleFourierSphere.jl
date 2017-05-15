export laplace_sphere, laplace_sphere_inv

#TODO clean up the interfaces after screwing with FFT code, esp. dimension

# computes the laplacian of U over a sphere
# FIXME there's damping on the m0, neven spherical modes?
function laplace_sphere(U)
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

# solve ΔU = G
function laplace_sphere_inv(G)
    Gf = fftsphere(G)
    Uf = zeros(Gf) # same size, type

    M = convert(Int64, round(size(Gf,1)/2))
    Ms = [0:M-1 ; -M:-1]

    # setup differentiation matrices
    # m=0
    J = size(Gf,2)
    D, A = DAzero( (J,J) )

    # some hands on surgery to fix constant term.
    # TODO this constant term is trickier than expected,
    # since m0, neven modes have a non-zero (0,0) component
    # god, I hope this gets compiles into something better
    Uf[1, 2:end] = D[:, 2:end] \ A * Gf[1, :] 
    Uf[1,1] = 0

    # m odd
    for mi in 2:2:size(Gf,1)
        m = Ms[mi]
        DAodd!(D, A, m)
        Uf[mi, :] = D \ A * Gf[mi, :]
    end

    # m even
    for mi in 3:2:size(Gf,1)
        m = Ms[mi]
        DAeven!(D, A, m)
        Uf[mi, :] = D \ A * Gf[mi, :]
    end

    ifftsphere(Uf)
end

# each of these has two forms, one of which accepts a result buffer
# TODO consider using a sparse matrix set up. just call sparse(A), sparse(D)
function DAodd(dims, M)
    D, A = zeros(dims), zeros(dims)
    DAodd!(D, A, M)

    D, A
end

function DAodd!(D, A, M)
    @assert isodd(M)
    @assert size(D) == size(A)

    # apparently, inner loop should be row in Julia
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

function DAeven(dims, M)
    D, A = zeros(dims), zeros(dims)
    DAeven!(D, A, M)

    D, A
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

function DAzero(dims)
    D, A = zeros(dims), zeros(dims)
    DAzero!(D, A)
    D, A
end

# These are giving me trouble. I need to sit down and do algebra.
function DAzero!(D, A)
    M = 0 # by definition
    @assert size(D) == size(A)

    # @inbounds
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
