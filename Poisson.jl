export invertPoisson

# solve ΔU = G
function invertPoisson(G)
    X, Y = size(G)

    Gf = fftsphere(G)
    Uf = zeros(Gf) # same size, type


    # setup differentiation matrices
    # m=0
    D, A = DAzero( size(Gf) )

    # some hands on surgery to fix constant term.
    # god, I hope this gets compiles into something better
    Uf[1, 2:end] = D[:, 2:end] \ A * Gf[1, :] 
    Uf[1,1] = 0

    # m odd
    for mi in 2:2:size(Gf,1)
        m = mi-1
        DAodd!(D, A, m)
        Uf[mi, :] = D \ A * Gf[mi, :]
    end

    # m even
    for mi in 3:2:size(Gf,1)
        m = mi-1
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

    # two more equalities given in Cheong that bewilder me
    A[1,3] = -1/4
    D[1,3] = 1/2
end
