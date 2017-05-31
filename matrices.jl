# returns the matrix corresponding to the sinφdφ
# operator in frequency space

# this guy operates on whole enchilada
function dλ(Uf)
    M, Ms = zonal_modes(Uf)
    1.0im*spdiagm(Ms)
end

function sinφdφ_zero!(D)
    @assert size(D,1) == size(D,2)
    N = size(D,1)

    for j = 0:N-1, i = 0:N-1
        D[i+1,j+1] = j==i+1 ?  -(i+1)/2:
                     j==i-1 ?   (i-1)/2:
                     0
    end

    D
end

function sinφdφ_odd!(D)
    @assert size(D,1) == size(D,2)
    N = size(D,1)

    for j = 1:N, i = 1:N
        D[i,j] =  j==i+1 ?  -(i+1)/2:
                  j==i-1 ?   (i-1)/2:
                  0
    end

    D
end

function sinφdφ_even!(D)
    @assert size(D,1) == size(D,2)
    N = size(D,1)

    for j = 1:N, i = 1:N
        D[i,j] =  j==i+1 ?  -i/2:
                  j==i-1 ?   i/2:
                  0
    end

    D
end

# sine squared matrices
function sin2_zero!(A)
    for ji in 1:size(A,2), ii in 1:size(A,1)
        i, j = ii-1, ji-1

        A[ii,ji] =
            i==j   ? 1/2 :
            j==i-2 ? -1/4 :
            j==i+2 ? -1/4 :
            0
    end

    A[2,2] = 1/4
    A[3,1] = -1/2
    A[1,3] = -1/4

    A
end

function sin2_odd!(A)
    for j in 1:size(A,2), i in 1:size(A,1)
        A[i,j] =
            i==j   ? 1/2 :
            j==i-2 ? -1/4 :
            j==i+2 ? -1/4 :
            0
    end

    A[1,1] = 3/4;

    A
end

function sin2_even!(A)
    for j in 1:size(A,2), i in 1:size(A,1)
        A[i,j] =
            i==j   ?  1/2 :
            j==i-2 ? -1/4 :
            j==i+2 ? -1/4 :
            0
    end

    A[1,1] = 3/4;

    A
end

# sin2 laplacian matrices, coupled with sin2
function DAzero!(D,A)
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

        A[ii,ji] =
            i==j   ? 1/2 :
            j==i-2 ? -1/4 :
            j==i+2 ? -1/4 :
            0
    end

    D[1,1] = 0 # from comment in Cheong
    D[3,1] = 0 # ditto

    # from comment in cheong, §2.4
    A[2,2] = 1/4
    A[3,1] = -1/2

    # two more equalities given in Cheong
    A[1,3] = -1/4
    D[1,3] = 1/2

    D, A
end

function DAodd!(D, A, M)
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

    D, A
end

function DAeven!(D, A, M)
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

    D, A
end
