module DiffMatrices

# following Cheong

# tranposes? check matrix dimensions

function DAodd(N, M)
    D = [
          i==j   ? -(i^2 + 2*M^2)/2 :
          j==i-2 ? (i-1)*(i-2)/4 :
          j==i+2 ? (i+1)*(i+2)/4 :
          0
          for i in 1:N, j in 1:N
        ]

    A = [
          i==j   ? 1/2 :
          j==i-2 ? -1/4 :
          j==i+2 ? -1/4 :
          0
          for i in 1:N, j in 1:N
        ]

    A[1,1] = 3/4;
    # TODO figure out what the hell Cheong means by 'for odd n'
    # pretty sure he decomposes it into two true tridiagonal matrices,
    # but I don't think it's necessary.
    # consider using a sparse matrix set up. just call sparse(A), sparse(D)

    return D, A
end

function DAeven(N, M)
    D = [
          i==j   ? -(i^2 + 2*M^2)/2 :
          j==i-2 ? i*(i-1)/4 :
          j==i+2 ? i*(i+1)/4 :
          0
          for i in 1:N, j in 1:N
        ]

    A = [
          i==j   ? 1/2 :
          j==i-2 ? -1/4 :
          j==i+2 ? -1/4 :
          0
          for i in 1:N, j in 1:N
        ]

    A[1,1] = 3/4;

    return D, A
end

function DAzero(N)
    M = 0 # by definition

    # TODO, check ranges for i,j here, cuz they are different
    D = [
          i==j   ? -(i^2 + 2*M^2)/2 :
          j==i-2 ? i*(i-1)/4 :
          j==i+2 ? i*(i+1)/4 :
          0
          for i in 0:N-1, j in 0:N-1
        ]

    D[1,1] = 0 # from comment in Cheong
    D[3,1] = 0 # ditto

    A = [
          i==j   ? 1/2 :
          j==i-2 ? -1/4 :
          j==i+2 ? -1/4 :
          0
          for i in 1:N, j in 1:N
        ]

    A[1,1] = 1/4
    A[2,1] = -1/2 # actually A[2,1], since starting from zero on evens.

    # two more equalities given in Cheong that bewilder me

    return D, A
end

end
