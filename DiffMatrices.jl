export DAodd, DAeven, DAzero

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

# These are giving me trouble. I need to sit down and do algebra.
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

    #D[1,1] = 0 # from comment in Cheong, but we need non-zero element otherwise singular
    D[3,1] = 0 # ditto

    A = [
          i==j   ? 1/2 :
          j==i-2 ? -1/4 :
          j==i+2 ? -1/4 :
          0
          for i in 1:N, j in 1:N
        ]

    # my own guesses CRAP
    # based off of sums of columns and rows
    D[1,1] = 0
    D[3,1] = 1/2

    # remember, n=0,1,2,3,...
    # A[2,2] = 1/4
    # A[3,1] = -1/2 

    # two more equalities given in Cheong that bewilder me
    # A[1,2] = -1/4
    # D[1,2] = 1/2

    # TODO as of right now, D is singular
    D[1,1] = 1.0 # cheat, reset after solve

    return D, A
end
