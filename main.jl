# G(λ,θ) is real scalar field, equally spaced
# over lat-long interior grid.
# returns F, which is similar.
function poisson(G::Array{Float64})
    Gc = fftsphere(G)
    L,Y,Z = size(Gc,1), size(G,2), size(Gc, 2)

    #solution matrix
    Fc = Array{Complex128}(L,Z)

    #poisson
    # TODO check this with fine-tooth comb
    A = [i==j ? 1/2 : abs(i-j)==2 ? -1/4 :0 for i in 1:L, j in 1:Z]
    A[1,1] = 3/4 # periodicity of g 
    for mi in 2:L
        m = mi-1 #TODO replace this with shifted index lookup
        if isodd(m)
            D = [ i  == j ? -(2*m^2+j^2)/2 :
                  i-j==-2 ? (j-1)*(j-2)/4 :
                  i-j== 2 ? (j+1)*(j+2)/4 :
                 0
                 for i in 1:L, j in 1:Z]
        else
            D = [  i  == j ? -(2*m^2+j^2)/2 :
                  i-j==-2 ? j*(j-1)/4 :
                  i-j== 2 ? j*(j+1)/4 :
                 0
                 for i in 1:L, j in 1:Z]
        end

        Fc[mi] = D\(A*G[mi])
    end

    F = ifftsphere(Fc)
end

function dcost(X)
    N = size(X,1)
    θs = 
    [sum(X .* cos(2*pi*l*((1:N)-0.5)/N) ) for l in (-div(N,2)+1):div(N,2)]/N
end

function dsint(X)
    N = size(X,1)
    [sum(X .* sin(2*pi*l*((1:N)-0.5)/N) ) for l in (-div(N,2)+1):div(N,2)]/N
end

function dsinsint(X)
    N = size(X,1)
    [sum(X .* sin(2*pi*((1:N)-0.5)/N) .* sin(2*pi*l*((1:N)-0.5)/N) ) for l in (-div(N,2)+1):div(N,2)]/N
end

function idcost(X)
    N = size(X,1)
    [sum(X .* cos(θ * ((-div(N,2)+1):div(N,2) )) ) for θ in 2*pi*((1:N)-0.5)/N]
end

function idsint(X)
    N = size(X,1)
    [sum(X .* sin(θ * ((-div(N,2)+1):div(N,2) )) ) for θ in 2*pi*((1:N)-0.5)/N]
end

function idsinsint(X)
    N = size(X,1)
    [sin(θ) * sum(X .* sin(θ * ((-div(N,2)+1):div(N,2) )) ) for θ in 2*pi*((1:N)-0.5)/N]
end

function fftsphere(G::Array{Float64})
    Gc = complete2torus(G)
    # longitudinal transform: 1 means FFT wrt to longitude dim
    Gc = fft(Gc, 1)

    L,Y,Z = size(Gc,1), size(G,2), size(Gc, 2)
    dλ, dθ = 2*π/L, 2*π/Z
    off = -dθ/2 # θ offset for interior grid

    Gf = Array{Complex128}(size(Gc))
    #  manually calc coefficients
    ms = [(0:div(L,2))... ((-div(L,2)+1):-1)...]
    for mi in 1:L
          m = ms[mi]
          if m == 0
              Gf[mi,:] = dcost(Gc[mi, :])
          elseif iseven(m) 
              Gf[mi,:] = dsinsint(Gc[mi, :])
          else
              Gf[mi,:] = dsint(Gc[mi, :])
          end
    end

    # ifftshift(Gf, 2)
    Gf
end

function ifftsphere(Gc::Array{Complex128})
    L,Y,Z = size(Gc,1), div(size(Gc,2),2), size(Gc, 2)
    dλ, dθ = 2*π/L, 2*π/Z
    off = -dθ/2 # θ offset for interior grid

    Gf = Array{Complex128}(size(Gc))
    #  manually calc coefficients
    ms = [(0:div(L,2))... ((-div(L,2)+1):-1)...]
    for mi in 1:L
          m = ms[mi]
          if m == 0
              Gf[mi,:] = idcost(Gc[mi, :])
          elseif iseven(m) 
              Gf[mi,:] = idsinsint(Gc[mi, :])
          else
              Gf[mi,:] = idsint(Gc[mi, :])
          end
    end

    Gf = ifft(Gf,1)
    Gf = real(Gf)
    
    #TODO return lower half
    Gf[:, 1:Y]
end

# returns the periodic completion of F
# to 2-torus, make it complex while we're at it
function complete2torus(F::Array{Float64})
    L, Y = size(F,1), size(F,2)
    Fc = Array{Complex128}(L, 2*Y)

    # identity on first half
    for i in 1:L
        for j in 1:Y
            Fc[i,j] = F[i,j]
        end
    end

    # shifted on second
    for i in 1:L
        for j in Y+1:2*Y
            Fc[i,j] = Fc[modu(-1+i+div(L,2),L), modu(1-j,Y)]
        end
    end

    Fc
end

function modu(x, m)
    x <= 0 ? m+(x%m) : (x%m)+1
end

# returns the m-l fourier mode on the sphere,
# m = λ wavenumber
# l = θ wavenumber
function fouriermode(m, l)
    if m==0 
        return (λ, θ) -> cos(l*θ)
    elseif m%2 == 0
        return (λ, θ) -> sin(θ)*sin(l*θ)*exp(1.0im*m*λ)
    else
        return (λ, θ) -> sin(l*θ)*exp(1.0im*m*λ)
    end
end

function sphericalmode(m,l)
    return (λ, θ) -> legendre(m,l, cos(θ))*exp(1.0im*m*λ)
end

# rotate around x-axis by α TODO, flipped orientation
function rotX(λ, θ, α)
    (
        atan2( sin(θ)*cos(λ),(sin(θ)*sin(λ)*cos(α) + cos(θ)*sin(α)) ), #phase shift here, TODO
        acos(-sin(θ)*sin(λ)*sin(α)+cos(θ)cos(α)),
    )
end

function fouriermodeα(m,l,α)
    Q = fouriermode(m,l)
    # atan2(x,y) == atan(y/x)
    (λ, θ) -> Q(rotX(λ, θ, α)...)
end

function sphericalmodeα(m,l,α)
    return function (λ, θ)
        λ, θ = rotX(λ, θ, α)
        legendre(m,l,cos(θ))*exp(1.0im*m*λ)
    end
end

using GSL
function legendre(m,l,x)
    sf_legendre_sphPlm(l,m,x)
end

# returns a real scalar field F on long-lat grid,
# with F(i,j) = f(2*pi*i/L, pi*j/M)
function fillsphere(L, M, f)
    #0.5== interior grid
    [f(2*pi*i/L, pi*(j+0.5)/M) for i in 0:(L-1), j in 0:(M-1)]
end

# utilities for plotting. crikey.
function λsphere(L, M)
    [2*pi*i/L for i in 0:(L-1), j in 0:(M-1)]
end

function θsphere(L, M)
    [pi*(j+0.5)/M for i in 0:(L-1), j in 0:(M-1)]
end

# a là Muraki
using PyPlot
function sphereplot(G)
    clf()
    L, M = size(G,1), size(G,2)
    contourf(λsphere(L,M), θsphere(L,M), G, cmap=ColorMap("inferno"))
end

function demo01(m, l, output_filename)
    clf()
    L, M = 40, 20
    Q = fouriermode(m, l)

    contourf(λsphere(L,M), θsphere(L,M), real(fillsphere(L,M,Q)), cmap=ColorMap("inferno"))
end

function demo02(m, l, α, output_filename)
    clf()
    L, M = 40, 20
    Q = fouriermodeα(m, l, α)

    contourf(λsphere(L,M), θsphere(L,M), real(fillsphere(L,M,Q)), cmap=ColorMap("inferno"))
end

function ckfftinv(m,l)
    L, M = 20,10
    # contourf(λsphere(L,M), θsphere(L,M), G, cmap=ColorMap("inferno"))
    # contourf(λsphere(L,M), θsphere(L,M), ifftsphere(fftsphere(G)), cmap=ColorMap("inferno"))
    # for m in 1:4
    #     for l in 1:4
            Q = fouriermode(m, l)
            G = real(fillsphere(L,M,Q))
            println(abs(G-ifftsphere(fftsphere(G))))
        # end
    # end
    G
end

function run()
    L, M = 16,8
    @show G = real(fillsphere(L,M,fouriermode(2,1)))
    println()
    @show round(real(complete2torus(G))*100)
end

run()
