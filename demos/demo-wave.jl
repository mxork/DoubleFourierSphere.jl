import DoubleFourierSphere; dfs = DoubleFourierSphere
using Plots
using JLD

Plots.gr()


G = zeros(Complex128, 128, 64)
G[13:33, 23:43] = [ e^(-0.03*x^2 - 0.03*y^2) for x in -10:10, y in -10:10]

G *= 3

G0 = copy(G)
G1 = copy(G)
G0f = dfs.fft_sphere(G0)
G1f = dfs.fft_sphere(G1)

# dt arbitrarily chosen
Wave! = dfs.plan_wave_spectral!(G0f, 6e-3)

anim = Animation()

for i= 1:1000
    G1f, G0f = Wave!(G1f, G0f)
    if i%10 == 0
        heatmap( real(dfs.ift_sphere(G1f))', color_limits=(-1,1), colorbar=:none)
        Plots.frame(anim)
    end
end

gif(anim, "wave.gif", fps=10)
