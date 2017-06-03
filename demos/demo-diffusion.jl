import DoubleFourierSphere; dfs = DoubleFourierSphere
using Plots
using JLD

Plots.gr()

G = zeros(Complex128, 128, 64)
G[13:33, 43:63] = [ e^(-0.03*x^2 - 0.03*y^2) for x in -10:10, y in -10:10]
G0 = copy(G)
Gf = dfs.fft_sphere(G)
Gf *= 3

# dt arbitrarily chosen
Diffuse! = dfs.plan_diffuse_spectral!(Gf, 1e-4)

anim = Animation()

for i= 1:1000
    Diffuse!(Gf)
    if i%10 == 0
		Plots.heatmap(real( dfs.ift_sphere(Gf))', color_limits=(0,1), colorbar=:none)
        Plots.frame(anim)
    end
end

gif(anim, "diffusion.gif", fps=10)
