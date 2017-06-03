import DoubleFourierSphere; dfs = DoubleFourierSphere
using Plots
using JLD

Plots.gr()

@load "../../demos/G.jld"
@load "../../demos/tilted-stream.jld" # called "KK"

G0 = copy(G)
G0 *= 3
G0f = dfs.fft_sphere(G0)

# import stream function, take gradient
Vy, Vx = dfs.gradient(KK)
Vx *= -1

# dt arbitrarily chosen
dt = 5e-3
Advect! = dfs.plan_advection_spectral!(G0f)
Diffuse! = dfs.plan_diffuse_spectral!(G0f, dt*2*pi/128/16)

anim = Animation()

for i= 1:1000
    Advect!(G0f, Vx, Vy, dt); Diffuse!(G0f)
    if i%10 == 0
        println("tick")
        heatmap( real(dfs.ift_sphere(G0f))', color_limits=(0,1), colorbar=:none)
        Plots.frame(anim)
    end
end

gif(anim, "advect.gif", fps=10)
