module DoubleFourierSphere
	include("fft_sphere.jl")
	include("laplace.jl")
	include("gradient.jl")
	include("diffusion.jl")
	include("advection.jl")
	include("wave.jl")
	include("matrices.jl")
	include("vorticity.jl")
	include("helper.jl")
end
