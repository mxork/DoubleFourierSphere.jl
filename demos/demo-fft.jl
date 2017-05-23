using DoubleFourierSphere

M = 4
N = 3

Gλ, Gφ = spheregrids(64,32)
U = fouriermode(M,N)(Gλ, Gφ)
Uf = fft_sphere(U)

# display

import Plots
Plots.gr(size=(1200,400))

Plots.plot( plot_sphere(real(U)), plot_frequency(abs(Uf)) )

Plots.savefig("output-images/output.png")
