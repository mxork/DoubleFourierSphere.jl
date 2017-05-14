# Status

need to recode the fft/ifft routine so that the fourier modes
line up with the computed coefficients, or justify to myself that
the mismatch comes from somewhere else.

then we can address the damping of the Neven, M0 modes in the laplace routine
(I suspect is related).

*ANSWER*: there is a mismatch between muraki's coefficients and Cheong's 
(see §2-8, constants b, c), which are getting in the works when we setup
the linear systems. probably can just hack in the appropriate constants
in the fft routines (or reimplement and use naive as weather rod).
