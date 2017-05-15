# Status

Shadowed the fast sphere transform with some naive ones; they
go plenty fast for the size of data I'm working on. Fixed 
some dimensional issues where I had assumed real input.
FFT/IFFT tests are passing, as is forward laplace. Inverse
laplace is breaking, but only off by a constant term. This
is because I set the constant-mode coefficient to zero during
the inversion process, when I should be setting it to the appropriate
power level of the spherical harmonic (it is non-zero, somewhat
unexpectedly). 

Gonna look up the formula, plop it in, test and should be g2go. Can
either revisit the fast transforms and hack out the details (worthwhile?),
Can also add some '_into' versions of the solvers to speed up the iteration
we'll be doing in the next section:

*Timestepping*.
