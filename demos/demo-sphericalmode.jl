using Helper

print("Longitudinal wavenumber m: ")
M = parse(Int, readline())

print("Latitudinal wavenumber n: ")
N = parse(Int, readline())

sphereplot(sphericalmode(M,N)( spheregrids(128, 64)... ))
