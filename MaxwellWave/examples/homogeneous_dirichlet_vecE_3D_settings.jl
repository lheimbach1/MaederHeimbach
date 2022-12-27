# physics
lx = 10
ly = lx
lz = lx
# nonphysical constants for proof of concept
mu0 = 1
epsilon0 = 1
c2 = 1 / mu0 / epsilon0
lambda = 1
k0 = 2 * pi / lambda
w0 = k0 * sqrt(c2)
p0 = zeros(3)

# linearly polarized in y direction
p0[2] = 1
#pulse shape 
sigma2 = [1 1 1]

# numerics
nx = 40
ny = nx
nz = nx
nt = 1000
nvis = 10
