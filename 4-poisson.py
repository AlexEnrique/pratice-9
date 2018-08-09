from numpy import empty,zeros,max
from pylab import imshow,gray,show

# Constants
M = 100          # Grid square on a side
target = 1e-5    # Target accuracy
rho0 = 1         # Charge density (C/m**2)
len = 14         # Length of square charge (cm)
pX, pY = 70, 50  # Positive charge center coord
nX, nY = 30, 50  # Negative charge center coord

# Create arrays to hold potential values
phi = zeros([M+1,M+1], float)
phiprime = empty([M+1,M+1], float)

# Create array to hold values of charge density
rho = zeros([M+1,M+1], float)
rho[int(pX - len / 2):int(pX + len / 2), int(pY - len / 2):int(pY + len / 2)] = rho0
rho[int(nX - len / 2):int(nX + len / 2), int(nY - len / 2):int(nY + len / 2)] = -rho0

# Main loop
delta = 1.0
while (delta > target):
    # Calculate new values of the potential
    for i in range(M+1):
        for j in range(M+1):
                if i == 0 or i == M or j == 0 or j == M:
                    phiprime[i,j] = phi[i,j]
                else:
                    phiprime[i,j] = (phi[i+1,j] + phi[i-1, j] \
                                    + phi[i, j+1] + phi[i,j-1])/4 \
                                    + rho[i,j] / (4)

    delta = max(abs(phi-phiprime))

    phi, phiprime = phiprime, phi

imshow(phi)
gray()
show()
