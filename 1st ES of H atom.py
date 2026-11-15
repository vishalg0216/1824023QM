import numpy as np
import matplotlib.pyplot as plt
import math

N = int(input("Enter number of grid points: "))
rmax = float(input("Enter maximum radius in Angstrom: "))

e = 3.795
hbar = 1973
m = 0.511e6

alpha = hbar**2 / (2*m)

r = np.linspace(0, rmax, N+2)
dr = r[1] - r[0]
x = r[1:-1]

V = -e**2 / x

H = np.zeros((N, N))

for i in range(N):
    H[i, i] = 2*alpha/dr**2 + V[i]

    if i > 0:
        H[i, i-1] = -alpha/dr**2

    if i < N-1:
        H[i, i+1] = -alpha/dr**2

energy, wave = np.linalg.eigh(H)

u = np.zeros(N+2)
u[1:-1] = wave[:, 1]

u = u / math.sqrt(np.trapezoid(u**2, r))

print("First Excited State Energy =", round(energy[1], 3), "eV")
print("Analytical Energy =", -3.4, "eV")

plt.plot(r, u, label="First Excited State")
plt.xlabel("Radius (Angstrom)")
plt.ylabel("Radial Wavefunction")
plt.title("Hydrogen Atom First Excited State")
plt.legend()
plt.grid()
plt.show()
