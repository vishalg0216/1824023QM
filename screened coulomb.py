import numpy as np
import matplotlib.pyplot as plt
import math

N = int(input("Enter number of grid points: "))
rmax = float(input("Enter maximum radius in Angstrom: "))
a = float(input("Enter screening length a (3, 5 or 7): "))

e = 3.795
hbar = 1973
m = 0.511e6

alpha = hbar**2 / (2*m)

r = np.linspace(0, rmax, N+2)
dr = r[1] - r[0]
x = r[1:-1]

V = -e**2 * np.exp(-x/a) / x

H = np.zeros((N, N))

for i in range(N):
    H[i, i] = 2*alpha/dr**2 + V[i]

    if i > 0:
        H[i, i-1] = -alpha/dr**2

    if i < N-1:
        H[i, i+1] = -alpha/dr**2

energy, wave = np.linalg.eigh(H)

print("\nScreened Coulomb Potential")
print("Screening length =", a, "Angstrom")
print("Ground State Energy =", round(energy[0], 3), "eV")

u = np.zeros(N+2)
u[1:-1] = wave[:, 0]

u = u / math.sqrt(np.trapezoid(u**2, r))

plt.figure(figsize=(9,5))
plt.plot(r, u)
plt.xlabel("Radius (Angstrom)")
plt.ylabel("Radial Wavefunction u(r)")
plt.title("Screened Coulomb Potential")
plt.grid()
plt.show()
