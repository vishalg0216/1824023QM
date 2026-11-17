import numpy as np
import matplotlib.pyplot as plt
import math

m = 940 * 10**6
D = 0.755501
alpha = 1.44
r0 = 0.131349
h = 1973

r = np.linspace(0.05, 0.8, 800)
dr = r[1] - r[0]

rp = (r - r0)/r0

V = D*(np.exp(-2*alpha*rp) - np.exp(-alpha*rp))

n = len(r) - 2

A = np.zeros((n, n))

for i in range(n):
    A[i, i] = 2
    if i > 0:
        A[i, i-1] = -1
    if i < n-1:
        A[i, i+1] = -1

H = (h**2/(2*m*dr**2))*A

for i in range(n):
    H[i, i] += V[i+1]

E, psi = np.linalg.eigh(H)

wave = psi[:, 0]

wave = wave/math.sqrt(np.sum(wave**2)*dr)

print("Ground state energy =", round(E[0], 6), "eV")

plt.plot(r[1:-1], wave)
plt.xlabel("Internuclear distance (Angstrom)")
plt.ylabel("Wave Function")
plt.title("Morse Potential Ground State")
plt.grid()
plt.show()