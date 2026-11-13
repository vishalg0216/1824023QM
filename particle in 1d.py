
import numpy as np
import matplotlib.pyplot as plt
import math

L = float(input("Enter box length: "))
n = int(input("Enter quantum number: "))

N = 100
dx = L / (N + 1)

x = np.linspace(dx, L-dx, N)

H = (2*np.eye(N) - np.diag(np.ones(N-1), 1) - np.diag(np.ones(N-1), -1)) / (2*dx**2)

E, V = np.linalg.eigh(H)

psi = V[:, n-1] / np.sqrt(dx)

psi_exact = np.sqrt(2/L) * np.sin(n*math.pi*x/L)

if np.sum(psi*psi_exact) < 0:
    psi = -psi

E_exact = n**2 * math.pi**2 / (2*L**2)

print("Numerical Energy =", E[n-1])
print("Analytical Energy =", E_exact)

plt.plot(x, psi, label="Numerical")
plt.plot(x, psi_exact, "--", label="Analytical")
plt.xlabel("Position (x)")
plt.ylabel("Wave Function")
plt.title("Particle in a 1D Infinite Potential Box")
plt.legend()
plt.grid()
plt.show()
