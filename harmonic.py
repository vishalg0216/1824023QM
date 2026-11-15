
import numpy as np
import matplotlib.pyplot as plt
import math

k = float(input("Enter force constant k: "))
N = int(input("Enter number of grid points: "))
xmax = int(input("Enter maximum position: "))
states = int(input("Enter number of states: "))

x = np.linspace(-xmax, xmax, N + 2)
dx = x[1] - x[0]

V = 0.5 * k * x**2

H = np.zeros((N, N))

for i in range(N):
    H[i, i] = 2 / dx**2 + V[i + 1]

    if i > 0:
        H[i, i-1] = -1 / dx**2

    if i < N-1:
        H[i, i+1] = -1 / dx**2

energy, wavefunctions = np.linalg.eigh(H)

print("\nEnergy Eigenvalues")

for i in range(states):
    E = energy[i]
    E_analytical = (i + 0.5) * math.sqrt(2 * k)

    print("State", i + 1)
    print("Numerical Energy =", round(E, 5))
    print("Analytical Energy =", round(E_analytical, 5))

plt.figure(figsize=(9, 5))
plt.plot(x, V)
plt.xlabel("Position")
plt.ylabel("Potential Energy")
plt.title("Harmonic Potential")
plt.grid()
plt.show()

plt.figure(figsize=(9, 5))

for i in range(states):
    psi = np.zeros(N + 2)
    psi[1:-1] = wavefunctions[:, i]
    norm = math.sqrt(np.trapezoid(psi**2, x))
    psi = psi / norm

    plt.plot(x, psi, label=f"State {i+1}")

plt.xlabel("Position")
plt.ylabel("Wavefunction")
plt.title("Wavefunctions")
plt.legend()
plt.grid()
plt.show()

plt.figure(figsize=(9, 5))

for i in range(states):
    psi = np.zeros(N + 2)
    psi[1:-1] = wavefunctions[:, i]
    norm = math.sqrt(np.trapezoid(psi**2, x))
    psi = psi / norm

    plt.plot(x, psi**2, label=f"State {i+1}")

plt.xlabel("Position")
plt.ylabel("Probability Density")
plt.title("Probability Density")
plt.legend()
plt.grid()
plt.show()