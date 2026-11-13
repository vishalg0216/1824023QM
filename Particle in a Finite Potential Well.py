
import numpy as np
import matplotlib.pyplot as plt
import math

V0 = float(input("Enter potential depth: "))
a = float(input("Enter well width: "))
n = int(input("Enter quantum state number: "))

L = 5*a
N = 300
dx = L/(N+1)

x = np.linspace(-L/2, L/2, N)

V = np.zeros(N)
V[np.abs(x) <= a/2] = -V0

H = (2*np.eye(N) - np.diag(np.ones(N-1), 1) - np.diag(np.ones(N-1), -1))/(2*dx**2)
H = H + np.diag(V)

E, psi = np.linalg.eigh(H)

bound = np.where(E < 0)[0]

if n < 1 or n > len(bound):
    print("Invalid state number or no bound state")
else:
    k = bound[n-1]
    wave = psi[:, k]/math.sqrt(dx)

    if np.sum(wave) < 0:
        wave = -wave

    print("Energy Eigenvalue =", E[k])

    plt.plot(x, wave, label="Wave Function")
    plt.xlabel("Position (x)")
    plt.ylabel("Wave Function")
    plt.title("Particle in a Finite Potential Well")
    plt.grid()
    plt.legend()
    plt.show()

    plt.plot(x, V)
    plt.xlabel("Position (x)")
    plt.ylabel("Potential Energy")
    plt.title("Finite Potential Well")
    plt.grid()
    plt.show()
