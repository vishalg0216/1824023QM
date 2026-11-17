import numpy as np
import matplotlib.pyplot as plt
import math

m = 940
k = 100
h = 197.3

r = np.linspace(0.01, 5, 500)
dr = r[1] - r[0]

b_values = [0, 10, 30]

for b in b_values:

    V = 0.5*k*r**2 + (b/3)*r**3

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

    wave = wave / math.sqrt(np.sum(wave**2)*dr)

    print("b =", b)
    print("Ground state energy =", round(E[0], 3), "MeV")

    plt.plot(r[1:-1], wave, label="b = "+str(b))

plt.xlabel("r (fm)")
plt.ylabel("Wave Function")
plt.title("Anharmonic Oscillator Ground State")
plt.legend()
plt.grid()
plt.show()