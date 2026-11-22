
import numpy as np
import matplotlib.pyplot as plt
import math

hbar = 1.054e-34
m = 9.1e-31
e = 1.602e-19


a = 0.3e-9
V0 = 10*e


z0 = a*math.sqrt(2*m*V0)/hbar


def root(f, x1, x2):
    for i in range(80):
        xm = (x1+x2)/2

        if f(x1)*f(xm) < 0:
            x2 = xm
        else:
            x1 = xm

    return (x1+x2)/2

states = []

for n in range(int(z0/math.pi)+2):

    
    x1 = n*math.pi + 0.00001
    x2 = min(n*math.pi+math.pi/2-0.00001, z0-0.00001)

    if x2 > x1:
        f = lambda z: z*math.tan(z)-math.sqrt(z0**2-z**2)

        if f(x1)*f(x2) < 0:
            states.append((root(f, x1, x2), "even"))


    x1 = (n+0.5)*math.pi+0.00001
    x2 = min((n+1)*math.pi-0.00001, z0-0.00001)

    if x2 > x1:
        f = lambda z: -z/math.tan(z)-math.sqrt(z0**2-z**2)

        if f(x1)*f(x2) < 0:
            states.append((root(f, x1, x2), "odd"))


states.sort()


x = np.linspace(-4*a, 4*a, 1000)


for i, (z, parity) in enumerate(states):

    k = z/a
    E = hbar**2*k**2/(2*m)-V0
    alpha = math.sqrt(-2*m*E)/hbar

    if parity == "even":
        psi = np.where(abs(x) <= a,
                       np.cos(k*x),
                       np.cos(z)*np.exp(-alpha*(abs(x)-a)))
    else:
        psi = np.where(abs(x) <= a,
                       np.sin(k*x),
                       np.sign(x)*np.sin(z)*np.exp(-alpha*(abs(x)-a)))

    
    dx = x[1]-x[0]
    psi = psi/math.sqrt(np.sum(psi**2)*dx)

    print("Energy", i+1, "=", round(E/e, 4), "eV")

    plt.plot(x*1e9, psi, label="State "+str(i+1))

plt.xlabel("Position (nm)")
plt.ylabel("Normalized Wavefunction")
plt.title("Particle in Finite Potential Well")
plt.legend()
plt.grid()
plt.show()
