
import numpy as np
import matplotlib.pyplot as plt
import math


l = 2
m = 1

m = abs(m)


theta = np.linspace(0, math.pi, 100)
phi = np.linspace(0, 2*math.pi, 100)

TH, PH = np.meshgrid(theta, phi)

x = np.cos(TH)


P = np.ones_like(x)

for i in range(m):
    P = -(2*i+1)*np.sqrt(1-x**2)*P

if l > m:
    P1 = x*(2*m+1)*P

    for i in range(m+2, l+1):
        P2 = ((2*i-1)*x*P1-(i+m-1)*P)/(i-m)
        P = P1
        P1 = P2

    P = P1


N = math.sqrt((2*l+1)/(4*math.pi) *
              math.factorial(l-m)/math.factorial(l+m))


Y = N * P
probability = Y**2


X = np.sin(TH)*np.cos(PH)
Yc = np.sin(TH)*np.sin(PH)
Z = np.cos(TH)


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

surface = ax.plot_surface(
    X, Yc, Z,
    facecolors=plt.cm.jet(probability/probability.max())
)

plt.title("Spherical Harmonics Probability Density")
plt.show()
