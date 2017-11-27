import scipy.integrate as integrate
from numpy import arcsin, cos, sin, sqrt
import scipy.special as special
import math
import numpy as np


def cyl_calc_noncoaxial_w(R1, R2, Z1, Z2, Z3, Z4, N1, N2, I1, I2, p):
    coefficient = 4*10**(-4)*N1/(Z2-Z1)*N2/(Z4-Z3)*R1*R2*I1*I2
    [fx, error] = integrate.quad(xdir, 0, math.pi, args=(R1, R2, h, p))
    
    # [fz, error] = integrate.quad(zdir, 0, math.pi, args=(R1, R2, h, p))
    Fx = -fx * coefficient # [mN]


def xdir(t, r, R, h, p):
    X = sqrt(r**2+R**2-2*r*R*cos(t))
    f = np.zeros(4)
    g = np.zeros(4)
    for i in range(0,4):
        f[i] = sqrt((p+X)**2+h[i]**2)
        g[i] = sqrt((p-X)**2+h[i]**2)
    m = 1-g**2/f**2

    KK = special.ellipk(m)
    EE = special.ellipe(m)
    F2 = special.ellipkinc(arcsin(np.divide(h, g)), 1-m)
    E2 = special.ellipeinc(arcsin(np.divide(h, g)), 1-m)

    Ta = f*EE
    Tb = (p**2 - X**2) * KK / f
    Tc = np.sign(p-X) * h * (F2 * (EE-KK) + KK*E2 - 1)
    Td = -math.pi/2 * h
    T = cos(t)/p * (Ta + Tb + Tc + Td)

    gx = -T[0]+T[1]+T[2]-T[3]
    return gx


N1 = N2 = 100
I1 = I2 = 1

Z1 = 0
Z2 = 4

Z3 = 4
Z4 = 6

p = 1
R1 = 1
R2 = 0.5

const = 4*10**(-4)*N1/(Z2-Z1)*N2/(Z4-Z3)*R1*R2*I1*I2

h = np.array([Z4-Z2, Z3-Z2, Z4-Z1, Z3-Z1])
cyl_calc_noncoaxial_w(R1, R2, Z1, Z2, Z3, Z4, N1, N2, I1, I2, p)