import scipy.integrate as integrate
from numpy import arcsin, cos, sin, sqrt
import scipy.special as special
import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# The following equations were implemented according
# Dr. Will Andersen (Univ of Adelaide) matlab code
# and article: Forces Between Thin Coils with Parallel Axes Using Bessel Functions
# IEEE Transactions on Magnetism
# Author: Conway, John T.; University of Agder, Engineering and Science

mu_0 =  4*math.pi*10**(-7) # [ H/m ]

# Exact Solutions for the Magnetic Fields of Axisymmetric Solenoids and Current 
# Distributions





def get_eliptic_param(h, g, m):
    KK = special.ellipk(m)
    EE = special.ellipe(m)
    F2 = special.ellipkinc(arcsin(np.divide(h, g)), 1-m)
    E2 = special.ellipeinc(arcsin(np.divide(h, g)), 1-m)
    return KK, EE, F2, E2


def get_l2_norm(a, b, t):
    return sqrt(a**2 + b**2 - 2*a*b*cos(t))


def cyl_calc_noncoaxial_w(R1, R2, X1, X2, X3, X4, N1, N2, I1, I2, p):
    coefficient = mu_0*N1/(X2-X1)*N2/(X4-X3)*R1*R2*I1*I2
    h = np.array([X4-X2, X3-X2, X4-X1, X3-X1])
    [fx, error] = integrate.quad(force_xdir, 0, math.pi, args=(R1, R2, h, p))
    [fz, error] = integrate.quad(force_zdir, 0, math.pi, args=(R1, R2, h, p))
    Fx = -fx * coefficient # [mN]
    print(f'Fx = {Fx}')
    Fz = -fz * coefficient # [mN]
    print(f'Fz = {Fz}')
    return Fx, Fz
    #return  Fz



def force_zdir(t, r, R, h, p):
    X = get_l2_norm(p, R, t)
    f = np.zeros(4)
    g = np.zeros(4)
    for i in range(0, 4):
        f[i] = sqrt((r+X)**2+h[i]**2)
        g[i] = sqrt((r-X)**2+h[i]**2)
    m = 1-g**2/f**2

    KK, EE, F2, E2 = get_eliptic_param(h, g, m)

    Ta = h*f*(EE-KK)
    Tb = -h*KK*(r-X)**2 / f
    Tc = np.abs(r**2 - X**2)*(F2 * (EE-KK) + KK*E2 - 1)
    Td = 4/math.pi * (r**2 + X**2 - np.abs(r**2 - X**2))/2
    T = (R - p*cos(t)) /(2*r*X**2)* (Ta + Tb + Tc + Td)

    gz = -T[0]+T[1]+T[2]-T[3]
    return gz


def force_xdir(t, r, R, h, p):
    X = get_l2_norm(r, R, t)
    f = np.zeros(4)
    g = np.zeros(4)
    for i in range(0,4):
        f[i] = sqrt((p+X)**2+h[i]**2)
        g[i] = sqrt((p-X)**2+h[i]**2)
    m = 1-g**2/f**2
    KK, EE, F2, E2 = get_eliptic_param(h, g, m)

    Ta = f*EE
    Tb = (p**2 - X**2) * KK / f
    Tc = np.sign(p-X) * h * (F2 * (EE-KK) + KK*E2 - 1)
    Td = -math.pi/2 * h
    if (p == 0):
        T = 0
    else:
        T = cos(t)/p * (Ta + Tb + Tc + Td)
    # print(T)
    gx = -T[0]+T[1]+T[2]-T[3]
    return gx





# R1 = 1
# R2 = 0.5
# length = 0.3
# X1 = 0
# X2 = X1 + length

# distance = 0.15
# X3 = distance + X2
# X4 = X3 + length
# I1 = I2 = 1
# N1 = 100
# N2 = 100

# print (f'{X1} {X2} and {X3} {X4}')
# p = 0


# resistivity_copper = 1.68 * 10**(-8) # [ohm meters]
# density_copper = 8940 # [kg/m^3]
# L_per_turn = 2*np.pi*R1
# L_total = N1 * L_per_turn # [m]
# voltage = 10 # [V]
# print(f'L total {L_total}')

# # R = voltage/I1; #V

# # A_wire = resistivity_copper * L_total/R;
# A_wire = 0.823*10**(-6)

# R = resistivity_copper * L_total/A_wire
# V = R*I1

# weight = L_total * A_wire * density_copper
# print(f'Weight {weight}')
# print(f'V {V}')

# fz = cyl_calc_noncoaxial_w(R1, R2, X1, X2, X3, X4, N1, N2, I1, I2, p)
# print(f'Force Z {fz}')



