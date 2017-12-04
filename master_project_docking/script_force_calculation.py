from oct2py import octave
import numpy as np
import math
import matplotlib.pyplot as plt
from oct2py import Struct
# import fcn_non_axial_cylinders as fcyl

# Optimize for force and weight over Number of spires, current, diameter and length


N = float(200)                        # number of spirals in the coil
current = 2.0                         # [A] Current in the coil
diameter = 0.5                        # [m] Diameter of the coil
voltage = 24                          # [V]
length_cylinder = 0.2                 # [m]

radius_cylinder = diameter/2          # [m] enclosed by the loop
area = math.pi*radius_cylinder**2     # [m^2] 

mu_0 = 4*math.pi*10**(-7)             # [H/m] magnetic constant
mu_fixed_magnet = N * current * area  # [H/m] magnetic permeability 1st magnet
mu_float_magnet = N * current * area  # [H/m] magnetic permeability 1st magnet

resistivity_copper = 1.68 * 10**(-8)  # [ohm meters]
density_copper = 8960                 # [kg/m^3]
L_per_turn = 2*math.pi*radius_cylinder
L_total = N * L_per_turn


R = voltage/current                    # [ohm] resistance of the wire
A_wire = resistivity_copper * L_total/R

weight = L_total * A_wire * density_copper #[kg]

print(f'total weight: {weight} kg')
print(f'wire length: {L_total} m')
print(f'wire surf: {A_wire*10**6} mm2')


# Force calculation using standard far-field model
# https://en.wikipedia.org/wiki/Force_between_magnets

x = np.linspace(0.01, 0.5, 1000)
force_far_field = -3*mu_0*mu_fixed_magnet*mu_float_magnet/(2*math.pi*x**4)

# Force calculation using 
# A public framework for calculating the forces between magnets and multipole arrays of magnets
# Author Will Andersen (Univ of Adelaide)

## Magnet definition

coil_fixed = Struct()
coil_float = Struct()

coil_fixed.type = 'cylinder'
coil_float.type = 'cylinder'

coil_fixed.dim = np.array([radius_cylinder, length_cylinder])
coil_float.dim = np.array([radius_cylinder, length_cylinder])

coil_fixed.turns = N
coil_float.turns = N

coil_fixed.magdir = np.array([0, 0, 1.0]) # z
coil_float.magdir = np.array([0, 0, 1.0]) # z

coil_fixed.current = current
coil_float.current = current

displ = np.linspace(0, 1, 400)
direction = np.array([[0], [0], [1.0]])
displ_range = direction*displ




p = 0
X1 = 0
X2 = X1 + length_cylinder
X3 = displ + length_cylinder
X4 = displ + 2*length_cylinder

force = np.zeros(len(displ))

f_xyz_coaxial = octave.magnetforces(coil_fixed, coil_float, displ_range)

# for i in range(0, len(displ)):
#     force[i] = fcyl.cyl_calc_noncoaxial_w(radius_cylinder, radius_cylinder, X1, X2, X3[i], X4[i], N, N, current, current, p)

# plt.plot(x, force_far_field)
# print(force)
plt.plot(displ, f_xyz_coaxial[2], '*')
plt.plot(displ, force, 'g')




plt.ylabel('Force Z[N]')
plt.xlabel('dist [m]')
plt.show()

