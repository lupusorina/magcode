import unittest  
import numpy as np
import fcn_non_axial_cylinders as fnx
import math

class FooTestCase(unittest.TestCase):
    def test_force_table_6(self):
        N1 = N2 = 100
        I1 = I2 = 1

        Z1 = 0
        Z2 = 4

        Z3 = 4
        Z4 = 6

        R1 = 1
        R2 = 0.5

        p = np.array([0.1, 0.4999, 1.0, 1.5001, 2.0, 5.0])
        true_FX = np.array([-0.03159192174752516, -0.2001494249235013, -0.3554204990676609, -0.1689290127079929, -0.06024079635233924, -0.00003576808160726246])
        true_FZ = np.array([-0.5440315334495091, -0.5478220750465542, -0.2166979779748779, 0.04601077650686824, 0.03407448362790638, 0.00497123544052176])
        h = np.array([Z4-Z2, Z3-Z2, Z4-Z1, Z3-Z1])
        fx = np.zeros(4)
        fz = np.zeros(4)
        for i in range(0, 4):
            fx[i], fz[i] = fnx.cyl_calc_noncoaxial_w(R1, R2, Z1, Z2, Z3, Z4, N1, N2, I1, I2, p[i])
            self.assertAlmostEqual(fx[i], true_FX[i])
            self.assertAlmostEqual(fz[i], true_FZ[i])

    def test_force_table_7(self):
        N1 = N2 = 100
        I1 = I2 = 1
        c = 1
        Z1 = 0
        Z2 = 4

        Z3 = 1 + c
        Z4 = 3 + c

        R1 = 1
        R2 = 0.5

        p = np.array([0.25, 0.5, 1.6, 1.8, 2.0])
        true_FX = np.array([0.08803276352092297, 0.2119335033836623, 0.1531385579896640, 0.1091886537529803, 0.08267048939348621])
        true_FZ = np.array([-0.5106532345800797, -0.5158240831644936, 0.05732233766284198, 0.04935774531856688,0.04212255772204719])

        h = np.array([Z4-Z2, Z3-Z2, Z4-Z1, Z3-Z1])
        fx = np.zeros(4)
        fz = np.zeros(4)
        for i in range(0, 4):
            fx[i], fz[i] = fnx.cyl_calc_noncoaxial_w(R1, R2, Z1, Z2, Z3, Z4, N1, N2, I1, I2, p[i])
            print(fx[i])

            self.assertAlmostEqual(fx[i], true_FX[i])
            self.assertAlmostEqual(fz[i], true_FZ[i])

    def test_get_l2_norm(self):
        self.assertAlmostEqual(fnx.get_l2_norm(1,1, np.pi/2), np.sqrt(2))

    # def test_lambda_function(self):
    #     # calculated in matlab
    #     E = 1.53075
    #     K = 1.61244
    #     F2 = 0.90145
    #     E2 = 0.71781
    #     beta = 0.8
    #     k = 0.1
    #     lambda_test = 2/math.pi*(E*F2 + K*E2 - K*F2)
    #     self.assertAlmostEqual(lambda_test, fnx.lambda_function(beta, k), 4)

    # def test_beta_k(self):
    #     R = 0.25
    #     r = 0.3
    #     z = 0.2
    #     Z2 = 0.1
    #     #self.assertAlmostEqual(1.325817, fnx.beta_k(r, R, z)[0], 3)
    #     self.assertAlmostEqual(0.9349019, fnx.beta_k(r, R, np.abs(Z2-z))[1], 3)







