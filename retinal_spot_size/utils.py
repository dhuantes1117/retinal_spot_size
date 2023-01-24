import numpy as np
from scipy.constants import pi

def complex_matmul(q_hat, M):
        A = M[0, 0]
        B = M[0, 1]
        C = M[1, 0]
        D = M[1, 1]
        NUM = A * q_hat + B
        DEN = C * q_hat + D
        return NUM / DEN

def q(z, _z_R):
	return z + 1j*_z_R

def propagate_q(q_0, System_Matrix):
	return complex_matmul(q_0, System_Matrix)

def w_z(lambda_, q):
	Imag = np.imag(1 / q)
	waist_sq = -lambda_ / Imag / pi
	return np.sqrt(waist_sq)

def z_R(lambda_, div):
	return lambda_ / pi / div**2

def sellmeier(n_0, A, B, _lambda):
	return n_0
