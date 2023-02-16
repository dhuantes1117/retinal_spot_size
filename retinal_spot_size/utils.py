from retinal_spot_size import *
from .eye_properties import *

def complex_matmul(q_hat, M):
	A = M[0, 0]
	B = Q_(M[0, 1], 'm')
	C = Q_(M[1, 0], 'm^-1')
	D = M[1, 1]
	assert(q_hat.units == ureg.meter)
	NUM = A * q_hat + B
	DEN = C * q_hat + D
	q_prime = NUM / DEN
	assert(q_prime.units == ureg.meter)
	return q_prime

@ureg.wraps(ureg.meter, [ureg.meter, ureg.meter, None])
def q_hat(z, _z_R, n=1):
	return z + 1j*_z_R / n

def propagate_q_hat(q_hat, System_Matrix):
	assert(q_hat.units == ureg.meter)
	return complex_matmul(q_hat, System_Matrix)

@ureg.wraps(ureg.meter, [ureg.meter, ureg.meter])
def current_beam_rad(lambda_, q_hat):
	Imag = np.imag(1 / q_hat)
	waist_sq = -lambda_ / pi / Imag
	return np.sqrt(waist_sq)

@ureg.wraps(ureg.meter, [ureg.meter, ureg.rad])
def z_R(lambda_, div):
	return lambda_ / pi / div**2

@ureg.wraps(ureg.rad, [ureg.meter, ureg.meter])
def ang_div(_lambda, waist):
	return _lambda / pi / waist

@ureg.wraps(None, [None, None, None, ureg.um**2, ureg.um**2, ureg.um**2, ureg.um])
def sellmeier_true(B1, B2, B3, C1, C2, C3, _lambda):
	T0 = 1
	T1 = B1 * _lambda**2 / (_lambda**2 - C1) 
	T2 = B2 * _lambda**2 / (_lambda**2 - C2)
	T3 = B3 * _lambda**2 / (_lambda**2 - C3)
	return np.sqrt(T0 + T1 + T2 + T3)

def sellmeier(B, C, _lambda):
	return sellmeier_true(*B, *C, _lambda)

def sellmeier_vitreous(_lambda):
	return sellmeier_true(*B_COEFFICIENTS_VITREOUS, *C_COEFFICIENTS_VITREOUS, _lambda)
