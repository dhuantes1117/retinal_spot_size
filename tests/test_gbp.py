from retinal_spot_size import *
from retinal_spot_size import ABCD
from retinal_spot_size import eye_ABCD
from retinal_spot_size.utils import *
from pytest import *
import numpy.random as rand

def test_construct_thin_lens_matrix():
	mat = ABCD.thin_lens(2 * ureg.meter)
	assert mat[0,0] == 1
	assert mat[0,1] == 0
	assert mat[1,0] == -0.5
	assert mat[1,1] == 1
	
def test_construct_free_space_matrix():
	mat = ABCD.free_space(6 * ureg.meter, 7)
	assert mat[0,0] == 1
	assert mat[0,1] == approx(6/7)
	assert mat[1,0] == 0
	assert mat[1,1] == 1

def test_construct_curved_interface_matrix():
	mat = ABCD.curved_interface(1 * ureg.meter, 2, 3)
	assert mat[0,0] == 1
	assert mat[0,1] == 0
	assert mat[1,0] == 1
	assert mat[1,1] == 1

def test_matmul():
	z_0 	= Q_(rand.rand(), 'm')
	z_R_0 	= Q_(rand.rand(), 'm')
	q_0 = q_hat(z_0, z_R_0)
	q_1 = complex_matmul(q_0, np.identity(2))
	assert(q_0.units == ureg.meter)
	assert(q_0 == q_1)


def test_free_space():
	L = 10 * rand.rand() * ureg.meter
	n = rand.rand()
	z_0 	= Q_(rand.rand(), 'm')
	z_R_0 	= Q_(rand.rand(), 'm')
	q_0 = q_hat(z_0, z_R_0)
	q_1 = propagate_q_hat(q_0, ABCD.free_space(L, n))
	q_2 = propagate_q_hat(q_1, ABCD.free_space(-L, n))
	assert(q_0.units == q_2.units)
	assert(q_0.magnitude == approx(q_2.magnitude))

# update to involve every matrix
def test_optical_train():
	f = 10 * rand.rand() * ureg.mm
	n = 1 #+ rand.rand()
	n2 = 2
	R = 10 * ureg.cm
	free_space_focal_length = ABCD.free_space(f, n)
	thin_lens_mat		= ABCD.thin_lens(f)
	curved_int		= ABCD.curved_interface(R, n, n2)
	OpticalSystem = thin_lens_mat @ free_space_focal_length @ curved_int
	assert(approx(np.linalg.det(OpticalSystem)) == 1)

def test_eye_ABCD():
	assert(approx(np.linalg.det(eye_ABCD.eye_equiv_ABCD(532 * ureg.nm))) == 1)

def test_thin_lens():
	f = 10 * rand.rand() * ureg.mm
	n = 1 #+ rand.rand()
	z_0 	= 0 * ureg.meter
	z_R_0 	= Q_(rand.rand(), 'mm')
	free_space_focal_length = ABCD.free_space(f, n)
	thin_lens_mat		= ABCD.thin_lens(f)
	#OpticalSystem = free_space_focal_length @ thin_lens_mat @ free_space_focal_length
	OpticalSystem = np.matmul(np.matmul(ABCD.free_space(f, n), ABCD.thin_lens(f)), ABCD.free_space(f, n))
	q_0 = q_hat(z_0, z_R_0)
	q_1 = propagate_q_hat(q_0, OpticalSystem)
	assert(np.real(q_0.magnitude) == approx(np.real(q_1.magnitude), abs=1e-6))
	
def test_curved_dielectric():
	# Lensmaker equation for same radius of curvature is given by
	# (1 / f) = (2 / R) - (n - 1)^2 * d / (n * R^2)
	n = 1.2 + rand.rand() * 0.5
	R = (15 + 10 * rand.rand()) * ureg.mm
	d = 1e-5 * rand.rand() * ureg.mm
	f = (((n - 1) * ((2 / R) - (d * (n - 1) /  n / R**2)))**-1)
	front_curve 	= ABCD.curved_interface(-R, 1, n)
	distance 	= ABCD.free_space(d, n)
	back_curve 	= ABCD.curved_interface(R, n, 1)
	thin_lens_approx = np.matmul(np.matmul(back_curve, distance), front_curve)
	thin_lens_approx = back_curve @ distance @ front_curve
	thin_lens_real = ABCD.thin_lens(f)
	assert(thin_lens_approx[0, 0] == approx(thin_lens_real[0, 0], rel=1e-5))
	assert(thin_lens_approx[0, 1] == approx(thin_lens_real[0, 1], abs=1e-5))
	assert(thin_lens_approx[1, 0] == approx(thin_lens_real[1, 0], rel=1e-5))
	assert(thin_lens_approx[1, 1] == approx(thin_lens_real[1, 1], abs=1e-5))

# maybe cauchy coefficients
def test_sellmeier_coefficients():
	B = np.array([1, 1, 1])
	C = np.array([3, 3, 3]) * ureg.um**2
	_lambda_1 = 0 * ureg.um
	_lambda_2 = 2 * ureg.um
	assert(sellmeier(B, C, _lambda_1) == 1)
	assert(sellmeier(B, C, _lambda_2) == approx(np.sqrt(13)))

def test_sellmeier_vitreous():
	_lambda_1 = 1400 * ureg.nanometer
	_lambda_2 = 400 * ureg.nanometer
	n1 = sellmeier_vitreous(_lambda_1)
	n2 = sellmeier_vitreous(_lambda_2)
	assert(n1 == approx(1.31326094))
	assert(n2 == approx(1.351098373))
