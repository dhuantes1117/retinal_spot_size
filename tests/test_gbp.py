from retinal_spot_size import ABCD
from retinal_spot_size.utils import *
from pytest import *
import numpy as np
import numpy.random as rand

def test_construct_thin_lens_matrix():
	mat = ABCD.thin_lens(2)
	assert mat[0,0] == 1
	assert mat[0,1] == 0
	assert mat[1,0] == -0.5
	assert mat[1,1] == 1
	
def test_construct_free_space_matrix():
	mat = ABCD.free_space(6, 7)
	assert mat[0,0] == 1
	assert mat[0,1] == approx(6/7)
	assert mat[1,0] == 0
	assert mat[1,1] == 1

def test_construct_curved_interface_matrix():
	mat = ABCD.curved_interface(1, 2, 3)
	assert mat[0,0] == 1
	assert mat[0,1] == 0
	assert mat[1,0] == -1 / 3
	assert mat[1,1] == 2 / 3

def test_matmul():
	q_0 = q(rand.rand(), rand.rand())
	q_1 = complex_matmul(q_0, np.identity(2))
	assert(q_0 == q_1)

def test_free_space():
	L = 10 * rand.rand()
	n = rand.rand()
	q_0 = q(rand.rand(), rand.rand())
	q_1 = propagate_q(q_0, ABCD.free_space(L, n))
	q_2 = propagate_q(q_1, ABCD.free_space(-L, n))
	assert(q_0 == approx(q_2))

# update to involve every matrix
def test_optical_train():
	f = 10 * rand.rand()
	n = 1 #+ rand.rand()
	free_space_focal_length = ABCD.free_space(f, n)
	thin_lens_mat		= ABCD.thin_lens(f)
	OpticalSystem = np.matmul(np.matmul(ABCD.free_space(f, n), ABCD.thin_lens(f)), ABCD.free_space(f, n))
	assert(approx(np.linalg.det(OpticalSystem)) == 1)

def test_thin_lens():
	f = 10 * rand.rand()
	n = 1 #+ rand.rand()
	free_space_focal_length = ABCD.free_space(f, n)
	thin_lens_mat		= ABCD.thin_lens(f)
	#OpticalSystem = free_space_focal_length @ thin_lens_mat @ free_space_focal_length
	OpticalSystem = np.matmul(np.matmul(ABCD.free_space(f, n), ABCD.thin_lens(f)), ABCD.free_space(f, n))
	q_0 = q(0, rand.rand())
	q_1 = propagate_q(q_0, OpticalSystem)
	assert(np.real(q_0) == approx(np.real(q_1), abs=1e-6))
	
def test_curved_dielectric():
	# Lensmaker equation for same radius of curvature is given by
	# (1 / f) = (2 / R) - (n - 1)^2 * d / (n * R^2)
	n = 1.2 + rand.rand() * 0.5
	R = 15 + 10 * rand.rand()
	d = 1e-5 * rand.rand()
	f = ((n - 1) * ((2 / R) - (d * (n - 1) /  n / R**2)))**-1
	front_curve = ABCD.curved_interface(R, 1, n)
	distance = ABCD.free_space(d, n)
	back_curve = ABCD.curved_interface(-R, n, 1)
	thin_lens_approx = np.matmul(np.matmul(back_curve, distance), front_curve)
	thin_lens_approx = back_curve @ distance @ front_curve
	thin_lens_real = ABCD.thin_lens(f)
	assert(thin_lens_approx[0, 0] == approx(thin_lens_real[0, 0], rel=1e-5))
	assert(thin_lens_approx[0, 1] == approx(thin_lens_real[0, 1], abs=1e-5))
	assert(thin_lens_approx[1, 0] == approx(thin_lens_real[1, 0], rel=1e-5))
	assert(thin_lens_approx[1, 1] == approx(thin_lens_real[1, 1], abs=1e-5))

