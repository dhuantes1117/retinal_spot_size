from retinal_spot_size import *
from .eye_properties import *

def complex_matmul(q_hat, M):
	"""
	preforms the calculation of the change in reduced complex beam parameter
	following propagation through an optical system represented by ABCD matrix 'M'.

	Parameters
	-----------
	q_hat: 	pint.Quantity
			The reduced complex beam parameter of the incoming beam, expressed in units of length.
	M: 		numpy.ndarray
			Unitless ABCD matrix representing the optical system of interest.

	Returns
	----------
	pint.Quantity
		The reduced complex beam parameter of the outgoing beam, expressed in units of length.

	"""
	A = M[0, 0]
	B = Q_(M[0, 1], 'm')
	C = Q_(M[1, 0], 'm^-1')
	D = M[1, 1]
	assert(q_hat.check('[length]'))
	NUM = A * q_hat + B
	DEN = C * q_hat + D
	q_prime = NUM / DEN
	assert(q_prime.check('[length]'))
	return q_prime

@ureg.wraps(ureg.meter, (ureg.meter, ureg.meter, None))
def q_hat(z, _z_R, n=1):
	"""
	function for generating reduced complex beam parameters
	
	Parameters
	-----------
	z: 		pint.Quantity
			The displacement from the beam waist in the propagation direction,
			expressed in units of length
	_z_R: 	pint.Quantity
			The rayleigh range for the beam, expressed in units of length
	n:		pint.Quantity, default=1
			The refractive index of the media the beam is in as a unitless
			Quantity OR a float (default value indicates propagation in air, or
			vaccuum)
			
	Returns
	----------
	pint.Quantity
		The reduced complex beam parameter of the beam at the present location
	"""
	return z + 1j*_z_R / n

def propagate_q_hat(q_hat, System_Matrix):
	"""
	function for generating the reduced complex beam parameter immediately
	following propagation through an optical system.

	Parameters
	-----------
	q_hat: 			pint.Quantity 
					The reduced complex beam parameter immediately preceding
					entering the optical system, expressed as a pint Quantity
					with units of length
	System_Matrix: 	numpy.ndarray
					The ABCD matrix corresponding to the physical parameters of the optical system

	Returns
	----------
	pint.Quantity 
		The reduced complex beam parameter of the beam exiting the optical
		system expressed as a pint Quantity with units of length

	"""
	assert(q_hat.check('[length]'))
	return complex_matmul(q_hat, System_Matrix)

# IS THIS 1 / e or 1 / e**2
@ureg.wraps(ureg.meter, (ureg.meter, ureg.meter))
def current_beam_rad(lambda_, q_hat):
	"""
	calculates and returns the value for the beam radius at the current
	position based on wavelength and the complex beam parameter
	
	Parameters
	-----------
	lambda_:	pint.Quantity
		wavelenth of the beam IN VACCUUM, expressed as a pint Quantity with units of length. Important for not double accounting for refractive index
	q_hat:	pint.Quantity
		the reduced complex beam parameter at the present location expressed as a pint Quantity with units of length
	
	Returns
	----------
	pint.Quantity
		the radius of the beam at the present location, expressed as a pint Quantity with units of length 
	
	"""
	Imag = np.imag(1 / q_hat)
	waist_sq = -lambda_ / pi / Imag
	return np.sqrt(waist_sq)

# IS THIS 1 / e or 1 / e**2
@ureg.wraps(ureg.meter, (None, ureg.meter))
def current_beam_roc(n, q_hat):
	"""
	calculates the current beam radius of curvature
	
	Parameters
	-----------
	n:	float
		refractive index of the medium the beam is immersed in
		
	q_hat:	Quantity
		The complex beam parameter at the current location expressed as a pint
		Quantity with units of length
	
	Returns
	----------
	Quantity
		The beam's current radius of curvature expressed as a pint Quantity with units of length
	
	"""
	real_ = np.real(n / q_hat)
	return 1 / real_

# 1 / e**2
@ureg.wraps(ureg.meter, (ureg.meter, ureg.rad))
def z_R(lambda_, div):
	"""
	calculates the rayleigh range based on wavelength and angular divergence
	
	Parameters
	-----------
	lambda_:	Quantity
		wavelenth of the beam, expressed as a pint Quantity with
		units of length. 
	div:	Quantity
		angular divergence of the beam, expressed as a pint Quantity with units of radians (anything angular really)
	
	Returns
	----------
	Quantity
		The Rayleigh range for the beam expressed as a pint Quantity with units
		of length	

	"""
	return lambda_ / pi / div**2

# 1 / e**2
@ureg.wraps(ureg.rad, (ureg.meter, ureg.meter))
def ang_div(_lambda, waist):
	"""
	calculates the theoretical angular divergence of a beam with a given
	wavelength and beam waist
	
	Parameters
	-----------
	_lambda:	Quantity
		wavelenth of the beam, expressed as a pint Quantity with units of
		length. 
	waist:	Quantity
		beam waist radius, expressed as a pint Quantity with units of length
			
	Returns
	----------
	Quantity
		the theoretical beam angular divergence, expressed as a pint Quantity with units of angular-ity (radians)
			
	"""
	return _lambda / pi / waist

@ureg.wraps(None, (None, None, None, ureg.um**2, ureg.um**2, ureg.um**2, ureg.um))
def sellmeier_true(B1, B2, B3, C1, C2, C3, _lambda):
	"""
		Sellmeier_true calculates the refractive index of a given wavelength in
		a given material with sellmeier coefficients given by the first 6
		parameters. This function is 'true' because it checks for units, and
		should not be called directly in favor of
		'retinal_spot_size.utils.sellmeier' and
		'retinal_spot_size.utils.sellmeier_vitreous'
	
	Parameters
	-----------
	B1:	float
		B1 Sellmeier coefficient, unitless
	B2:	float
		B2 Sellmeier coefficient, unitless
	B3:	float
		B3 Sellmeier coefficient, unitless
	C1:	Quantity
		C1 Sellmeier coefficient, expressed as a pint Quantity with units of
		square micron, or length squared in general
	C2:	Quantity
		C2 Sellmeier coefficient, expressed as a pint Quantity with units of
		square micron, or length squared in general
	C3:	Quantity
		C3 Sellmeier coefficient, expressed as a pint Quantity with units of
		square micron, or length squared in general
	_lambda:	Quantity
		wavelenth of the beam, expressed as a pint Quantity with units of
		length. 
		
	
	Returns
	----------
	float
		The refractive index of the given material accounting for chromatic abberation
	
	"""
	T0 = 1
	T1 = B1 * _lambda**2 / (_lambda**2 - C1) 
	T2 = B2 * _lambda**2 / (_lambda**2 - C2)
	T3 = B3 * _lambda**2 / (_lambda**2 - C3)
	return np.sqrt(T0 + T1 + T2 + T3)

def sellmeier(B, C, _lambda):
	"""
	Function wrapper for 'retinal_spot_size.utils.sellmeier_true' that takes
	the B and C coefficients as an array-like input and returns the refractive
	index at the provided wavelength
	
	Parameters
	-----------
	B:	numpy.ndarray
		a 1D numpy array whose first 3 entries are the first three sellmeier B
		coefficients
	C:	numpy.ndarray
		a 1D numpy array whose first 3 entries are the first three sellmeier C
		coefficients, expressed as a pint Quantity-s with units of length
		squared. All units will be converted to square micron before
		calculations are preformed
	_lambda:	Quantity
		wavelenth of the beam, expressed as a pint Quantity with units of
		length. 
		
	
	Returns
	----------
	float
		The refractive index of the given wavelength accounting for chromatic
		abberation with the Sellmeier equation
	
	"""
	return sellmeier_true(*B, *C, _lambda)

def sellmeier_vitreous(_lambda):
	"""
	Function wrapper for 'retinal_spot_size.utils.sellmeier_true' that takes
	the B and C coefficients for the refractive index from the data collected
	by Liou et al. and Fernandez et al. in Vincelette et al. for the sellmeier
	coefficients for the human eye and returns the refractive index at the
	provided wavelength. Name is a misnomer, as the fit they do is for the eye
	as a whole, not just the vitreous
	
	Parameters
	-----------
	_lambda:	Quantity
		wavelenth of the beam, expressed as a pint Quantity with units of
		length. 
	
	Returns
	----------
	float
		The net refractive index of the given wavelength inside the reduced eye
		model, accounting for chromatic abberation of the human eye	
	"""
	return sellmeier_true(*B_COEFFICIENTS_VITREOUS, *C_COEFFICIENTS_VITREOUS, _lambda)
