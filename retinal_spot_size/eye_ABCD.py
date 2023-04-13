from retinal_spot_size import *
from .utils import *
from .eye_properties import *
from .ABCD import *

def reduced_eye(_lambda, R):
	"""
	should return an ABCD matrix representing the reduced eye model with
	chromatic aberration, and corneal radius of curvature of R
		
	Parameters
	-----------
	_lambda:	Quantity
		The wavelength of the light, expressed as a pint Quantity with units of
		length
	R:	Quantity
		The corneal radius of cornea, expressed as a pint Quantity with units
		of length
	
	Returns
	----------
	numpy.ndarray
		ABCD matrix representing the air-cornea interface and traversal to
		retina
	
	"""
	M_front_cornea  = curved_interface(-1*R, N_AIR, sellmeier_vitreous(_lambda))
	M_cornea 	= free_space(CORNEAL_THICKNESS, sellmeier_vitreous(_lambda))
	M_back_cornea 	= curved_interface(-1*CORNEA_RAD_BACK, sellmeier_vitreous(_lambda), sellmeier_vitreous(_lambda))
	M_to_lens 	= free_space(CORNEA_TO_LENS_DISTANCE, sellmeier_vitreous(_lambda)) # AQEUOUS
	M_lens_front	= curved_interface(-1*LENS_RAD_FRONT, sellmeier_vitreous(_lambda), sellmeier_vitreous(_lambda))
	M_lens_prop	= free_space(LENS_THICKNESS, sellmeier_vitreous(_lambda))
	M_lens_back	= curved_interface(LENS_RAD_BACK, sellmeier_vitreous(_lambda), sellmeier_vitreous(_lambda))
	M_lens 		= np.matmul(M_lens_back, np.matmul(M_lens_prop, M_lens_front))
	M_vitreous 	= free_space(LENS_TO_RETINA_DISTANCE, sellmeier_vitreous(_lambda))
	M_vitreous 	= free_space(CORNEAL_THICKNESS + CORNEA_TO_LENS_DISTANCE + LENS_THICKNESS + LENS_TO_RETINA_DISTANCE, sellmeier_vitreous(_lambda))
	#print(f"Eye_Length: {CORNEAL_THICKNESS + CORNEA_TO_LENS_DISTANCE + LENS_THICKNESS + LENS_TO_RETINA_DISTANCE}")
	M_vitreous 	= free_space(2.44 * ureg.cm, sellmeier_vitreous(_lambda))

	optical_system = [M_front_cornea, M_cornea, M_back_cornea, M_to_lens, M_lens_front, M_lens_prop, M_lens_back, M_vitreous]
	optical_system = [M_front_cornea, M_vitreous]

	ABCD = np.identity(2)
	for interface in optical_system:
		ABCD = np.matmul(interface, ABCD)
	return ABCD

def eye_equiv_ABCD(_lambda):
	"""
	[eye_equiv_ABCD description]
	
	Parameters
	-----------
	_lambda:	Type
		Description
	
	Returns
	----------
	Type
		Description
	
	"""
	#return just_cornea(_lambda)
	M_front_cornea  = curved_interface(-1*CORNEA_RAD_FRONT, N_AIR, N_CORNEA)
	M_cornea 	= free_space(CORNEAL_THICKNESS, N_CORNEA)
	M_back_cornea 	= curved_interface(-1*CORNEA_RAD_BACK, N_CORNEA, sellmeier_vitreous(_lambda))
	M_to_lens 	= free_space(CORNEA_TO_LENS_DISTANCE, sellmeier_vitreous(_lambda)) # AQEUOUS
	M_lens_front	= curved_interface(-1*LENS_RAD_FRONT, sellmeier_vitreous(_lambda), N_LENS)
	M_lens_prop	= free_space(LENS_THICKNESS, N_LENS)
	M_lens_back	= curved_interface(LENS_RAD_BACK, N_LENS, sellmeier_vitreous(_lambda))
	M_lens 		= np.matmul(M_lens_back, np.matmul(M_lens_prop, M_lens_front))
	M_vitreous 	= free_space(LENS_TO_RETINA_DISTANCE, sellmeier_vitreous(_lambda))

	# M_front_cornea  = curved_interface(-1*CORNEA_RAD_FRONT, N_AIR, sellmeier_vitreous(_lambda))
	# M_vitreous 	= free_space(EYE_LENGTH/2, sellmeier_vitreous(_lambda))

	optical_system = [M_front_cornea, M_cornea, M_back_cornea, M_to_lens, M_lens_front, M_lens_prop, M_lens_back, M_vitreous]
	# optical_system = [M_front_cornea, M_cornea, M_back_cornea, M_to_lens, M_lens, M_vitreous]
	# optical_system = [M_front_cornea, M_vitreous]

	ABCD = np.identity(2)
	for interface in optical_system:
		ABCD = np.matmul(interface, ABCD)
	return ABCD


def eye_lens_ABCD(_lambda):
	"""
	[eye_lens_ABCD description]
	
	Parameters
	-----------
	_lambda:	Type
		Description
	
	Returns
	----------
	Type
		Description
	
	"""
	M_lens_front	= curved_interface(-1*LENS_RAD_FRONT, sellmeier_vitreous(_lambda), N_LENS)
	M_lens_prop	= free_space(LENS_THICKNESS, N_LENS)
	M_lens_back	= curved_interface(LENS_RAD_BACK, N_LENS, sellmeier_vitreous(_lambda))

	optical_system = [M_lens_front, M_lens_prop, M_lens_back]

	ABCD = np.identity(2)
	for interface in optical_system:
		ABCD = np.matmul(interface, ABCD)
	return ABCD

