from retinal_spot_size import *
from .utils import *

# return the ABCD Matrix for a thin lens w focal length f 
@ureg.wraps(None, [ureg.meter])
def thin_lens(f):
		"""
		function to generate ABCD matrix representing a thin lens.

		Parameters
		-----------
		f: 	Quantity 
			Focal length of thin lens expressed as a Pint Quantity with dimensions of length.

		Returns
		----------
		numpy ndarray
			ABCD matrix representing a thin lens with focal length 'f' with C in units of m^-1.

		Note: Function is wrapped, so arguments are converted to expected units and stripped of their units before calculations.

		"""
        return np.array([[1    , 0],
                        [-1 / f, 1]])

@ureg.wraps(None, [ureg.meter, None])
def free_space(L, n_0):
		"""
		function to generate ABCD matrix representing free space propagation.

		Parameters
		-----------
		L: 		Quantity 
				Length of propagation expressed as a Pint Quantity with dimensions of length.
		n_0:	Quantity 
				Refractive index of media expressed as a dimensionless Pint Quantity OR float.

		Returns
		----------
		numpy ndarray
			ABCD matrix representing propagation through a media of refractive index 'n_0' over a distance 'L'.

		Note: Function is wrapped, so arguments are converted to expected units and stripped of their units before calculations.

		"""
        return np.array([[1, L / n_0],
                         [0,       1]])

@ureg.wraps(None, [ureg.meter, None, None])
def curved_interface(R, n_1, n_2):
		"""
		function to generate ABCD matrix representing free space propagation.

		Parameters
		-----------
		R: 		Quantity 
				Radius of curvature of interface expressed as a Pint Quantity
				with dimensions of length. Follows sign convention of Siegmans, AKA for
				propagation from left to right ) has 'R'>0 and ( has 'R'<0.
		n_1:	Quantity 
				Refractive index of media the light exits expressed as a dimensionless Pint Quantity OR float.
		n_2:	Quantity 
				Refractive index of media the light enters expressed as a dimensionless Pint Quantity OR float.

		Returns
		----------
		numpy ndarray
			ABCD matrix representing propagation through an interface of radius of curvature 'R' from media of refractive index 'n_1' to 'n_2'.

		Note: Function is wrapped, so arguments are converted to expected units and stripped of their units before calculations.

		"""
        return np.array([[               	1, 	  0],
                         [(n_2 - n_1) / R 	 ,	  1]])



