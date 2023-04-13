from retinal_spot_size import *
from .utils import *
# NEED TO SOURCE THESE CONSTANTS. Some from Sliney, some from paper
N_AIR 	   = 1    
N_CORNEA   = 1.376
N_VITREOUS = 1.337
N_LENS 	   = 1.4

CORNEA_RAD_FRONT = 7.98  * ureg.mm
CORNEA_RAD_BACK  = 6.22  * ureg.mm
LENS_RAD_FRONT   = 10.20 * ureg.mm
LENS_RAD_BACK    = 6.17  * ureg.mm
FOCAL_LENGTH     = 22.89 * ureg.mm

CORNEAL_THICKNESS 	= 1.15  * ureg.mm
CORNEA_TO_LENS_DISTANCE = 2.39  * ureg.mm # AQUEOUS
LENS_TO_RETINA_DISTANCE = 17.15 * ureg.mm 
EYE_LENGTH 		= 24.75 * ureg.mm 
LENS_THICKNESS 		= 4.06  * ureg.mm

VITREOUS_LENGTH 	= 14.75 * ureg.mm
B_COEFFICIENTS_VITREOUS = [7.516e-1, -4.484e-3, -1.503 * 10] # unitless
C_COEFFICIENTS_VITREOUS = [1.641e-2 * ureg.um**2, 8.596e-2 * ureg.um**2, -1.028e3 * ureg.um**2] # probably um**2
