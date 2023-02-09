from retinal_spot_size import *
from retinal_spot_size.ABCD import *
from retinal_spot_size.eye_properties import *
from retinal_spot_size.utils import *
import matplotlib.pyplot as plt


M_front_cornea  = curved_interface(-1*CORNEA_RAD_FRONT, N_AIR, N_CORNEA)
M_cornea 	= free_space(CORNEAL_THICKNESS, N_CORNEA)
M_back_cornea 	= curved_interface(-1*CORNEA_RAD_BACK, N_CORNEA, N_VITREOUS)
M_to_lens 	= free_space(CORNEA_TO_LENS_DISTANCE, N_VITREOUS) # AQEUOUS
M_lens_front	= curved_interface(-1*LENS_RAD_FRONT, N_VITREOUS, N_LENS)
M_lens_prop	= free_space(LENS_THICKNESS, N_LENS)
M_lens_back	= curved_interface(LENS_RAD_BACK, N_LENS, N_VITREOUS)
M_lens 		= np.matmul(M_lens_front, np.matmul(M_lens_prop, M_lens_back))
M_lens 		= thin_lens(FOCAL_LENGTH)
M_vitreous 	= free_space(CORNEA_TO_LENS_DISTANCE, N_VITREOUS)

optical_system = [M_front_cornea, M_cornea, M_back_cornea, M_to_lens, M_lens, M_vitreous]

eye_equiv_ABCD = np.identity(2)
for interface in optical_system[::-1]:
	eye_equiv_ABCD = np.matmul(interface, eye_equiv_ABCD)


_lambda = Q_(532e-9, 'm')
_ang_div = 10-3 * ureg.rad


def test_prop():
	q_0 = q_hat(0 * ureg.meter, z_R(_lambda, _ang_div))
	q_1 = complex_matmul(q_0, eye_equiv_ABCD)

