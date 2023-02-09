from retinal_spot_size import *
from .utils import *
# return the ABCD Matrix for a thin lens w focal length f 
@ureg.wraps(None, [ureg.meter])
def thin_lens(f):
        return np.array([[1    , 0],
                        [-1 / f, 1]])

@ureg.wraps(None, [ureg.meter, None])
def free_space(L, n_0):
        return np.array([[1, L / n_0],
                         [0,       1]])

@ureg.wraps(None, [ureg.meter, None, None])
def curved_interface(R, n_1, n_2):
        return np.array([[               	1, 	  0],
                         [(n_2 - n_1) / R 	 ,	  1]])



