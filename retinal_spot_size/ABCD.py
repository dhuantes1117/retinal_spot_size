import numpy as np
from .utils import *

# return the ABCD Matrix for a thin lens w focal length f 
def thin_lens(f):
        return np.array([[1    , 0],
                        [-1 / f, 1]])

def free_space(L, n_0):
        return np.array([[1, L / n_0],
                         [0,       1]])

def curved_interface(R, n_1, n_2):
        return np.array([[               	1, 	  0],
                         [(n_1 - n_2) / R / n_2 , n_1 / n_2]])

def propagate_q(q_0, System_Matrix):
        return complexMatmul(q_0, System_Matrix)

