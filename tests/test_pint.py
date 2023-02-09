from retinal_spot_size import *
from retinal_spot_size import ABCD
from retinal_spot_size.utils import *
from pytest import *
import numpy.random as rand

def test_pint():
	A = (1 + 2j) * ureg.meter
	assert(A.units == ureg.meter) 
