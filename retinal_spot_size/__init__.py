from scipy.constants import pi
import numpy as np
# does this cause recursion in utils file? no clue
# yup, it did
import pint

ureg = pint.UnitRegistry()
Q_ = ureg.Quantity

