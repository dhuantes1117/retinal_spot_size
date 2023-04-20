import numpy as np
import csv
from retinal_spot_size import *
from .utils import *
from .eye_ABCD import *
from .eye_properties import *
import typer
import pathlib
import json
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline as interp

"""
Pseudo:
Goal: a list of the error in retinal spot size caused by slight shifts in each parameter
Questions: What does 'error in retinal spot size' mean?
	for now, put into function
Questions: What parameters are we shifting
	definitely z (by varying z we vary the waist location)
	Radius of curvture in reduced eye model
	Refractive index distribution
	maybe eye_length but I'd have to add a parameter to reduced_eye()
Questions: How do I run the sim without changing one thing
	make a function!


Algorithm:
	run calculation with nothing varied and store baseline distribution
	run calculation with M = 1.1 and store baseline distribution
	for each parameter we want to vary
		run calculation and then store (dictionary?)
	// One thing of interest would be plotting, but not of substance for numerical calculations
	for each parameter we varied
		error(varied_results, true_results) and store in dictionary(?)
	error(M=1.2_results, true_results) and store in dictionary(?)
	print neatly resultant errors and see if M has a significant effect (2 error bars)
	
"""
def percent_error_calc():
	N = 100
	_lambda_lower_bound = Q_(545, 'nm')
	_lambda_upper_bound = Q_(645, 'nm')
	_lambda_range = np.linspace(_lambda_lower_bound, _lambda_upper_bound, N)

	z0 			= Q_(0, 'm')
	eye_ROC 	= Q_(6.1, 'mm')
	eye_length 	= Q_(2.4, 'cm')
	M = 1
	
	Delta_z0 			= Q_(10, 'mm')
	Delta_eye_ROC 		= eye_ROC * 0.01
	Delta_eye_length = eye_length * 0.01

	Nominal_Rad_M1 = generate_radius(z0, eye_ROC, eye_length, 0, M)
	Nominal_Rad_M1.2 = generate_radius(z0, eye_ROC, eye_length, 0, 1.2)
	
	err_dict = {"z0": Delta_z0, "eye_ROC": Delta_eye_ROC, "eye_length": Delta_eye_length, "n_perc": 0.01}
	adj_matrix = np.diag(list(err_dict))
	
	for label, adj_vector in zip(err_dict.keys(), adj_matrix):
		e1, e2, e3, e4 = adj_vector
		r_percent = generate_radius(z0 + e1, eye_ROC + e2, eye_length + e3, e4, M)
		with open(f"{label}.csv", 'w') as f
			writer = csv.writer(f, delimiter=' ')
			writer.writerows(zip(_lambda_range, Nominal_Rad_M1, r_percent))

	

def calc_err(A, B):
	return 0.1
	pass

def generate_radius(z0 = Q_(0, 'm'), eye_ROC = Q_(6.1, 'mm'), eye_length, n_perc, M=1): # OR error in B, C
	"""
	[rss_figure_8 description]
	
	Parameters
	-----------
	:	Type
		Description
	
	Returns
	----------
	Type
		Description
	
	"""
	# Wavelength range discretization
	N = 100
	_lambda_lower_bound = Q_(545, 'nm')
	_lambda_upper_bound = Q_(645, 'nm')
	_lambda_range = np.linspace(_lambda_lower_bound, _lambda_upper_bound, N)

	# Initialize physical constants using pint
	beam_waist = Q_(3.24 / 2 / M, 'mm')

	# Initialize lists for beam radius and radius of curvature
	radius_list = []
	roc_list = []
	for _lambda in _lambda_range:
		# Physical parameters that depend on lambda
		_ang_div 		 = ang_div(_lambda, beam_waist)
		rayleigh_range 	 = z_R(_lambda, _ang_div)
		reduced_eye_ABCD = reduced_eye(_lambda, eye_ROC, eye_length)
		# Initialize reduced complex beam parameter
		q_0 = q_hat(z0, rayleigh_range)
		# Propogate reduced complex beam parameter
		q_1 = propagate_q_hat(q_0, reduced_eye_ABCD)
		# Find radius of curvature and radius
		retinal_roc = current_beam_roc(sellmeier_vitreous(_lambda), q_1)
		retinal_rad = current_beam_rad(_lambda, q_1)
		retinal_rad.ito("um")
		# Plotting goes easier with non-Quantities so units are stripped, and rad is scaled up
		radius_list.append(M * retinal_rad.magnitude)
		roc_list.append(retinal_roc.magnitude)
	radius_list = np.array(radius_list)
	roc_list = np.array(roc_list)
	return np.array([_lambda_range, radius_list])

	plt.plot(_lambda_range.magnitude, 2 * radius_list, label=f"M$^2$ = {M**2:.3}")		
	plt.title("Retinal Diameter Variation with Wavelength")#, waist size %f %s" %(beam_waist.magnitude, beam_waist.units))
	plt.xlabel(f"wavelength, ({_lambda_range.units})")
	plt.ylabel(f"retinal diameter, ({retinal_rad.units})")
	plt.legend()
	plt.show()
