from retinal_spot_size import *
from .utils import *
from .eye_ABCD import *
from .eye_properties import *
import typer
import pathlib
import json
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline as interp


app = typer.Typer()

@app.command()
def pds():
	print(sellmeier_vitreous.__doc__)

@app.command()
def rss_figure_8():
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
	z0 = Q_(0, 'm')
	for M in [m for m in np.linspace(1, 1.5, 3)]:
			beam_waist = Q_(3.24 / 2 / M, 'mm')

			# Initialize lists for beam radius and radius of curvature
			radius_list = []
			roc_list = []
			for _lambda in _lambda_range:
				# Physical parameters that depend on lambda
				_ang_div 		= ang_div(_lambda, beam_waist)
				rayleigh_range 	= z_R(_lambda, _ang_div)
				reduced_eye_ABCD = reduced_eye(_lambda, 6.1 * ureg.mm)
				# Initialize reduced complex beam parameter
				q_0 = q_hat(z0, rayleigh_range)
				# Propogate reduced complex beam parameter
				q_1 = propagate_q_hat(q_0, reduced_eye_ABCD)
				# Find radius of curvature and radius
				retinal_roc = current_beam_roc(sellmeier_vitreous(_lambda), q_1)
				retinal_rad = current_beam_rad(_lambda, q_1)
				retinal_rad.ito("um")
				# Plotting goes easier with non-Quantities so units are stripped
				radius_list.append(M*retinal_rad.magnitude)
				roc_list.append(retinal_roc.magnitude)
			radius_list = np.array(radius_list)
			roc_list = np.array(roc_list)
			# Code is generating a reflected, shifted version, these numbers
			# qualitatively fit what is happening, so the question is: What is
			# causing our function to have this issue?
			plt.plot(_lambda_range.magnitude, 2 * radius_list, label=f"M$^2$ = {M**2:.3}")		
			#plt.plot(_lambda_range.magnitude, 2 * roc_list)		
	plt.title("Retinal Radius Variation with Wavelength")#, waist size %f %s" %(beam_waist.magnitude, beam_waist.units))
	plt.xlabel(f"wavelength, ({_lambda_range.units})")
	plt.ylabel(f"retinal radius, ({retinal_rad.units})")
	plt.legend()
	plt.show()


@app.command()
def rf8():
	"""
	[rf8 description]
	
	Parameters
	-----------
	:	Type
		Description
	
	Returns
	----------
	Type
		Description
	
	"""
	import pudb
	#pu.db
	rss_figure_8()	
	

@app.command()
def spot_size_z_range():
	"""
	[spot_size_z_range description]
	
	Parameters
	-----------
	:	Type
		Description
	
	Returns
	----------
	Type
		Description
	
	"""
	_lambda = Q_(700, 'nm')
	beam_waist = Q_(100, 'um')
	z0_range = np.linspace(-50, 50, 100) * ureg.cm
	div = ang_div(_lambda, beam_waist)
	radius_list = np.zeros(len(z0_range))
	for i, z0 in enumerate(z0_range):
		q_0 = q_hat(z0, z_R(_lambda, div))
		q_1 = complex_matmul(q_0, reduced_eye(_lambda, 6.1 * ureg.mm))
		beam_waist = current_beam_rad(_lambda, q_0) 
		retinal_rad = current_beam_rad(_lambda, q_1)
		retinal_rad.ito('um')
		radius_list[i] = retinal_rad.magnitude
	plt.plot(z0_range, radius_list)		
	plt.title("retinal radius as function of waist location, waist size %f (%s)" %(beam_waist.magnitude, beam_waist.units))
	plt.show()


@app.command()
def ssz():
	"""
	[ssz description]
	
	Parameters
	-----------
	:	Type
		Description
	
	Returns
	----------
	Type
		Description
	
	"""
	spot_size_z_range()	

@app.command()
def spot_size_w_range():
	"""
	[spot_size_w_range description]
	
	Parameters
	-----------
	:	Type
		Description
	
	Returns
	----------
	Type
		Description
	
	"""
	_lambda = Q_(700, 'nm')
	z0 = 0 * ureg.mm 
	w0_range = np.linspace(10, 1000, 100) * ureg.um
	radius_list = np.zeros(len(w0_range))
	for i, beam_waist in enumerate(w0_range):
		div = ang_div(_lambda, beam_waist)
		q_0 = q_hat(z0, z_R(_lambda, div))
		q_1 = complex_matmul(q_0, reduced_eye(_lambda, 6.1 * ureg.mm))
		retinal_rad = current_beam_rad(_lambda, q_1)
		retinal_rad.ito('um')
		radius_list[i] = retinal_rad.magnitude
	plt.plot(w0_range, radius_list)		
	plt.title("retinal radius as function of waist size, z0 %f (%s)" %(z0.magnitude, z0.units))
	plt.xlabel("waist size, (%s)" %(beam_waist.units))
	plt.ylabel("retinal radius, (%s)" %(retinal_rad.units))
	plt.show()


@app.command()
def ssw():
	"""
	[ssw description]
	
	Parameters
	-----------
	:	Type
		Description
	
	Returns
	----------
	Type
		Description
	
	"""
	spot_size_w_range()	

@app.command()
def display_sizes():
	"""
	[display_sizes description]
	
	Parameters
	-----------
	:	Type
		Description
	
	Returns
	----------
	Type
		Description
	
	"""
	_lambda = 700 * ureg.nm
	M_lens = reduced_eye(_lambda, 6.1 * ureg.mm)
	print(M_lens)
	print("f equiv for lens: %f" %(1 / M_lens[1][0]))
	M_front_cornea  = curved_interface(-1*CORNEA_RAD_FRONT, N_AIR, N_CORNEA)
	M_cornea 	= free_space(CORNEAL_THICKNESS, N_CORNEA)
	M_back_cornea 	= curved_interface(-1*CORNEA_RAD_BACK, N_CORNEA, sellmeier_vitreous(_lambda))
	M_to_lens 	= free_space(CORNEA_TO_LENS_DISTANCE, sellmeier_vitreous(_lambda)) # AQEUOUS
	M_lens_front	= curved_interface(-1*LENS_RAD_FRONT, sellmeier_vitreous(_lambda), N_LENS)
	M_lens_prop	= free_space(LENS_THICKNESS, N_LENS)
	M_lens_back	= curved_interface(LENS_RAD_BACK, N_LENS, sellmeier_vitreous(_lambda))
	M_lens 		= np.matmul(M_lens_back, np.matmul(M_lens_prop, M_lens_front))
	M_vitreous 	= free_space(LENS_TO_RETINA_DISTANCE, sellmeier_vitreous(_lambda))

	optical_system = [M_front_cornea, M_cornea, M_back_cornea, M_to_lens, M_lens_front, M_lens_prop, M_lens_back, M_vitreous]

	ABCD = np.identity(2)
	z0 = Q_(0, 'mm')
	beam_waist = Q_(1, 'mm')
	_ang_div = ang_div(_lambda, beam_waist)
	q_0 = q_hat(z0, z_R(_lambda, _ang_div))
	waist = current_beam_rad(_lambda, q_0)
	waist.ito("um")
	radius_list = [waist.magnitude]
	for interface in optical_system[::-1]:
		q_1 = propagate_q_hat(q_0, interface)
		retinal_rad = current_beam_rad(_lambda, q_1)
		retinal_rad.ito('um')
		radius_list.append(retinal_rad.magnitude)
	
	radius_list = np.array(radius_list)
	plt.plot(np.linspace(0, 10, len(radius_list)), radius_list)		
	plt.title("retinal radius as function of wavelength, waist size %f %s" %(beam_waist.magnitude, beam_waist.units))
	plt.xlabel("interface")
	plt.ylabel("retinal radius, (%s)" %(retinal_rad.units))
	plt.show()

@app.command()
def ds():
	"""
	[ds description]
	
	Parameters
	-----------
	:	Type
		Description
	
	Returns
	----------
	Type
		Description
	
	"""
	display_sizes()

@app.command()
def propagation_progress():
	"""
	[propagation_progress description]
	
	Parameters
	-----------
	:	Type
		Description
	
	Returns
	----------
	Type
		Description
	
	"""
	_lambda = Q_(700, 'nm')
	beam_waist = Q_(51, 'um')
	Delta_z = Q_(1, 'mm')
	z0 = Q_(0, 'mm')
	L = Q_(20, 'mm')
	s = L / 2
	f = (pi**2 * beam_waist**4 / 2 / s / _lambda**2) + s / 2
	f = pi * beam_waist**2 /  _lambda
	f.ito('mm')
	L = f * 2
	f.ito('m')
	# f = -f
	# print(f / Q_(260, 'm'))
	# f = Q_(260, 'm')
	# f = Q_(260, 'm')
	#f = L * 0.5
	N = L / Delta_z
	N.ito("")
	N = int(N.magnitude)
	print(N)
	distance = np.linspace(-L.magnitude, L.magnitude, N)
	div = ang_div(_lambda, beam_waist)
	q_0 = q_hat(-L / 2, z_R(_lambda, div))
	n_water = 1.33
	ABCD_forward 		= free_space(Delta_z, 1)
	ABCD_forward_water 	= free_space(Delta_z, n_water)
	ABCD_thin_lens 	= thin_lens(f)
	ABCD_air_water 	= curved_interface(-Q_(6, 'mm'), 1, n_water)
	radii = []
	q_i = q_0
	print("focal length: %s" %(f))
	print("           s: %s" %(s))
	for i in range(N // 2):
		q_i = propagate_q_hat(q_i, ABCD_forward)
		beam_radius_i = current_beam_rad(_lambda, q_i)
		radii.append(beam_radius_i.magnitude)
	# q_i = propagate_q_hat(q_i, ABCD_thin_lens)
	q_i = propagate_q_hat(q_i, ABCD_forward)
	beam_radius_i = current_beam_rad(_lambda, q_i)
	radii.append(beam_radius_i.magnitude)
	for i in range(N // 2):
		q_i = propagate_q_hat(q_i, ABCD_forward)
		beam_radius_i = current_beam_rad(_lambda, q_i)
		radii.append(beam_radius_i.magnitude)
	
	radii = np.array(radii)
	plt.plot(distance, radii)
	plt.plot(distance, -(radii))
	plt.title("")
	plt.xlabel("distane (%s)" %L.units)
	plt.ylabel("radius (%s)" %beam_radius_i.units)
	plt.show()
	print("initial complex beam parameter: %s" %(q_0))
	print("  final complex beam parameter: %s" %(q_i))

@app.command()
def pp():
	"""
	[pp description]
	
	Parameters
	-----------
	:	Type
		Description
	
	Returns
	----------
	Type
		Description
	
	"""
	propagation_progress()

@app.command()
def wsq():
	"""
	[wsq description]
	
	Parameters
	-----------
	:	Type
		Description
	
	Returns
	----------
	Type
		Description
	
	"""
	for M, color in zip([1, 1.2], ['orange', 'blue']):
		_lambda = Q_(700, 'nm')
		beam_waist = Q_(53 / M, 'um')
		Delta_z = Q_(1, 'mm')
		z0 = Q_(0, 'mm')
		L = Q_(20, 'mm')
		s = L / 2
		f = (pi**2 * beam_waist**4 / 2 / s / _lambda**2) + s / 2
		f = pi * beam_waist**2 /  _lambda
		f.ito('mm')
		L = f * 2 * M**2
		f.ito('m')
		N = L / Delta_z
		N.ito("")
		N = int(N.magnitude)
		distance = np.linspace(0, L.magnitude, N)
		div = ang_div(_lambda, beam_waist)
		q_0 = q_hat(-L / 2, z_R(_lambda, div))
		q_0 = q_hat(0 * ureg.mm, z_R(_lambda, div))
		n_water = 1.33
		ABCD_forward 		= free_space(Delta_z, 1)
		ABCD_forward_water 	= free_space(Delta_z, n_water)
		ABCD_thin_lens 	= thin_lens(f)
		ABCD_air_water 	= curved_interface(-Q_(6, 'mm'), 1, n_water)
		radii = []
		q_i = q_0
		print("focal length: %s" %(f))
		print("           s: %s" %(s))
		for i in range(N // 2):
			q_i = propagate_q_hat(q_i, ABCD_forward)
			beam_radius_i = current_beam_rad(_lambda, q_i)
			radii.append(beam_radius_i.magnitude)
		# q_i = propagate_q_hat(q_i, ABCD_thin_lens)
		q_i = propagate_q_hat(q_i, ABCD_forward)
		beam_radius_i = current_beam_rad(_lambda, q_i)
		radii.append(beam_radius_i.magnitude)
		for i in range(N // 2):
			q_i = propagate_q_hat(q_i, ABCD_forward)
			beam_radius_i = current_beam_rad(_lambda, q_i)
			radii.append(beam_radius_i.magnitude)
		
		radii = np.array(radii) * 1000
		plt.plot(distance, radii*M, color, label=f"M={M*1.0:.2}")
		plt.plot(distance, -(radii*M), color)
		plt.plot(-1 * distance[::-1], radii[::-1]*M, color)
		plt.plot(-1 * distance[::-1], -(radii[::-1]*M), color)
	plt.title("Effect of Beam Quality Factor")
	plt.legend()
	L.ito("m")
	plt.xlabel("distance (%s)" %L.units)
	plt.ylabel("radius (%s)" %beam_radius_i.units)
	plt.show()
	print("initial complex beam parameter: %s" %(q_0))
	print("  final complex beam parameter: %s" %(q_i))

@app.command()
def w():
	"""
	[w(): description]
	
	Parameters
	-----------
	:	Type
		Description
	
	Returns
	----------
	Type
		Description
	
	"""
	wsq()

@app.command()
def sellmeier_figure_7():
	"""
	[sellmeier_figure_7 description]
	
	Parameters
	-----------
	:	Type
		Description
	
	Returns
	----------
	Type
		Description
	
	"""
	lambda_ = np.linspace(400, 1400, 1000) * ureg.nm
	n	= sellmeier_vitreous(lambda_)
	plt.plot(lambda_, n)
	plt.xlabel(f"$\lambda$ (Wavelength), ({lambda_[0].units})")
	plt.ylabel(f"n (refractive index)")
	plt.title (f"Refractive Index for Eye")
	plt.show()

@app.command()
def rf7():
	"""
	[rf7 description]
	
	Parameters
	-----------
	:	Type
		Description
	
	Returns
	----------
	Type
		Description
	
	"""
	sellmeier_figure_7()

@app.command()
def reduced_eye_ROC():
	"""
	[reduced_eye_ROC description]
	
	Parameters
	-----------
	:	Type
		Description
	
	Returns
	----------
	Type
		Description
	
	"""
	R_range = [6.1 * ureg.mm]
	R_range = np.linspace(6.16 * ureg.mm, 6.20 * ureg.mm, 5)
	spot_size_roc = []
	for cornea_roc in R_range:
		N = 50
		_lambda_lower_bound = Q_(8000, 'nm')
		_lambda_upper_bound = Q_(10400, 'nm')
		_lambda_lower_bound = Q_(400, 'nm')
		_lambda_upper_bound = Q_(1400, 'nm')
		_lambda_range = np.linspace(_lambda_lower_bound, _lambda_upper_bound, 100)
		z0 = Q_(0, 'm')
		beam_waist = Q_(4.24, 'mm')

		radius_list = []
		z_list = []
		for _lambda in _lambda_range:
			_ang_div = ang_div(_lambda, beam_waist)
			q_0 = q_hat(z0, z_R(_lambda, _ang_div))
			q_1 = propagate_q_hat(q_0, reduced_eye(_lambda, cornea_roc))
			#print(f"{_lambda} {q_1}")
			retinal_rad = current_beam_rad(_lambda / N_VITREOUS, q_1)
			retinal_rad.ito("um")
			radius_list.append(retinal_rad.magnitude)
		radius_list = np.array(radius_list)
		spot_size_roc.append(radius_list)
		# Code is generating a reflected, shifted version, these numbers
		# qualitatively fit what is happening, so the question is: What is
		# causing our function to have this issue?
		#plt.plot(-0.35 * (_lambda_range.magnitude - 11600), radius_list * 0.1)		
	for spot_size_i, R_i in zip(spot_size_roc, R_range):
		plt.plot(_lambda_range.magnitude, spot_size_i, label=f"{R_i}")		
	plt.title("retinal radius as function of wavelength, waist size %f %s" %(beam_waist.magnitude, beam_waist.units))
	plt.xlabel("wavelength, (%s)" %(_lambda_range.units))
	plt.ylabel("retinal radius, (%s)" %(retinal_rad.units))
	plt.legend()
	plt.show()

@app.command()
def rer():
	"""
	[rer description]
	
	Parameters
	-----------
	:	Type
		Description
	
	Returns
	----------
	Type
		Description
	
	"""
	reduced_eye_ROC()
