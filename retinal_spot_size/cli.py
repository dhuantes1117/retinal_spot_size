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
def rss(config_file : pathlib.Path):
	if not config_file.exists():
		print(f"Config file does not exits `{config_file}`")
		return 1


@app.command()
def rss_figure_8():
	N = 100
	_lambda_lower_bound = Q_(8000, 'nm')
	_lambda_upper_bound = Q_(10400, 'nm')
	_lambda_lower_bound = Q_(400, 'nm')
	_lambda_upper_bound = Q_(1400, 'nm')
	_lambda_range = np.linspace(_lambda_lower_bound, _lambda_upper_bound, 101)
	z0 = Q_(0, 'm')
	beam_waist = Q_(4.24 / 2, 'mm')

	print(eye_equiv_ABCD(Q_(700, 'nm')))
	radius_list = []
	roc_list = []
	z_list = []
	for _lambda in _lambda_range:
		_ang_div = ang_div(_lambda, beam_waist)
		q_0 = q_hat(z0, z_R(_lambda, _ang_div))
		q_1 = propagate_q_hat(q_0, reduced_eye(_lambda, 6.1 * ureg.mm))
		#print(f"{_lambda} {q_1}")
		retinal_roc = current_beam_roc(N_VITREOUS, q_1)
		retinal_rad = current_beam_rad(_lambda, q_1)
		retinal_rad.ito("um")
		radius_list.append(retinal_rad.magnitude)
		roc_list.append(retinal_roc.magnitude)
	radius_list = np.array(radius_list)
	# Code is generating a reflected, shifted version, these numbers
	# qualitatively fit what is happening, so the question is: What is
	# causing our function to have this issue?
	#plt.plot(-0.35 * (_lambda_range.magnitude - 11600), radius_list * 0.1)		
	plt.plot(_lambda_range.magnitude, 2 * radius_list)		
	# plt.plot(_lambda_range.magnitude, roc_list)		
	plt.title("retinal radius as function of wavelength, waist size %f %s" %(beam_waist.magnitude, beam_waist.units))
	plt.xlabel("wavelength, (%s)" %(_lambda_range.units))
	plt.ylabel("retinal radius, (%s)" %(retinal_rad.units))
	plt.show()
	print(interp(_lambda_range.magnitude, 2 * radius_list)(1300 * ureg.nm))


@app.command()
def rf8():
	import pudb
	#pu.db
	rss_figure_8()	
	

@app.command()
def spot_size_z_range():
	_lambda = Q_(700, 'nm')
	beam_waist = Q_(100, 'um')
	z0_range = np.linspace(-100, 100, 100) * ureg.mm
	div = ang_div(_lambda, beam_waist)
	radius_list = np.zeros(len(z0_range))
	for i, z0 in enumerate(z0_range):
		q_0 = q_hat(z0, z_R(_lambda, div))
		q_1 = complex_matmul(q_0, eye_equiv_ABCD(_lambda))
		beam_waist = current_beam_rad(_lambda, q_0) 
		retinal_rad = current_beam_rad(_lambda, q_1)
		retinal_rad.ito('um')
		radius_list[i] = retinal_rad.magnitude
	plt.plot(z0_range, radius_list)		
	plt.title("retinal radius as function of waist location, waist size %f (%s)" %(beam_waist.magnitude, beam_waist.units))
	plt.show()


@app.command()
def ssz():
	spot_size_z_range()	

@app.command()
def spot_size_w_range():
	_lambda = Q_(700, 'nm')
	z0 = 0 * ureg.mm 
	w0_range = np.linspace(10, 1000, 100) * ureg.um
	radius_list = np.zeros(len(w0_range))
	for i, beam_waist in enumerate(w0_range):
		div = ang_div(_lambda, beam_waist)
		q_0 = q_hat(z0, z_R(_lambda, div))
		q_1 = complex_matmul(q_0, eye_equiv_ABCD(_lambda))
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
	spot_size_w_range()	

@app.command()
def display_sizes():
	_lambda = 700 * ureg.nm
	M_lens = eye_lens_ABCD(_lambda)
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
	print(waist.magnitude)
	radius_list = [waist.magnitude]
	for interface in optical_system[::-1]:
		q_1 = propagate_q_hat(q_0, interface)
		retinal_rad = current_beam_rad(_lambda, q_1)
		retinal_rad.ito('um')
		radius_list.append(retinal_rad.magnitude)
	print(radius_list)
	
	radius_list = np.array(radius_list)
	plt.plot(np.linspace(0, 10, len(radius_list)), radius_list)		
	plt.title("retinal radius as function of wavelength, waist size %f %s" %(beam_waist.magnitude, beam_waist.units))
	plt.xlabel("interface")
	plt.ylabel("retinal radius, (%s)" %(retinal_rad.units))
	plt.show()

@app.command()
def ds():
	display_sizes()

@app.command()
def propagation_progress():
	_lambda = Q_(700, 'nm')
	beam_waist = Q_(100, 'um')
	Delta_z = Q_(1, 'mm')
	z0 = Q_(0, 'mm')
	L = Q_(2, 'mm')
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
	distance = np.linspace(0, L.magnitude, N)
	div = ang_div(_lambda, beam_waist)
	q_0 = q_hat(z0, z_R(_lambda, div))
	n_water = 1.33
	ABCD_forward 	= free_space(Delta_z, 1)
	ABCD_forward_water = free_space(Delta_z, n_water)
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
	#q_i = propagate_q_hat(q_i, ABCD_thin_lens)
	q_i = propagate_q_hat(q_i, ABCD_air_water)
	beam_radius_i = current_beam_rad(_lambda, q_i)
	radii.append(beam_radius_i.magnitude)
	for i in range(N // 2):
		q_i = propagate_q_hat(q_i, ABCD_forward_water)
		beam_radius_i = current_beam_rad(_lambda, q_i)
		radii.append(beam_radius_i.magnitude)
	
	radii = np.array(radii)
	plt.plot(distance, radii)
	plt.plot(distance, -(radii))
	plt.title("")
	plt.xlabel("distacne (%s)" %L.units)
	plt.ylabel("radius (%s)" %beam_radius_i.units)
	plt.show()
	print("initial complex beam parameter: %s" %(q_0))
	print("  final complex beam parameter: %s" %(q_i))

@app.command()
def pp():
	propagation_progress()

@app.command()
def sellmeier_figure_7():
	lambda_ = np.linspace(400, 1400, 1000) * ureg.nm
	n	= sellmeier_vitreous(lambda_)
	plt.plot(lambda_, n)
	plt.show()

@app.command()
def rf7():
	sellmeier_figure_7()

@app.command()
def reduced_eye_ROC():
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

		print(eye_equiv_ABCD(Q_(700, 'nm')))
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
	reduced_eye_ROC()
