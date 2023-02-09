from retinal_spot_size import *
from .utils import *
from .eye_ABCD import *
import typer
import pathlib
import json
import matplotlib.pyplot as plt


app = typer.Typer()

@app.command()
def rss(config_file : pathlib.Path):
	if not config_file.exists():
		print(f"Config file does not exits `{config_file}`")
		return 1


@app.command()
def rss_figure_7():
	#if not config_file.exists():
	if False:
		print(f"Config file does not exits `{config_file}`")
		return 1
	# Handle in arguments later
	N = 100
	_lambda_lower_bound = Q_(400, 'nm')
	_lambda_upper_bound = Q_(1400, 'nm')
	_lambda_range = np.linspace(_lambda_lower_bound, _lambda_upper_bound, 100)

	z0 = Q_(0, 'mm')
	beam_waist = Q_(3, 'mm')
	radius_list = []
	for _lambda in _lambda_range:
		_ang_div = ang_div(_lambda, beam_waist)
		q_0 = q_hat(z0, z_R(_lambda, _ang_div))
		q_1 = complex_matmul(q_0, eye_equiv_ABCD(_lambda))
		retinal_rad = current_beam_rad(_lambda, q_1)
		retinal_rad.ito('um')
		print(retinal_rad)
		radius_list.append(retinal_rad.magnitude)
	radius_list = np.array(radius_list)
	plt.plot(_lambda_range.magnitude, radius_list)		
	plt.title("retinal radius as function of wavelength, waist size %f %s" %(beam_waist.magnitude, beam_waist.units))
	plt.xlabel("wavelength, (%s)" %(_lambda_range.units))
	plt.ylabel("retinal radius, (%s)" %(retinal_rad.units))
	plt.show()
	

@app.command()
def spot_size_z_range():
	lambda_range = np.array([1315, 1318, 1319, 1330, 1338, 1356]) * ureg.nm
	ang_div_range = np.linspace(4, 5, 1) * ureg.mrad
	_lambda = lambda_range[0]
	z0_range = np.linspace(-10, 10, 50) * ureg.mm
	for div in ang_div_range:
		radius_list = np.zeros(len(z0_range))
		for i, z0 in enumerate(z0_range):
			q_0 = q_hat(z0, z_R(_lambda, div))
			q_1 = complex_matmul(q_0, eye_equiv_ABCD(_lambda))
			beam_waist = current_beam_rad(_lambda, q_0) 
			retinal_rad = current_beam_rad(_lambda, q_1 / N_VITREOUS)
			radius_list[i] = retinal_rad.magnitude
		plt.plot(z0_range, radius_list)		
		plt.title("retinal radius as function of waist location, waist size %f" %beam_waist.magnitude)
		plt.show()

@app.command()
def ssz():
	spot_size_z_range()	
