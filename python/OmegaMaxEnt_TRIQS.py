###################################################################################
#
# TRIQS interface for the analytic continuation program OmegaMaxEnt
#
# Copyright (C) Simons Foundation
#
# Author: Dominic Bergeron (dominic.bergeron@usherbrooke.ca)
#
# TRIQS is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# TRIQS. If not, see <http://www.gnu.org/licenses/>.
#
###################################################################################

from OmegaMaxEnt_parameters import *
from math import pi
import numpy as np
import subprocess as sp
from pytriqs.gf import *
from pytriqs.archive import HDFArchive as HA
import os
from os import path
from collections import Iterable

tol_Gi_tau = 1e-8

OME_cmd = "OmegaMaxEnt"

# parameter file generated by compute_GfReFreq() that can also be modified during the calculation
params_file = "OmegaMaxEnt_input_params.dat"
other_params_file = "OmegaMaxEnt_other_params.dat"
# params_file_copy = "OmegaMaxEnt_input_params_copy.dat"

# template file for OmegaMaxEnt_input_params.dat
# template_params_file = "OmegaMaxEnt_input_params_template.dat"

file_name = "G.dat"
error_file_name = "error_G.dat"
FT_G_file_name = "Fourier_transformed_data/Fourier_transform_G_ascii.dat"
result_file_name="OmegaMaxEnt_final_result/real_frequency_Green_function.dat"

# top_section_line = "PARAMETERS IN THIS SECTION ARE SET BY compute_GfReFreq(). DO NOT MODIFY THEM. PARAMETERS IN THE OTHER SECTIONS CAN BE MODIFIED."

# #parameters set exclusively by the function compute_scalar_GfReFreq()
# data_str = "data file:"
# boson_str = "bosonic data (yes/[no]):"
# time_str = "imaginary time data (yes/[no]):"
# temp_str = "temperature (in energy units, k_B=1):"
# grid_params_str = "output real frequency grid parameters (w_min dw w_max):"
#
# #parameters set by compute_scalar_GfReFreq() if passed to the function, otherwise they can be optionally set in file OmegaMaxEnt_input_params.dat
# err_str = "error file:"
# inter_m_str = "interactive mode ([yes]/no):"
# SW_str="spectral function width:"
# SC_str="spectral function center:"
# dw_str="real frequency step:"
# nu_grid_str="use non uniform grid in main spectral range (yes/[no]):"

# def compute_GfReFreq(G, ERR=None, grid_params=[], name="$G^R$", interactive_mode=True, save_figures_data=True,
#					 save_G=True, comp_grid_params=[], non_uniform_grid=False, inv_sym=False, mu=1, nu=1):
def compute_GfReFreq(G, **kwa):
	"""
	Compute a GfReFreq object or a BlockGf containing GfReFreq objects from a Matsubara input Green function using the
	program OmegaMaxEnt. For more details, see the OmegaMaxEnt-TRIQS interface documentation and the OmegaMaxEnt
	user guide at https://www.physique.usherbrooke.ca/MaxEnt/index.php/User_Guide.

	Most frequently used parameters:
	--------------------------------
	G:	Gf, GfImFreq, GfImTime or BlockGf object.
		Input Matsubara Green function.

	ERR:	Optional numpy array.
		Standard deviation if G is scalar or has a single element.
		ERR must have the same shape as the G.data.
		For a non-diagonal covariance, see the interface documentation or the OmegaMaxEnt User Guide.

	grid_params:	Optional list of the form [omega_min, omega_step, omega_max].
			Defines the real frequency grid of the output Green function. If empty, the output grid is set by
			OmegaMaxEnt.

	name:	Optional string. Default: "$G^R$"
		Name parameter of the returned GfReFreq object

	interactive_mode:	Optional boolean. Default: True
				Turns off the interactive mode of OmegaMaxEnt if set to False.

	save_figures_data:	Optional boolean. Default: True
				Tells OmegaMaxEnt not to save figure files if set to False.
				Set to False if G is a matrix or a BlockGf

	save_G:		Optional boolean.
			By default, the result is save in hdf5 format in file G_Re_Freq.h5.

	comp_grid_params:	Optional list of the form [freq_step] or [freq_step, SW] or [freq_step, SW, SC].
				Grid parameters used in the computation.
				freq_step is the frequency step used in the main spectral region, namely, the part of the grid where
				most of the spectral weight is located.
				SW is the width of the main spectral region (typically between 2 and 4 standard deviations).
				SC is the center of the main spectral region.
				freq_step and SW are ignored if not positive.

	non_uniform_grid:	Optional boolean. Tells OmegaMaxEnt to use a non-uniform grid in the main spectral region for
				the computation. This will accelerate the calculation if the spectrum has a peak at zero frequency with
				a width much smaller than the total width of the spectrum.

	inv_sym:	Optional boolean.
			If G is a matrix or a BlockGf, set inv_sym to True if G[i,j]=G[j,i]. This simplifies the calculation of the
			off diagonal elements.

	mu, nu:		Optional parameters involved in the calculation of off-diagonal elements of matrix-valued Green
			functions. See appendix C of the OmegaMaxEnt user guide for more details.

	The other parameters are defined in OmegaMaxEnt_parameters.py and are descibed in the OmegaMaxEnt user guide
	(https://www.physique.usherbrooke.ca/MaxEnt/index.php/User_Guide).
	"""

	if not isinstance(G, Gf) and not isinstance(G, GfImFreq) and not isinstance(G, GfImTime) and not isinstance(G, BlockGf):
		print "compute_GfReFreq(): input type " + str(G.__class__) + " not accepted. Only objects of types Gf, GfImFreq, GfImTime or BlockGf are accepted."
		return None

	print kwa

	save_G = True
	if 'save_G' in kwa:
		if isinstance(kwa['save_G'], bool):
			save_G = kwa['save_G']
		else:
			print "compute_GfReFreq() warning: save_G parameter must be boolean"
	name = "$G^R$"
	if 'name' in kwa:
		name = kwa['name']
	grid_params = []
	if 'grid_params' in kwa:
		grid_params = kwa['grid_params']

	if not isinstance(G, BlockGf):
		if len(G.target_shape)==2:
			ind = G.indices[0]
			if G.target_shape[0]==1 and G.target_shape[1]==1 and G.indices[0]==G.indices[1]:
				Gtmp=compute_scalar_GfReFreq(G[ind[0],ind[0]], **kwa)
				if not isinstance(Gtmp, GfReFreq):
					print "continuation failed"
					return None
				n_freq = len(Gtmp.mesh)
				GR=GfReFreq(indices=ind, window=(Gtmp.mesh.omega_min,Gtmp.mesh.omega_max), n_points=n_freq, name=name)
				GR[ind[0],ind[0]] = Gtmp
			elif G.indices[0]==G.indices[1]:
				GR = compute_matrix_GfReFreq(G, **kwa)
				if not isinstance(GR, GfReFreq):
					print "continuation failed"
					return None
			else:
				print "compute_GfReFreq() only treats Green functions with the same indices along both dimensions"
				return None
		elif not len(G.target_shape):
			GR = compute_scalar_GfReFreq(G, **kwa)
			if not isinstance(GR, GfReFreq):
				print "continuation failed"
				return None
		else:
			print "compute_GfReFreq() only treats matrix or scalar Green functions"
			return None
	else: #BlockGf
		kwa.update(dict(name=''))
		kwa.update(dict(save_G=False))
		list_G = []
		for bl,Gbl in G:
			Gtmp = compute_GfReFreq(Gbl, **kwa)
			if not isinstance(Gtmp, GfReFreq):
				return None
			list_G.append(Gtmp)
			if len(grid_params) != 3:
				n_freq = len(Gtmp.mesh)
				step = (Gtmp.mesh.omega_max - Gtmp.mesh.omega_min) / (n_freq - 1)
				grid_params = [Gtmp.mesh.omega_min, step, Gtmp.mesh.omega_max]
				kwa.update(dict(grid_params=grid_params))
		GR = BlockGf(name_list = list(G.indices), block_list = list_G, name=name)


	if save_G:
		with HA("G_Re_Freq.h5", 'w') as A:
			A['G'] = GR

	print "continuation done"

	return GR


#def compute_matrix_GfReFreq(G, grid_params=[], inv_sym=False, mu=1, nu=1, interactive_mode=True, save_figures_data=False, name="$G^R$", comp_grid_params=[], non_uniform_grid=False, save_G=False):
def compute_matrix_GfReFreq(G, **kwa):
	"""
	Used by compute_GfReFreq() to compute a matrix-valued GfReFreq from a matrix-valued Matsubara function G.
	"""

	if not isinstance(G,Gf) and not isinstance(G,GfImFreq) and not isinstance(G,GfImTime):
		print "compute_matrix_GfReFreq(): input type " + str(G.__class__) + " not accepted. Only objects of types Gf, GfImFreq or GfImTime are accepted."
		return None

	if len(G.target_shape)!=2:
		print "compute_matrix_GfReFreq(): the Green function must be a matrix"
		return None

	if G.target_shape[0]<2 or G.target_shape[0]<2 or G.target_shape[0]!=G.target_shape[1]:
		print "compute_matrix_GfReFreq(): the Green function must be a square matrix of dimension at least 2"
		return None

	if G.indices.data[0]!=G.indices.data[1]:
		print "compute_matrix_GfReFreq(): row and column indices must the same"
		return None

	grid_params = []
	if 'grid_params' in kwa:
		if len(kwa['grid_params'])==3:
			grid_params = kwa['grid_params']
		else:
			print "compute_matrix_GfReFreq() warning: 'grid_params' parameter must contain three elements"
	inv_sym = False
	if 'inv_sym' in kwa:
		if isinstance(kwa['inv_sym'],bool):
			inv_sym=kwa['inv_sym']
		else:
			print "compute_matrix_GfReFreq() warning: 'inv_sym' parameter must be boolean"
	mu = 1
	if 'mu' in kwa:
		mu=kwa['mu']
	nu = 1
	if 'nu' in kwa:
		nu=kwa['nu']
	name = "$G^R$"
	if 'name' in kwa:
		name=kwa['name']
		kwa['name']=''
	save_G = False
	if 'save_G' in kwa:
		save_G=kwa['save_G']

	#remove all parameters that do not make sense for matrix-valued Green's function
	if 'ERR' in kwa:
		print "compute_matrix_GfReFreq() warning: 'ERR' parameter is applicable only to scalar Green's functions. Parameter discarded."
		del kwa['ERR']
	if 'norm' in kwa:
		print "compute_matrix_GfReFreq() warning: 'norm' parameter is applicable only to scalar Green's functions. Parameter discarded."
		del kwa['norm']
	if 'M1' in kwa:
		print "compute_matrix_GfReFreq() warning: 'M1' parameter is applicable only to scalar Green's functions. Parameter discarded."
		del kwa['M1']
	if 'M2' in kwa:
		print "compute_matrix_GfReFreq() warning: 'M2' parameter is applicable only to scalar Green's functions. Parameter discarded."
		del kwa['M2']
	if 'M3' in kwa:
		print "compute_matrix_GfReFreq() warning: 'M3' parameter is applicable only to scalar Green's functions. Parameter discarded."
		del kwa['M3']
	if 'error_file' in kwa:
		print "compute_matrix_GfReFreq() warning: 'error_file' parameter is applicable only to scalar Green's functions. Parameter discarded."
		del kwa['error_file']
	if 'cov_re_re' in kwa:
		print "compute_matrix_GfReFreq() warning: 'cov_re_re' parameter is applicable only to scalar Green's functions. Parameter discarded."
		del kwa['cov_re_re']
	if 'cov_im_im' in kwa:
		print "compute_matrix_GfReFreq() warning: 'cov_im_im' parameter is applicable only to scalar Green's functions. Parameter discarded."
		del kwa['cov_im_im']
	if 'cov_re_im' in kwa:
		print "compute_matrix_GfReFreq() warning: 'cov_re_im' parameter is applicable only to scalar Green's functions. Parameter discarded."
		del kwa['cov_re_im']

	# if not path.exists(params_file):
	# 	create_params_file(False)

	ind =G.indices[0]

	Gtmp=compute_scalar_GfReFreq(G[ind[0],ind[0]], **kwa)

	if not isinstance(Gtmp, GfReFreq):
		return None

	n_freq = len(Gtmp.mesh)

	if len(grid_params)!=3:
		step=(Gtmp.mesh.omega_max-Gtmp.mesh.omega_min)/(n_freq-1)
		grid_params = [Gtmp.mesh.omega_min, step, Gtmp.mesh.omega_max]
		kwa.update(dict(grid_params=grid_params))

	GM=GfReFreq(indices=G.indices, window = (Gtmp.mesh.omega_min, Gtmp.mesh.omega_max), n_points = n_freq, name=name)

	GM[ind[0], ind[0]] = Gtmp

	print "G[" + ind[0] + "," + ind[0] + "] computed"
	if save_G:
		with HA("G_Re_Freq_" + ind[0] + "_" + ind[0] + ".h5", 'w') as A:
			A['G'] = GM[ind[0], ind[0]]

	for l in ind[1:]:
		Gtmp=compute_scalar_GfReFreq(G[l, l], **kwa)
		if not isinstance(Gtmp, GfReFreq):
			return None
		GM[l, l]=Gtmp
		print "G[" + l + "," + l + "] computed"
		if save_G:
			with HA("G_Re_Freq_" + l + "_" + l + ".h5", 'w') as A:
				A['G'] = GM[l, l]

	if not inv_sym:
		for l in range(len(ind)):
			for m in range(l+1,len(ind)):
				GO=G[ind[l],ind[l]]+mu*G[ind[l],ind[m]]+mu*G[ind[m],ind[l]]+mu*mu*G[ind[m],ind[m]]
				Gtmp = compute_scalar_GfReFreq(GO, **kwa)
				if not isinstance(Gtmp, GfReFreq):
					return None
				GOR = Gtmp
				GP=G[ind[l],ind[l]]-1j*nu*G[ind[l],ind[m]]+1j*nu*G[ind[m],ind[l]]+nu*nu*G[ind[m],ind[m]]
				Gtmp = compute_scalar_GfReFreq(GP, **kwa)
				if not isinstance(Gtmp, GfReFreq):
					return None
				GPR = Gtmp
				R=GOR-GM[ind[l],ind[l]]-mu*mu*GM[ind[m],ind[m]]
				S=GPR-GM[ind[l],ind[l]]-nu*nu*GM[ind[m], ind[m]]
				GM[ind[l],ind[m]]=(R/mu+1j*S/nu)/2
				GM[ind[m],ind[l]]=(R/mu-1j*S/nu)/2
				print "G[" + ind[l] + "," + ind[m] + "] computed"
				print "G[" + ind[m] + "," + ind[l] + "] computed"
				if save_G:
					with HA("G_Re_Freq_" + ind[l] + "_" + ind[m] + ".h5", 'w') as A:
						A['G'] = GM[ind[l], ind[m]]
					with HA("G_Re_Freq_" + ind[m] + "_" + ind[l] + ".h5", 'w') as A:
						A['G'] = GM[ind[m], ind[l]]
	else:
		for l in range(len(ind)):
			for m in range(l+1,len(ind)):
				GO=G[ind[l],ind[l]]+2*mu*G[ind[l],ind[m]]+mu*mu*G[ind[m],ind[m]]
				Gtmp = compute_scalar_GfReFreq(GO, **kwa)
				if not isinstance(Gtmp, GfReFreq):
					return None
				GOR = Gtmp
				GM[ind[l],ind[m]]=(GOR-GM[ind[l],ind[l]]-mu*mu*GM[ind[m],ind[m]])/(2*mu)
				GM[ind[m],ind[l]]=GM[ind[l], ind[m]]
				print "G[" + ind[l] + "," + ind[m] + "] computed"
				if save_G:
					with HA("G_Re_Freq_"+ind[l]+"_"+ind[m]+".h5", 'w') as A:
						A['G'] = GM[ind[l],ind[m]]

	#print "matrix continuation done"

	return GM


#def compute_scalar_GfReFreq(G, ERR=None, grid_params=[], name="$G^R$", interactive_mode=True, save_figures_data=True, comp_grid_params=[], non_uniform_grid=False):
def compute_scalar_GfReFreq(G, **kwa):
	"""
	Used by compute_GfReFreq() and compute_matrix_GfReFreq() to compute a scalar GfReFreq object from a scalar Matsubara function G.
	"""

	if not isinstance(G,Gf) and not isinstance(G,GfImFreq) and not isinstance(G,GfImTime):
		print "compute_scalar_GfReFreq(): input type " + str(G.__class__) + " not accepted. Only objects of types Gf, GfImFreq or GfImTime are accepted."
		return None

	ERR = None
	if 'ERR' in kwa:
		ERR = kwa['ERR']
	name = "$G^R$"
	if 'name' in kwa:
		name = kwa['name']
	comp_grid_params = []
	if 'comp_grid_params' in kwa:
		comp_grid_params = kwa['comp_grid_params']
		del kwa['comp_grid_params']
		comp_grid_params_dict={}
		if len(comp_grid_params) > 0 and comp_grid_params[0]>0:
			comp_grid_params_dict.update(dict(freq_step=comp_grid_params[0]))
		if len(comp_grid_params) > 1 and comp_grid_params[1]>0:
			comp_grid_params_dict.update(dict(SW=comp_grid_params[1]))
		if len(comp_grid_params) > 2:
			comp_grid_params_dict.update(dict(SC=comp_grid_params[2]))
		if len(comp_grid_params_dict) > 0:
			kwa.update(comp_grid_params_dict)
	save_figures_data = True
	if 'save_figures_data' in kwa:
		save_figures_data = kwa['save_figures_data']

	interactive_mode = True
	if 'interactive_mode' in kwa:
		interactive_mode = kwa['interactive_mode']

	if len(G.target_shape):
		print "compute_scalar_GfReFreq(): the Green function must be scalar"

	cmd = [OME_cmd]

	if not save_figures_data:
		cmd = cmd+ ["-np"]

	if not interactive_mode:
		cmd = cmd + ["-ni"]

	# if not path.exists(params_file):
	# 	create_params_file(False)
	
	im_t=isinstance(G.mesh,MeshImTime)
#	if im_t:
#		print "imaginary time Green function provided to compute_scalar_GfReFreq()"

	bosonic=(G.mesh.statistic=="Boson")
#	if bosonic:
#		print "bosonic Green function provided to compute_scalar_GfReFreq()"

	error_provided = isinstance(ERR, np.ndarray)

	if not path.exists(other_params_file):
		pf = open(other_params_file, "w")
		for key, val in OmegaMaxEnt_other_params.iteritems():
			v=Other_params_default_values[key]
			if key in kwa:
				v=kwa[key]
			str_tmp=val + " " + str(v)+'\n'
			pf.write(str_tmp)
		pf.close()

	# pf=open(params_file,"r")
	# str_tmp=pf.read()
	# cpf=open(params_file_copy,"w")
	# cpf.write(str_tmp)
	# pf.close()
	# cpf.close()
	#
	# cpf=open(params_file_copy,"r")
	# pf = open(params_file, "w")
	# for str_tmp in cpf:
	# 	if str_tmp[0:len(data_str)]==data_str:
	# 		str_tmp = data_str + file_name + '\n'
	# 	elif str_tmp[0:len(boson_str)]==boson_str:
	# 		if bosonic:
	# 			str_tmp = boson_str + "yes" + '\n'
	# 		else:
	# 			str_tmp = boson_str + '\n'
	# 	elif str_tmp[0:len(time_str)]==time_str:
	# 		if im_t:
	# 			str_tmp = time_str + "yes" + '\n'
	# 		else:
	# 			str_tmp = time_str + '\n'
	# 	elif  str_tmp[0:len(grid_params_str)]==grid_params_str:
	# 		if len(grid_params)==3:
	# 			str_tmp=grid_params_str+str(grid_params[0])+" "+str(grid_params[1])+" "+str(grid_params[2])+'\n'
	# 		else:
	# 			str_tmp = grid_params_str+'\n'
	# 	elif error_provided and str_tmp[0:len(err_str)]==err_str:
	# 		str_tmp=err_str+error_file_name+'\n'
	# 	elif len(comp_grid_params)>0 and str_tmp[0:len(dw_str)]==dw_str and comp_grid_params[0]>0:
	# 		str_tmp =dw_str+str(comp_grid_params[0])+'\n'
	# 	elif len(comp_grid_params)>1 and str_tmp[0:len(SW_str)]==SW_str and comp_grid_params[1]>0:
	# 		str_tmp =SW_str+str(comp_grid_params[1])+'\n'
	# 	elif len(comp_grid_params)>2 and str_tmp[0:len(SC_str)]==SC_str:
	# 		str_tmp =SC_str+str(comp_grid_params[2])+'\n'
	# 	elif non_uniform_grid and str_tmp[0:len(nu_grid_str)]==nu_grid_str:
	# 		str_tmp = nu_grid_str + "yes" + '\n'
	# 	elif not interactive_mode and str_tmp[0:len(inter_m_str)]==inter_m_str:
	# 		str_tmp=inter_m_str+"no"+'\n'
	# 	pf.write(str_tmp)
	# cpf.close()
	# pf.close()
	#
	# os.remove(params_file_copy)

	Gr = G.data.real
	Gi = G.data.imag

	if not im_t:
		iwn=np.array([[w.value for w in G.mesh]])
		wn=iwn.imag.T
		Gr=np.array([Gr])
		Gr=Gr.transpose()
		Gi=np.array([Gi])
		Gi=Gi.transpose()
		data_array=np.concatenate((wn,Gr,Gi),axis=1)
	else:
		if abs(Gi).max()/abs(Gr).max()>tol_Gi_tau:
			print "compute_scalar_GfReFreq(): warning, only the real part of imaginary time data is used"
		tau=np.array([[t.value for t in G.mesh]])
		tau=tau.T
		Gr = np.array([Gr])
		Gr = Gr.transpose()
		data_array=np.concatenate((tau,Gr),axis=1)

	n_points=data_array.shape[0]

	np.savetxt(file_name,data_array)

	if error_provided:
		dim_ERR=np.array(ERR.shape)
		if dim_ERR.max()!=n_points:
			print "compute_scalar_GfReFreq(): provided error array does not have the same size as the data."
			return None
		if ERR.shape[1]>ERR.shape[0]:
			ERR=ERR.transpose()
		if not im_t:
			error_array=np.concatenate((wn,ERR),axis=1)
		else:
			error_array=np.concatenate((tau,ERR),axis=1)
		np.savetxt(error_file_name,error_array)

	pf = open(params_file, "w")
	str_tmp = data_str + file_name + '\n'
	pf.write(str_tmp)
	if bosonic:
		str_tmp = boson_str + "yes" + '\n'
		pf.write(str_tmp)
	if im_t:
		str_tmp = time_str + "yes" + '\n'
		pf.write(str_tmp)
	if error_provided:
		str_tmp = err_str + error_file_name + '\n'
		pf.write(str_tmp)

	for key, val in kwa.iteritems():
		if key in OmegaMaxEnt_input_params:
			if isinstance(val, bool):
				val_str = "no"
				if val:
					val_str = "yes"
			elif isinstance(val, Iterable):
				val_str = ""
				for elem in val:
					val_str = val_str + " " + str(elem)
			else:
				val_str = str(val)
			str_tmp = OmegaMaxEnt_input_params[key] + val_str + '\n'
			pf.write(str_tmp)
	pf.close()

	# call OmegaMaxEnt
	rval=sp.call(cmd)

	if rval:
		return None

	if im_t:
		save_Fourier_transform_G_hdf5()
	
	#retrieve the real frequency Green function
	result_file=open(result_file_name,"r")
	G_Re_w_data=np.loadtxt(result_file)
	result_file.close()

	if os.path.exists(file_name):
		os.remove(file_name)
	if error_provided:
		if os.path.exists(error_file_name):
			os.remove(error_file_name)

	GR_omega=GfReFreq(target_shape=(),window = (G_Re_w_data[0,0], G_Re_w_data[-1,0]), n_points = G_Re_w_data.shape[0], name = name)

	GR_omega.data.real = G_Re_w_data[:, 1]
	GR_omega.data.imag = G_Re_w_data[:, 2]

	return GR_omega

def create_params_file(overwrite=True):
	"""
	create the parameter files OmegaMaxEnt_input_params.dat and OmegaMaxEnt_other_params.dat used by OmegaMaxEnt.

	parameters:
	----------
	overwrite:	optional boolean.
			If True, existing parameter files will be destroyed.
	"""
	cmd = [OME_cmd, "-ni"]

	if not path.exists(params_file) or overwrite:
		if not path.exists(template_params_file):
			if path.exists(params_file):
				os.remove(params_file)
			sp.call(cmd)
		tf=open(template_params_file,"r")
		uf=open(params_file,"w")
		uf.write(top_section_line+"\n")
		uf.write(data_str+"\n")
		uf.write(boson_str + "\n")
		uf.write(time_str + "\n")
		uf.write(temp_str + "\n")
		uf.write(grid_params_str+"\n")
		uf.write("\n")
		for line in tf:
			str_tmp=line[0:-1]
			if str_tmp!=data_str and str_tmp!=boson_str and str_tmp!=time_str and str_tmp!=temp_str and str_tmp!=grid_params_str:
				uf.write(line)
		tf.close()
		uf.close()
		os.remove(template_params_file)

def display_figures():
	"""
	Display the figures showing the result after compute_GfReFreq() was called with save_figures_data=True (default). The figures are the ones displayed
	if interactive_mode=True (default).
	"""
	figs_ind=np.int_(np.loadtxt("figs_ind.dat"))
	if not isinstance(figs_ind,np.ndarray):
		figs_ind=np.array([figs_ind])
	if figs_ind.shape[0]>0:
		for i in figs_ind:
			sname="OmegaMaxEnt_figs_"+str(i)+".py"
			fig_file=open(sname,"r")
			fig_cmd=fig_file.read()
			exec(fig_cmd)

def save_Fourier_transform_G_hdf5():
	"""
	Called by compute_scalar_GfReFreq() to save the Fourier transform of a scalar GfImTime object as a GfImFreq in hdf5 format
	"""
	data_file = open(FT_G_file_name, "r")
	G_data = np.loadtxt(data_file)
	data_file.close()

	statistic="Fermion"
	bosonic=False
	if not G_data[0,0]:
		statistic="Boson"
		bosonic=True

	ind0 = 0
	N = G_data.shape[0] - 1
	if bosonic:
		N = G_data.shape[0]
		ind0 = 1

	if bosonic:
		beta = 2*pi/G_data[1,0]
	else:
		beta = pi/G_data[0, 0]

	Gwn = GfImFreq(target_shape=(), n_points=N, beta=beta, statistic=statistic)

	Gr_p = G_data[:N, 1]
	Gi_p = G_data[:N, 2]

	Gi_n = -Gi_p[ind0:]
	Gi_n = np.flipud(Gi_n)
	Gr_n = Gr_p[ind0:]
	Gr_n = np.flipud(Gr_n)

	Gr = np.concatenate((Gr_n, Gr_p), 0)
	Gi = np.concatenate((Gi_n, Gi_p), 0)

	Gwn.data.real = Gr
	Gwn.data.imag = Gi

	A = HA("G_im_freq.h5", "w")
	A['G'] = Gwn