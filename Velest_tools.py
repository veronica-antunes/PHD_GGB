#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Veronica Antunes 07/11/2018 modified 03/12/2019
# AIM: Set of tools that can be used with velest output files

import sys, os, re
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from shutil import copyfile

def best_vpvs(infile='invers.out', outfile_vel='input2.mod', vm_path = 'all_vm_models', xlim=[0,10]) :
	#'AIM: reads the best RMS value from invers.out file and writes the velocity model file
	#'input.mod in a way that can be used in the next VELEST inversion'
	#invers='invers.out'
	#outfile_name='input2.mod'
	rms_list=[]
	invers_file = open(infile, 'r')

	for line in invers_file.readlines():
		if re.search('RMS RESIDUAL', line, re.I):
			rms_list.append(float(line.split(' ')[-1][0:-1]))

	min_rms=min(rms_list)

	invers_file_read = open(infile, 'r').read()

	datvar=re.split(str(min_rms), invers_file_read)[0][-62:-53]
	sqrd=re.split(str(min_rms), invers_file_read)[0][-29:-20]

	try:
		VM_P=re.split(str(min_rms), invers_file_read)[1].split(
			' Velocity model   1')[1].split(
			'Calculation of average velocity starts at layer')[0].split('\n')

		VM_S=re.split(str(min_rms), invers_file_read)[1].split(
			' Velocity model   2')[1].split(
			'Calculation of average velocity starts at layer')[0].split('\n')
	except:
		print ('problem: check script')
	# CHECK IF THE PROBLEM CAN BE SOLVED WITH THE BELOW OPTIONS. NOT TESTED!
	#	VM=re.split(str(min_rms), invers_file_read)[1].split(
	#		'nlay   top ..... bottom     velocity   NHYP NREF %len  
	#		NHIT xy-km  z-km  RFLX')[1].split(
	#		'Total nr of events was')[0].split('\n')
	#	vmlist=[]
	#	vmlist.append('velocity model')
	#	for i in range(len(VM)-2) :
	#		if i <= 1:
	#			vmlist.append(VM[i])
	#		if i > 1 :
	#			vel_value=float(VM[i][8:13])
	#			if vel_value != 0 :
	#				vmlist.append('   '+VM[i][29:34]+'  0.000   
	#					'+VM[i][8:13]+'       ')
	#	vmlist.append(' ')
	#	vmlist.append(' ')
	#	VM=vmlist

	g=open( 'input.mod' , 'r').readlines()
	f = open( outfile_vel , 'w')
	f.write(' Model: input.mod\n')
	nlayer=str(len(VM_P)-3)
	f.write(' '+nlayer+'        vel,depth,vdamp,phase(f5.2,5x,f7.2,2x,f7.3,'+
		'3x,a1)                    \n')
	VP=[]
	VS=[]
	Z=[]

	for l in range(len(VM_P)-2) :
		if l == 1 :
			vel=round(float(VM_P[l][3:8]),2)
			dep=round(float(VM_P[l][18:23]),2)
			damp=g[l+1][19:26]

			f.write(' '+'{:04.2f}'.format(vel)+'       '+'{:05.2f}'.format(dep)+'  '+
				damp+'           P-VELOCITY MODEL                           \n')
			Z.append(dep)
			VP.append(vel)
			VP.append(vel)

		if l > 1 :
			vel=round(float(VM_P[l][3:8]),2)
			dep=round(float(VM_P[l][18:23]),2)
			damp=g[l+1][19:26]
			f.write(' '+'{:04.2f}'.format(vel)+'       '+'{:05.2f}'.format(dep)+'  '+
				damp+'                                                      \n')
			Z.append(dep)
			Z.append(dep)
			VP.append(vel)
			VP.append(vel)

	Z.append(Z[-1]+20)
	f.write(' '+str(nlayer)+'\n')

	for l in range(len(VM_S)-2) :
		if l == 1 :
			vel=round(float(VM_S[l][3:8]),2)
			dep=round(float(VM_S[l][18:23]),2)
			damp=g[l+2+int(nlayer)][19:26]
			f.write(' '+'{:04.2f}'.format(vel)+'       '+'{:05.2f}'.format(dep)+'  '+
				damp+'           S-VELOCITY MODEL                           \n')
			VS.append(vel)
			VS.append(vel)

		if l > 1 :
			vel=round(float(VM_S[l][3:8]),2)
			dep=round(float(VM_S[l][18:23]),2)
			damp=g[l+2+int(nlayer)][19:26]
			f.write(' '+'{:04.2f}'.format(vel)+'       '+'{:05.2f}'.format(dep)+'  '+
				damp+'                                                      \n')
			VS.append(vel)
			VS.append(vel)

	f.close()

	#vm_path = 'all_vm_models'
	if not os.path.exists(vm_path):
	    os.makedirs(vm_path)

	copyfile('input2.mod', vm_path+'/'+str('{:1.6f}'.format(min_rms))[2::]+'_'+outfile_vel)

	fig=plt.figure()
	fig.suptitle('DATVAR: ' + datvar + ' sqrd: ' + sqrd +' RMS: '+str(min_rms) , fontsize=14)
	plt.plot(VP,Z, 'b', label='Vp')
	plt.plot(VS,Z, 'r', label='Vs')
	plt.legend()
	plt.xlim(xlim[0],xlim[1])
	plt.ylim(Z[0],Z[-1])
	plt.gca().invert_yaxis()
	return plt.savefig('input2.png', bbox_inches='tight')
	return plt.savefig(vm_path+'/'+str(min_rms)[2::]+'_'+'input.png', 
		bbox_inches='tight')

def km2deg(deg) :
	km=deg*111.133
	return km

def get_time(t) :
	from obspy import UTCDateTime
	time=UTCDateTime(year=int(t[0:4]),month=int(t[5:7]), day=int(t[7:9]), 
		hour=int(t[10:12]), minute=int(t[12:14])) + float(t[15:19])
	return time

def plot_diff(array, name, units):
	x=np.arange(1,len(array)+1,1)
	colors = np.where(array<0,'b', np.where(array>0,'r','grey'))
	mean=round(np.mean(array), 2)
	plt.scatter(x, array, c=colors, s=10, label='mean: '+str(mean))
	plt.axhline(y=0, color='grey', linewidth=0.3)
	plt.vlines(x, ymin=0, ymax=array, color=colors, linewidth=0.5)
	plt.legend(fontsize=18)
	plt.xticks(fontsize=18)
	plt.yticks(fontsize=18)
	plt.title(name, fontsize=20)
	plt.ylabel('difference in '+name+' ['+units+']', fontsize=20)
	plt.xlabel('earthquake #', fontsize=20)
	plt.savefig(name+'_diff.png', bbox_inches='tight')	
	plt.close()

def read_veloutdif(velfile='velout.dif') :
	a=open(velfile,'r').read()
	eve_n=len(re.split('depth:', a))
	depth_list=[]
	time_list=[]
	lat_list=[]
	lon_list=[]

	for i in range(eve_n):
		if i > 0 :
			depth_list.append(float(re.split('depth:', a)[i][0:7]))
			lat_list.append(km2deg(float(re.split('latitude:', a)[i][0:9])))
			lon_list.append(km2deg(float(re.split('longitude:', a)[i][0:10])))
			t1=get_time((re.split('\n\n', a))[i].split('\n')[0][1:21])
			t2=get_time((re.split('\n\n', a))[i].split('\n')[1][1:21])
			time_list.append(float(t2-t1))

	depth=np.array(depth_list)
	time=np.array(time_list)
	lat=np.array(lat_list)
	lon=np.array(lon_list)

	return time, lat, lon, depth

def plot_veloutdif(velfile='velout.dif') :
	#AIM: plots diff depth, time, lat and lon
	time, lat, lon, depth = read_veloutdif(velfile=velfile)
		
	plot_diff(depth, 'Depth','km')
	plot_diff(time, 'Time','s')
	plot_diff(lat, 'Latitude','km')
	plot_diff(lon, 'Longitude','km')

def quality_solutions(invers_path) :
	#AIM: checks velest solutions: plots rms vs nruns
	#scenarios_path=glob('P_*')
	from glob import glob

	#for s in scenarios_path:
	runs=sorted(glob(invers_path), key=os.path.getmtime)
	min_rms=[]

	for r in runs:
		try:
			#invers=(r+'/'+'invers.out')
			invers_file = open(r, 'r')
			rms_list=[]
			for line in invers_file.readlines():
				if re.search('RMS RESIDUAL', line, re.I):
					rms_list.append(float(line.split(' ')[-1][0:-1]))
		except:
			print ('no invers file for .../' +r)
		min_rms.append(min(rms_list))

	fig=plt.figure()
	#fig.suptitle( '' , fontsize=14)
	plt.plot(min_rms, c='steelblue', linewidth=1.8)
	y=min(min_rms)
	x=min_rms.index(min(min_rms))
	plt.plot(x, y, marker='o', markersize=5, color="red", 
		label='best solution \n'+str(min(min_rms)))
	plt.legend(loc='upper right', numpoints=1, fontsize=18)
	plt.axhline(y=y, color='r', linestyle='--')
	plt.xlabel('Number of runs', fontsize=20)
	plt.ylabel('rms misfit', fontsize=20)
	plt.ylim(0,max(min_rms)+0.2)
	plt.xticks(fontsize=18)
	plt.yticks(fontsize=18)
	plt.savefig('plot_run_rms.png', bbox_inches='tight')
	plt.close()

def get_bestvm(models_path) :
	#AIM: gets best vm for sta_corr
	#models_path='P_*/all_vm_models/*'
	vm_list=glob(models_path)
	best_vm=sorted([x.split('/')[-1] for x in vm_list])[0]
	copyfile(glob(models_path+best_vm)[0], 'input_FINAL.mod')

def plot_P_vms(initial_models, min1D_model, final_models=None, all_models='*.mod', 
	xlim=[2,13], ylim=[-5,60], outfile='all_vms.png') :

	from glob import glob
	all_models=glob(all_models)
	vm_models=all_models
	if final_models:
		final_models=glob(final_models)
		vm_models+=final_models
	initial_models=glob(initial_models)
	vm_models+=initial_models
	vm_models+=glob(min1D_model)

	vm_models = sorted(list(dict.fromkeys(vm_models)))

	VP={}
	Z={}

	color_list=['darkred', 'mediumblue', 'green', 'purple', 'darkblue', 
		'violet', 'yellow', 'brown', 'pink']

	fig=plt.figure()
	fig.suptitle('Velocity models' , fontsize=22)
	plt.figure(figsize=(6,9))
	plt.xlim(xlim[0],xlim[1])
	plt.ylim(ylim[0],ylim[1])
	plt.gca().invert_yaxis()

	count=0
	for vm_file in vm_models :
		vm = open(vm_file, 'r').read()
		nlayers=str(int(re.split('\n', vm)[1][0:5]))
		vm_list=(re.split('\n '+nlayers+'\n', vm)[0]).split(
			'\n '+nlayers+'  ')[1].split('\n')[1::]

		VP[vm_file] = []
		Z[vm_file] = []

		for l in range(len(vm_list)) :
			vel=float(vm_list[l][1:6])
			dep=float(vm_list[l][10:17])

			if l == 0:
				Z[vm_file].append(dep)
				VP[vm_file].append(vel)
				VP[vm_file].append(vel)
			if l > 0:
				Z[vm_file].append(dep)
				Z[vm_file].append(dep)
				VP[vm_file].append(vel)
				VP[vm_file].append(vel)

		Z[vm_file].append(Z[vm_file][-1]+50)

		if vm_file in all_models and vm_file not in initial_models : #intermediate models
			if vm_file == all_models[0] :
				plt.plot(VP[vm_file],Z[vm_file], '-', color='grey', linewidth=0.5,
					label='Best '+str(len(all_models))+' models')
			plt.plot(VP[vm_file],Z[vm_file], '-', color='grey', linewidth=0.5)

		if vm_file in initial_models: #initial and final models
			plt.plot(VP[vm_file],Z[vm_file], '--', color=color_list[count], 
				linewidth=2, label=vm_file.split('.mod')[0][2::], alpha=0.8)
			count=count+1

		if final_models and vm_file in final_models: #initial and final models
			plt.plot(VP[vm_file],Z[vm_file], '-', color=color_list[count], 
				linewidth=2, label=vm_file.split('.mod')[0][2::], alpha=0.8)
			count=count+1

		if vm_file in min1D_model : #min 1D
			plt.plot(VP[vm_file],Z[vm_file], '-', color='red', linewidth=4, 
				label=vm_file.split('.mod')[0], alpha=0.8)

	plt.grid()
	plt.legend(loc='upper right', fontsize=14)
	plt.xticks(fontsize=16)
	plt.yticks(fontsize=16)
	plt.xlabel('Velocity [km/s]', fontsize=18)
	plt.ylabel('Depth [km]', fontsize=18)

	plt.savefig(outfile, bbox_inches='tight')
	plt.close()

