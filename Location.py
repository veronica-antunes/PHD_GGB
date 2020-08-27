
###########################################################################################
# Veronica Antunes 03/07/2017                                                             #
# Last update 16/2020                                                                     #
# veronica.antunes@unige.ch                                                               #
# Aim: Automatically locate Nirano signals using a sliding window                         #
###########################################################################################

from __future__ import division # to allow division using /

import matplotlib.pyplot as plt
import numpy as np
from obspy import read
from glob import glob


def read_stafile(sta_file):

	from math import sin, cos
	from obspy.geodetics import gps2dist_azimuth

	stat_file = open(sta_file)		#read the file with the stat locations
	stat_lines = stat_file.readlines()

	x_stat=[]
	y_stat=[]
	stat=[]

	for line in stat_lines :

		if line != '\n':

			fields = line.strip().split()

			if fields[0] == str('REF') :

				ref_lat=float(fields[1])
				ref_lon=float(fields[2])

			if fields[0] != str('REF') and fields[0] != str('station'):
				try:
					sta_lat=float(fields[1])
					sta_lon=float(fields[2])

				except:
					raise ValueError('ERROR: Incorrect file format!' +
						'Check file: '+sta_file)


				stat.append(fields[0])

				try: 
					dist, az, baz = gps2dist_azimuth(sta_lat, sta_lon, 
						ref_lat, ref_lon)
				except :
					raise ValueError('ERROR: Missing REF station in file : '+sta_file)

				hipotenusa=dist
				angle=baz

				alpha=90-angle

				x=cos(np.deg2rad(alpha))*hipotenusa
				y=sin(np.deg2rad(alpha))*hipotenusa

				x_stat.append(x)
				y_stat.append(y)

	return stat, x_stat, y_stat

def step1_grid(sta_file, x_start, x_end, y_start, y_end, step, outfile=None) :

	# Create the grid: np.mgrid creates automatically the grid ##########
	#station_file='stations_file.txt' #should be in the same format as the "template"
	# Grid ####
	#x_start=0
	#x_end=1000
	#y_start=0
	#y_end=800
	#step=1

	from numpy import mgrid

	print ('Creating grid...')

	plt.close('all')

	glob=np.mgrid[x_start:x_end+step:step, y_start:y_end+step:step]


	x_glob=glob[0]
	y_glob=glob[1]

	ax = plt.subplot(111)

	plt.plot(x_glob, y_glob, '.', color='0.75') 
				     #color: give a number between 0-1 to have gray shadows

	ax.set_xlim([0-step, x_end+step])
	ax.set_ylim([0-step, y_end+step])

	print ('Converting coordinates...')
	stat, x_stat, y_stat = read_stafile(sta_file)

	stat_struct=np.array( (stat, x_stat, y_stat) ).T

	plt.plot(x_stat, y_stat, '^', color='b', markersize=14, markeredgecolor='darkgrey') 

	plt.xlabel('distance [m]', size=20)
	plt.ylabel('distance [m]', size=20)


	plt.xticks(fontsize=18)
	plt.yticks(fontsize=18)

	if outfile:
		plt.savefig(outfile, bbox_inches='tight' )
	#else :
	#	plt.savefig('grid.png', bbox_inches='tight' )


	np.save( 'stat_struct', stat_struct)
	np.save( 'mesh_struct', glob)

	print ('Step 1 Done!')


def step2_crosscorr(wavs_file, delta_time_w) :

	from scipy.signal import hilbert, chirp
	from obspy.signal.cross_correlation import xcorr_pick_correction, xcorr, xcorr_max
	from scipy.signal import savgol_filter
	from numpy import array, sign, zeros
	from scipy.interpolate import interp1d

	## INPUT VALUES
	#wavs_file='filtered_wavs.mseed'	#waveforms with all Z channels
	#delta_time_w=1				#delta time window in (seconds) 
	#slide_w=5

	# variables for delta window ###
	wavs=read(wavs_file)
	samp_rate=wavs[0].stats.sampling_rate		#sample rate in Hz
	nstat=len(wavs)

	##### Do small windows for cross_correlation
	#nt=int(time_duration*samp_rate) 	#number of time sample, legnth of signal
						# vector #this is for 20 min of data
	nt=wavs[0].stats.npts
	dt=int(delta_time_w*samp_rate) 		#number of time samples, short window

	#To check if the trace size is lower than the last point
	check=range(0,nt,dt)[-1]+dt
	if check > nt :
		nt=range(0,nt,dt)[-1]
		print 'trace is too long... resizing for', nt, 'npts'


	n_stat=int(nstat) 		#number of stations

	#max_dt=int((dt/2)+1)		# max shift in samples (for obspy cross-correlation)
	max_dt=dt #(dt*2)-1		# max shift in samples (for python cross-correlation full)
	n_pair=int(n_stat*(n_stat-1)/2)	#number of stat pairs
	n_segments=int(nt/dt)		#number of segments
	#nslide=int(dt/4)		# sliding window for cc. in nsamples obspy method


	#t=linspace(0, nt, nt) #in case we want to plot. put in x axis

	pair=[]
	si=[]
	sj=[]

	#xcorr_struct=(dt, n_segments, n_pair) ex: xcorr_struct=(251, 240, 21)
	xcorr_struct=np.zeros( (max_dt,n_segments,n_pair) ) #size fct, n segm, n stat pairs

	env_struct=np.zeros( (max_dt,n_segments,n_pair) ) #size fct, n segm, n stat pairs

	senv_struct=np.zeros( (max_dt,n_segments,n_pair) ) #size fct, n segm, n stat pairs

	npair=-1	# to start with 0

	for i in range(0,n_stat,1):
		for j in range(i+1,n_stat,1):

			namestat1=wavs[i].stats.station
			namestat2=wavs[j].stats.station

			print  'getting cc and env for pair...', str(namestat1), str(namestat2)
			npair=npair+1

			#tr1=stat[i-1]	#syntetics
			#tr2=stat[j-1]	#syntetics

			tr1=wavs.select(station=namestat1) 
			if tr1 :
				tr1=tr1[0].data
			else:
				tr1=np.zeros((nt))


			tr2=wavs.select(station=namestat2)
			if tr2 :
				tr2=tr2[0].data
			else:
				tr2=np.zeros((nt))

			si.append( namestat1 )
			sj.append( namestat2 )
			pair.append(npair)

			nsegm=-1	# to start with 0

			for ti in range(0,nt,dt): #int(nt)

				nsegm=nsegm+1
				tf=int(ti+dt)

				sg1=tr1[ti:tf]
				sg2=tr2[ti:tf]
				sg1=np.reshape(sg1, (dt))
				sg2=np.reshape(sg2, (dt))


				### Trace interpolation
				#max_sg1=sg1.max()
				#max_sg2=sg2.max()
				#i_sg1=sg1/max_sg1
				#i_sg2=sg2/max_sg2

				#cross-correlation using obspy. normalized cc function
				#index,value1,fct=xcorr(sg2, sg1, nslide, full_xcorr=True)
				#xcorr_struct[ : , nsegm , npair] = fct

				#cross-correlation using python
				fct=np.correlate(sg2, sg1, mode='same')
				xcorr_struct[ : , nsegm , npair] = fct

				#envelope
				analytic_signal = hilbert(fct)
				env = np.abs(analytic_signal)
				env_struct[ : , nsegm , npair] = env

				#smoothing the envelope
				#senv=interpolate.interp1d(env, kind="linear")
				#senv=savgol_filter(env, 75, 6,mode='constant')
				#env_struct[ : , nsegm , npair] = senv

				#s=env
				#q_u = zeros(s.shape)
				#u_x = [0,]
				#u_y = [s[0],]
				#for k in xrange(1,len(s)-1):
				#	if (sign(s[k]-s[k-1])==1) and (sign(s[k]-s[k+1])==1):
				#		u_x.append(k)
				#		u_y.append(s[k])
				#u_x.append(len(s)-1)
				#u_y.append(s[-1])
				#u_p = interp1d(u_x,u_y, kind = 'cubic',
					#bounds_error = False, fill_value=0.0)
				#for k in xrange(0,len(s)):
				#	q_u[k] = u_p(k)
				#senv=q_u
				#senv_struct[ : , nsegm , npair] = senv


	##a=np.array((pair, sj, si)).T

	np.save( 'pair', np.array((pair, si, sj)).T )
	np.save( 'xcorr_struct', xcorr_struct)
	np.save( 'env_struct', env_struct)
	#np.save( 'senv_struct', senv_struct)
	print 'Step 2 Done!'


def step3_difftime(vel1, vel2, step) :

	import math as m
	import os

	#input values
	#vel_string=np.arange(380,385, 5)
	vel_string=np.arange(vel1, vel2, step)

	#import files that we built in step1 and step2
	xcorr_struct=np.load('xcorr_struct.npy')	# cross correlation
	env_struct=np.load('env_struct.npy')		# envelope	
	mesh_struct=np.load('mesh_struct.npy')		# xy grid in meters

	stat_list=np.load('stat_struct.npy')		# station coordinates in meters
	pair_list=np.load('pair.npy')			# pair and stations	

	#remove old stuff
	remove=glob('A_array_*')
	remove+=glob('B_array_*')
	remove+=glob('C_array_*')
	remove+=glob('D_array_*')
	remove+=glob('E_array_*')

	for r in remove:
		os.remove(r)

	## Do calculations for each point of the grid

	pair=len(pair_list)
	gx=mesh_struct.shape[1]	#len(mesh_struct[0,:,0])
	gy=mesh_struct.shape[2]	#len(mesh_struct[0,0,:])


	for v in range(len(vel_string)) :
		vel=vel_string[v]
		print ('calculating for velocity: ' + str(vel) + ' ...')
		A=np.zeros( (gx,gy,pair) ) #size pair, grid x, grid y ex: (21, 17, 21)
		for p in range(len(pair_list)) :

			npair, s1, s2=pair_list[p]

			print 'getting dt for station pair', p, '... ', s1, s2

			#find the raw and column corresponding to certain station
			s1n=np.where(stat_list==s1)[0][0]
			s2n=np.where(stat_list==s2)[0][0]

			#station pair coordinates
			st1, x1, y1=stat_list[s1n]
			st2, x2, y2=stat_list[s2n]


			for i in range(0, gx):
				for j in range(0, gy):

					xi, yi= mesh_struct[:, i, j]

					d1=m.sqrt( (float(xi)-float(x1))**2 + 
						(float(yi)-float(y1))**2 )

					d2=m.sqrt( (float(xi)-float(x2))**2 + 
						(float(yi)-float(y2))**2 )

					t1=d1/vel
					t2=d2/vel

					dt=t2-t1

					A[ i , j , p] = dt

		np.save( 'A_array_'+str("%04d" % vel), A)

	print 'Step 3 Done!'



def step4_envvalue(wavs_file, outfile=None, point=None) :

	from numpy import isnan

	mesh_struct=np.load('mesh_struct.npy')	# array with step1 results: xy grid in meters
	xcorr_struct=np.load('xcorr_struct.npy')# array with step2 results: cross correlation
	env_struct=np.load('env_struct.npy')	# array with step2 results: envelope	
	#dt_vector=np.load('A_array.npy')	# array with step3 results: dt


	if point :
		i1,j1,p1,w1 = point
		#print point, i1, j1, p1, w1

	wavs=read(wavs_file)
	samp_rate=wavs[0].stats.sampling_rate	#sample rate in Hz

	#dt=int(delta_time_w*samp_rate) 	#number of time samples, short window
	#dt=(env_struct.shape[0]*2)-2
	dt=((env_struct.shape[0]))

	nt=wavs[0].stats.npts

	nslide=int(dt/2)
	#nslide=int(dt)

	#dt_vector=np.load('A_array_*.npy')

	#import array lists
	stat_list=np.load('stat_struct.npy')		# station coordinates in meters
	pair_list=np.load('pair.npy')			# pair and stations	


	# Getting array sizes
	gp=len(pair_list)		#number of pairs
	gx=len(mesh_struct[0,:,0])	#number of x grid elements
	gy=len(mesh_struct[0,0,:])	#number of y grid elements
	gw=len(env_struct[0,:,0])	#number of time window segments

	#Create 4D matrix to store results

	vel_file=sorted(glob('A_array_*'))
	#vel_fileB=sorted(glob('B_array_*'))

	for v in range(len(vel_file)) :
		#vel=vel_file[v][8:11]
		vel=vel_file[v].split('_')[-1].split('.npy')[0]

		print ('calculating for velocity: ' + str(vel) + ' ...')

		dt_vector=np.load(vel_file[v])

		B=np.zeros( (gx, gy, gw, gp) ) #size grid x, grid y, window size, pair
		D=np.zeros( (gx, gy, gw, gp) ) #size grid x, grid y, window size, 
							#pair to do normalization

		## Calculations for each pair (gp), each point of the grid (gx, gy), 
			#each window segment (gw) ###
		for p in range(gp): #(gp) :

			npair, s1, s2=pair_list[p]

			print 'getting yenv for dt, station pair', p, '... ', s1, s2

			for i in range(gx): #(gx):
				for j in range(gy): #(gy):

		#			xi, yi= mesh_struct[:, i, j]

					for w in range(gw): #(gw):

						dtstat = dt_vector[i, j, p]	
					# time differece betwen station pair p (seconds)

						dsstat= float(dtstat*samp_rate)	
					# time differece betwen station pair p (samples)

						if dsstat < ((nslide-1)*-1) or dsstat > (nslide-1) :
							#print ('dt id too large, 
								#0 value will be used instead')
							yenv=0
						else :

							xenv=nslide+dsstat
							# The middle point '0' is actually 
								#nslide(125)

							env=env_struct[:, w, p]

							#Without interpolation
			#				yenv=env[int(xenv)]
				

							#Interpolation to get the value of yenv
							xenv_1=int(xenv)
							yenv_1=env[xenv_1]

							xenv_2=int(xenv+1)

	
							yenv_2=env[xenv_2]

							xenv_decimal=xenv-xenv_1
							yenv_decimal=(yenv_2-yenv_1)*xenv_decimal

							yenv=yenv_1+yenv_decimal

						if isnan(yenv) :	#in the case yenv is nan 
									#(when channels are missing)

							yenv=0

						B[ i , j, w , p] = yenv

		## To test if necessary ###################################################
						if point :
							if p==p1 and i==i1 and j==j1 and w==w1 :
							### Plot function
		#					print dtstat, dsstat, xenv, yenv
								print 'plotting'
								plt.close('all')
								new_xlabel=np.linspace(-nslide/samp_rate, 
									nslide/samp_rate, len(env))
								plt.plot(new_xlabel,env)
								new_xenv=(xenv-nslide)/samp_rate

								plt.axvline(x=new_xenv, 
									color='g',
									linestyle='-',
									label='xenv='+
									str(round(new_xenv,3)))
								plt.axhline(y=yenv, 
									color='r',
									linestyle='-' ,
									label='yenv='+
									str(round(yenv,5)))

								plt.legend(loc="upper left", fontsize=16)
								plt.title('Grid point ('+ str(i)+','+
									str(j)+	'), pair '+str(p) + 
									', window ' + str(w), size=20)

								plt.xlim(-nslide/samp_rate,
									nslide/samp_rate)

								plt.xlabel('time [s]', size=20)
								plt.ylabel('Env-CC value', size=20)

								plt.xticks(np.linspace(-nslide/samp_rate, 
									nslide/samp_rate, 5), size=18)
								plt.yticks(size=18)
							#plt.show()

								if outfile:
									plt.savefig(outfile,
										bbox_inches='tight')
								else:
									plt.savefig(str(i)+ '-' + 
										str(j)+ '_' + str(p)+ 
										'_' + str(w) + '.png',
										bbox_inches='tight' )
									
								plt.close()

	#########################################################################
		####### For normalization #################################################
		
	#		#B=np.load(vel_fileB[v])		# array with step3 results: dt without
	#		print 'Normalizing amplitude values for pair', p , '...'

	#		max_value=np.amax(B[ :, :, :, p])

	#		for i in range(gx): #(gx):
	#			for j in range(gy): #(gy):
	#				for w in range(gw): #(gw):

	#					amplt_value=B[i, j, w, p]
			
	#					norm_value=float(amplt_value/max_value)

	#					D[ i , j, w, p ] = norm_value

	#	np.save( 'D_array_'+str("%03d" % int(vel) ), D)
		###########################################################################

		np.save( 'B_array_'+str("%04d" % int(vel) ), B)

	print 'Step 4 Done!'


def step5_stackall(wavs_file) :

	print 'importing files... '

	#import files that we built in step1 step2 and step3
	mesh_struct=np.load('mesh_struct.npy')	# array with step1 results: xy grid in meters
	xcorr_struct=np.load('xcorr_struct.npy')# array with step2 results: cross correlation
	env_struct=np.load('env_struct.npy')	# array with step2 results: envelope	

	#dt_vector=np.load('A_array.npy')	# array with step3 results: dt

	#yenv_vector=np.load('B_array.npy')	# array with step3 results: dt without normalization
	#yenv_vector=np.load('D_array.npy')	# array with step3 results: dt with normalization

	#read wavs
	wavs=read(wavs_file)

	samp_rate=wavs[0].stats.sampling_rate		#sample rate in Hz

	nt=wavs[0].stats.npts

	#dt=((env_struct.shape[0])+1)/2
	#dt=int(delta_time_w*samp_rate) 	#number of time samples, short window
	#dt=int(delta_time_w*samp_rate) 	#number of time samples, short window
	#nslide=int(dt/4)


	#import array lists
	stat_list=np.load('stat_struct.npy')		# station coordinates in meters
	pair_list=np.load('pair.npy')			# pair and stations	

	# Getting array sizes
	gp=len(pair_list)		#number of pairs
	gx=len(mesh_struct[0,:,0])	#number of x grid elements
	gy=len(mesh_struct[0,0,:])	#number of y grid elements
	gw=len(env_struct[0,:,0])	#number of time window segments


	vel_fileB=sorted(glob('B_array_*'))

	for v in range(len(vel_fileB)) :

		#vel=vel_fileB[v][8:11]
		vel=vel_fileB[v].split('_')[-1].split('.npy')[0]
		print ('calculating for velocity: ' + str(vel) + ' ...')

		yenv_vector=np.load(vel_fileB[v])	# array with step3 results: dt without normalization


		#Create 3D matrix to store results
		C=np.zeros( (gx, gy, gw) ) #size grid x, grid y, window segments
		E=np.zeros( (gx, gy, gw) ) #for normalization 2


		## Get envelope for all pairs for each point of the grid (gx, gy), 
			# each window segment (gw) ###
		print 'Get absolute value for amplitude envelope ...'

		for i in range(gx): #(gx):
			for j in range(gy): #(gy):
				for w in range(gw): #(gw):

					ampl=yenv_vector[i, j, w, :]
			
					amplt_s=sum(ampl)

					C[ i , j, w ] = amplt_s

		np.save( 'C_array_'+str("%04d" % int(vel) ), C)

		# For normalization
		print 'Normalizing amplitude values ...'
		max_value=np.amax(C[:,:,:])

		for i in range(gx): #(gx):
			for j in range(gy): #(gy):
				for w in range(gw): #(gw):

					amplt_value=C[i, j, w]

					if max_value== 0:

						norm_value=0

					else :
						norm_value=float(amplt_value/max_value)

					E[ i , j, w ] = norm_value


		np.save( 'E_array_'+str("%04d" % int(vel)), E)

	print ('Step 5 Done!')

def step6_plotall(name_folder, wavs_file, ntotal_stat=None, normalized=False, source_loc_file = None,
	source_loc = None) :

	import matplotlib
	from matplotlib import dates, ticker, cm
	from matplotlib.colors import LinearSegmentedColormap
	from matplotlib.dates import date2num, num2date, DateFormatter
	from matplotlib.ticker import AutoMinorLocator, FormatStrFormatter
	import matplotlib.image as mpimg
	#from decimal import Decimal
	from obspy.imaging.cm import pqlx, obspy_divergent
	import os
	#from os import system

	#INPUT PARAMETERS
	#path to save files
	#save_path='new/180828-1'
	#name_folder='d2-3stat-'
	#wavs_file='filtered_wavs.mseed'			#waveforms with all Z channels

	if not ntotal_stat :
		#number_stat=False takes the nstat from trace. True takes the given ntotal_stat value
		number_stat = False
		ntotal_stat=ntotal_stat

	if ntotal_stat :
		number_stat = True


	colours=['b','forestgreen','r','magenta','BlueViolet','SaddleBrown', 'grey', 'darkorange', 'y',
		'lightgreen','cyan','lightpink','darkgrey','k','g','purple','brown','orange','olive',
		'darkred','darkblue','darkgreen']
	colours2=['b','forestgreen','r','magenta','BlueViolet','SaddleBrown', 'grey', 'darkorange', 'y',
		'lightgreen','cyan','lightpink','darkgrey','k','g','purple','brown','orange','olive',
		'darkred','darkblue','darkgreen']

	C=np.load(glob('E_array_*.npy')[0]) #with normalization1, with normalization2

	print 'importing files... '

	wavs=read(wavs_file)

	#stat_net
	#station_name=[]
	#for sta_name in wavs:
	#	station_name.append(sta_name.stats.station)

	nstat=len(wavs)					#number of station
	samp_rate=wavs[0].stats.sampling_rate		#sample rate in Hz

	t0=wavs[0].stats.starttime
	t1=wavs[0].stats.endtime			
	date=t0.strftime(format='%Y-%m-%d %H:%M:%S')

	nt=wavs[0].stats.npts
	#nt=int(time_duration*samp_rate) #number of time sample, legnth of signal vector 
	#this is for 20 min of data
	env_struct=np.load('env_struct.npy')
	#dt=int(delta_time_w*samp_rate) 	#number of time samples, short window
	dt=(env_struct.shape[0])
	
	dtw=float(dt /samp_rate)

	ts=int((nt/samp_rate)/1)  #time ticks in seconds for xplot UTC

	#import matlab files
	#synthetic_data = loadmat(mat_file) 	#, squeeze_me=True, struct_as_record=False)
	#synthetic_data = synthetic_data['signals']
	#stat=np.hsplit(synthetic_data, 7) 	#split by 7 stations

	#A=np.load('A_array.npy')
	#B=np.load('B_array.npy')


	stat_net=np.load('stat_struct.npy')
	mesh=np.load('mesh_struct.npy')
	pair_list=np.load('pair.npy')

	matplotlib.rcParams.update({'font.size': 14})

	# Getting array sizes

	gw=len(C[0,0,:])	#number of time window segments

	#x_mesh=mesh[0]
	#y_mesh=mesh[1]
	#plt.plot(x_mesh, y_mesh, '.', color='0.75') 


	x_stat=stat_net[:,1]
	y_stat=stat_net[:,2]
	stat_list=stat_net[:,0]

	#plt.plot(x_stat, y_stat, '^', color='b') 

	#ax.set_xlim([0-50, 1000+50])
	#ax.set_ylim([0-50, 800+50])

	#region limits, takes automatically the values from mesh matrix
	step=mesh[0][1][0]
	corner=0
	xstart=mesh[0][0][0]
	xend=mesh[0][-1][-1]
	ystart=mesh[1][0][0]
	yend=mesh[1][-1][-1]

	#try to make movie
	#FFMpegWriter = animation.writers #['ImageMagickFileWriter']
	#metadata = dict(title='Movie Test', artist='Matplotlib',
	#               comment='Movie support!')
	#writer = FFMpegWriter(fps=15, metadata=metadata)
	#files = []

	print 'plotting waveforms... '

	#fig = plt.figure() #[20,18]15.15   figsize=(10,8)
	#ax1.plot(x, y)
	#other ideas to improve the graphs
	#fig, ax2 = plt.subplots(nrows=1, ncols=2)
	#ax2=plt.subplots(nrows=len(stat), sharex=True)

	###########################################################
	#fig, ax1 = plt.subplots(nrows=len(stat), sharex=True)
	#for s in range(len(stat)) :
	#	ax1[s].plot(stat[s])
	###########################################################

	#count=-1
	#for s in range(len(stat)) :
	#	count=count+2
	#	plt.subplot(13,2,count)
	#	ax1=plt.plot(stat[s])

	#To check if the trace size is lower than the last point
	#check=range(0,nt,dt)[-1]+dt
	#if check > nt :
	#	nt=range(0,nt,dt)[-1]
	#	print 'trace is too long. resizing for', nt, 'npts'

	if number_stat == False:
		ntotal_stat=len(wavs)
		stat_net=[]
		for tr in wavs:
			stat_net.append(str(tr.stats.station))
		#print stat_net
		#stat_net=sorted(stat_net)

	vel_fileC=sorted(glob('E_array_*'))

	for v in range(len(vel_fileC)) :

		#vel=vel_fileC[v][8:11]
		vel=vel_fileC[v].split('_')[-1].split('.npy')[0]
		print ('doing map for velocity: '+ vel +' m/s')

		C=np.load(vel_fileC[v])	# array with step3 results: dt without normalization


		path=name_folder+'-'+str(dtw)+'_'+str(int(vel))
		path=path[:-len(path.split('/')[-1])]
		if not os.path.exists(path):
			try:
				os.makedirs(path)
			except OSError as exc: # Guard against race condition
				if exc.errno != errno.EEXIST:
					raise
		path=name_folder+'-'+str(dtw)+'_'+str(int(vel))

		print ( 'plotting location maps in: ' + path )


		for w in range(gw): #(20,21,1) : #(gw)

			tw=range(0,nt+1,dt)[w]

			plt.close('all')
			fig=plt.figure() #[20,18]15.15 figsize=(10,8) #figsize=(10,6)
			fig.suptitle(str(date)[0:10], fontsize=16, fontweight='bold')
			#fig.subplots_adjust(hspace=.08)	
			#plot waveforms
			count=-1
			for s in range(int(ntotal_stat)) :
				count=count+1
		#		ax1=fig.add_subplot(7,2,count) #VERY GOOD TOO
				size_plot=int(62)
				ax1=plt.subplot2grid((ntotal_stat, size_plot), 
					(count, 0), colspan=11)
				ax2=plt.subplot2grid((ntotal_stat, size_plot), 
					(count, 12), colspan=8)

				if number_stat == True :		
					stat=wavs.select( station=stat_net[s][0] )

				else:
					stat=wavs.select( station=stat_net[s] )


				#plot x axis with time
				start =date2num(t0.datetime)
				end = date2num(t1.datetime)

				#if count==0 :
				#	ax1.set_title(str(date)[0:10], fontsize=10)


				if stat :
					stat=stat[0].data

				else:
					stat=np.zeros((nt))
					ax1.yaxis.set_ticks(np.linspace(-1, 1, 3 ))

				#plot x axis with time
				t = np.linspace(start, end, nt) #nt+1)
				if t.size == stat.size:
					t=t
				else:
					print ('changing t size ')
					t = np.linspace(start, end, stat.size)

				asd = np.linspace(0, stat.size/samp_rate, stat.size )

				ax1.plot_date(t, stat, '-' ,linewidth=0.15, color=colours[s])
				ax2.plot(asd, stat, linewidth=0.3, color=colours[s])
				#ax2.plot(stat[tw:tw+dt], linewidth=0.3, color=colours[s])



				#date_tw=np.linspace(start,end,gw+1)[w]
				#date_tw2=np.linspace(start,end,gw+1)[w+1]

				date_tw=date2num((t0+(tw/samp_rate)).datetime)
				date_tw2=date2num((t0+(tw+dt)/samp_rate).datetime)
				ax1.axvline(x=date_tw, color='r', alpha=0.5) 
				ax1.axvline(x=date_tw2, color='r', alpha=0.5)

				if normalized :
					starty, endy = [-1,1]

				else:
					starty, endy = ax1.get_ylim() #get trace standard starty, endy 
					#starty = "{:.0E}".format(Decimal(starty))
					#endy = '%.0e' % Decimal(endy)

				#ax1.yaxis.set_ticks(np.linspace(starty, endy, 3 ))
				ax1.yaxis.set_ticks([0,endy])

				#ax2.yaxis.tick_right()

				#startx, endx = ax1.get_xlim()
				#ax1.xaxis.set_ticks(np.arange(startx, endx+0.001, ( abs(startx
										#)+abs(endx) )/4 ))
#				startx2, endx2 = ax2.get_xlim()

				ax2.set_xlim([tw/samp_rate, (tw+dt)/samp_rate])

				
				ax2.xaxis.set_ticks(np.linspace( tw/samp_rate , 
					(tw+dt)/samp_rate, 2 ) )


	#			ax2.set_ylim([tw,tw+dt])

				#set on this to get y scale diferent
				#starty = np.amax(stat[tw:tw+dt])
				#endy = np.amin(stat[tw:tw+dt])
				#ax2.set_ylim([starty,endy])

				ax2.yaxis.set_ticks(np.linspace(starty, endy, 3 ))

	#			ax2.yaxis.set_ticks(np.linspace(6e-9, -6e-9, 3 ))

				if endy < 0.001 :

					ax1.yaxis.set_major_formatter(FormatStrFormatter('%.0E'))
				#ax1.ticklabel_format(style='sci', axis='y') #scilimits=(-5,-20)) 
				#, useMathText=True, useLocale=True)
				#ax1.ticklabel_format(style='sci', scilimits=(0,0), axis='x')
				#ax2.ticklabel_format(style='sci', scilimits=(0,0), axis='y')

				#change size and position of 1e in the plot
				#ax1.get_yaxis().get_offset_text().set_position((0.1,-0.9))
				#ax2.get_yaxis().get_offset_text().set_position((1,0))

				#change size and size of 1e in the plot
				#ax2.get_yaxis().get_offset_text().set_size(6)
				#ax1.get_yaxis().get_offset_text().set_size(6)

				ax1.tick_params(labelsize=12, axis='y')
				#ax2.tick_params(labelsize=12, axis='y')
				ax1.tick_params(labelsize=12, axis='x')
				ax2.tick_params(labelsize=12, axis='x')

				#ax1.get_yaxis().get_offset_text().set_rotation(0)

				#hide 1e in the plot and add it as axis lable
				ax1.yaxis.offsetText.set_visible(False)

				#print '%.0e' % Decimal(starty),'%.0e' % Decimal(endy)

				#offset = "{:.0E}".format(Decimal(starty))
				#offset = str(offset)[3:6]

				#ax1.yaxis.set_label_text( '1e '+offset, fontsize=7)

				ax2.yaxis.offsetText.set_visible(False)
				ax2.axes.get_yaxis().set_ticklabels([])
				#offset2 = ax2.yaxis.get_major_formatter().get_offset()
				#ax2.yaxis.set_label_text('ampl (' + offset2 + ')')




				# TO get proper ticks automatically
				a=float(end-start)/4
				b=(end-start)
				c=float(a*ts)/b

			

				ax1.xaxis.set_major_formatter( DateFormatter('%M:%S') )
				ax1.xaxis.set_major_locator(dates.SecondLocator(interval=int(c)))
				ax1.xaxis.set_minor_locator(AutoMinorLocator(5))

				xticks = ax1.xaxis.get_major_ticks()

				if len(xticks) == 4 : 

					xticks[1].label1.set_visible(False)
					xticks[-1].label1.set_visible(False)

				if len(xticks) == 5 : 
					xticks[0].label1.set_visible(False)
					xticks[2].label1.set_visible(False)
					xticks[-1].label1.set_visible(False)


				if s < ntotal_stat-1 :
					ax1.axes.get_xaxis().set_ticklabels([]) #to remove the labels
										#on top figures
					ax2.axes.get_xaxis().set_ticklabels([])


				#ax2.axes.get_yaxis().set_ticklabels([])

				if s ==  ntotal_stat-1 :
					#ax1.set_position((.1, .3, .3, .1)) 
					# to make a bit of room for extra text
					ax1.set_xlabel("time [MM:SS]\nafter "+str(date)[11::]+ 
						' UTC', fontsize=12)

					ax2.set_xlabel("[s]", fontsize=12)

		#	plt.show()

			############ END WAVS PLOT ##################################

			img=C[:,:,w]	#results for each time window

		#	ax2=fig.add_subplot(1,2,2) #VERY GOOD TOO
			ax2 = plt.subplot2grid((ntotal_stat, size_plot), (0, 29) , rowspan=50,
				 colspan=40)

			plt.scatter(x_stat, y_stat, color='white', marker='v',  edgecolor='black',
				linewidth='0.8', s=80)
			c2=0
			for each_stat in stat_net :

				loc=np.where(stat_list == each_stat)[0][0]
				#loc=[pos for pos, char in enumerate(stat_list) if char == each_stat][0]
				plt.scatter(float(x_stat[loc]), float(y_stat[loc]),
					color=colours2[c2], marker='v',
					edgecolor='white', linewidth=0.8, s=80) #, s=80) 
				c2=c2+1

#			for c2 in range(len(x_stat)):
#				plt.scatter(float(x_stat[c2]), float(y_stat[c2]), 
#					color=colours2[c2], marker='v',
#					edgecolor='white', linewidth=0.8, s=80) #, s=80) 

			if source_loc:

				shot, x_shot, y_shot=read_stafile(source_loc_file)

				loc=[pos for pos, char in enumerate(shot) if char == source_loc][0]
				plt.scatter(float(x_shot[loc]), float(y_shot[loc]),
					color='darkorange', marker=(5, 1), #darkorange
					edgecolor='k', 
					linewidth=1.0, s=100) #, s=80) 

			vmax=1
			cmap = LinearSegmentedColormap.from_list('mycmap', 
				[(0 / vmax, 'darkblue'),
				(0.50 / vmax, 'blue'), 
				(0.65 / vmax, 'cyan'),
				(0.73 / vmax, 'lightgreen'),
				(0.81 / vmax, 'yellow'),
				(0.86 / vmax, 'orange'),
				(0.94 / vmax, 'red'),
				(0.99 / vmax, 'darkred'),
				(1.00 / vmax, 'black')])


			img2=plt.imshow( img.T, 
				cmap=cmap, clim=(0, 1), #'jet'
				interpolation='none',	#interpolation: nearest or bicubic
				origin='lower',
				extent=[xstart-corner,xend+corner,ystart-corner,yend+corner] ) 
				#extent=[600,800,400,600] ) 
		##		aspect='auto' or 'equal'
		##		clim=(-0.25, 0.25)

			#plt.plot(x_mesh, y_mesh, '.', color='0.75')

			colorticks=np.linspace(0,1,6)

			cbar=plt.colorbar(orientation='horizontal', ticks=colorticks) #, shrink=0.7)
			cbar.ax.set_xlabel('Source likelyhood')
				#'stacked back-projections \n normalized')
				#rotation=180)

			if xend <= 1000 :

				ax2.set_xlim([xstart-step, xend+step])
				ax2.set_ylim([ystart-step, yend+step])
				#ax2.set_xlim([450, 850])
				#ax2.set_ylim([350, 650])

				ax2.set_xlabel('distance [m]', fontsize=14)
				ax2.set_ylabel('distance [m]', fontsize=14)

			if xend > 1000 :

				ticks_x = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x/1000.0))
				ax2.xaxis.set_major_formatter(ticks_x)

				ticks_y = ticker.FuncFormatter(lambda y, pos: '{0:g}'.format(y/1000.0))
				ax2.yaxis.set_major_formatter(ticks_y)

				ax2.set_xlabel('distance [km]', fontsize=14)
				ax2.set_ylabel('distance [km]', fontsize=14)


			ax2.tick_params(labelsize=14)


	#		ax2.set_xlim([500, 800])
	#		ax2.set_ylim([400, 600])

			ax2.set_title('Window '+str( "%04d" % w)+'  vel: '+str(int(vel))+'m/s\n',
				fontsize=15)

		#	plt.tight_layout()


			#print (path+'/w_'+str( "%03d" % w)+'_.png')

			pdf=plt

			pdf.savefig(path+'_w_'+str( "%04d" % w)+'_'+
				str(int(vel))+'.pdf', bbox_inches='tight' )

			plt.savefig(path+'_w_'+str( "%04d" % w)+'_'+
				str(int(vel))+'.png', bbox_inches='tight' )
		#	plt.show()
			plt.close()
		#	writer.grab_frame()

	
		#print 'making movie'

	#os.chdir('/home/veronica/Dropbox/scripts/z.nir_locations/test_step5' #+str(pfiltered_wavs.mseedath))

	#	os.system("ffmpeg -framerate 10 -pattern_type glob -i 'test_step5/*plot_*.png' -c:v 
	#libx264 -r 30 -pix_fmt yuv420p test_step5/nirano_test.mp4") 
	#framerate gives the number of frames per second

	#	os.rename('test_step5/nirano_test.mp4', 'test_step5/nirano_test_v'+str(int(vel))+'.mp4' )

	#os.chdir('/home/veronica/Dropbox/scripts/z.nir_locations')

	#os.system("mencoder /home/veronica/Dropbox/scripts/z.nir_locations/test_step5/plot*.png 
	#-mf type=png:fps=10 -ovc lavc -lavcopts vcodec=wmv2 -oac copy -o nirano_movie2.mpg") 



	#try to make movie
	#from glob import glob
	#fname=glob('test_step5/*')
	#ani = animation.FuncAnimation(fig, updatefig, interval=50, blit=True)
	#ani=animation.ArtistAnimation(fig, img, interval=20, blit=True, repeat_delay=0)
	#ani.save('MovWave.mpeg', writer="ffmpeg")
	#writer.save('demo.mp4') #,writer=writer,dpi=100)
	#writer.saving(fig, "writer_test.mp4", 100)

	print ('Step 6 Done!')


def help_location() :

	print ('step1_grid(sta_file, x_start, x_end, y_start, y_end, step)')

	print ('step2_crosscorr(wavs_file, delta_time_w)')

	print ('step3_difftime(vel1, vel2, step)')

	print ('step4_envvalue(wavs_file)')

	print ('step5_stackall(wavs_file)')

	print ('step6_plotall(name_folder, wavs_file, ntotal_stat=None, normalized=False), source_loc=None')

