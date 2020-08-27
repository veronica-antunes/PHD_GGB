# -*- coding: utf-8 -*-

import sys, math, os
import numpy as np
import matplotlib.pyplot as plt
from glob import glob
from obspy import read, UTCDateTime
from math import log10
from scipy import interpolate


def assign_weight(tw_noise, tw_signal, rea_db, wav_db, plot_uncertainty = False, plot_fullwaveform = False,
	exclude_events = False, n_pick = 6, tw_interval = 0.7) :

	if not os.path.exists("figures"):
		os.makedirs("figures")

	sfile_path=sorted(glob(rea_db+'*/*/*/*'))
	wav_path=glob(wav_db)[0]
	count=0

	average_SNR=[]
	eve_nst_P=[]
	eve_nst_S=[]

	stat_snr = {} #creates a dictionary with data from each seismic station
	stat_rms = {}
	stat_weight = {}
	stat_Pdiff = {}

	n_P_pick=[]
	n_S_pick=[]

	weight_picking_list=[]
	pick_difference_list=[]
	uncertainty_list=[]
	weight_list=[]
	snr_list=[]

	for event in sfile_path :

		s_file_path = event

		sfile = open(s_file_path,'r').readlines()  #r+ allows read and write

		stat_pick_P=[]
		stat_pick_S=[]
		weight_picking_check=[]
		waveform=[]
		for lines in sfile:
			phase=[]
			if lines[-2] == '6' :
				waveform.append(lines.split(' ')[1])

			if lines[-2] == ' ':
				if lines[10:11] == 'P' :
					stat_pick_P.append(lines[1:6]+', '+
						lines[18:20]+', '+ lines[20:22] +', '+
						lines[22:28]+', '+ lines[71:76] + ', ' +
						lines[63:70]+', '+ lines[1:15])

				if lines[10:11] == 'S' :
					stat_pick_S.append(lines[1:6]+', '+
						lines[18:20]+', '+ lines[20:22] +', '+
						lines[22:28]+', '+ lines[71:76]+ ', ' +
						lines[63:70]+', '+ lines[1:15])

		if stat_pick_P : # or stat_pick_S
			wav=read(wav_path+'/*/*/'+waveform[0])

			if len(waveform) > 1 :
				for w in range(1, len(waveform)) :
					wav+=read(wav_path+'/*/*/'+waveform[w])

			wav.detrend(type='demean')
			wav.taper(max_percentage=0.005, type='hann', side='both') #max_lenght=None
			wav.filter('bandpass', freqmin=1, freqmax=15, corners=4) #, zerophase='False'

			y=wav[0].stats.starttime.year
			m=wav[0].stats.starttime.month
			d=wav[0].stats.starttime.day

			SNR_eve=[]

			for info in stat_pick_P :

				stat=info.split(',')[0].split(' ')[0]

				#get picking
				hh=int(info.split(',')[1])
				mm=int(info.split(',')[2])
				ss=float(info.split(',')[3])
				t=UTCDateTime(y,m,d,hh,mm,ss)

				pick='   '+'{:2}'.format(hh)+'{:2}'.format(mm)+' '+str(ss)
				info_pick=info.split(',')[6]

				# for P pick
				tr_z=wav.select(station=stat, channel='??Z')
				tr_z.trim(starttime=t-tw_noise, endtime=t+tw_signal)

				try :
					wav_pick=(tr_z[0].data)
				except :
					raise ValueError('ERROR: Missing waveform! check '+ event)

				samp_rate=tr_z[0].stats.sampling_rate
				picking_npts=(samp_rate*tw_noise)

				# compute SNR
				noise=wav_pick[0:int(picking_npts-tw_interval*samp_rate)]
				signal=wav_pick[int(picking_npts):-1]

				A_rms_noise=np.sqrt(np.mean(noise**2))
				A_rms_signal=np.sqrt(np.mean(signal**2))

				SNR=(A_rms_signal/A_rms_noise)**2

				#SNR_db=SNR
				SNR_db=10*log10(SNR)

				snr_list.append(SNR_db)

				# plot
				if plot_pickings == True :

					if not os.path.exists("figures/plot_picks"):
						os.makedirs("figures/plot_picks")

					fig=plt.figure()
					ax1=fig.add_subplot(1,1,1)
					plt.plot(wav_pick)
					ax1.axvline(x=picking_npts, color='r')
					ax1.set_title(str(stat) + '  SNR: ' + str(SNR_db), fontsize=16)
					plt.savefig('figures/plot_picks/'+str(SNR_db)+'_'+str(t)[0:-8]+
						'_'+stat+'.png')
					plt.close('all')

				if SNR_db > 0 : #and SNR_db < 20 :
					#choose big interval in order to be sure your picking will be
					#t1 the picking would be too early
					#t2 the picking would be too late
					#we will work around this interval				
					t1=int(picking_npts-(tw_interval*samp_rate))
					t2=int(picking_npts+(tw_interval*samp_rate))
					pick_interval=wav_pick[t1:t2]

					#get the max of the noise. our picking must be above the noise
					max_noise=max(abs(noise))*1.5
					min_signal=max(abs(noise))*1.05 

					#tL will be the position point which will be our latest 
					#pick of uncertainty.
					#For sure after this point everything will be signal
					#try because one of the values may be empty. this way we 
					#avoid the script to stop

					try:
						tL_positive=np.where(pick_interval > max_noise)[0][0]
						tA_positive=np.where(pick_interval > min_signal)[0][0]
					except:
						tL_positive=np.NaN
						tA_positive=np.NaN
					try:
						tL_negative=np.where(pick_interval < -max_noise)[0][0]
						tA_negative=np.where(pick_interval < -min_signal)[0][0]
					except:
						tL_negative=np.NaN
						tA_negative=np.NaN

					try: 
						tL_x=int(np.nanmin((tL_positive, tL_negative), axis=0))
						tA_x=int(np.nanmin((tA_positive, tA_negative), axis=0))

						tL_y=pick_interval[tL_x]
						tA_y=pick_interval[tA_x]

						#tE will be the point which tangent is equal to 0
						x_vector=np.arange(0,len(pick_interval),1)

						#get tangent of point in a recursive way
						tck = interpolate.splrep(x_vector, pick_interval)

						x0=tL_x-1
						slope_list=[]
						x0_list=[]
						slope_control=1
						slope=1

						while slope != 0 :
							y0 = interpolate.splev(x0,tck)
							dydx = interpolate.splev(x0,tck,der=1)
							tngnt = lambda x: dydx*x + (y0-dydx*x0)

							#get tangent slope
							tan_x1=20
							tan_y1=tngnt(tan_x1)
							tan_x2=100
							tan_y2=tngnt(tan_x2)
							slope=int((tan_y2 - tan_y1) / (tan_x2 - tan_x1))
							slope_factor=max(signal)*0.01

							#print (slope, slope_factor)
							if slope>-slope_factor and slope<slope_factor :
								slope=0
							slope_list.append(slope)
							x0_list.append(x0)

							x0=x0-1

						#print 'slope '+str(slope)
						myarray1 = np.asarray(slope_list)

						min_list= min(abs(myarray1))
						try:
							index=slope_list.index(min_list)
						except:
							index=slope_list.index(-min_list)

						#slope_list[index]
						tE_x=x0_list[index]

						#tE_x=x0+1 
						tE_y=pick_interval[tE_x] 

						tM=len(pick_interval)/2 #manual picking
						tA=((tL_x-tE_x)/2)+tE_x
						pick_difference=(tM-tA)/samp_rate

						tmin=min(tL_x, tE_x, tM, tA)
						tmax=max(tL_x, tE_x, tM, tA)

						uncertainty=(tmax-tmin)/samp_rate

						pick_difference_list.append(pick_difference)
						uncertainty_list.append(uncertainty)

						if SNR_db > 25 and uncertainty > 0.4 :
							print ('recheck picking for '+ stat + 
								' ' + event.split('/')[-1])

						#plot uncertainty for control
						if plot_uncertainty == True : # and uncertainty < 0.4:
							if not os.path.exists("figures/uncertainty"):
								os.makedirs("figures/uncertainty")

							fig=plt.figure()
							ax1=fig.add_subplot(1,1,1)
							ax1.set_title(str(t)[0:-8]+'\n'+str(stat) + 
								'  SNR: ' + str(round(SNR_db,2)) +
								' unc: ' + str(round(uncertainty,3)), 
								fontsize=16)
							plt.plot(x_vector, pick_interval) #label='slope:'

							plt.legend()

							if uncertainty != np.NaN :
								shift_factor=0
								plt.plot(tL_x+shift_factor, tL_y, "og")
								plt.plot(tE_x+shift_factor, tE_y, "og")
								ax1.axvline(x=tA+shift_factor, color='g')

							ax1.axvline(x=tM+shift_factor, color='r')

							plt.savefig('figures/uncertainty/'+str(t)[0:-8]+
								'_'+stat+'.png')
							plt.close('all')
					except:
						weight_picking = 5
						uncertainty = np.NaN
						#print event, stat, '1.5*max of noise 
						#condition not reached'

				# gives a weight taking into consideration the SNR_db. 
				# weight from 0 (good) to 5 (bad)
				# -- SNR --
				if SNR_db < 12 :
					weight_picking_factor = 2 
				if SNR_db >= 12 and SNR_db < 20:
					weight_picking_factor = 1
				if SNR_db >= 20:
					weight_picking_factor = 0

				# -- uncertaity --
				if uncertainty <= 0.05 :
					weight_picking=0 + weight_picking_factor

				if uncertainty > 0.05 and uncertainty <= 0.1 :
					weight_picking=1 + weight_picking_factor

				if uncertainty > 0.1 and uncertainty <= 0.2 :
					weight_picking=2 + weight_picking_factor

				if uncertainty > 0.2 and uncertainty <= 0.4 :
					weight_picking=3 + weight_picking_factor

				if uncertainty > 0.4 :
					weight_picking = 4 + weight_picking_factor

				if weight_picking > 5 :
					weight_picking = 5

				if SNR_db <= 0 : #and SNR_db < 20 :
					weight_picking = 5

				st_rms=float(info[30:35])

				if st_rms > 0.8 or st_rms < -0.8 :
					weight_picking = 5

				weight_list.append(weight_picking)

				#print SNR_db, uncertainty, weight_picking, stat
				#print event, weight_picking
				if plot_fullwaveform == True:

					if not os.path.exists("figures/weight"):
						os.makedirs("figures/weight")

					fig=plt.figure()
					ax1=fig.add_subplot(1,1,1)
					ax1.set_title(str(t)[0:-8]+'\n'+str(stat) + '  SNR: ' + 
						str(round(SNR_db,2)) + ' unc: ' + 
						str(round(uncertainty,3)) +
						 ' rms: ' + str(st_rms), fontsize=16)

					plt.plot(wav_pick,'k')
					try: 
						plt.plot(np.arange(t1,t2,1), pick_interval,'b', 
							label=weight_picking)
					except:
						plt.plot(np.arange(t1,t1+len(pick_interval),1), 
							pick_interval ,'b', label=weight_picking)

					plt.plot(np.arange(0,len(noise),1), noise,'r')
					plt.legend()

					shift_factor=(samp_rate*tw_noise)-(samp_rate*tw_interval)

					if math.isnan(uncertainty) == False :
						plt.plot(tL_x+shift_factor, tL_y, "og")
						plt.plot(tE_x+shift_factor, tE_y, "og")
						ax1.axvline(x=tA+shift_factor, color='g')

					ax1.axvline(x=tM+shift_factor, color='r')
					#print tM+shift_factor, samp_rate

					plt.savefig('figures/weight/'+str(weight_picking)+'-'+
						str(t)[0:-8]+'_'+stat+'.png')
					plt.close('all')

				#statistics: difference between original picking and repick (s)
				weight_picking_list.append(weight_picking_list)

				#Save the weight value in the sfile
				#new_info_pick = info_pick[0:14]+str(weight_picking)+info_pick[15::]
				new_info_pick = info_pick[0:14]+str(weight_picking)
				with open(s_file_path, "r+") as fin:
					with open(s_file_path, "r+") as fout:
						for line in fin:
							fout.write(line.replace(info_pick, new_info_pick))

				#Create dictionary associated with each station to save SNR
				if stat in stat_snr :
					stat_snr[stat].append(SNR_db)
				else:
					stat_snr[stat] = []
					stat_snr[stat].append(SNR_db)

				#Create dictionary associated with each station to save RMS
				if stat in stat_rms :
					try:
						stat_rms[stat].append(float(info.split(',')[5][0:7]))
					except:
						pass
				else:
					try:
						stat_rms[stat] = []
						stat_rms[stat].append(float(info.split(',')[5][0:7]))
					except:
						pass


				if weight_picking <= 4 :
					weight_picking_check.append(weight_picking)

	#	print weight_picking_check
		n_P_pick.append(len(stat_pick_P))
		n_S_pick.append(len(stat_pick_S))

	#	if len(stat_pick_P) < 6 :
		if exclude_events == True:
			if not os.path.exists("figures/excluded_events"):
				os.makedirs("figures/excluded_events")

			if len(weight_picking_check) < n_pick :
				os.rename(event, 'figures/excluded_events/'+event.split('/')[-1])

		plt.close('all')

		count=count+1

	for stat in stat_snr:
		if not os.path.exists("figures/stations"):
			os.makedirs("figures/stations")

		rms=stat_rms[stat]
		snr=stat_snr[stat]

		stat_npicks=len(stat_rms[stat])

		fig=plt.figure()
		fig.suptitle(str(stat) + ' pickings: ' + str(stat_npicks), fontsize=16)

		ax1 = fig.add_subplot(2,1,1)
		ax1.set_title('SNR   ', fontsize=16, x=0.95, y=1)
		nbins = np.arange(0, max(snr)+1, 1)

		snr = [x for x in snr if str(x) != 'nan'] #cleans all nan
		plt.hist(snr, bins=nbins, color='blue')

		ax2 = fig.add_subplot(2,1,2)
		ax2.set_title('RMS   ', fontsize=16, x=0.95, y=1)
		nbins = np.arange(min(rms)-1, max(rms)+1, 0.1)
		plt.hist(rms, bins=nbins, color='green')
		fig.subplots_adjust(hspace=.2)
		plt.savefig('figures/stations/'+stat+'.png')
		plt.close('all')

	plt.close('all')

	np.save('nP_event.npy', np.array(n_P_pick)) 
	np.save('SNR_pickings.npy', np.array(snr_list))
	np.save('uncertainty.npy', np.array(uncertainty_list))
	np.save('wheight.npy', np.array(weight_list))
	np.save('pick_difference.npy', np.array(pick_difference_list))
	 
	fs=16

	fig=plt.figure()
	fig.suptitle('Pickings Quality', fontsize=fs)

	ax1 = fig.add_subplot(5,1,1)
	ax1.set_title('nP/event', fontsize=fs, x=0.87, y=0.6)
	plt.hist(n_P_pick, bins=np.arange(0,50,5), color='lightblue')
	start, end = ax1.get_ylim()
	ax1.set_yticks(np.linspace(start, end+1, 3))
	ax1.set_xticks(ax1.get_xticks()[1:])
	ax1.tick_params(labelsize=fs)

	ax2 = fig.add_subplot(5,1,2)
	ax2.set_title('SNR [db]', fontsize=fs, x=0.87, y=0.6)
	snr_list = [x for x in snr_list if str(x) != 'nan'] 

	plt.hist(snr_list, bins=np.arange(int(min(snr_list)), int(max(snr_list)+1), 5), color='forestgreen')
	start, end = ax2.get_ylim()
	ax2.set_yticks(np.linspace(start, end+20, 3))
	ax2.set_xticks(ax2.get_xticks()[1:])
	ax2.tick_params(labelsize=fs)

	ax3 = fig.add_subplot(5,1,3)
	ax3.set_title('Unc. [s]', fontsize=fs, x=0.87, y=0.6)
	plt.hist(uncertainty_list, bins=np.arange(0,1,0.05), color='orange')
	start, end = ax3.get_ylim()
	ax3.set_yticks(np.linspace(start, end+10, 3))
	ax3.set_xticks(ax3.get_xticks()[1:])
	ax3.tick_params(labelsize=fs)

	ax4 = fig.add_subplot(5,1,4)
	ax4.set_title('Wheight', fontsize=fs, x=0.87, y=0.6)
	plt.hist(weight_list, bins=(np.arange(9) - 0.5), color='c')
	start, end = ax4.get_ylim()
	ax4.set_yticks(np.linspace(start, end+10, 3))
	ax4.set_xticks(ax4.get_xticks()[1:])
	ax4.tick_params(labelsize=fs)

	ax5 = fig.add_subplot(5,1,5)
	ax5.set_title('tM-tA [s]', fontsize=fs, x=0.87, y=0.6)
	plt.hist(pick_difference_list, bins=np.arange(-1,1,0.05), color='salmon')
	start, end = ax5.get_ylim()
	ax5.set_yticks(np.linspace(start, end+10, 3))
	ax5.set_xticks(ax5.get_xticks()[1:])
	ax5.tick_params(labelsize=fs)

	fig.subplots_adjust(hspace=.3)
	plt.savefig('figures/quality_of_pickings2.png')

