# -*- coding: utf-8 -*-

import sys, os, math, re
import numpy as np
import matplotlib.pyplot as plt
from glob import glob
from obspy import read, UTCDateTime, Stream
from obspy.io.xseed import Parser
from math import log10
from obspy.geodetics import gps2dist_azimuth
from datetime import datetime
import matplotlib.ticker as mtick

def determ_magn(s_file, t_window, wav_path, dataless_path, distance_type='hypocentral', 
	plot_wavs = False, print_stat_mag = False ) :
	#distance_type='hypocentral'#epicentral or hypocentral
	outfile_name=s_file.split('.')[0]+'_'+str(t_window)+'s_M.out'

	# parameters for M formula: ML=a*log10(ampl)+b*log10(dist)+c*Dist+d. Two formulas:
	distance=15
	distance2=60

	#dist <= distance km :
	A, B, C, D = [1, 1.5, 0, 0.45-log10(2.8)]
	# distance < dist < distance2 km :
	E, F, G, H = [1, 0, 0.0180, 1.77+0.1]
	#dist > distance2 km :
	I, J, K, L = [1, 0, 0.0038, 2.62+0.1]

	dataless=glob(dataless_path)

	#dataless_ug='/home/veronica/ownCloud/PSD/dataless_all.dseed'
	#dataless_permanent='/home/veronica/ownCloud/PSD/dataless_pstations_good.dseed'
	#True or False to plot the amplitudes for E and N channels
	#plot_wavs = False
	#saves the magnitudes for each station
	#print_stat_mag = False

	#wood-anderson stantard seismograph instrument response according to SED
	paz_wa = {'sensitivity': 2800, 'zeros': [0j], 'gain': 1,
		  'poles': [-6.2832 - 4.7124j, -6.2832 + 4.7124j]}

	#parser_ug = Parser(dataless_ug)
	#parser_pstat = Parser(dataless_permanent)

	#sfile_path=glob(rea_path)

	#wav_path=glob('/home/veronica/SEISAN/WAV/'+wav_database+'*')
	#wav_path=glob(wavs_db_path+)
	#count=0
	magnitude=[]

	g=open(outfile_name,'w')

	f=open(s_file,'r').read()	
	all_sfile=re.split('                                  ' +
		'                                              \n', f)

	for sfile in all_sfile :

		station_pickings=[]
		waveforms=[]
		for lines in sfile.split('\n'):

			if len(lines) > 0 and lines[-1] == '6' :
				waveforms.append(lines.split(' ')[1])

			if len(lines) > 0 and lines[-1] == ' ':
				if lines[10:11] == 'S' :
					if  lines[14:15]==' ' or lines[14:15]<4:
						station_pickings.append(lines[1:6]+', '+
							lines[18:20]+', '+ lines[20:22] +', '+
							lines[22:28]+', '+ lines[71:76])

			if len(lines) > 0 and lines[-1] == '1' and lines[23:30] != '       ' :
				event=UTCDateTime(year=int(lines[1:5]),
					month=int(lines[6:8]),
					day=int(lines[8:10]),
					hour=int(lines[11:13]),
					minute=int(lines[13:15]))+float(lines[16:20])

				type_eve=lines[21:23]
				lat=float(lines[23:30])
				lon=float(lines[31:38])
				dep=float(lines[38:43])
				rms=float(lines[51:55])

		if station_pickings:
			try:
				try:
					wav=read(wav_path+'/'+
						str(event.year)+'/'+
						"%02d" %event.month+
						'/'+waveforms[0].split('.')[0]
						+'*.mseed')
					#if len(wav) > 0 :
					#	for w in range(1,len(wav)):
					#		wav+=read(wav_path+'/'+
					#			str(event.year)+'/'+
					#			"%02d" %event.month+
					#			'/'+waveforms[w])
				except:
					wav=read(wav_path+'/'+waveforms[0].split('.')[0]
						+'*.mseed')
					#if len(wav) > 0 :
					#	for w in range(1,len(wav)):
					#		wav+=read(wav_path+'/'+waveforms[w])

				wav.detrend(type='demean')
				wav.taper(max_percentage=0.005, type='hann', side='both')

				y=wav[0].stats.starttime.year
				m=wav[0].stats.starttime.month
				d=wav[0].stats.starttime.day

			#	print 'Magnitude calculation for ', waveform		
				stat_mag=[]
				stat_info=[]
			except:
				print('no waveform found!')
				wav=read()
				wav.clear()
				pass

			if wav :

				for info in station_pickings :

					stat=info.split(',')[0].split(' ')[0]

					hh=int(info.split(',')[1])
					mm=int(info.split(',')[2])
					ss=float(info.split(',')[3])

					#print info
					t=UTCDateTime(y,m,d,hh,mm)+ss

					tr_e=wav.select(station=stat, channel='??E')
					tr_n=wav.select(station=stat, channel='??N')

					if not (tr_e and tr_n) :
						tr_e=wav.select(station=stat, channel='??3')
						tr_n=wav.select(station=stat, channel='??2')

					for each_dataless in dataless:
						try:
							parser=Parser(each_dataless)

							paz_stat_e = parser.get_paz(
								tr_e[0].id, datetime=t)
							paz_stat_n = parser.get_paz(
								tr_e[0].id, datetime=t)

							tr_e[0].stats.latitude=parser.get_coordinates(
								tr_e[0].id, datetime=t).get('latitude')
							tr_e[0].stats.longitude=parser.get_coordinates(
								tr_e[0].id, datetime=t).get('longitude')

							tr_n[0].stats.latitude=parser.get_coordinates(
								tr_n[0].id,datetime=t).get('latitude')
							tr_n[0].stats.longitude=parser.get_coordinates(
								tr_n[0].id,datetime=t).get('longitude')
						except:
							pass

					if not (tr_e and tr_n):
						print ('Missing waveform for '+ stat)

					if tr_n and not tr_n[0].stats.longitude :
						print  'Missing instrument response for ' + stat
						paz_stat_e = False
						paz_stat_n = False
						#print 'Unknown problems with station ' + stat

					if paz_stat_e and paz_stat_n and (tr_e or tr_n) :
						try :
							sta_lat = tr_e[0].stats.latitude
							sta_lon = tr_e[0].stats.longitude
						except:
							sta_lat = tr_n[0].stats.latitude
							sta_lon = tr_n[0].stats.longitude

						if len(tr_e) > 1 :
							#print 'Merging traces for '+ stat
							#tr_e.plot()
							tr_e.merge(fill_value='interpolate')
							if max(abs(tr_e[0].data)) == 0 :
								tr_e=Stream(traces=[wav.select(
									station=stat, 
									channel='??E')[1]])

						tr_e.simulate(paz_remove=paz_stat_e, 
							paz_simulate=paz_wa,
							water_level=10)
						tr_e.trim(starttime=t-1, endtime=t+t_window)

						if tr_e:
							amplt_e=max(abs(tr_e[0].data))
							#sta_lat = tr_e[0].stats.latitude
							#sta_lon = tr_e[0].stats.longitude

						if not tr_e:
							amplt_e=np.nan

						if len(tr_n) > 1:
							#print 'Merging traces for '+ stat
							#tr_n.plot()
							tr_n.merge(fill_value='interpolate')
							#tr_n.plot()
							if max(abs(tr_n[0].data)) == 0 :
								tr_n=Stream(traces=[wav.select(
									station=stat, 
									channel='??N')[1]])

						tr_n.simulate(paz_remove=paz_stat_n, 
							paz_simulate=paz_wa,
							water_level=10)
						tr_n.trim(starttime=t-1, endtime=t+t_window)

						if tr_n:
							amplt_n=max(abs(tr_n[0].data))
							#sta_lat = tr_n[0].stats.latitude
							#sta_lon = tr_n[0].stats.longitude
						if not tr_n:
							amplt_n=np.nan

						#ampl=max(amplt_e, amplt_n)
						ampl=np.nanmax([amplt_e,amplt_n])
						#dist=float(info.split(',')[4]) 
						#distance taken from s-file
						#calculation of distance between station and event

						event_lat = lat
						event_lon = lon
						event_dep = dep
						dist, az, baz = gps2dist_azimuth(event_lat, event_lon, 
							sta_lat, sta_lon)
						dist = dist / 1000
						#print 'epi_dist: ', epi_dist, ' dist: ', dist

						if distance_type=='epicentral':
							dist=dist

						elif distance_type=='hypocentral':
							dist=math.sqrt((float(dist)**2)+
								(float(event_dep)**2))

						else:
							raise ValueError("ERROR: Wrong type "+
								"of distance! "+ distance_type +
								"' selected. Needs to be"+ 
								"'epicentral' or 'hypocentral'")

						if dist <= distance :
							#ml = (A*log10(ampl * 1000) + B*log10(dist) +
							#	C*dist + D)	
							ml = np.NaN
							print ('Distance lower than 15km for station '+ 
								stat + ' ' +str(round(dist,1)) )

						if dist > distance and dist <= distance2 :
							ml = (E*log10(ampl * 1000) + F*log10(dist) + 
								G*dist + H)

						if dist > distance2:
							ml = (I*log10(ampl * 1000) + J*log10(dist) +
								K*dist + L)

						## print(ml), stat
						stat_mag.append(ml)

						if print_stat_mag == True :
							#f.write( str(stat) + '  ' + str(ml) + '  ' +
							#	str(dist) + '  ' + str(ampl*1000) +'\n')
							stat_info.append( str(stat) + ' ' + str(ml) +
								' '+ str(dist) + ' ' + str(ampl*1000) )
							#separate list by space and convert to array
							#a=np.array([x.split(' ') for x in stat_info])
							#b=np.asarray(a)

						#plot for control wavs
						if plot_wavs == True :
							if not os.path.exists("figures"):
								os.makedirs("figures")

							fig=plt.figure()
							fig.suptitle('Amplitude (A) mm - time window '+ 
								str(t_window+1) + 's' +
								'\n' + str(event)[0:16] + '  ' + stat+ 
								' - ML ' + str(round(ml,1)) +' - D '
								+ str(round(dist,2)))

							ax1 = fig.add_subplot(4,1,1)
							ax1.plot(tr_e[0].data, 'g')	
							ax1.set_title(' A     E', x=0.92, y=0.6)
							ax1.yaxis.set_ticks(np.linspace(
								min(tr_e[0].data), 
								max(tr_e[0].data), 3))
							ax1.axvline(amplt_e, color='b')
							ax1.set_xlim(0,len(tr_e[0].data))
							ax1.axes.get_xaxis().set_ticklabels([])
							abc=np.where(abs(tr_e[0].data) == amplt_e)
							ax1.axvline(float(abc[0]), color='b', alpha=0.5)
							ax1.yaxis.set_major_formatter(
								mtick.FormatStrFormatter('%.2e'))

							ax2 = fig.add_subplot(4,1,2)
							ax2.plot(abs(tr_e[0].data), 'g', 
								label='max '+'{:.2E}'.format(amplt_e))
							ax2.legend(loc="upper left", fontsize=11)
							ax2.set_title('|A|    E', x=0.92, y=0.6)
							ax2.yaxis.set_ticks(np.linspace(0, ampl, 3))
							ax2.set_xlim(0,len(tr_e[0].data))
							ax2.axvline(float(abc[0]), color='b', alpha=0.5)
								#, scatter options
								#marker='o',  edgecolor='black',
								#linewidth='0.8', s=80)
							ax2.axes.get_xaxis().set_ticklabels([])
							ax2.set_ylim(0,ampl)
							ax2.yaxis.set_major_formatter(
								mtick.FormatStrFormatter('%.2e'))

							ax3 = fig.add_subplot(4,1,3)
							ax3.plot(tr_n[0].data, 'r')	
							ax3.set_title(' A     N', x=0.92, y=0.6)
							ax3.yaxis.set_ticks(np.linspace(
								min(tr_n[0].data), 
								max(tr_n[0].data), 3))
							ax3.set_xlim(0,len(tr_n[0].data))
							ax3.axes.get_xaxis().set_ticklabels([])
							abc=np.where(abs(tr_n[0].data) == amplt_n)
							ax3.axvline(float(abc[0]), color='b', alpha=0.5)
							ax3.yaxis.set_major_formatter(
								mtick.FormatStrFormatter('%.2e'))

							ax4 = fig.add_subplot(4,1,4)
							ax4.plot(abs(tr_n[0].data), 'r', 
								label='max '+'{:.2E}'.format(amplt_n))
							ax4.legend(loc="upper left", fontsize=11)
							ax4.set_title('|A|    N', x=0.92, y=0.6)
							ax4.yaxis.set_ticks(np.linspace(0, ampl, 3))
							ax4.set_xlim(0,len(tr_n[0].data))
							ax4.set_xlabel('npts')
							ax4.axvline(float(abc[0]), color='b', alpha=0.5)
							ax4.set_ylim(0,ampl)
							ax4.yaxis.set_major_formatter(
								mtick.FormatStrFormatter('%.2e'))

							plt.savefig('figures/'+str(event)[0:16] + '_' + 
								str(round(dist,2)) + '_' + stat+'.png',
								bbox_inches='tight')
							plt.close()

				eve_mag=round(np.nanmedian(stat_mag),1)	
				magnitude.append(eve_mag)

				print ( str(event) + ' ' + str(type_eve) + ' ' + str(eve_mag) + ' ' 
					+ str(len(stat_mag)) + ' ' +str(lat) + ' '+ str(lon) 
					+ ' ' + str(dep) )

				g.write( str(event) + ' ' + str(type_eve) + ' ' + str(eve_mag) + ' ' 
					+ '{:02.0f}'.format(len(stat_mag)) + ' ' + '{:01.3f}'.format(lat) 
					+ ' ' + '{:01.3f}'.format(lon) + ' '+ '{:02.2f}'.format(dep)
					+' ' + str(rms) + ' \n' )

			#		f.write( ' ' + '{:%Y   %m   %d   %H   %M   %S.%f}'.format(event.
			#			datetime)[0:-4]+'     '+'{:07.3f}'.format(lat[count])+
			#			'     '+ '{:07.3f}'.format(lon[count]) + 
			#			'     ' + '{:07.3f}'.format
			#			(dep[count]) + '   ' + str(eve_mag)+'\n' )

				if print_stat_mag :
					stat_info=sorted(stat_info, key = lambda x: x.split(' ')[1])
					for stat_info_lines in stat_info: 
						f.write(stat_info_lines+'\n')
					g.write('--------------------------------------------------------' 
						+ '\n')

		#
		#count=count+1
	g.close()
