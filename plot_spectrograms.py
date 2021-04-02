##!/usr/bin/env python2
# -*- coding: utf-8 -*-

from obspy.core import UTCDateTime, read, Trace, Stream
from obspy.imaging.spectrogram import spectrogram
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import *
from pylab import figure, axes, pie, title, show
from matplotlib.dates import date2num, num2date, DateFormatter
from matplotlib import dates
import datetime, matplotlib
from obspy.io.xseed import Parser
from obspy.signal.tf_misfit import cwt
from glob import glob
from matplotlib.ticker import AutoMinorLocator


def wav_spectrogram(waveforms, fmi1, fma1, fmi2, fma2, clip, yscale='linear') :

	files=glob(waveforms)

	matplotlib.rcParams.update({'font.size': 16})

	for f in range(len(files)) :
		x=read(files[f])

		t0=x[0].stats.starttime
		t1=x[0].stats.endtime

		zx=x.select(component='Z')

		zx = zx.slice(starttime=t0, endtime=t1)

		for i in range(len(zx)) :
			b=zx[i]
			a=b.copy()
			station=a.stats.station

			ac=a.copy()

			ac.detrend(type='demean')
			ac.taper(max_percentage=0.005, type='hann', side='both') #max_lenght=None
			ac.filter('bandpass', freqmin=fmi1, freqmax=fma1, corners=2, zerophase='True')

			af=a.copy()
			af.detrend(type='demean')
			af.taper(max_percentage=0.005, type='hann', side='both') #max_lenght=None
			af.filter('bandpass', freqmin=fmi2, freqmax=fma2, corners=2, zerophase='True')

			t0=a.stats.starttime
			t1=a.stats.endtime
			stat=a.stats.station

			#plot original filtered and spectrogram
			fig2 = plt.figure(figsize=(10,8)) #[20,18]15.15
			date=t0.strftime(format='%Y-%m-%d %H:%M:%S')

			ax1 = fig2.add_axes([0.1, 0.79, 0.8, 0.12]) #[left bottom width height]
			ax4 = fig2.add_axes([0.1, 0.67, 0.8, 0.12], sharex=ax1) #[left bottom width height]
			ax2 = fig2.add_axes([0.1, 0.04, 0.8, 0.52]) #, sharex=ax1)
			ax3 = fig2.add_axes([0.92, 0.04, 0.03, 0.52])
			#ax1.set_title(str(stat) + '\n' + date + '\n', fontsize=16, fontweight='bold')

			#insert dates
			start =date2num(t0.datetime)
			end = date2num(t1.datetime)

			#make time vector
			t = np.linspace(start, end, a.stats.npts)

			#plot waveform original (top subfigure)
			ax1.plot_date(t, ac.data, 'g', label=station+': '+
				str(fmi1) + '-' + str(fma1) + 'Hz')

			ax1.xaxis.set_major_locator(dates.MinuteLocator(interval=1))

			ax1.xaxis.set_major_formatter( dates.DateFormatter('%H:%M:%S') )
			ax1.xaxis.set_minor_locator(dates.SecondLocator(interval=2))

			ax1.set_ylabel("Amplitude", fontsize=18)

			plt.setp(ax1.get_xticklabels(),visible=0) #hide tick text but shows tick marks
			ax1.legend(loc="upper right", prop={'size': 16, 'weight':'bold'}, 
				framealpha=1.0, bbox_to_anchor=(1.0, 1.15), frameon=False)

			#max_Y=max(abs(ac.data))
			#ax1.set_yticks(np.linspace(-max_Y, max_Y, 3))

			#plot filtered wav (2nd figure)
			ax4.plot_date(t, af.data, 'r', label=station+ ': ' +
				str(fmi2) + '-' + str(fma2) + 'Hz' )
			ax1.locator_params(nbins=3)
			ax4.locator_params(nbins=3)
			ax4.set_xlabel("time [HH:MM:SS]", fontsize=18)
			ax4.yaxis.set_ticks_position('right')
			ax4.yaxis.set_label_position('right')
			ax4.set_ylabel("Amplitude", fontsize=18)
			ax4.legend(loc="upper right", prop={'size': 16, 'weight':'bold'}, 
				framealpha=1.0, bbox_to_anchor=(1.0, 1.15), frameon=False) #increase size
			ax4.set_xticks(np.linspace(start, end, 6)[:-1])
			#max_Y=max(abs(af.data))
			#ax4.set_yticks(np.linspace(-max_Y, max_Y, 3))
			###### plot spectrogram (bottom subfigure)
	
			name=str(t0)[0:10] + '_'   + stat

			fig2 = a.spectrogram(log=True, show=False ,  axes=ax2, clip=[0,clip],
				cmap=plt.cm.jet)      		

			mappable = ax2.collections[0] #for logaritmic scale on spectrogram function
			#mappable = ax2.images[0]      #for linear scale on spectrogram function
			plt.colorbar(mappable=mappable, cax=ax3)

			ax2.set_ylabel("Frequency [Hz]", fontsize=18)
			
			ax2.set_yscale(yscale) #linear #log
			ax2.tick_params(axis='x', color='white')
			ax2.tick_params(axis='y',which='minor',color='white')
			ax2.set_xlabel('time after ' +str(date)+' [s]', fontsize=18)
			ax2.set_title("All Spectrogram" , fontsize=18)

			#subplots_adjust(top=10,bottom=5)
			plt.savefig( str(t0)[0:16] + '_' + str(stat) + '_spectrogram_'+ yscale +'.png',
				bbox_inches='tight')

			plt.clf()
			plt.close('all')



def all_spectrogram(waveform, maxv, fx=5, fy=4, fs=11, outfile='spectrogram.png') :

	from matplotlib.pyplot import figure

	figure(num=None, figsize=(fx, fy), dpi=300)

	matplotlib.rcParams.update({'font.size': fs})
	number=0.71

	zx=read(waveform)
	zx=zx.select(channel='*Z')
	#zx=zx.trim(starttime=zx[0].stats.starttime+(19*60),endtime=zx[0].stats.starttime+(23*60))

	t0=zx[0].stats.starttime

	fig=plt.figure()
	fig.subplots_adjust(hspace=.25)

	for n in range(len(zx)):
		stat=zx[n].stats.station

		#print 'Doing envelope for: ', stat

		ax = fig.add_subplot(len(zx),1,n+1)
		ax = zx[n].spectrogram(show=False,
			axes=ax,  cmap=plt.cm.jet,
			v0=0, v1=maxv)
		#ax.set_xlim([ti,tf])

		ax.set_title(stat, fontsize=fs, x=0.93, y=0.70, weight='bold', color='white')
		ax.set_yticks(np.arange(0, zx[n].stats.sampling_rate/2,zx[n].stats.sampling_rate/5))

		if n < len(zx)-1 :
			ax.axes.get_xaxis().set_ticklabels([]) #to remove the labels on top figures
		if n ==  len(zx)-1 :
			ax.set_xlabel('time [s] after '+str(t0)[0:-8], fontsize=fs)

		subplots_adjust(left=0.1, bottom=None, right=0.87, top=0.95, wspace=None, hspace=0.05)
		number=number-0.205

	ax2 = fig.add_axes([0.89, 0.18, 0.03, 0.65])
	mappable = ax.images[0] 
	cb=plt.colorbar(mappable=mappable, cax=ax2, orientation='vertical')
			#labels=np.round(np.arange(0, maxv,maxv/3.0),1)
			#cb.set_ticks(labels)
	cb.ax.set_ylabel('Amplitudes')
	labels=np.linspace(0, maxv,7)
	cb.set_ticks(labels)

	fig.text(0.00, 0.5, 'Frequency [Hz]', va='center', rotation='vertical', fontsize=fs)
	time=str(t0)[0:-8].split(':')
	name=time[0]+'-'+time[1]+'-'+time[2]

	if maxv < 0.1 :
		cb.formatter.set_powerlimits((0, 0))
		cb.ax.yaxis.set_offset_position('left')                         
		cb.update_ticks()
		labels=np.arange(0, maxv,maxv/5)
		cb.set_ticks(labels)

	fig.set_figheight(fy)
	fig.set_figwidth(fx)

	plt.savefig( outfile, bbox_inches='tight')
	plt.close('all')
