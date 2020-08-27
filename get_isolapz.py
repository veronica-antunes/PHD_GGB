#!/usr/bin/env python2
# -*- coding: utf-8 -*-
# Verónica Antunes 06/2019 last modified 08/2020
# Automatically writes input files for ISOLA:
# 	- Converts waveforms to sac
# 	- Gets isola-pz file from dataless or online depository
#	- Gives precendence to dataless
#	- Takes instrument response of the event date
# veronica.antunes@unige.ch

from obspy.io.xseed import Parser
from obspy import UTCDateTime, read
from glob import glob
import os, shutil, sys
from decimal import Decimal

def get_isolapz(waveforms, dataless = None, download=False, server='GFZ', station_type='B'):
	waveforms=glob(waveforms)

	for waveform in waveforms :
		wav=read(waveform)
		otime=wav[0].stats.starttime
		try :
			path_name=waveform.split('.')[0]

		except :
			raise ValueError('ERROR: No format detected! Write waveform as '+
				'"waveform_name.format"')
		try:
			os.makedirs(path_name)
			os.makedirs(path_name+'/pzfiles')
			os.makedirs(path_name+'/SAC')
		except:
			pass

		for w in wav:

			### converts and writes waveform to SAC for isola waveform input ###
			sta_id=w.id
			station=w.stats.station
			channel=w.stats.channel
			network=w.stats.network
			w.write(path_name+'/SAC/'+sta_id+'.SAC', format='SAC')	

			### converts instrument response and write pz-file for isola instrument input ###
			if dataless :
				for d in dataless :
					try:
						parser_dataless=Parser(d)
						a=p.get_coordinates(sta_id, datetime=otime)
					except:
						pass

			if download :
				try:
					from obspy.clients.arclink.client import Client
					client = Client(user=server)

					parser_download=client.get_paz(network, station, '', 
						channel, otime)
					a=p.get_coordinates(sta_id, datetime=otime)
				except:
					pass

			if parser_dataless :
				seism_gain=parser_dataless.get_paz(sta_id, 
					datetime=otime)['seismometer_gain']
				digi_gain=parser_dataless.get_paz(sta_id, 
					datetime=otime)['digitizer_gain']
				sensitivity=parser_dataless.get_paz(sta_id, 
					datetime=otime)['sensitivity']
				count=1/sensitivity
				A0=parser_dataless.get_paz(sta_id, 
					datetime=otime)['gain']
				npoles=len(parser_dataless.get_paz(sta_id, 
					datetime=otime)['poles'])
				nzeros=len(parser_dataless.get_paz(sta_id, 
					datetime=otime)['zeros'])
				zeros=(parser_dataless.get_paz(sta_id, 
					datetime=otime)['zeros'])
				poles=(parser_dataless.get_paz(sta_id, 
					datetime=otime)['poles'])

			elif parser_download :
				sensitivity=parser_download['sensitivity']
				count=1/sensitivity
				try:
					seism_gain=float(parser_download['name'].split('=')[1])
					digi_gain=sensitivity/seism_gain
				except:
					seism_gain = 'unknown'
					digi_gain = 'unknown'
				A0=parser_download['gain']
				npoles=len(parser_download['poles'])
				nzeros=len(parser_download['zeros'])
				zeros=(parser_download['zeros'])
				poles=(parser_download['poles'])

			if sensitivity:
				f = open(path_name+'/pzfiles/'+station + station_type + 
					channel[1::]+'.pz','w')
				f.write('A0\n')
				f.write(str(A0)+'\n')
				f.write('count-->m/sec\n')
				f.write(str(count)+'\n')
				f.write('zeros\n')
				f.write(str(nzeros)+'\n')

				for z in zeros:
					f.write('%.6e' % Decimal(z.real) +'     ' +
						'%.6e' % Decimal(z.imag)+'\n')
				f.write('poles\n')
				f.write(str(npoles)+'\n')

				for p in poles:
					f.write('%.6e' % Decimal(p.real) +'    ' +
						'%.6e' % Decimal(p.imag)+'\n')
				f.write('Info:  '+otime.strftime('%d-%m-%Y')+'  '+station+
					'  Digi sens  ' + 
					str(digi_gain) + 
					'  Seism sens   '+
					str(seism_gain))
				f.close()

			else:
				print ('WARNING: missing instrument response for station '+ station)
				continue


