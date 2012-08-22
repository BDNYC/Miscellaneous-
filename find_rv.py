#!/usr/bin/env python
# encoding: utf-8
"""
find_rv.py

Created by Vivienne Baldassare on 2012-03-30.
Copyright (c) 2012 __MyCompanyName__. All rights reserved.

Description: 
	This code finds the radial velocity of a target when supplied with data for the target and data for a standard object
	whose radial velocity is known.

Usage:
	Note: Data used should already be corrected for heliocentric velocity.
	
	Inputs:
		wv_obj, fx_obj, and sig_obj are arrays containing data for the the wavelength, flux, and flux uncertainty of the target.
		wv_std, fx_std, and sig_std are arrays containing data for the the wavelength, flux, and flux uncertainty of the standard.
		rv_std is the radial velocity of the standard.
		rv_std_err is the uncertainty in the radial velocity of the standard.
		obj_name and std_name are strings containing the names of the target and standard.  These are used in the production of plots.

	Example:
		>>> import find_rv
		>>> find_rv.radial_velocity(wv_obj,fx_obj,sig_obj,wv_std,fx_std,sig_std,rv_std,rv_std_err,obj_name,std_name)

"""

from array import array
import asciitable
import math
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import numpy
import pyfits
import pylab
import random
import scipy
from scipy.stats import norm
from scipy import constants
from scipy import interpolate
import scipy.optimize as op
import scipy.ndimage
import sys
import pdb
import pyspeckit

def radial_velocity(wv_obj,fx_obj,sig_obj,wv_std,fx_std,sig_std,rv_std,rv_std_err,obj_name,std_name):

# Find where standard and object overlap ---------------

	wv_min = max([min(wv_std),min(wv_obj)])
	wv_max = min([max(wv_std),max(wv_obj)])


# Overlap wv_obj and wv_std arrays.  Where they do not overlap, flux is set to 1 ---------------- 
	length=len(wv_obj)
	# For standard
	a_param = wv_std > wv_min
	b_param = wv_std < wv_max  
	i=0
	j=0
	while i < length:
		if a_param[i] == False:
			fx_std[i]=1
			i=i+1
		else:
			i=i+1
	while j < length:
		if b_param[j] == False:
			fx_std[j]=1
			j=j+1
		else:
			j=j+1	
	n_pix_std = len(wv_std)
	# For object
	a_param = wv_obj > wv_min
	b_param = wv_obj < wv_max
	i=0
	j=0
	while i < length:
		if a_param[i] == False:
			fx_obj[i]=1
			i=i+1
		else:
			i=i+1
	while j< length: 
		if b_param[j] == False:
			fx_obj[j]=1
			j=j+1
		else:
			j=j+1
	n_pix_obj = len(wv_obj)



# Creates ln standard wavelength array ---------------------------------
	min_wv_std = min(wv_std)
	max_wv_std = max(wv_std)
	acoef_std = (n_pix_std -1)/(math.log(max_wv_std) - math.log(min_wv_std))
	bcoef_std = (n_pix_std) - (acoef_std * math.log(max_wv_std))

	arr = numpy.arange(n_pix_std)+1
	wv_ln_std = numpy.exp((arr - bcoef_std)/acoef_std)
	
	
# Interpolate data onto same ln wavelength scale -------------------------------

	fx_interp_std = numpy.interp(wv_ln_std, wv_std, fx_std) 
	fx_interp_obj = numpy.interp(wv_ln_std, wv_obj, fx_obj)


# Rebin Data ----------------------------

	wv_arr_std=numpy.asarray(wv_ln_std,dtype=float)
	fx_arr_obj=numpy.asarray(fx_interp_obj,dtype=float)
	fx_arr_std=numpy.asarray(fx_interp_std,dtype=float)
	sig_arr_obj=numpy.asarray(sig_obj,dtype=float)
	sig_arr_std=numpy.asarray(sig_std,dtype=float)
	
	wv_ln_rebin_std=scipy.ndimage.interpolation.zoom(wv_arr_std,10)		#data rebinned by factor of 10
	fx_rebin_obj=scipy.ndimage.interpolation.zoom(fx_arr_obj,10)
	fx_rebin_std=scipy.ndimage.interpolation.zoom(fx_arr_std,10)
	sig_rebin_obj=scipy.ndimage.interpolation.zoom(sig_arr_obj,10)
	sig_rebin_std=scipy.ndimage.interpolation.zoom(sig_arr_std,10)

	
# Plot object and standard so you can clearly see that shift exists --------------------------------
	plt.figure(1)
	plt.plot(wv_ln_rebin_std,fx_rebin_obj,'r')
	plt.plot(wv_ln_rebin_std,fx_rebin_std,'b')
	v=[1.545,1.570,0,2]
	plt.axis(v)	
	

# Cross correlation loop -------------------------------- 
	pix_shift=[]		#initialize array for pixel shift values
	l = 0

	for l in range(0,500):
	
	# GETTING ARRAYS READY FOR CROSS CORRELATION
		
		# Randomize noise:
		# create gaussian distribution of random numbers b/t 1 and -1, multiply err by numbers, add numbers to flux
		fx_temp_obj=[None]*len(fx_rebin_obj)
		fx_temp_std=[None]*len(fx_rebin_std)
		rand_dist=[None]*len(fx_rebin_std)
		rand_dist2=[None]*len(fx_rebin_std)
		rand_dist=[random.gauss(0,.34) for i in rand_dist]
		rand_dist2=[random.gauss(0,.34) for i in rand_dist2]
		rand_dist=numpy.array(rand_dist)
		rand_dist2=numpy.array(rand_dist2)
		fx_temp_obj=numpy.array(fx_temp_obj)
		fx_temp_std=numpy.array(fx_temp_std)
		fx_temp_obj = fx_rebin_obj + (sig_rebin_obj * rand_dist)
		fx_temp_std = fx_rebin_std + (sig_rebin_std * rand_dist2)
		
		# Find std dev and mean of flux data
		mean_obj=fx_temp_obj.mean()
		mean_std=fx_temp_std.mean()
		stddev_obj=fx_temp_obj.std()
		stddev_std=fx_temp_std.std()
		
		# Regularize data (subtract mean, divide by std dev)
		fx_reg_temp_obj = fx_temp_obj-mean_obj
		fx_reg_temp_obj = fx_reg_temp_obj/stddev_obj
		fx_reg_temp_std = fx_temp_std-mean_std
		fx_reg_temp_std = fx_reg_temp_std/stddev_std
	

	# CROSS CORRELATION 

		# what you're comparing: obj flux and std flux
		y1=fx_reg_temp_obj
		y2=fx_reg_temp_std
		
		# compute the cross-correlation between y1 and y2
		ycorr = scipy.correlate(y1, y2, mode='full')
		ycorr1=ycorr[9750:10750]	#isolate section of array with gaussian

		length=len(ycorr1)
		xcorr=range(length)	#create x axis values
		#print xcorr		
		
		def chi2(p):	#define gaussian function for fitting
			sig2=p[2] ** 2
			m = (p[0] * numpy.exp(-0.5 * (xcorr - p[1]) ** 2 / sig2)) + p[3]
			return (ycorr1 - m)	
		
		amp = 6000	# guess some values
		mean = 300
		sig = 100
		sky = 1000	
		
		amp, mean, sig, sky = op.leastsq(chi2, [amp, mean, sig, sky])[0]
		
		
		#print 'amp=',amp,' mu=',mean, ' sig=',sig, ' sky=',sky
		
		print_num=l%50		#prints data every 100 fits
		if print_num == 0:
			print 'amp=',amp,' mu=',mean, ' sig=',sig, ' sky=',sky
		
		mean1=mean+9750	#add 9750 because I cut array down to just include gaussian

		ycorr_length=len(ycorr)
		pix_shift_val=(ycorr_length/2) - mean1

		pix_shift.append(pix_shift_val)
		
		l=l+1	

# End cross correlation loop --------------------------------- 
	

	#print len(ycorr)	
	#print pix_shift
	#pix_shift=numpy.array(pix_shift)	
	(mu,sigma)=norm.fit(pix_shift)	# get mean and std dev of array of pixel shift values
	print mu,sigma	
	
	my_gauss=[None]*len(xcorr)
	i=0
	while i < len(xcorr):	#creating an array based on values determined by gaussian fit
		sig2=sig ** 2	
		my_gauss[i] = (amp * (numpy.exp(-0.5 * ((xcorr[i] - mean) ** 2) / sig2))) + sky
		i=i+1
	
# Apply shift to arrays -------------------------------- 
	
	fx_rebin_list_obj=fx_rebin_obj.tolist()
	fx_rebin_list_std=fx_rebin_std.tolist()
	print 'mu=',mu
	if mu < 0:
		val= abs(mu)	# so we can shift properly	
		i=0
		while i < val:
			del fx_rebin_list_obj[0]
			fx_rebin_list_obj.append(1)
			i=i+1
		print 'mu is negative'		
	elif mu >= 0:
		val=mu
		i=0
		while i < val:
			del fx_rebin_list_std[0]	
			fx_rebin_list_std.append(1)
			i=i+1
	print 'mu=',mu		

# Create plots --------------------------------- 
	
	fig=plt.figure(l+1, figsize=(10,10))
	plt.plot([1,2,3])
	
	#Plots target and standard with shift applied
	plt.subplot(311)
	plt.plot(wv_ln_rebin_std, fx_rebin_list_obj, 'red')
	plt.plot(wv_ln_rebin_std, fx_rebin_list_std, 'blue')
	plt.xlabel('wavelength (microns)')
	plt.ylabel('normalized flux')
	target = 'Target: %s' %(obj_name)
	standard = 'Standard: %s' %(std_name)
	plt.annotate(target,xy=(.6,.9),xycoords='axes fraction',xytext=(.6,.9),textcoords='axes fraction',color='red') 
	plt.annotate(standard,xy=(.6,.8),xycoords='axes fraction',xytext=(.6,.8),textcoords='axes fraction',color='blue') 
	#plt.subplots_adjust(hspace=.5)
	
	#Plots example of gaussian fit to cross correlation function
	plt.subplot(312)
	plt.plot(xcorr, ycorr1, 'k.')
	plt.plot(xcorr, my_gauss, 'r--', linewidth=2)
	plt.xlabel('example of fit to cross correlation function')
	
	#print pix_shift

# Transform pixel shift to shift in radial velocity -------------------------------- 
	
	vshift=.426*mu
	err=.426*sigma
	print "vshift=",vshift
	
	rv_obj=rv_std-vshift
	
	rv_err=err+rv_std_err
	print "rv_obj=",rv_obj, "+/-", rv_err, ' km/s'
	
	rv_obj_round=round(rv_obj,4)
	err_round=round(rv_err,4)
	
	pix_shift_conv1=[.426*i for i in pix_shift]
	rv_arr=[(rv_std-i) for i in pix_shift_conv1]
	
	
# Plot histogram of pixel shift values -------------------------------- 
	plt.subplot(313)
	n, bins, patches=plt.hist(rv_arr,normed=1.0,facecolor='green',align='mid') 
	#Plot best fit gaussian over histogram
	y=mlab.normpdf(bins,rv_obj,err)
	plt.plot(bins,y,'r--',linewidth=2)
	plt.xlabel('radial velocity of target')
	plt.ylabel('frequency (normalized)')
	rad='RV = %s +/- %s' %(rv_obj_round,err_round)
	plt.annotate(rad,xy=(.6,.9),xycoords='axes fraction',xytext=(.65,.9),textcoords='axes fraction',color='black')
	plt.subplots_adjust(hspace=.4)

	figname='rv_%s.pdf' %(obj_name)
	plt.savefig(figname)
	
	#plt.figure(l+1)
	#plt.hist(pix_shift)
	
#END RADIAL VELOCITY FUNCTION -----------------------------------
