#!/usr/bin/python
# -*- coding: UTF-8 -*-
# 2D random aperture generator
#
# Assumes a guassian distribution with specified corelation lengths
#
#
# ap_gen(nx,ny,lx,ly,hurst,mean,stdev)
#
#

import numpy as np, matplotlib.pyplot as plot

def ap_gen(model = None,fracture = None):
	nx = model.nx
	ny = model.ny
	
	xx = np.arange(-nx/2,nx/2-1,step=1)
	yy = np.arange(-ny/2,ny/2-1,step=1)
	xx = xx*fracture.lx/nx
	yy = yy*fracture.ly/ny
	[x,y] = np.meshgrid(xx,yy)
	r = (x**2 + y**2) * 2 * np.pi
	r = np.fft.fftshift(r)
	
	# generate the power spectrum:
	f = 2.0*np.sqrt(1/(1+r)**(1+fracture.hurst)/(nx*ny))
	p = np.pi*(2*np.random.rand(nx-1,ny-1))

	re = f*np.cos(p)
	im = f*np.sin(p)
	
	print(np.fft.ifft2(np.complex(re,im)))
	b = np.real(np.fft.ifft2(np.complex(re,im)))
	b = (b-np.mean(b[:]))/np.std(b[:])*fracture.sigma+fracture.mu
	
	# set all b values less than eps == eps
	eps=1e-8
	b[b<1e-8]=1e-8
	return b	


if __name__ == '__main__':
	print(__name__)
	class model(object):
		def __init__(self):
			self.nx = 100
			self.ny = 100
	class fracture(object):
		def __init__(self):
			self.lx = 10
			self.ly = 10
			self.mu = 150
			self.sigma = 20
			self.hurst = 0.8
	
	
	b = ap_gen(model = model(), fracture = fracture())
