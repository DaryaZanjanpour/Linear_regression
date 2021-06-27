# -*- coding: utf-8 -*-
"""
Created on Tue Nov  5 14:00:09 2019

@author: Z
"""

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from astropy.io import fits
from astropy.utils.data import download_file
from astropy.io import fits
from PIL import Image
from PIL import ImageOps




#---------------------------------------------
#Neon data set
path2data = 'Neon.fit'
hdu = fits.open(path2data)
hdu[0].header
data = hdu[0].data
print(data)



Peaks = np.loadtxt('DemoPeaks.txt', skiprows=0, delimiter='\t')
YPeaks = np.loadtxt('DemoYPeaks.txt', skiprows=0, delimiter='\t')

fig, ax = plt.subplots(figsize=(15, 15))
ax.imshow(data, origin='lower')
plt.show()



fig, ax = plt.subplots(figsize=(15, 15))



peaks_neon, _ = find_peaks(data[300,::-1], height=1280)

print(peaks_neon)
flipdata= data[::,::-1]
print(flipdata[300,peaks_neon])

ax.plot(data[300,::-1])
ax.set_xlabel('pixel')
ax.set_ylabel('Intensity')
plt.scatter(peaks_neon, flipdata[300,peaks_neon], color='purple')
plt.show()


#Wavelenght table 

table = np.loadtxt('DemoTable.txt', skiprows=0, delimiter='\t')
plt.figure(10,figsize=(15,15))
plt.scatter(table[:,1],table[:,0])
plt.xlabel('wavelenght'), plt.ylabel('Intensity')
plt.show()


#Linear fitting by website codes
#Test least squares fitting by simulating some data.

nx= 11 #Number of data poimts 
m = 1.0 #gradiant 
c= 0.0 #Intercept

x = peaks_neon # pixels 
y = table[:,1]  # dependent variable
#Generate Gaussian errors 
sigma = 1.0     # Measurement error 
np.random.seed(1)    # init random no. generator 
errors = sigma*np.random.randn(nx)  # Gaussian distributed errors 
ye = y + errors 
plt.figure(6,figsize=(15,15))                    # Add the noise
plt.plot(x,ye,'o',label='data')
plt.xlabel('x')
plt.ylabel('y')
# Construct the matrices 
ma = np.array([ [np.sum(x**2), np.sum(x)],[np.sum(x), nx ] ]  )
mc = np.array([ [np.sum(x*ye)],[np.sum(ye)]]) 
# Compute the gradient and intercept
mai = np.linalg.inv(ma) 
print('Test matrix inversion gives identity',np.dot(mai,ma))
md = np.dot(mai,mc)     # matrix multiply is dot
# Overplot the best fit
mfit = md[0,0] 
cfit = md[1,0] 
plt.plot(x, mfit*x + cfit)
plt.axis('scaled') 
plt.text(5,15,'m = {:.3f}\nc = {:.3f}'.format(mfit,cfit))
plt.savefig('lsq2.png') 
plt.show()

#Plot the calibrated version of each spectrum
m=0.225
c=487.315
#-----------