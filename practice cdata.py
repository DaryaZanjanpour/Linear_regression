# -*- coding: utf-8 -*-
"""
Created on Wed Nov  6 11:09:32 2019

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
#Flat data
 #flat frames avereging for 120 ex time data
 
path2data = 'flat1.fit'
hdu = fits.open(path2data)
hdu[0].header
data1 = hdu[0].data
#print(data1)

path2data = 'flat2.fit'
hdu = fits.open(path2data)
hdu[0].header
data2 = hdu[0].data
#print(data2)

path2data = 'flat3.fit'
hdu = fits.open(path2data)
hdu[0].header
data3 = hdu[0].data
#print(data3)

path2data = 'flat4.fit'
hdu = fits.open(path2data)
hdu[0].header
data4 = hdu[0].data
#print(data4)

path2data = 'flat5.fit'
hdu = fits.open(path2data)
hdu[0].header
data5 = hdu[0].data
#print(data5)

Avgflat = (data1+data2+data3+data4+data5)/5
    
#print("Avg flat Data=", Avgflat[300,:])
#-------------------------------------------
#dark frames avereging for 120 ex time data

path2data = 'dark1_120.fit'
hdu = fits.open(path2data)
hdu[0].header
data1 = hdu[0].data
#print(data1)

path2data = 'dark2_120.fit'
hdu = fits.open(path2data)
hdu[0].header
data2 = hdu[0].data
#print(data2)

path2data = 'dark3_120.fit'
hdu = fits.open(path2data)
hdu[0].header
data3 = hdu[0].data
#print(data3)

path2data = 'dark4_120.fit'
hdu = fits.open(path2data)
hdu[0].header
data4 = hdu[0].data
#print(data4)

path2data = 'dark5_120.fit'
hdu = fits.open(path2data)
hdu[0].header
data5 = hdu[0].data
#print(data5)

AvgDark = (data1+data2+data3+data4+data5)/5
    
#print("Avg Dark Data=", AvgDark[300,:])

#---------------------------------------------
#planck function 
pixel = np.loadtxt(fname = 'pixel.txt')
Planck = np.loadtxt(fname = 'empty.txt')
h= 6.62e-34
k= 1.38e-23
clight= 3.0e8
m=0.225
c=487.315
WL = (pixel[:,0])*m+c
v= clight/WL


def Planck(T):
    
    return (2*h*v**3/clight**2)*(1/np.exp(h*v/k*T)-1)

    


#fitiing














#Neon data set
path2data = 'Neon.fit'
hdu = fits.open(path2data)
hdu[0].header
data = hdu[0].data
print(data)

cdata= (-1)*((data-AvgDark)/(Avgflat - AvgDark))*Planck(3200) #correcteddata #3200k for lamp
print(cdata)


Peaks = np.loadtxt('DemoPeaks.txt', skiprows=0, delimiter='\t')
YPeaks = np.loadtxt('DemoYPeaks.txt', skiprows=0, delimiter='\t')

#fig, ax = plt.subplots(figsize=(15, 15))
#ax.imshow(cdata, origin='lower')
#plt.show()



fig, ax = plt.subplots(figsize=(15, 15))



flipdata= cdata[::,::-1]


peaks_neon, _ = find_peaks(flipdata[300,:], height=1.96e-35, distance=400)

print("peaks neon=",peaks_neon)

print(flipdata[300,peaks_neon])



ax.plot(cdata[300,::-1])
ax.set_xlabel('pixel')
ax.set_ylabel('Relative Intensity')
plt.scatter(peaks_neon, flipdata[300,peaks_neon], color='purple')
plt.show()