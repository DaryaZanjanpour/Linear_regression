# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 15:56:17 2019

@author: Z
"""

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from astropy.io import fits





m = 0.539
c= 19.428



from astropy.utils.data import download_file
from astropy.io import fits
image_file = download_file('http://www.astro.utoronto.ca/~astrolab/files/Data/Spectroscopy/Oct24/GroupC/CCD%20Image%2092%20-%20Vega.fit', cache=True )
hdu_list = fits.open(image_file)
hdu_list.info()
image_data = hdu_list[0].data

print("type_image_data=", type(image_data))
print(image_data.shape)
print("image_data=", image_data)
print("image_data 00=", image_data[300,:])

plt.figure(11,figsize=(15,15))

print("pixel=", value(image_data))
#plt.plot(image_data[300,:],image_data[:,0])
#plt.plot(image_data[:,0],image_data[:,1])
plt.show()


hdu_list.close()
image_data = fits.getdata(image_file)
print(type(image_data))
print(image_data.shape)
plt.imshow(image_data, cmap='gray')
plt.colorbar()

print('Min:', np.min(image_data))
print('Max:', np.max(image_data))
print('Mean:', np.mean(image_data))
print('Stdev:', np.std(image_data))
print(type(image_data.flatten()))


#--------------------------------------------
'''


wavelenght_vega = image_data[:,0]*m+c




plt.figure(20,figsize=(15,15))
plt.plot(wavelenght_vega, image_data[:,1])
plt.show()


'''