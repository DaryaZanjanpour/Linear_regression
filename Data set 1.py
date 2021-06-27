#Data set 1 

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from astropy.io import fits

#----------------------------------------------------------------
#importing Bulb Data
Bulb_data = np.loadtxt(fname = 'Group C - Bulb.txt', delimiter = '\t', skiprows = 17)

plt.title("Bulb")
plt.figure(1,figsize=(15,15))
plt.plot(Bulb_data[:,0], Bulb_data[:,1])
plt.xlabel('Pixel'), plt.ylabel('Intensity')
plt.show()



#------------------------------------------------------------------------
#importing Purple Data
Purple_data = np.loadtxt(fname = 'Group C - Purple.txt', delimiter = '\t', skiprows = 17)

plt.title("Neon")
plt.figure(2,figsize=(15,15))
plt.plot(Purple_data[:,0], Purple_data[:,1])
plt.xlabel('Pixel'), plt.ylabel('Intensity')
plt.show()
#--------------------------------------------------------------------------









#importing Sun Data
Sun_data = np.loadtxt(fname = 'S.txt', delimiter = '\t', skiprows = 17)

plt.title("Sun")
plt.figure(4,figsize=(15,15))
plt.plot(Sun_data[:,0], Sun_data[:,1])
plt.xlabel('Pixel'), plt.ylabel('Intensity')
plt.show()




#--------------------------------------------------------------------------
#importing Red Data
Red_data = np.loadtxt(fname = 'Group C - Red.txt', delimiter = '\t', skiprows = 17)

plt.title("Hydrogen")
plt.figure(3,figsize=(15,15))
plt.plot(Red_data[:,0], Red_data[:,1])
plt.xlabel('Pixel'), plt.ylabel('Intensity')
plt.show()


peaks_red, _ = find_peaks(Red_data[:,1], height=5000)

print("peaks_red=", peaks_red)


v = [585,594,607,609,614,626,633,638,640,650,667]

#____________________________________________________
#us_table


Tableus = np.loadtxt('Tableus.txt', skiprows=0, delimiter='\t')
print(Tableus)
plt.figure(5,figsize=(15,15))
plt.scatter(Tableus[:,1],Tableus[:,0])
plt.show()

peaks_Tableus, _ = find_peaks(Tableus[:,1], height=10)
print(peaks_Tableus)
new = np.array([1224, 1282, 1365, 1379, 1409, 1490, 1536, 1568, 1582, 1653, 1773])
plt.scatter(Tableus[:,1],v)
plt.show()
 


#___________________________________________
#matrix and calibration

ye= Tableus[:,1]
x = new
nx = 11 # Number of data points
m = 1.0 # Gradient  
c = 0.0 # Intercept


y = m * x + c
sigma = 1.0


plt.figure(6,figsize=(15,15))
plt.plot(x,ye,'o',label='data') 
plt.xlabel('x')
plt.ylabel('y') 
#plt.show()

ma = np.array([ [np.sum(x**2), np.sum(x)],[np.sum(x), nx ]]) 
mc = np.array([ [np.sum(x*ye)],[np.sum(ye)]]) 

mai = np.linalg.inv(ma)
print ('Test matrix inversion gives identity',np.dot(mai,ma))
md = np.dot(mai,mc)


mfit = md[0,0]
cfit = md[1,0]
plt.plot(x, mfit*x + cfit)

plt.axis('scaled')

plt.show()
print(mfit,cfit)


#----------------------------------
#difference regression 

plt.figure(16,figsize=(15,15))

plt.scatter(x,v-(mfit*x + cfit))
plt.axis('scaled')

plt.show()

#----------------------------------------------------------------
#using calibration components


m = 0.2
c= 318.058

wavelenght_Bulb_data = Bulb_data[:,0]*m+c


plt.figure(7,figsize=(15,15))
plt.plot(wavelenght_Bulb_data,Bulb_data[:,1] )
plt.show()


#------------------------------------------------------------

wavelenght_Sun_data = Sun_data[:,0]*m+c
plt.figure(8,figsize=(15,15))
plt.plot(wavelenght_Sun_data,Sun_data[:,1] )
plt.show()



#--------------------------------------------------------------------------


wavelenght_Purple_data = Purple_data[:,0]*m+c
plt.figure(9,figsize=(15,15))
plt.plot(wavelenght_Purple_data,Purple_data[:,1] )
plt.show()



#-----------------------------------------------------------------
#NEON

from astropy.utils.data import download_file
from astropy.io import fits
image_file = download_file('http://www.astro.utoronto.ca/~astrolab/files/Data/Spectroscopy/Oct24/GroupC/CCD%20Image%2091%20-%20Neon.fit', cache=True )
hdu_list = fits.open(image_file)
hdu_list.info()
image_data = hdu_list[0].data

print(type(image_data))
print(image_data.shape)


plt.figure(10,figsize=(15,15))
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
NBINS = 1000
histogram = plt.hist(image_data.flatten(), NBINS)


from matplotlib.colors import LogNorm
plt.imshow(image_data, cmap='gray', norm=LogNorm())
plt.show()


"""

same process that we did for 2D array

will happed here

"""
slices=(image_data[:,300])
print(slices)
peaks_image, _ = find_peaks(slices, height=1077)
print('peaks=', peaks_image)
#peaks_image_y = peaks_image[:,300]
#print('y=', peaks_image_y )
print('this=', slices[11], slices[14], slices[22], slices[27], slices[29], slices[33], slices[37], slices[40], slices[46], slices[50])
new2 = np.array([1078, 1083, 1083, 1085, 1086, 1090, 1101, 1102, 1152])



#___________________________________________
#tabelus
v2 = [585,594,607,609,614,638,640,650,667]

Tableus2 = np.loadtxt('Tableus2.txt', skiprows=0, delimiter='\t')
print(Tableus2)
plt.figure(17,figsize=(15,15))
plt.scatter(Tableus2[:,1],Tableus2[:,0])
plt.show()

peaks_Tableus2, _ = find_peaks(Tableus2[:,1], height=10)
print(peaks_Tableus2)
new2 = np.array([1090, 1078, 1102, 1083, 1086, 1101, 1083, 1152, 1085])
plt.scatter(Tableus2[:,1],v2)
plt.show()

#___________________________________________
#matrix and calibration

ye= Tableus2[:,1]
x = new2
nx = 9 # Number of data points
m = 0.539
c= 19.428


y = m * x + c
sigma = 1.0


plt.figure(18,figsize=(15,15))
plt.plot(x,ye,'o',label='data') 
plt.xlabel('x')
plt.ylabel('y') 
#plt.show()

ma = np.array([ [np.sum(x**2), np.sum(x)],[np.sum(x), nx ]]) 
mc = np.array([ [np.sum(x*ye)],[np.sum(ye)]]) 

mai = np.linalg.inv(ma)
print ('Test matrix inversion gives identity',np.dot(mai,ma))
md = np.dot(mai,mc)


mfit = md[0,0]
cfit = md[1,0]
plt.plot(x, mfit*x + cfit)

plt.axis('scaled')

plt.show()
print(mfit,cfit)


#----------------------------------
#difference regression 

plt.figure(19,figsize=(15,15))

plt.scatter(x,v2-(mfit*x + cfit))
plt.axis('scaled')

plt.show()








#_____________________________________________________________________
#vega
m = 0.539
c= 19.428



from astropy.utils.data import download_file
from astropy.io import fits
image_file = download_file('http://www.astro.utoronto.ca/~astrolab/files/Data/Spectroscopy/Oct24/GroupC/CCD%20Image%2092%20-%20Vega.fit', cache=True )
hdu_list = fits.open(image_file)
hdu_list.info()
image_data = hdu_list[0].data

print(type(image_data))
print(image_data.shape)


plt.figure(11,figsize=(15,15))


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
NBINS = 1000
histogram = plt.hist(image_data.flatten(), NBINS)


from matplotlib.colors import LogNorm
plt.imshow(image_data, cmap='gray', norm=LogNorm())
plt.show()




m = 0.2
c= 318.058

wavelenght_Bulb_data = Bulb_data[:,0]*m+c


plt.figure(7,figsize=(15,15))
plt.plot(wavelenght_Bulb_data,Bulb_data[:,1] )
plt.show()






#--------------------------------------------

"""

wavelenght_vega = image_data[300,0]*m+c




plt.figure(20,figsize=(15,15))
plt.plot(wavelenght_vega, image_data[:,1])
plt.show()
"""
#_____________________________________________________________________
#navi



from astropy.utils.data import download_file
from astropy.io import fits
image_file = download_file('http://www.astro.utoronto.ca/~astrolab/files/Data/Spectroscopy/Oct24/GroupC/CCD%20Image%2094%20-%20Navi.fit', cache=True )
hdu_list = fits.open(image_file)
hdu_list.info()
image_data = hdu_list[0].data

print(type(image_data))
print(image_data.shape)


plt.figure(12,figsize=(15,15))

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
NBINS = 1000
histogram = plt.hist(image_data.flatten(), NBINS)


from matplotlib.colors import LogNorm
plt.imshow(image_data, cmap='gray', norm=LogNorm())
plt.show()

#_____________________________________________________________________
#enif



from astropy.utils.data import download_file
from astropy.io import fits
image_file = download_file('http://www.astro.utoronto.ca/~astrolab/files/Data/Spectroscopy/Oct24/GroupC/CCD%20Image%2093%20-%20Enif.fit', cache=True )
hdu_list = fits.open(image_file)
hdu_list.info()
image_data = hdu_list[0].data

print(type(image_data))
print(image_data.shape)



plt.figure(13,figsize=(15,15))

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
NBINS = 1000
histogram = plt.hist(image_data.flatten(), NBINS)


from matplotlib.colors import LogNorm
plt.imshow(image_data, cmap='gray', norm=LogNorm())
plt.show()



#_____________________________________________________________________
#schea



from astropy.utils.data import download_file
from astropy.io import fits
image_file = download_file('http://www.astro.utoronto.ca/~astrolab/files/Data/Spectroscopy/Oct24/GroupC/CCD%20Image%2095%20-%20scheat.fit', cache=True )
hdu_list = fits.open(image_file)
hdu_list.info()
image_data = hdu_list[0].data

print(type(image_data))
print(image_data.shape)


plt.figure(14,figsize=(15,15))


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
NBINS = 1000
histogram = plt.hist(image_data.flatten(), NBINS)


from matplotlib.colors import LogNorm
plt.imshow(image_data, cmap='gray', norm=LogNorm())
plt.show()


#_____________________________________________________________________
#Nepton



from astropy.utils.data import download_file
from astropy.io import fits
image_file = download_file('http://www.astro.utoronto.ca/~astrolab/files/Data/Spectroscopy/Oct24/GroupC/CCD%20Image%2096-%20Neptune.fit', cache=True )
hdu_list = fits.open(image_file)
hdu_list.info()
image_data = hdu_list[0].data

print(type(image_data))
print(image_data.shape)


plt.figure(15,figsize=(15,15))

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
NBINS = 1000
histogram = plt.hist(image_data.flatten(), NBINS)


from matplotlib.colors import LogNorm
plt.imshow(image_data, cmap='gray', norm=LogNorm())
plt.show()

