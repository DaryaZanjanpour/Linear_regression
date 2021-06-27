#first part of the lab s to export the different data of indoor sources

#libraries

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from astropy.io import fits
from astropy.utils.data import download_file
from astropy.io import fits

#importing
#----------------------------------------------------------------
#importing Bulb Data
Bulb_data = np.loadtxt(fname = 'Group C - Bulb.txt', delimiter = '\t', skiprows = 17)
plt.figure(1,figsize=(15,15))
plt.plot(Bulb_data[:,0], Bulb_data[:,1])
plt.xlabel('Pixel'), plt.ylabel('Intensity')
plt.title('Light bulb')
plt.show()
#------------------------------------------------------------------------
#importing Purple Data
Purple_data = np.loadtxt(fname = 'Group C - Purple.txt', delimiter = '\t', skiprows = 17)
plt.figure(2,figsize=(15,15))
plt.plot(Purple_data[:,0], Purple_data[:,1])
plt.xlabel('Pixel'), plt.ylabel('Intensity')
plt.title('Neon')
plt.show()
#--------------------------------------------------------------------------
#importing Sun Data
Sun_data = np.loadtxt(fname = 'S.txt', delimiter = '\t', skiprows = 17)
plt.figure(4,figsize=(15,15))
plt.plot(Sun_data[:,0], Sun_data[:,1])
plt.xlabel('Pixel'), plt.ylabel('Intensity')
plt.title('Sun')
plt.show()
#--------------------------------------------------------------------------
#importing Red Data
Red_data = np.loadtxt(fname = 'Group C - Red.txt', delimiter = '\t', skiprows = 17)
plt.figure(3,figsize=(15,15))
plt.plot(Red_data[:,0], Red_data[:,1])
plt.xlabel('Pixel'), plt.ylabel('Intensity')
plt.title("Hydrogen")
plt.show()


#Calibration for Ocean optic Process

Ocean_table = np.loadtxt('Ocean.txt', skiprows=0, delimiter='\t')
plt.figure(5,figsize=(15,15))
plt.scatter(Ocean_table[:,0], Ocean_table[:,1])
plt.xlabel('Pixel'), plt.ylabel('Wavelenght')
plt.title("Ocean Optic Callibration table")
plt.show()



#Linear fitting by website codes
#Test least squares fitting by simulating some data.

nx= 26 #Number of data poimts 
m = 1.0 #gradiant 
c= 0.0 #Intercept

x = Ocean_table[:,0] # Independent variable  
y = Ocean_table[:,1]  # dependent variable
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
plt.savefig('lsq1.png') 
#difference regression 
plt.figure(7,figsize=(15,15))
plt.scatter(x,10*(y-(mfit*x + cfit))) #rescaling by factor 10, for better vision 
plt.axis('scaled')
plt.show()

#Plot the calibrated version of each spectrum
m=0.165
c=377.798

#Bulb-data
wavelenght_Bulb_data = Bulb_data[:,0]*m+c
plt.figure(8,figsize=(15,15))
plt.plot(wavelenght_Bulb_data,Bulb_data[:,1] )
plt.show()
#sun data
wavelenght_Sun_data = Sun_data[:,0]*m+c
plt.figure(9,figsize=(15,15))
plt.plot(wavelenght_Sun_data,Sun_data[:,1] )
plt.show()
#Neon data
wavelenght_Purple_data = Purple_data[:,0]*m+c
plt.figure(10,figsize=(15,15))
plt.plot(wavelenght_Purple_data,Purple_data[:,1] )
plt.show()
#Hydrogen data
wavelenght_Red_data = Red_data[:,0]*m+c
plt.figure(11,figsize=(15,15))
plt.plot(wavelenght_Red_data,Red_data[:,1] )
plt.show()

#part 2 
#telescope data 
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
#Neon data set
path2data = 'Neon.fit'
hdu = fits.open(path2data)
hdu[0].header
data = hdu[0].data
print(data)

Peaks = np.loadtxt('Peaks.txt', skiprows=0, delimiter='\t')
YPeaks = np.loadtxt('YPeaks.txt', skiprows=0, delimiter='\t')

fig, ax = plt.subplots(figsize=(15, 15))
ax.imshow(data, origin='lower')
plt.show()



fig, ax = plt.subplots(figsize=(15, 15))

ax.plot(data[300, :])
ax.set_xlabel('pixel')
ax.set_ylabel('Intensity')
plt.scatter(Peaks, YPeaks, color='purple')
plt.show()

peaks_neon, _ = find_peaks(data[300,:], height=1096)
print(peaks_neon)

#Wavelenght table 

table = np.loadtxt('Table.txt', skiprows=0, delimiter='\t')
plt.figure(10,figsize=(15,15))
plt.scatter(table[:,1],table[:,0])
plt.xlabel('wavelenght'), plt.ylabel('Intensity')
plt.show()


#Linear fitting by website codes
#Test least squares fitting by simulating some data.

nx= 20 #Number of data poimts 
m = 1.0 #gradiant 
c= 0.0 #Intercept

x = Peaks # pixels 
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

#difference regression 
plt.figure(7,figsize=(15,15))
plt.scatter(x,10*(y-(mfit*x + cfit))) #rescaling by factor 10, for better vision 
plt.axis('scaled')
plt.show()

#Plot the calibrated version of each spectrum
m=0.329
c=551.772

#Vega
path2data = 'Vega.fit'
hdu = fits.open(path2data)
hdu[0].header
data = hdu[0].data
print(data)

pixel = np.loadtxt(fname = 'pixel.txt')

plt.figure(16,figsize=(15,15))
plt.plot((pixel[:,0])*m+c, data[340,:], color='orange')
plt.plot((pixel[:,0])*m+c, (data[340,:]-AvgDark[340,:]//(Avgflat[340,:]-AvgDark[340,:])), color='purple')

plt.show()

#Enif
path2data = 'CCD Image 93 - Enif.fit'
hdu = fits.open(path2data)
hdu[0].header
data = hdu[0].data
print(data)

pixel = np.loadtxt(fname = 'pixel.txt')

plt.figure(16,figsize=(15,15))
plt.plot((pixel[:,0])*m+c, data[340,:], color='orange')
plt.plot((pixel[:,0])*m+c, (data[340,:]-AvgDark[340,:]//(Avgflat[340,:]-AvgDark[340,:])), color='purple')

plt.show()

#Navi
path2data = 'CCD Image 94 - Navi.fit'
hdu = fits.open(path2data)
hdu[0].header
data = hdu[0].data
print(data)

pixel = np.loadtxt(fname = 'pixel.txt')

plt.figure(16,figsize=(15,15))
plt.plot((pixel[:,0])*m+c, data[340,:], color='orange')
plt.plot((pixel[:,0])*m+c, (data[340,:]-AvgDark[340,:]//(Avgflat[340,:]-AvgDark[340,:])), color='purple')

plt.show()

#scheat
path2data = 'CCD Image 95 - scheat.fit'
hdu = fits.open(path2data)
hdu[0].header
data = hdu[0].data
print(data)

pixel = np.loadtxt(fname = 'pixel.txt')

plt.figure(16,figsize=(15,15))
plt.plot((pixel[:,0])*m+c, data[340,:], color='orange')
plt.plot((pixel[:,0])*m+c, (data[340,:]-AvgDark[340,:]//(Avgflat[340,:]-AvgDark[340,:])), color='purple')

plt.show()

#Neptune
path2data = 'CCD Image 96- Neptune.fit'
hdu = fits.open(path2data)
hdu[0].header
data = hdu[0].data
print(data)

pixel = np.loadtxt(fname = 'pixel.txt')

plt.figure(16,figsize=(15,15))
plt.plot((pixel[:,0])*m+c, data[340,:], color='orange')
plt.plot((pixel[:,0])*m+c, (data[340,:]-AvgDark[340,:]//(Avgflat[340,:]-AvgDark[340,:])), color='purple')

plt.show()