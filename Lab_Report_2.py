#Telescope_data

#Neon_calibration_data


from astropy import units as u
import numpy as np
from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from astropy.io import fits

#-------------------------------
 #flat frames avereging for 120 ex time data
 
 

 

path2data = 'flat1.fit'
hdu = fits.open(path2data)
hdu[0].header
data1 = hdu[0].data
print(data1)


path2data = 'flat2.fit'
hdu = fits.open(path2data)
hdu[0].header
data2 = hdu[0].data
print(data2)


path2data = 'flat3.fit'
hdu = fits.open(path2data)
hdu[0].header
data3 = hdu[0].data
print(data3)


path2data = 'flat4.fit'
hdu = fits.open(path2data)
hdu[0].header
data4 = hdu[0].data
print(data4)


path2data = 'flat5.fit'
hdu = fits.open(path2data)
hdu[0].header
data5 = hdu[0].data
print(data5)

Avgflat = (data1+data2+data3+data4+data5)/5
    
print("Avg flat Data=", Avgflat[300,:])
 #-------------------------------
 #dark frames avereging for 120 ex time data
 
 

 

path2data = 'dark1_120.fit'
hdu = fits.open(path2data)
hdu[0].header
data1 = hdu[0].data
print(data1)


path2data = 'dark2_120.fit'
hdu = fits.open(path2data)
hdu[0].header
data2 = hdu[0].data
print(data2)


path2data = 'dark3_120.fit'
hdu = fits.open(path2data)
hdu[0].header
data3 = hdu[0].data
print(data3)


path2data = 'dark4_120.fit'
hdu = fits.open(path2data)
hdu[0].header
data4 = hdu[0].data
print(data4)


path2data = 'dark5_120.fit'
hdu = fits.open(path2data)
hdu[0].header
data5 = hdu[0].data
print(data5)

AvgDark = (data1+data2+data3+data4+data5)/5
    
print("Avg Dark Data=", AvgDark[300,:])
#---------------------------------------------------


path2data = 'Neon.fit'
hdu = fits.open(path2data)
hdu[0].header
data = hdu[0].data
print(data)



fig, ax = plt.subplots(figsize=(15, 15))
ax.imshow(data, origin='lower')
plt.show()



fig, ax = plt.subplots(figsize=(15, 15))

ax.plot(data[300, :])
ax.set_xlabel('axis=1')
ax.set_ylabel('Intensity')
plt.show()



peaks_neon, _ = find_peaks(data[300,:], height=1096)

print(peaks_neon)

print("paek point as bellow")
print(data[300,2]-AvgDark[300,2]/(Avgflat[300,2]-AvgDark[300,2]))
print(data[300,6]-AvgDark[300,6]/(Avgflat[300,6]-AvgDark[300,6]))
print(data[300,30]-AvgDark[300,30]/(Avgflat[300,30]-AvgDark[300,30]))
print(data[300,38]-AvgDark[300,38]/(Avgflat[300,38]-AvgDark[300,38]))
print(data[300,44]-AvgDark[300,44]/(Avgflat[300,44]-AvgDark[300,44]))
print(data[300,62]-AvgDark[300,62]/(Avgflat[300,62]-AvgDark[300,62]))
print(data[300,112]-AvgDark[300,112]/(Avgflat[300,112]-AvgDark[300,112]))
print(data[300,121]-AvgDark[300,121]/(Avgflat[300,121]-AvgDark[300,121]))
print(data[300,140]-AvgDark[300,140]/(Avgflat[300,140]-AvgDark[300,140]))
print(data[300,155]-AvgDark[300,155]/(Avgflat[300,155]-AvgDark[300,155]))
print(data[300,161]-AvgDark[300,161]/(Avgflat[300,161]-AvgDark[300,161]))
print(data[300,185]-AvgDark[300,185]/(Avgflat[300,185]-AvgDark[300,185]))
print(data[300,190]-AvgDark[300,190]/(Avgflat[300,190]-AvgDark[300,190]))
print(data[300,201]-AvgDark[300,201]/(Avgflat[300,201]-AvgDark[300,201]))
print(data[300,217]-AvgDark[300,217]/(Avgflat[300,217]-AvgDark[300,217]))
print(data[300,229]-AvgDark[300,229]/(Avgflat[300,229]-AvgDark[300,229]))
print(data[300,241]-AvgDark[300,241]/(Avgflat[300,241]-AvgDark[300,241]))
print(data[300,246]-AvgDark[300,246]/(Avgflat[300,246]-AvgDark[300,246]))
print(data[300,257]-AvgDark[300,257]/(Avgflat[300,257]-AvgDark[300,257]))
print(data[300,261]-AvgDark[300,261]/(Avgflat[300,261]-AvgDark[300,261]))
print(data[300,292]-AvgDark[300,292]/(Avgflat[300,292]-AvgDark[300,292]))
print(data[300,306]-AvgDark[300,306]/(Avgflat[300,306]-AvgDark[300,306]))
print(data[300,314]-AvgDark[300,314]/(Avgflat[300,314]-AvgDark[300,314]))
print(data[300,755]-AvgDark[300,755]/(Avgflat[300,755]-AvgDark[300,755]))
print(data[300,759]-AvgDark[300,759]/(Avgflat[300,759]-AvgDark[300,759]))
print(data[300,763]-AvgDark[300,763]/(Avgflat[300,763]-AvgDark[300,763]))



#new = np.array( [114, 127, 154, 156, 193, 206.4, 231, 266, 270, 328, 333, 342, 343, 348, 498, 499.6, 511, 574, 602, 694, 717, 1067, 1652])

new =  np.loadtxt('new20.txt')

Tableus = np.loadtxt('Table.txt', skiprows=0, delimiter='\t')
print(Tableus)
plt.figure(5,figsize=(15,15))
plt.scatter(Tableus[:,0],Tableus[:,1])
plt.show()

peaks_Tableus, _ = find_peaks(Tableus[:,1], height=1100)
print(peaks_Tableus)

plt.scatter(new, Tableus[:,1])
plt.show()



ye= Tableus[:,1]
x = new
nx = 20 # Number of data points
m = 0 # Gradient  
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

plt.scatter(x,new-(mfit*x + cfit))
plt.axis('scaled')

plt.show()


#---------------------------------------------------
#vega
m=0.035
c=553.67
path2data = 'Vega.fit'
hdu = fits.open(path2data)
hdu[0].header
data = hdu[0].data
print(data)

pixel = np.loadtxt(fname = 'pixel.txt')

plt.figure(16,figsize=(15,15))
plt.plot((pixel[:,0])*m+c,data[340,:])

plt.show()


