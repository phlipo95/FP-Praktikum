import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import sem
from scipy.optimize import curve_fit
from scipy.signal import find_peaks_cwt
from uncertainties import ufloat
import scipy.constants as c
import math
from lmfit import minimize, Parameter, Model

def stab(r1,r2,L):
    return (1-L/r1)*(1-L/r2)

r1 = 1000000; r2 = 1.4;
if r1<r2 : a = r1
elif r1>r2: a = r2
else: a = 2 * r1
L = np.linspace(0,a,100)
plt.plot(L ,stab(r1, r2, L),label=r'$r_1 = \, \text{flat} , r_2 = 1.4 \, \text{m} $')


r1 = 1.4; r2 = 1.4;
if r1<r2 : a = r1
elif r1>r2: a = r2
else: a = 2*r1
L = np.linspace(0,a,100)
plt.plot(L ,stab(r1, r2, L),label=r'$r_1 = \, \text{flat}, r_2 = 1.4 \, \text{m}$')

plt.grid()
plt.xlim(xmin=0, xmax=2.8)
plt.ylim(ymax=1)
plt.xlabel('Resonatorlänge / m')
plt.ylabel('Stabilitätsparameter')
plt.legend(loc='best')
plt.savefig('build/Stabilisationsparameter.pdf')
plt.close()

#Intensität in Abhängikeit der Polarisation
winkel, intensi = np.loadtxt('Polarisation.txt', unpack=True)
winkel = np.radians(winkel)
def fsin(x, a, b, c, d):
    return (np.sin(b*x +c))*a +d 

params, cov = curve_fit(fsin, winkel, intensi,
        bounds=([80,1.8,0,80],[150,2.2,2*np.pi,120]))
linw = np.linspace(-0.5,max(winkel)+0.5,100)
plt.plot(winkel, intensi, 'x', label='Messwerte')
plt.plot(linw , fsin(linw, *params), label='Fit')
plt.xlim(xmin=-0.5, xmax=(max(winkel)+0.5))
plt.ylim((-0.5, max(intensi)+0.5))
plt.xlabel(r'Polarisationsrichtung / 2 $\pi$')
plt.ylabel('Intensität / µA' )
plt.legend(loc='best')
plt.savefig('build/Polar.pdf')
plt.close()

#Intensität der verschiedenen TEM Moden
A, TEM00, TEM10 = np.loadtxt('TEM-Moden.txt', unpack=True)

def gausian(x, a, v, my):
    return a*np.exp(-2*(x+my)**2/(v**2))
params, cov = curve_fit(gausian, A, TEM00, bounds=([3,5,-20],[6,10,-10]))
errors= np.sqrt(np.diag(cov))
print('Parameter des Gaußfits a', params[0], errors[0], 'v ', params[1],
        errors[1], 'my ', params[2], errors[2])

def gausian2( x, a, v, my):
    return a*8*(x-my)**2/v**2*np.exp(-2*(x-my)**2/v**2)
params2, cov2 = curve_fit(gausian2, A, TEM10,bounds=([0,0,12],[1,10,20])) 
errors2= np.sqrt(np.diag(cov2))
print('Parameter2 des Gaußfits a', params2[0], errors2[0], 'v ', params2[1],
        errors2[1], 'my ', params2[2], errors2[2])

x = np.linspace(0,max(A+1),100)
plt.plot(x, gausian(x, *params), label='Fit')
plt.plot(A, TEM00, 'x', label='Messwerte')
plt.xlim(xmin=0, xmax=max(A+1))
plt.xlabel('Abstand / mm')
plt.ylabel(r'$\text{TEM}_{00}$ I / µA')
plt.legend(loc='best')
plt.savefig('build/TEM00.pdf')
plt.close()

plt.plot(x, gausian2(x, *params2),label='Fit')
plt.plot(A, TEM10, 'x', label='Messwerte')
plt.xlim(xmin=0, xmax=max(A+1))
plt.xlabel('Auslenkung / mm')
plt.ylabel(r'$\text{TEM}_{10}$ I / µA')
plt.legend(loc='best')
plt.savefig('build/TEM10.pdf')
plt.close()

#Bestimmung der Gitterkonstanten in Abhängigkeit der Wellenlänge

def wl(a):
    return a/(10**(5)*1.741)
abstand = np.array([0.112,0.113,0.118,0.110,0.112,0.116])
Wl = np.array([np.mean(wl(abstand)), sem(wl(abstand))])
print('Berechnete Wellenlänge des HeNe Lasers: ',Wl)
