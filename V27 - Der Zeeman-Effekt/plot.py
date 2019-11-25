import matplotlib.pyplot as plt
import numpy as np
from uncertainties import ufloat
from uncertainties import unumpy as unp
from scipy.optimize import curve_fit

# Hysterese Kurve
# B-Feld in mT, Strom in A
Bhoch = np.array([4, 61, 119, 178, 245, 300, 374, 427, 488, 549, 600, 666, 732, 790, 835, 879, 935, 968, 999, 1030, 1057])
Ahoch = np.arange(0, 21)

Brunter = np.array([1057, 1036, 1004, 975, 935, 893, 849, 791, 734, 678, 600, 562, 496, 432, 371, 303, 244, 182, 129, 58, 7])
Arunter = np.arange(0, 21)[::-1]

# Fitfunktion
def f(I, a, b, c, d):
    return a * np.arctan(b * I + c) + d

print('Ansteigender Strom:')
params1, covariance1 = curve_fit(f, Ahoch, Bhoch)
for i in range(4):
    print('Parameter =', params1[i], np.sqrt(covariance1[i, i]))

print('Abfallender Strom:')
params2, covariance2 = curve_fit(f, Arunter, Brunter)
for i in range(4):
    print('Parameter =', params2[i], np.sqrt(covariance2[i, i]))

plt.plot(Ahoch, Bhoch, 'rx', label='Messwerte beim hochdrehen')
plt.plot(Ahoch, f(Ahoch, *params1), 'r-', linewidth=0.5, label='Fit')

plt.plot(Arunter, Brunter, 'gx', label='Messwerte beim runterdrehen')
plt.plot(Arunter, f(Arunter, *params2), 'g--', linewidth=0.5, label='Fit')

plt.legend(loc='best')
plt.xlabel(r'$I$ / A')
plt.ylabel(r'$B$ / mT')
plt.savefig('Bilder/Hysterese.pdf')
plt.close()



# Berechnung des Dispersionsgebietes
def Dispersionsgebiet(lam, n, d=0.004):
    return lam**2 / (2*d) * np.sqrt(1 / (n**2 - 1))

WellenlängeRot = 643.8 * 10**(-9)
nRot = 1.4567
WellenlängeBlau = 480.0 * 10**(-9)
nBlau = 1.4635

lamRot = Dispersionsgebiet(WellenlängeRot, nRot)
lamBlau = Dispersionsgebiet(WellenlängeBlau, nBlau)
print('--------------------------------------')
print('Dispersionsgebiet Rot =', lamRot)
print('Dispersionsgebiet Blau =', lamBlau)

# Berechnung von Delta S (großes Delta!)
def DS(A):
    return A[1:] - A[:10]

# Berechnung von delta S (kleines delta!)
def dS(A, B):
    return B - A

# Berechnung der Wellenlängen Verschiebung
def deltaLambda(dS, DS, lam):
    return 0.5 * lam * dS / DS

# Alles in Pixeln!!
# BoB = Blau ohne B-Feld
BoB = np.array([
np.mean([1600, 1598, 1596]),
np.mean([1734, 1728, 1730]),
np.mean([1866, 1866, 1866]),
np.mean([1996, 1998, 2006]),
np.mean([2132, 2128, 2126]),
np.mean([2260, 2258, 2260]),
np.mean([2384, 2384, 2384]),
np.mean([2504, 2508, 2510]),
np.mean([2632, 2628, 2628]),
np.mean([2750, 2752, 2752]),
np.mean([2878, 2870, 2872]), ])

# B6ASig = Blau mit B-Feld (6 Ampere), Sigma-Komponente
B6ASig1 = np.array([
np.mean([1568, 1564, 1558]),
np.mean([1708, 1704, 1700]),
np.mean([1840, 1836, 1832]),
np.mean([1966, 1968, 1966]),
np.mean([2102, 2102, 2098]),
np.mean([2228, 2230, 2226]),
np.mean([2358, 2352, 2356]),
np.mean([2486, 2478, 2482]),
np.mean([2602, 2602, 2606]),
np.mean([2722, 2726, 2726]), ])

B6ASig2 = np.array([
np.mean([1624, 1626, 1628]),
np.mean([1760, 1758, 1762]),
np.mean([1892, 1894, 1894]),
np.mean([2026, 2024, 2026]),
np.mean([2154, 2156, 2158]),
np.mean([2280, 2282, 2284]),
np.mean([2412, 2408, 2412]),
np.mean([2530, 2536, 2536]),
np.mean([2656, 2654, 2654]),
np.mean([2780, 2778, 2780]), ])

# B20APi = Blau mit B-Feld (20 Ampere), Pi-Komponente
B20APi1 = np.array([
np.mean([1604, 1600, 1600]),
np.mean([1744, 1744, 1748]),
np.mean([1884, 1882, 1882]),
np.mean([2020, 2022, 2026]),
np.mean([2154, 2158, 2162]),
np.mean([2284, 2290, 2292]),
np.mean([2420, 2426, 2430]),
np.mean([2548, 2554, 2556]),
np.mean([2678, 2682, 2684]),
np.mean([2804, 2806, 2814]), ])

B20APi2 = np.array([
np.mean([1660, 1664, 1666]),
np.mean([1814, 1810, 1808]),
np.mean([1946, 1952, 1948]),
np.mean([2092, 2086, 2086]),
np.mean([2220, 2220, 2218]),
np.mean([2348, 2354, 2354]),
np.mean([2480, 2478, 2480]),
np.mean([2612, 2608, 2610]),
np.mean([2736, 2742, 2734]),
np.mean([2866, 2858, 2862]), ])

print('--------------------------------------')
DSB = DS(BoB)
dS1 = dS(B6ASig1, B6ASig2)
dS2 = dS(B20APi1, B20APi2)
deltaLambda1 = deltaLambda(dS1, DSB, lamBlau)
deltaLambda1 = ufloat(np.mean(deltaLambda1), np.std(deltaLambda1))
deltaLambda2 = deltaLambda(dS2, DSB, lamBlau)
deltaLambda2 = ufloat(np.mean(deltaLambda2), np.std(deltaLambda2))
print('Ds =', DSB)
print('ds1 =', dS1)
#print('Verschiebung mit der blauen Linie(6A) =', deltaLambda1)
print('Mittelwert der Verschiebung =', deltaLambda1)
print('------')
print('ds2 =', dS2)
#print('Verschiebung mit der blauen Linie(20A) =', deltaLambda2)
print('Mittelwert der Verschiebung =', deltaLambda2)


# RoB = Rot ohne B-Feld
RoB = np.array([
np.mean([736, 724, 724]),
np.mean([963, 972, 981]),
np.mean([1218, 1212, 1206]),
np.mean([1433, 1439, 1442]),
np.mean([1666, 1654, 1654]),
np.mean([1866, 1869, 1878]),
np.mean([2084, 2075, 2072]),
np.mean([2272, 2272, 2290]),
np.mean([2478, 2472, 2463]),
np.mean([2648, 2663, 2666]),
np.mean([2857, 2851, 2842]), ])

# R9ASig = Rot mit B-Feld (9 Ampere), Sigma-Komponente
R9ASig1 = np.array([
np.mean([672, 662, 657]),
np.mean([907, 917, 927]),
np.mean([1148, 1153, 1163]),
np.mean([1374, 1384, 1388]),
np.mean([1600, 1604, 1614]),
np.mean([1815, 1815, 1825]),
np.mean([2022, 2026, 2026]),
np.mean([2227, 2230, 2236]),
np.mean([2418, 2430, 2436]),
np.mean([2612, 2618, 2624]), ])

R9ASig2 = np.array([
np.mean([780, 785, 795]),
np.mean([1035, 1020, 1015]),
np.mean([1271, 1261, 1256]),
np.mean([1492, 1492, 1482]),
np.mean([1712, 1707, 1703]),
np.mean([1919, 1914, 1909]),
np.mean([2125, 2121, 2121]),
np.mean([2324, 2327, 2321]),
np.mean([2527, 2518, 2512]),
np.mean([2709, 2703, 2696]), ])

print('--------------------------------------')
DSR = DS(RoB)
dSR = dS(R9ASig1, R9ASig2)
deltaLambdaR = deltaLambda(dSR, DSR, lamRot)
deltaLambdaR = ufloat(np.mean(deltaLambdaR), np.std(deltaLambdaR))
print('Ds =', DSR)
print('ds1 =', dSR)
#print('Verschiebung mit der roten Linie(9A) =', deltaLambdaR)
print('Mittelwert der Verschiebung =', deltaLambdaR)





# Bestimmung der Lande Faktoren
muB = 9.274 * 10**(-24)
h = 6.626 * 10**(-34)
c = 299792458

print('-------------')
print('rote Linie')
BR = ufloat(0.549, 0.05*0.549) *1.09
gR = (h * c * deltaLambdaR) / (WellenlängeRot**2 * muB * BR)
print('B =', BR)
print('g =', gR)
print('Abweichung:', (1 - gR) / 1)

print('-------------')
print('sigma: blaue Linie')
BB1 = ufloat(0.3, 0.05*0.3) *1.09
gB1 = (h * c * deltaLambda1) / (WellenlängeBlau**2 * muB * BB1)
print('B =', BB1)
print('g =', gB1)
print('Abweichung:', (1.75 - gB1) / 1.75)

print('-------------')
print('pi: blaue Linie')
BB2 = ufloat(1.057, 0.05*1.057) *1.09
gB2 = (h * c * deltaLambda2) / (WellenlängeBlau**2 * muB * BB2)
print('B =', BB2)
print('g =', gB2)
print('Abweichung:', (0.5 - gB2) / 0.5)

print('--------------------------------------')
