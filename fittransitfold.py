import numpy as np
import matplotlib.pyplot as plt
#from scipy.optimize import leastsq
from lmfit import minimize,Parameters,Parameter,report_fit
import batman
import math



def transitmin(params, time, flux):
    import batman
    parbat = batman.TransitParams()
    parbat.t0 = params['t0'].value                      #time of inferior conjunction
    parbat.per = params['per'].value                    #orbital period
    parbat.rp = params['rp'].value                     #planet radius (in units of stellar radii)
    parbat.a = params['a'].value                     #semi-major axis (in units of stellar radii)
    parbat.inc = params['inc'].value                     #orbital inclination (in degrees)
    parbat.ecc = params['ecc'].value                      #eccentricity
    parbat.w = params['w'].value                       #longitude of periastron (in degrees)
    parbat.u = [params['u1'].value, params['u2'].value]                #limb darkening coefficients
    parbat.limb_dark = "quadratic"       #limb darkening model


    m = batman.TransitModel(parbat, time, supersample_factor=21, exp_time=0.5/24.)    #initializes model
    model = m.light_curve(parbat)
    return model-flux



# Read light curve
file = 'ktwo211089792c04_FILlpd_LC.txt'
time, lc, lcerror = np.genfromtxt(file, unpack=True)

#plt.errorbar(time,lc,lcerror)
#plt.show()


#input parameters from the BLS
period = 3.2608374029105653
epoch = 57064.391890782324
depth = 0.013977265927714035
duration = 0.07909265190038392


# time fold and cut close to transits 3 transitdurations
phase0 = (time-epoch)/period
phase = phase0%1
phadur = duration*2.2/period

phase[np.where(phase>0.5)[0]] -= 1
cond = np.where(abs(phase) < phadur)[0]

time = phase[cond]*period
flux = lc[cond]

flux = flux[np.argsort(time)]
time = np.sort(time)

plt.plot(time,flux, '.')
plt.show()


#set up iniital Parameters
#estimate a/r
g = 6.67428e-11
gm = 1.32712440041e20
msun = gm/g
rsun = 6.95508e8
sdensity=1.06390
aoverr = (sdensity*(gm*(period*24.*3600.0)**2.)/(4.*math.pi**2.*rsun**3.))**(1./3.)

# estimate the limb darkening
teff = 5358.0
dteff = 38.0
feh = 0.16
dfeh = 0.03
logg = 4.540
dlogg = 0.012



params = batman.TransitParams()
params.t0 = 0.0                      #time of inferior conjunction
params.per = period                  #orbital period
params.rp = math.sqrt(depth)         #planet radius (in units of stellar radii)
params.a = aoverr                    #semi-major axis (in units of stellar radii)
params.inc = 88.                     #orbital inclination (in degrees)
params.ecc = 0.                      #eccentricity
params.w = 90.                       #longitude of periastron (in degrees)
params.u = [0.1, 0.3]                #limb darkening coefficients
params.limb_dark = "quadratic"       #limb darkening model





#calculate and plot the initial model
m = batman.TransitModel(params, time, supersample_factor=21, exp_time=0.5/24.)    #initializes model
model = m.light_curve(params)

plt.plot(time, flux, '.')
plt.plot(time,  model, color='r')
plt.xlabel("Time")
plt.ylabel("Relative flux")
plt.show()


# create a set of parameters
params=Parameters()
params.add('t0', value=0.0)
params.add('per', value=period, vary=False)
params.add('rp', value=math.sqrt(depth) )
params.add('a', value=aoverr)
params.add('inc', value=88.)
params.add('ecc', value=0.0, vary=False)
params.add('w', value=90., vary=False)
params.add('u1', value=0.1, vary=False)
params.add('u2', value=0.3, vary=False)



result = minimize(transitmin, params, args=(time,flux))

report_fit(result)


parbat = batman.TransitParams()
parbat.t0 = result.params['t0'].value                     #time of inferior conjunction
parbat.per = period                    #orbital period
parbat.rp = result.params['rp'].value                      #planet radius (in units of stellar radii)
parbat.a = result.params['a'].value                      #semi-major axis (in units of stellar radii)
parbat.inc = result.params['inc'].value                     #orbital inclination (in degrees)
parbat.ecc = 0.                      #eccentricity
parbat.w = 90.                       #longitude of periastron (in degrees)
parbat.u = [0.1, 0.3]                #limb darkening coefficients
parbat.limb_dark = "quadratic"       #limb darkening model



#calculate and plot the initial model
m = batman.TransitModel(parbat, time, supersample_factor=21, exp_time=0.5/24.)    #initializes model
fitmodel = m.light_curve(parbat)


plt.plot(time, flux, '.')
plt.plot(time, model, color='g')
plt.plot(time, fitmodel, color='r')
plt.xlabel("Time")
plt.ylabel("Relative flux")
plt.show()
