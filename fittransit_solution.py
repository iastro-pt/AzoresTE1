import numpy as np
import matplotlib.pyplot as plt
#from scipy.optimize import leastsq
from lmfit import minimize,Parameters,Parameter,report_fit
import batman
import math



def transitmin(params, time, flux, fluxerr):
    import batman
    parbat = batman.TransitParams()
    parbat.t0 = params['t0'].value                      #time of transit
    parbat.per = params['per'].value                    #orbital period
    parbat.rp = params['rp'].value                     #planet radius (in units of stellar radii)
    parbat.a = params['a'].value                     #semi-major axis (in units of stellar radii)
    parbat.inc = params['inc'].value                     #orbital inclination (in degrees)
    parbat.ecc = params['ecc'].value                      #eccentricity
    parbat.w = params['w'].value                       #longitude of periastron (in degrees)
    parbat.u = [params['u1'].value, params['u2'].value]                #limb darkening coefficients
    parbat.limb_dark = "quadratic"       #limb darkening model


    m = batman.TransitModel(parbat, time, supersample_factor = 21, exp_time = 0.5/24.)    #initializes model
    model = m.light_curve(parbat)
    return (model - flux)/fluxerr

#### #### #### #### #### #### #### #### ####
#first exercise plot light curve and model
####


"""Read light curve"""

file = '/Users/sbarros/Documents/work/python/TP_acores/data/ktwo211089792c04_FILlpd_LC.txt'


time, lc, lcerror = np.genfromtxt(file, unpack=True)

#plt.errorbar(time,lc,lcerror)
#plt.show()


#input parameters from the BLS
period = 3.2608374029105653
epoch = 57064.4244992
depth = 0.013977265927714035
duration = 0.07909265190038392




time = time-epoch
flux = lc

#set up iniital Parameters
#estimate a/r and rp/rs
g = 6.67428e-11
gm = 1.32712440041e20
msun = gm/g
rsun = 6.95508e8
sdensity = 1.0
aoverr = (sdensity*(gm*(period*24.*3600.0)**2.)/(4.*math.pi**2.*rsun**3.))**(1./3.)
rp = math.sqrt(depth)

# Assume  limb darkening
u1 = 0.4983
u2 = 0.2042


#set up the parameters for batman
params2 = batman.TransitParams()
params2.t0 = 0.0                       #time of inferior conjunction
params2.per = period                    #orbital period
params2.rp = rp                     #planet radius (in units of stellar radii)
params2.a = aoverr                       #semi-major axis (in units of stellar radii)
params2.inc = 89.                     #orbital inclination (in degrees)
params2.ecc = 0.                      #eccentricity
params2.w = 90.                       #longitude of periastron (in degrees)
params2.u = [u1, u2]                #limb darkening coefficients
params2.limb_dark = "quadratic"       #limb darkening model





#calculate and plot the initial model
m = batman.TransitModel(params2, time, supersample_factor = 21, exp_time = 0.5/24.)    #initializes model
model = m.light_curve(params2)

plt.plot(time, flux, '.')
plt.plot(time,  model, color='r')
plt.xlabel("Time")
plt.ylabel("Relative flux")
plt.show()


#### #### #### #### #### #### #### #### ####
#### Second exercise fit the transit
####


# create a set of parameters for the lmfit
# The ones that are fixed and will not be fitted have vary=false
params = Parameters()
params.add('t0', value = 0.0)
params.add('per', value = period, vary = True)
params.add('rp', value = rp)
params.add('a', value = aoverr)
params.add('inc', value = 88.)
params.add('ecc', value = 0.0, vary = False)
params.add('w', value = 90., vary = False)
params.add('u1', value = u1, vary = False)
params.add('u2', value = u2, vary = False)


# use lmfit to estimate the transit parameters
result = minimize(transitmin,params, args=(time,flux,lcerror ))

report_fit(result)

# recover the fitted parameters from the result
parbat = batman.TransitParams()
parbat.t0 = result.params['t0'].value                     #time of inferior conjunction
parbat.per = result.params['per'].value                    #orbital period
parbat.rp = result.params['rp'].value                      #planet radius (in units of stellar radii)
parbat.a = result.params['a'].value                      #semi-major axis (in units of stellar radii)
parbat.inc = result.params['inc'].value                     #orbital inclination (in degrees)
parbat.ecc = 0.                      #eccentricity
parbat.w = 90.                       #longitude of periastron (in degrees)
parbat.u = [u1, u2]                #limb darkening coefficients
parbat.limb_dark = "quadratic"       #limb darkening model

period = result.params['per'].value
epoch0 = result.params['t0'].value

#phase fold to plot
phase0 = (time )/period
phase = phase0%1
phadur = duration*2.2/period

phase[np.where(phase>0.5)[0]]-=1
cond = np.where(abs(phase) < phadur)[0]

timefold = phase[cond]*period
fluxfold = lc[cond]

fluxfold = fluxfold[np.argsort(timefold)]
timefold = np.sort(timefold)


#calculate and plot the two models initial model for the time folded dataset
m = batman.TransitModel(params2, timefold, supersample_factor = 21, exp_time = 0.5/24.)    #initializes model
model = m.light_curve(params2)


m = batman.TransitModel(parbat, timefold, supersample_factor = 21, exp_time = 0.5/24.)    #initializes model
fitmodel = m.light_curve(parbat)


plt.plot(timefold, fluxfold, '.')
plt.plot(timefold,  model, color='g')
plt.plot(timefold,  fitmodel, color='r')
plt.xlabel("Time")
plt.ylabel("Relative flux")
plt.show()


#### #### #### #### #### #### #### #### ####
#### Cut the transits and write them to a file
####

filename = 'cuttransits.txt'

fp = open(filename, 'w')
for i in cond:
    fp.write('%.10f\t%f\t%f\n'%(time[i]+epoch, flux[i], lcerror[i]))
fp.close()
