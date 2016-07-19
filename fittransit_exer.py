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
    parbat.per =                     #orbital period
    parbat.rp =                    #planet radius (in units of stellar radii)
    parbat.a =                     #semi-major axis (in units of stellar radii)
    parbat.inc =                      #orbital inclination (in degrees)
    parbat.ecc =                       #eccentricity
    parbat.w =                       #longitude of periastron (in degrees)
    parbat.u =                 #limb darkening coefficients
    parbat.limb_dark = "quadratic"       #limb darkening model

    #Initialise the model
    m = batman.TransitModel(parbat, time, supersample_factor = 21, exp_time = 0.5/24.)    #initializes model

    return (model - flux)/fluxerr

#### #### #### #### #### #### #### #### ####
#first exercise plot light curve and model
####


"""Read light curve"""

file =


time, lc, lcerror = np.genfromtxt(file, unpack=True)



#input parameters from the BLS
period =
epoch =
depth =
duration =




time = time-epoch
flux = lc

#set up iniital Parameters
#estimate a/r and rp/rs
g = 6.67428e-11
gm = 1.32712440041e20
msun = gm/g
rsun = 6.95508e8
sdensity = 1.0
aoverr =
rp =

# Assume  limb darkening
u1 = 0.4983
u2 = 0.2042


#set up the parameters for batman
params2 = batman.TransitParams()
params2.t0 = 0.0                       #time of inferior conjunction
params2.per = period                    #orbital period
params2.rp = rp                     #planet radius (in units of stellar radii)
params2.a = aoverr                       #semi-major axis (in units of stellar radii)
params2.inc =                     #orbital inclination (in degrees)
params2.ecc = 0.                      #eccentricity
params2.w = 90.                       #longitude of periastron (in degrees)
params2.u = [u1, u2]                #limb darkening coefficients
params2.limb_dark = "quadratic"       #limb darkening model





#calculate and plot the initial model
m = batman.TransitModel(params2, time, supersample_factor = 21, exp_time = 0.5/24.)    #initializes model
model = m.light_curve(params2)



#### #### #### #### #### #### #### #### ####
#### Second exercise fit the transit
####


# create a set of parameters for the lmfit
# The ones that are fixed and will not be fitted have vary=false
params = Parameters()
params.add('t0', value =)
params.add('per', value = period, vary = True)
params.add('rp', value = rp)
params.add('a', value = aoverr)
params.add('inc', value =  )
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
parbat.per =                    #orbital period
parbat.rp =                     #planet radius (in units of stellar radii)
parbat.a =                     #semi-major axis (in units of stellar radii)
parbat.inc = r                   #orbital inclination (in degrees)
parbat.ecc = 0.                      #eccentricity
parbat.w = 90.                       #longitude of periastron (in degrees)
parbat.u = [u1, u2]                #limb darkening coefficients
parbat.limb_dark = "quadratic"       #limb darkening model

period =


#phase fold to plot


phadur = duration*2.2/period




#calculate and plot the two models initial model for the time folded dataset
m = batman.TransitModel(params2, timefold, supersample_factor = 21, exp_time = 0.5/24.)    #initializes model
model = m.light_curve(params2)


m = batman.TransitModel(parbat, timefold, supersample_factor = 21, exp_time = 0.5/24.)    #initializes model
fitmodel = m.light_curve(parbat)





#### #### #### #### #### #### #### #### ####
#### Cut the transits and write them to a file
####

filename = 'cuttransits.txt'

fp = open(filename, 'w')
for i in  :
    fp.write('%.10f\t%f\t%f\n'%(time[i]+epoch, flux[i], lcerror[i]))
fp.close()
