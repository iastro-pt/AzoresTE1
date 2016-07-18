import bls
import numpy as np
import matplotlib.pyplot as plt



# Read light curve
file = 'ktwo211089792c04_FILlpd_LC.txt'
time, lc, lcerror = np.genfromtxt(file, unpack=True)





# Calculate BLS
#bls.eebls(time, lc, u, v, nf, fmin, df, nb, qmi, qma)

"""
Input variables
time, lc: N-dimensional data arrays
U, V N-dimensional empty work arrays
df = frequency step,
nb = number of bins for folded light curve,
qmi = minimum fractional transit duration,
qma = maximum transit duration,
fmin = minimum frequency
"""

#Define input variables
u = np.zeros_like(time)
v = np.zeros_like(time)

df =   #frequency step
nb =    #number of bins
qmi =   #minimum fractional transit duration
qma =   #maximum transit duration

# frequencies to test
fmin =   #searching a maximum period 70% of the full duration
fmax =  #searching a mininum period of 0.5 days
nf = (fmax-fmin)/df
# frequency array
f =

power, best_period, best_power, depth, qtran, in1, in2 = bls.eebls(time, lc, u, v, nf, fmin, df, nb, qmi, qma)

per = 1./f
duration = best_period*qtran
epoch0 =


# power is the nf-dimensional power spectrum array at the trial frequencies
# best_period is the best-fit period in the same units as time
# best_power is the power at best_period
# depth is the depth of the transit at best_period
# q is the fractional transit duration
# in1 is the bin index at the start of transit
# in2 is the bin index at the end of transit

per = 1./f
duration = best_period*qtran

#plot the periodgram
fig = plt.figure()
ax = fig.add_subplot(111)
ax.semilogx(per, power, 'k-', lw=2)
ax.set_xlabel('Period [days]')
ax.set_ylabel('BLS power')
plt.show()



# print the results for the best period



# phase fold the light curve on the best period
