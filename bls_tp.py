import bls
import numpy as np
import matplotlib.pyplot as plt

def BLS(time, lc, df, nb, qmi, qma):
    """
    time, lc: arrays
    df = frequency step,
    nb = number of bins,
    qmi = minimum fractional transit duration,
    qma = maximum transit duration,
    fmin = minimum frequency 
    """


    u = np.zeros_like(time)
    v = np.zeros_like(time)

    fmin = 1.3/(time[-1]-time[0])  #searching a maximum period 70% of the full duration
    fmax = 1/0.5 #searching a mininum period of 0.5 days
    nf = (fmax-fmin)/df

    power, best_period, best_power, depth, qtran, in1, in2 = bls.eebls(time, lc, u, v, nf, fmin, df, nb, qmi, qma)
    f = fmin + (np.arange(len(power)))*df
    per = 1./f
    duration = best_period*qtran

    #plt.plot(1/f, power)
    #plt.show()

    
    # power is the nf-dimensional power spectrum array at frequencies
    #  f = fmin + arange(nf) * df
    # best_period is the best-fit period in the same units as time
    # best_power is the power at best_period
    # depth is the depth of the transit at best_period
    # q is the fractional transit duration
    # in1 is the bin index at the start of transit
    # in2 is the bin index at the end of transit

    return power, per, best_period, best_power, depth, duration, in1, in2



# Read light curve
file = 'ktwo211089792c04_FILlpd_LC.txt'
time, lc, lcerror = np.genfromtxt(file, unpack=True)

fig = plt.figure(figsize=(12,6))
ax = fig.add_subplot(111)
ax.errorbar(time, lc, lcerror)
plt.show()


# Calculate BLS
# fmin = minimum frequency , 
# nf of number of frequencies, these are calculated from this and the

df = 0.001  #frequency step
nb = 200    #number of bins
qmi = 0.01  #minimum fractional transit duration
qma = 0.8  #maximum transit duration


power, per, best_period, best_power, depth, duration, in1, in2 = BLS(time, lc, df,  nb, qmi, qma)

print 'best period:', best_period
print 'depth, duration:', depth, duration
# print best_period, best_power, depth, duration
#3.2608374029105653, 0.002150274085459753, 0.013977265927714035, 0.07909265190038392

fig = plt.figure()
ax = fig.add_subplot(111)
ax.semilogx(per, power, 'k-', lw=2)
ax.set_xlabel('Period [days]')
ax.set_ylabel('BLS power')
plt.show()

# phase fold
epoch0 = time[0] + in1/float(nb)*best_period
phase0 = (time-epoch0)/best_period
phase = phase0%1

phase[np.where(phase>0.5)[0]]-=1

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(phase, lc, '.')
ax.set_xlabel('Phase ?')
plt.show()





